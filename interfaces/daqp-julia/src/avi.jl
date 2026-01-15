Base.@kwdef mutable struct AVISettings 
    iter_limit::Int = 1000
    inner_iter_limit::Int = 10000
    primal_tol::Float64 = 1e-6
    dual_tol::Float64 = 1e-12
    alpha::Float64 = 1.0
    beta::Float64 = 0.25
end

struct AVIWorkspace
    H::Matrix{Float64}
    f::Vector{Float64}
    A::Matrix{Float64}
    bu::Vector{Float64}
    bl::Vector{Float64}

    H_factor::LU{Float64}
    x_unc::Vector{Float64}

    H1pI::Matrix{Float64}
    H2mI::Matrix{Float64}
    H2pI::LU{Float64} # LU
    H2::Matrix{Float64}
    H2xpf::Vector{Float64}
    x::Vector{Float64}
    y::Vector{Float64}
    daqp_workspace::DAQPBase.Model
    kkt_buffer::Vector{Float64}

    settings::AVISettings
end

AVIWorkspace() = AVIWorkspace(zeros(0,0),zeros(0),zeros(0,0),zeros(0),zeros(0), 
                              lu(zeros(0,0)),zeros(0),
                              zeros(0,0),zeros(0,0),lu(zeros(0,0)),zeros(0,0), zeros(0),
                              zeros(0),zeros(0), Model(), zeros(0), AVISettings())


function setup_avi(H,f,A,bu,bl,sense;x0 = zeros(0), daqp_workspace=nothing, settings=AVISettings())
    n,m = length(f),length(bu)
    bl = isempty(bl) ? fill(-1e30,m) : bl
    sense = isempty(sense) ? zeros(Cint,m) : sense

    H1 = (H+H') / 4
    H2 = H1+(H-H') / 2
    H1pI,H2mI,H2pI = H1+I, H2-I, lu!(H2+I)

    #H_factor = lu(H)
    H_factor = lu(zeros(0,0))

    x = isempty(x0) ? fill(NaN,n) : copy(x0)
    y,H2xpf = zeros(n),zeros(n)

    kkt_buffer = Float64[]
    sizehint!(kkt_buffer,4*(n+1)*(n+1))

    daqp_workspace = isnothing(daqp_workspace) ? Model() : daqp_workspace;
    DAQPBase.settings(daqp_workspace, Dict(:iter_limit => settings.inner_iter_limit, 
                                           :primal_tol => settings.primal_tol, 
                                           :dual_tol => settings.dual_tol))
    exitflag,tsetup = DAQPBase.setup(daqp_workspace, H1pI, H2xpf, A, bu, bl, sense)

    if m > size(A,1)
        A = [I(n)[1:m-size(A,1),:];A] # To handle simple bounds
    end

    return exitflag, AVIWorkspace(H,f,A,bu,bl,
                                  H_factor,zeros(n),
                                  H1pI,H2mI,H2pI,H2,H2xpf,
                                  x,y,daqp_workspace,kkt_buffer,
                                  settings)
end

function solve(ws::AVIWorkspace)
    exitflag = 0
    λstar, ASstar = zeros(0), Int[]
    n = length(ws.f)
    ϵp, ϵd = ws.settings.primal_tol, ws.settings.dual_tol
    α,β = ws.settings.alpha, ws.settings.beta
    tot_iter,outer_iter,nkkt = 0,0,0
    isnan(ws.x[1]) && (ws.x .= -ws.H\ws.f) # Use unconstrained optimum as x0
    @inbounds for k in 1:ws.settings.iter_limit
        ws.H2xpf .= ws.f
        mul!(ws.H2xpf, ws.H2mI, ws.x, 1.0, 1.0)
        DAQPBase.update(ws.daqp_workspace, nothing, ws.H2xpf, nothing, nothing,nothing)
        ws.y[:], _, exitflag, info = DAQPBase.solve(ws.daqp_workspace)
        exitflag < 0 && break 
        tot_iter += info.iterations

        # Same AS -> Check if KKT conditions are satisfied
        if info.iterations == 1
            nkkt += 1
            ASu,ASl = _get_AS(info.λ,ϵd) 
            #xt,λ,AS = _solve_kkt_schur(ws,ASu,ASl)
            xt,λ,AS = _solve_kkt(ws,ASu,ASl)
            nu = length(ASu)
            if all(λ[i] > -ϵd for i in 1:nu) && all(λ[i] < ϵd for i in nu+1:length(λ))
                Ax = ws.A*xt
                if all(Ax-ws.bu .<= ϵp) && all(Ax-ws.bl .>=-ϵp)
                    ws.x .= xt
                    λstar,ASstar = copy(λ),AS
                    exitflag, outer_iter = 1,k
                    break
                end
            end
            axpby!(β,xt,1.0-β,ws.x)
        else
            axpby!(1.0-α,ws.x,α,ws.y)
            mul!(ws.y, ws.H2, ws.x, 1.0, 1.0)
            ldiv!(ws.x, ws.H2pI, ws.y)
        end
        k == ws.settings.iter_limit && (exitflag = -4)
    end
    info = (status=flag2status[exitflag], AS=ASstar, outer_iterations=outer_iter, iterations=tot_iter, nkkt=nkkt)
    return ws.x,λstar,exitflag,info
end

function _get_AS(λ,dual_tol)
    ASu,ASl = Int[],Int[]
    for i in 1:length(λ)
        if λ[i] > dual_tol 
            push!(ASu,i)
        elseif λ[i] < -dual_tol
            push!(ASl,i)
        end
    end
    return ASu,ASl
end

# TODO kkt_schur is faster in isolation, but leads to slower total execution 
#function _solve_kkt_schur(ws::AVIWorkspace,ASu,ASl)
#    # Prepare some helpers
#    AS = [ASu;ASl]
#    n,nAS = length(ws.f),length(AS)
#
#    # Juggle some buffers
#    resize!(ws.kkt_buffer,nAS*(nAS+n+1)+n) # Is this necessary?
#    invH_At = reshape(view(ws.kkt_buffer,1:nAS*n),(n,nAS))
#    offset = nAS*n 
#    S = reshape(view(ws.kkt_buffer,offset+1:offset+nAS^2),(nAS,nAS))
#    offset +=nAS^2
#    λ = view(ws.kkt_buffer,offset+1:offset+nAS)
#    offset +=nAS
#    x = view(ws.kkt_buffer,offset+1:offset+n)
#
#    AAS = ws.A[AS,:]
#
#    # Schur complement 
#    @inbounds for (i,id) in enumerate(AS)
#        @views invH_At[:,i] .= ws.A[id,:]
#    end
#    ldiv!(ws.H_factor,invH_At)
#    mul!(S, AAS, invH_At)
#    S_fact = lu!(S) 
#
#    # Solve for λ
#    @views λ .= AAS * ws.x_unc
#    @views λ[1:length(ASu)] .-= ws.bu[ASu]
#    @views λ[length(ASu)+1:end] .-= ws.bl[ASl]
#    ldiv!(S_fact,λ)
#
#    # Solve for x
#    x .= ws.x_unc .- invH_At * λ
#
#    return x,λ,AS 
#end
function _solve_kkt(ws,ASu,ASl)
    AS = [ASu;ASl]
    n = length(ws.f)
    nkkt = n+length(AS)
    resize!(ws.kkt_buffer,nkkt^2)
    K = reshape(view(ws.kkt_buffer,1:nkkt^2),(nkkt,nkkt))
    @views K[1:n,1:n] = ws.H
    @views K[n+1:end,1:n] = ws.A[AS,:] 
    @views K[1:n,n+1:end] = ws.A[AS,:]'
    K[n+1:end,n+1:end] .= 0
    z = K\[-ws.f;ws.bu[ASu];ws.bl[ASl]]
    x = @view z[1:n]
    λ = @view z[n+1:end]
    return x,λ,AS
end

"""
    x, λ, info = DAQPBase.solve_avi(H,f,A,b)
    x, λ, info = DAQPBase.solve_avi(H,f,A,bupper,blower)

finds a solution `x` to the affine variational inequality 

```
x ∈ {x: blower ≤ Ax ≤ bupper} such that: 
(Hx+f)'(x-y) ≥ 0 ∀y ∈ {y: blower ≤ Ay ≤ bupper}
```

# Input 
* `H`  			- linear term in objective function, (`n x n`)-matrix 
* `f`  			- linear term in objective function, n-vector 
* `A`  			- constraint normals, (`m x n`)-matrix
* `bupper` 		- upper bounds for constraints, `m`-vector
* `blower` 		- lower bounds for constraints, `m`-vector (default: -Inf)

# Output
* `x` 		- solution provided by solver
* `λ` 	  	    - dual variables 
* `info` 		- tuple containing extra information from the solver such as status, active set and number of iterations.

"""
function solve_avi(H::AbstractMatrix,f::AbstractVector,A::AbstractMatrix,bupper::AbstractVector,blower::AbstractVector=Float64[], sense=Cint[]; x0 = zeros(0), settings=AVISettings())
    exitflag, ws = setup_avi(H,f,A,bupper,blower,sense;x0,settings)
    return solve(ws)
end
