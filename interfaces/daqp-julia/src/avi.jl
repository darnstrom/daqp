Base.@kwdef mutable struct AVISettings 
    iter_limit::Int = 1000
    inner_iter_limit::Int = 10000
    primal_tol::Float64 = 1e-6
    dual_tol::Float64 = 1e-12
    alpha::Float64 = 1.0
    beta::Float64 = 0.25
    regularization::Float64 = 0.5 
    max_terminate_counter::Int = 5
    min_terminate_counter::Int = 5
end

struct AVIWorkspace
    H::Matrix{Float64}
    f::Vector{Float64}
    At::Matrix{Float64}
    bu::Vector{Float64}
    bl::Vector{Float64}
    sense::Vector{Cint}
    norm_factors::Vector{Cdouble}

    rho_soft::Float64

    H_factor::LU{Float64, Matrix{Float64}, Vector{Int64}}
    x_unc::Vector{Float64}

    H1pI::Matrix{Float64}
    H2pI::LU{Float64, Matrix{Float64}, Vector{Int64}}
    Hsym::Matrix{Float64}
    xtemp::Vector{Float64}
    Hx::Vector{Float64}
    x::Vector{Float64}
    y::Vector{Float64}
    daqp_workspace::DAQPBase.Model
    kkt_buffer::Vector{Float64}

    settings::AVISettings
end

AVIWorkspace() = AVIWorkspace(zeros(0,0),zeros(0),zeros(0,0),zeros(0),zeros(0),zeros(Cint,0),zeros(0),
                              0,lu(zeros(0,0)),zeros(0),
                              zeros(0,0),lu(zeros(0,0)),zeros(0,0), zeros(0), zeros(0),
                              zeros(0),zeros(0), Model(), zeros(0), AVISettings())


function setup_avi(H,f,A,bu,bl,sense;x0 = zeros(0), rho_soft = 1e6, daqp_workspace=nothing, settings=AVISettings())
    n,m = length(f),length(bu)
    bl = isempty(bl) ? fill(-1e30,m) : bl
    sense = isempty(sense) ? zeros(Cint,m) : sense

    Hsym = (H+H')./2
    H1pI = Hsym+settings.regularization*I
    H2pI = lu!(H+settings.regularization*I)

    #H_factor = lu(H)
    H_factor = lu(zeros(0,0))

    x = isempty(x0) ? fill(NaN,n) : copy(x0)
    y,Hx,xtemp = zeros(n),zeros(n),zeros(n)

    kkt_buffer = zeros(4*(n+1)*(n+1))
    #sizehint!(kkt_buffer,4*(n+1)*(n+1))

    daqp_workspace = isnothing(daqp_workspace) ? Model() : daqp_workspace;
    DAQPBase.settings(daqp_workspace, Dict(:iter_limit => settings.inner_iter_limit, 
                                           :primal_tol => settings.primal_tol, 
                                           :dual_tol => settings.dual_tol))
    exitflag,tsetup = DAQPBase.setup(daqp_workspace, H1pI, xtemp, A, bu, bl, sense)

    workspace = unsafe_load(daqp_workspace.work);
    norm_factors = unsafe_wrap(Vector{Cdouble}, workspace.scaling, m; own = false)


    if m > size(A,1)
        #A = [I(n)[1:m-size(A,1),:];A] # To handle simple bounds
        At = [I(n)[:,1:m-size(A,1)] A'] # To handle simple bounds
    else
        At = A'[:,:]
    end

    return exitflag, AVIWorkspace(H,f,At,bu,bl,sense,norm_factors,
                                  rho_soft,H_factor,zeros(n),
                                  H1pI,H2pI,Hsym,xtemp,Hx,
                                  x,y,daqp_workspace,kkt_buffer,
                                  settings)
end

function _is_optimal(λ,xt,λfull,ws,ASu,ASl)
    nu = length(ASu)
    ϵp, ϵd = ws.settings.primal_tol, ws.settings.dual_tol
    @inbounds for i in 1:nu
        λ[i] < -ϵd && return false
    end
    @inbounds for i in nu+1:length(λ)
        λ[i] > ϵd && return false
    end
    @inbounds for i in 1:length(ws.bu)
        if λfull[i] == 0.0 # Only test inactive constraints
            Axi = dot(view(ws.At,:,i),xt)
            Axi-ws.bu[i] > ϵp && return false
            Axi-ws.bl[i] < -ϵp && return false
        end
    end
    return true
end

function solve(ws::AVIWorkspace)
    exitflag = 0
    λstar, ASstar = zeros(0), Int[]
    n = length(ws.f)
    ϵp, ϵd = ws.settings.primal_tol, ws.settings.dual_tol
    α,β = ws.settings.alpha, ws.settings.beta
    tot_iter,outer_iter,nkkt = 0,0,0
    counter,terminate_limit = 0,ws.settings.min_terminate_counter
    #isnan(ws.x[1]) && (ws.x .= -ws.H\ws.f) # Use unconstrained optimum as x0
    isnan(ws.x[1]) && (ws.x .= zeros(n)) # Use unconstrained optimum as x0
    @inbounds for k in 1:ws.settings.iter_limit
        mul!(ws.Hx,ws.H,ws.x)
        ws.xtemp .= ws.f .+ ws.Hx
        mul!(ws.xtemp,ws.H1pI,ws.x,-1.0,1.0)

        DAQPBase.update(ws.daqp_workspace, nothing, ws.xtemp, nothing, nothing,nothing)
        ws.y[:], _, exitflag, info = DAQPBase.solve(ws.daqp_workspace)
        exitflag < 0 && break 
        tot_iter += info.iterations

        # Same AS -> Check if KKT conditions are satisfied
        if info.iterations == 1 
            if (counter += 1) == terminate_limit 
                nkkt += 1
                ASu,ASl = _get_AS(info.λ,ϵd) 
                xt,λ,AS = _solve_kkt(ws,ASu,ASl)
                if _is_optimal(λ,xt,info.λ,ws,ASu,ASl)
                    ws.x .= xt
                    λstar,ASstar = copy(λ),AS
                    exitflag, outer_iter = 1,k
                    break
                end
                ws.x .= xt
                continue
            end
        else
            counter = 0
        end

        #axpby!(1.0-α,ws.x,α,ws.y)
        ws.xtemp .= ws.y
        ws.y .-= ws.x
        mul!(ws.xtemp, ws.Hsym, ws.y, 0.5, ws.settings.regularization)
        ws.xtemp .+= ws.Hx
        ldiv!(ws.x, ws.H2pI, ws.xtemp)
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

function _solve_kkt(ws,ASu,ASl)
    AS = Int[ASu;ASl]
    n = length(ws.f)
    nkkt = n+length(AS)
    resize!(ws.kkt_buffer,nkkt^2)
    K = reshape(view(ws.kkt_buffer,1:nkkt^2),(nkkt,nkkt))
    @views K[1:n,1:n] = ws.H
    @views K[n+1:end,1:n] = ws.At[:,AS]'
    @views K[1:n,n+1:end] = ws.At[:,AS]
    K[n+1:end,n+1:end] .= 0
    for (i,id) in enumerate(AS)
        if ws.sense[id] & DAQPBase.SOFT != 0
            K[n+i,n+i] = -1/((ws.norm_factors[id]^2)*ws.rho_soft)
        end
    end
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
