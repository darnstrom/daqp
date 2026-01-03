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


function setup_avi(H,f,A,bu,bl;x0 = zeros(0), daqp_workspace=nothing, settings=AVISettings())
    n,m = length(f),length(bu)
    bl = isempty(bl) ? fill(-1e30,m) : bl

    H1 = (H+H') / 4
    H2 = H1+(H-H') / 2
    H1pI,H2mI,H2pI = H1+I, H2-I, lu!(H2+I)

    x = isempty(x0) ? zeros(n) : copy(x0)
    y,H2xpf = zeros(n),zeros(n)

    kkt_buffer = Float64[]
    sizehint!(kkt_buffer,4*(n+1)*(n+1))

    daqp_workspace = isnothing(daqp_workspace) ? Model() : daqp_workspace;
    DAQPBase.settings(daqp_workspace, Dict(:iter_limit => settings.inner_iter_limit, 
                                           :primal_tol => settings.primal_tol, 
                                           :dual_tol => settings.dual_tol))
    exitflag,tsetup = DAQPBase.setup(daqp_workspace, H1pI, H2xpf, A, bu, bl)

    return exitflag, AVIWorkspace(H,f,A,bu,bl,
                                  H1pI,H2mI,H2pI,H2,H2xpf,
                                  x,y,daqp_workspace,kkt_buffer,
                                  settings)
end

function solve(ws::AVIWorkspace)
    exitflag = -4
    λstar, ASstar = zeros(0), Int[]
    n = length(ws.f)
    ϵp, ϵd = ws.settings.primal_tol, ws.settings.dual_tol
    α,β = ws.settings.alpha, ws.settings.beta
    tot_iter,outer_iter = 0,0
    @inbounds for k in 1:ws.settings.iter_limit
        ws.H2xpf .= ws.f
        mul!(ws.H2xpf, ws.H2mI, ws.x, 1.0, 1.0)
        DAQPBase.update(ws.daqp_workspace, nothing, ws.H2xpf, nothing, nothing,nothing)
        ws.y[:], _, exitflag, info = DAQPBase.solve(ws.daqp_workspace)
        exitflag < 0 && break 
        tot_iter += info.iterations

        # Same AS -> Check if KKT conditions are satisfied
        if info.iterations == 1
            ASu,ASl = _get_AS(info.λ,ϵd) 
            z,AS = _solve_kkt(ws.H,ws.f,ws.A,ws.bu,ws.bl,ASu,ASl,ws.kkt_buffer)
            xt = @view z[1:n]
            λ = @view z[n+1:end]
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
    end
    info = (status=flag2status[exitflag], AS=ASstar, outer_iterations=outer_iter, iterations=tot_iter)
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

function _solve_kkt(H,f,A,bu,bl,ASu,ASl,kkt_buffer)
    AS = [ASu;ASl]
    n = length(f)
    nkkt = n+length(AS)
    resize!(kkt_buffer,nkkt^2)
    K = reshape(view(kkt_buffer,1:nkkt^2),(nkkt,nkkt))
    @views K[1:n,1:n] = H
    @views K[n+1:end,1:n] = A[AS,:] 
    @views K[1:n,n+1:end] = A[AS,:]'
    K[n+1:end,n+1:end] .= 0
    return K\[-f;bu[ASu];bl[ASl]], AS
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
function solve_avi(H::AbstractMatrix,f::AbstractVector,A::AbstractMatrix,bupper::AbstractVector,blower::AbstractVector=Float64[]; x0 = zeros(0), settings=AVISettings())
    exitflag, ws = setup_avi(H,f,A,bupper,blower;x0=zeros(0),settings)
    return solve(ws)
end
