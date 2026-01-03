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
function solve_avi(H::AbstractMatrix,f::AbstractVector,A::AbstractMatrix,bupper::AbstractVector,blower::AbstractVector=Float64[];
        alpha=1.0,beta=0.25,primal_tol=1e-6, dual_tol = 1e-12, x0 = zeros(0), iter_limit=1000, inner_iter_limit=10000)

    # Solves VI(Hx +h, C) where C={x| blower <= Ax <= bupper} using Douglas-Rachford
    n,m = length(f),length(bupper)
    blower = isempty(blower) ? fill(-1e30,m) : blower

    H1 = (H+H') / 4
    H2 = H1+(H-H') / 2
    H1pI,H2mI,H2pI = H1+I, H2-I, lu!(H2+I)

    x = isempty(x0) ? zeros(n) : copy(x0)
    y,H2xpf = zeros(n),zeros(n),zeros(n)

    backstep = DAQPBase.Model()
    DAQPBase.settings(backstep, Dict(:iter_limit =>inner_iter_limit, :primal_tol=>primal_tol, :dual_tol=>dual_tol))
    exitflag,tsetup = DAQPBase.setup(backstep, H1pI, H2xpf, A, bupper, blower)
    exitflag < 0 && return nothing,nothing, (status=:Infeasible, exitflag=exitflag, outer_iterations=0, iterations=0) 

    tot_iter = 0

    kkt_buffer = Float64[]
    sizehint!(kkt_buffer,4*(n+1)*(n+1))

    @inbounds for k in 1:iter_limit
        H2xpf .= f
        mul!(H2xpf, H2mI, x, 1.0, 1.0)
        DAQPBase.update(backstep, nothing, H2xpf, nothing, nothing,nothing)

        y[:], _, exitflag, info = DAQPBase.solve(backstep)
        exitflag < 0 && return nothing,nothing, (status=:Infeasible, exitflag=exitflag, outer_iterations=k, iterations=tot_iter) 
        tot_iter += info.iterations

        # Same AS -> Check if KKT conditions are satisfied
        if info.iterations == 1
            ASu,ASl = _get_AS(info.λ,dual_tol) 
            z,AS = _solve_kkt(H,f,A,bupper,blower,ASu,ASl,kkt_buffer)
            xt = @view z[1:n]
            λ = @view z[n+1:end]
            if all(λ[i] > -dual_tol for i in 1:length(ASu)) && all(λ[i] < dual_tol for i in length(ASu)+1:length(λ))
                Ax = A*xt
                if all(Ax-bupper .<= primal_tol) && all(Ax-blower .>=-primal_tol)
                    return copy(xt),copy(λ),(AS=AS, status=:Solved, outer_iterations=k,iterations=tot_iter)
                end
            end
            axpby!(beta,xt,1.0-beta,x)
        else
            axpby!(1.0-alpha,x,alpha,y)
            mul!(y, H2, x, 1.0, 1.0)
            ldiv!(x, H2pI, y)
        end
    end
    return x,nothing, (status=:MaximumIterations, outer_iterations=iter_limit, iterations=tot_iter)
end

