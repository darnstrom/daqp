"""
# Example calls
    xstar, fval, exitflag, info = DAQPBase.quadprog(H,f,A,bupper)
    xstar, fval, exitflag, info = DAQPBase.quadprog(H,f,A,bupper,blower,sense) 

finds the solution `xstar` to the quadratic program

```
 minimize	   0.5 x' H x + f' x
subject to   blower <= A x <= bupper
```
If `bupper` and `blower` have more elements than rows of `A`, the first
elements are interpreted as simple bounds. For example:

```
    A = [7.0 8.0]
    blower = [-4.0; -5.0; -6.0]
    bupper = [ 1.0;  2.0;  3.0]
```
is interpreted as

```
        -4.0 <= x₁ <= 1.0
        -5.0 <= x₂ <= 2.0
    -6.0 <= 7 x₁ + 8 x₂ <= 3.0

```

# Input 
* `H`           - cost matrix
* `f`           - cost vector
* `A`           - linear constraint matrix
* `buppe`       - upper bounds for constraints
* `blower`      - lower bounds for constraints (default: -Inf)
* `sense`       - constraint types, as a vector of Cints (default: 0). Example types:
  * `0 ` : inequality
  * `1 ` : active inequality (used as warm start)
  * `5 ` : equality
  * `8 ` : soft (allowed to be violated if necessary)
  * `16` : binary (either upper or lower bound should hold with equality)

# Output
* `xstar`       - solution
* `fval`        - objective function value for `xstar`. 
* `exitflag`    - flag from solver (>0 success, <0 failure) 
* `info`        - tuple containing profiling information from the solver. 

"""
function quadprog(H::Matrix{Float64},f::Vector{Float64}, 
        A::Matrix{Float64},bupper::Vector{Float64},blower::Vector{Float64}=Float64[],sense::Vector{Cint}=Cint[];A_rowmaj=false,settings=nothing)
    return quadprog(QPj(H,f,A,bupper,blower,sense;A_rowmaj);settings)
end
function quadprog(qpj::QPj;settings=nothing)
    # TODO: check validity of dimensions
    # Setup QP
    qp = QPc(qpj);

    # Setup output struct
    xstar = zeros(Float64,qp.n); 
    lam= zeros(Float64,qp.m); 
    result= Ref(DAQPResult(xstar,lam));
    ptr_settings = settings isa DAQPSettings ? Ref(settings) : Ptr{DAQPBase.DAQPSettings}(C_NULL)


    ccall((:daqp_quadprog, DAQPBase.libdaqp), Nothing,
          (Ref{DAQPBase.DAQPResult},Ref{DAQPBase.QPc},Ref{DAQPBase.DAQPSettings}), 
          result,Ref(qp),ptr_settings)

    info = (x = xstar, λ=lam, fval=result[].fval,
            exitflag=result[].exitflag,
            status = DAQPBase.flag2status[result[].exitflag],
            solve_time = result[].solve_time,
            setup_time = result[].setup_time,
            iterations= result[].iter, nodes = result[].nodes)
    return xstar,result[].fval,result[].exitflag,info
end

"""
    xstar, fval, exitflag, info = DAQPBase.linprog(f,A,bupper)
    xstar, fval, exitflag, info = DAQPBase.linprog(f,A,bupper,blower,sense)

finds the solution `xstar` to the linear program

```
min_x	f' x
subject to 
    blower[1:ms]	<= x[1:ms] <= bupper[1:ms]
    blower[ms+1:m]  <= A*x 	   <= bupper[ms+1:m],
```
where `m = length(bupper)` and `ms = m-size(A,2)`.

# Input 
* `f`  			- linear term in objective function, n-vector 
* `A`  			- constraint normals, (`(m-ms) x n`)-matrix
* `bupper` 		- upper bounds for constraints, `m`-vector
* `blower` 		- lower bounds for constraints, `m`-vector (default: -Inf)
* `sense` 		- constraint types,  `m`-vector of Cints (default: 0). Example types:
  * `0` : inequality
  * `5` : equality

# Output
* `xstar` 		- solution provided by solver
* `fval` 		- objective function value for `xstar`. 
* `exitflag` 	- flag from solver (>0 success, <0 failure) 
* `info` 		- tuple containing profiling information from the solver. 

"""
function linprog(f::Vector{Float64}, 
        A::Matrix{Float64},bupper::Vector{Float64},blower::Vector{Float64}=Float64[],sense::Vector{Cint}=Cint[];A_rowmaj=false)
    d = DAQPBase.Model() 
    DAQPBase.setup(d,QPj(zeros(0,0),f,A,bupper,blower,sense;A_rowmaj))
    return DAQPBase.solve(d);
end
"""
    d = DAQPBase.Model() 
creates an empty optimization model `d`. 


The following functions acts on such models: 

* `setup(d,H,f,A,bupper,blower,sense)`: setup a QP problem (see `DAQPBase.quadprog` for problem details)
* `solve(d)`: solve a populated model
* `update(d,H,f,A,bupper,blower,sense)`: update an existing model 
* `dict = DAQPBase.settings(d)`: return a Dictionary with the current settings for the model `d` 
* `settings(d,dict)`: update the settings for `d` with the Dictionary `dict` 
"""
mutable struct Model 
    work::Ptr{DAQPBase.Workspace}
    qpj::QPj
    qpc::QPc
    qpc_ptr::Ptr{DAQPBase.QPc}
    has_model::Bool
    x::Vector{Float64}
    λ::Vector{Float64}
    function Model()
        # Setup initial model
        work = Libc.calloc(1,sizeof(DAQPBase.Workspace))
        daqp= new(Ptr{DAQPBase.Workspace}(work))
        daqp.qpc_ptr = Libc.calloc(1,sizeof(DAQPBase.QPc))
        ccall((:allocate_daqp_settings,DAQPBase.libdaqp),Nothing,(Ptr{DAQPBase.Workspace},),work)
        finalizer(DAQPBase.delete!, daqp)
        daqp.has_model=false
        return daqp 
    end
end


function delete!(daqp::DAQPBase.Model)
    if(daqp.work != C_NULL)
        ccall((:free_daqp_workspace,DAQPBase.libdaqp),Nothing,(Ptr{DAQPBase.Workspace},),daqp.work)
        ccall((:free_daqp_ldp,DAQPBase.libdaqp),Nothing,(Ptr{DAQPBase.Workspace},),daqp.work)
        Libc.free(daqp.work);
        daqp.work = C_NULL
        Libc.free(daqp.qpc_ptr);
        daqp.qpc_ptr = C_NULL
    end
end

function setup(daqp::DAQPBase.Model, qp::DAQPBase.QPj)
    daqp.qpj = qp
    daqp.qpc = DAQPBase.QPc(daqp.qpj)
    old_settings = settings(daqp); # in case setup fails
    unsafe_store!(daqp.qpc_ptr,daqp.qpc)
    setup_time = Cdouble(0);

    if(isempty(qp.H) && !isempty(qp.f))# LP
        # ensure their is no binary constraint
        @assert(!any((qp.sense.&BINARY).==BINARY),
                "DAQP requires the objective to be strictly convex to support binary variables")
        # ensure proximal-point iterations are used for LPs
        (old_settings.eps_prox == 0) && settings(daqp,Dict(:eps_prox=>1))
    end

    exitflag = ccall((:setup_daqp,DAQPBase.libdaqp),Cint,(Ptr{DAQPBase.QPc}, Ptr{DAQPBase.Workspace}, Ptr{Cdouble}), daqp.qpc_ptr, daqp.work, Ref{Cdouble}(setup_time))
    if(exitflag < 0)
        # XXX: if setup fails DAQP currently clears settings
        ccall((:allocate_daqp_settings,DAQPBase.libdaqp),Nothing,(Ptr{DAQPBase.Workspace},),daqp.work)
        settings(daqp,old_settings)
    else
        daqp.has_model = true
        daqp.x = Vector{Float64}(undef, qp.n)
        daqp.λ = Vector{Float64}(undef, qp.m)
    end

    return exitflag, setup_time
end

function setup(daqp::DAQPBase.Model, H::Matrix{Cdouble},f::Vector{Cdouble},
        A::Matrix{Cdouble},bupper::Vector{Cdouble},blower::Vector{Cdouble}=Cdouble[],
        sense::Vector{Cint}=Cint[];A_rowmaj=false,break_points = Cint[])
    return setup(daqp,QPj(H,f,A,bupper,blower,sense;A_rowmaj,break_points))
end

function solve(daqp::DAQPBase.Model)
    if(!daqp.has_model) return  zeros(0), NaN, -10, [] end
    result= Ref(DAQPResult(daqp.x,daqp.λ));

    exitflag=ccall((:daqp_solve, DAQPBase.libdaqp), Cint,
                   (Ref{DAQPBase.DAQPResult},Ref{DAQPBase.Workspace}), 
                   result,daqp.work)

    info = (x = daqp.x, λ=daqp.λ, fval=result[].fval,
            exitflag=result[].exitflag,
            status = DAQPBase.flag2status[result[].exitflag],
            solve_time = result[].solve_time,
            setup_time = result[].setup_time,
            iterations= result[].iter, nodes = result[].nodes)
    return copy(daqp.x),result[].fval,result[].exitflag,info
end


function settings(p::Ptr{DAQPBase.Workspace})
    workspace = unsafe_load(p);
    if(workspace.settings != C_NULL)
        return unsafe_load(workspace.settings)
    end
end

settings(daqp::DAQPBase.Model) = settings(daqp.work)

function settings(p::Ptr{DAQPBase.Workspace}, new_settings::DAQPBase.DAQPSettings)
    workspace = unsafe_load(p);
    if(workspace.settings != C_NULL)
        unsafe_store!(workspace.settings,new_settings)
    end
    return new_settings
end

settings(daqp::DAQPBase.Model, new_settings::DAQPBase.DAQPSettings) = settings(daqp.work,new_settings)

function settings(daqp::DAQPBase.Model,changes::Dict{Symbol,<:Any})
    return settings(daqp.work,changes)
end
function settings(p::Ptr{DAQPBase.Workspace},changes::Dict{Symbol,<:Any})
    workspace = unsafe_load(p);
    if(workspace.settings == C_NULL) return end
    settings = unsafe_load(workspace.settings)
    new = [haskey(changes,f) ? changes[f] : getfield(settings,f)
           for f in fieldnames(DAQPBase.DAQPSettings)];
    new_settings = DAQPBase.DAQPSettings(new...)
    unsafe_store!(workspace.settings,new_settings);
    return new_settings;
end

function update(daqp::DAQPBase.Model, H,f,A,bupper,blower,sense=nothing,break_points=nothing) 
    update_mask = Cint(0);
    work = unsafe_load(daqp.work);
    if(!isnothing(H) && work.n == size(H,1) && work.n == size(H,2))
        daqp.qpj.H.=H
        update_mask +=1
    end
    if(!isnothing(A) && size(A,1)==(work.m-work.ms) && size(A,2)==work.n)
        daqp.qpj.A.=A'
        update_mask+=2
    end

    if(!isnothing(f) && length(f)==work.n)
        daqp.qpj.f.=f
        update_mask+=4
    end

    if(!isnothing(bupper) && !isnothing(blower) &&
       length(bupper)==work.m && length(blower)==work.m)
        daqp.qpj.bupper.=bupper
        daqp.qpj.blower.=blower
        update_mask+=8
    end

    if(!isnothing(sense) && length(sense)== work.m)
        daqp.qpj.sense .= sense
        update_mask+=16
    end

    if(!isnothing(break_points) && length(break_points)== work.nh)
        daqp.qpj.break_points .= break_points
        update_mask+=32
    end
    daqp.qpc = QPc(daqp.qpj);

    exitflag = ccall((:update_ldp,DAQPBase.libdaqp),Cint,(Cint,Ptr{DAQPBase.Workspace},Ptr{DAQPBase.QPc}), 
                     update_mask, daqp.work,Ref(daqp.qpc));
end

function reset(p::Ptr{DAQPBase.Workspace})
    ccall((:deactivate_constraints,libdaqp),Cvoid,(Ptr{DAQPBase.Workspace},),p);
    ccall((:reset_daqp_workspace,libdaqp),Cvoid,(Ptr{DAQPBase.Workspace},),p);
end

function reset(d::DAQPBase.Model)
    reset(d.work)
end

using Downloads
function codegen(d::DAQPBase.Model; fname="daqp_workspace", dir="codegen", src=false)
    @assert(d.has_model, "setup the model before code generation")

    dir[end] != '/' && (dir*="/") ## Make sure it is correct directory path
    isdir(dir) || mkdir(dir)

    reset(d) # Make sure workspace is cleared

    exitflag = ccall((:render_daqp_workspace, libdaqp),Cvoid,
                     (Ptr{DAQPBase.Workspace},Cstring,Cstring,), d.work,fname,dir);
    if src
        cfiles = ["daqp.c","auxiliary.c","factorization.c"]
        hfiles = ["daqp.h","auxiliary.h","factorization.h","constants.h", "types.h"]
 
        # Append BnB code if there are any binary variables
        if any((d.qpj.sense.&BINARY).==BINARY)
            push!(cfiles,"bnb.c")
            push!(hfiles,"bnb.h")
        end
        if d.qpj.nh>1
            push!(cfiles,"hierarchical.c")
            push!(hfiles,"hierarchical.h")
        end

        # Copy source files from GitHub
        for f in cfiles
            Downloads.download("https://raw.githubusercontent.com/darnstrom/daqp/master/src/"*f, dir*f)
        end
        for f in hfiles
            Downloads.download("https://raw.githubusercontent.com/darnstrom/daqp/master/include/"*f, dir*f)
        end
    end
end

function setup_c_workspace(n)::Ptr{DAQPBase.Workspace}
    p = Libc.calloc(1,sizeof(DAQPBase.Workspace));
    ccall((:allocate_daqp_workspace,libdaqp), Cvoid, (Ptr{Cvoid},Cint,Cint),p, n, 0);
    ccall((:allocate_daqp_settings,libdaqp), Cvoid, (Ptr{Cvoid},),p);
    return p
end

function free_c_workspace(p::Ptr{DAQPBase.Workspace})
    ccall((:free_daqp_workspace,libdaqp), Cvoid, (Ptr{Cvoid},),p)
    Libc.free(p)
end

function init_c_workspace_ldp(p::Ptr{DAQPBase.Workspace},A::Matrix{Cdouble},bupper::Vector{Cdouble},blower::Vector{Cdouble},sense::Vector{Cint}; max_radius=nothing) 
    # Set fval_bound to maximal radius for early termination
    if(!isnothing(max_radius))
        d_work = unsafe_load(Ptr{DAQPBase.Workspace}(p));
        unsafe_store!(Ptr{Cdouble}(d_work.settings+fieldoffset(DAQPSettings,8)),max_radius);
    end

    # Pass pointers for LDP to DAQP
    unsafe_store!(Ptr{Ptr{Cdouble}}(p+fieldoffset(DAQPBase.Workspace,5)),pointer(A))
    unsafe_store!(Ptr{Ptr{Cdouble}}(p+fieldoffset(DAQPBase.Workspace,6)),pointer(bupper))
    unsafe_store!(Ptr{Ptr{Cdouble}}(p+fieldoffset(DAQPBase.Workspace,7)),pointer(blower))
    unsafe_store!(Ptr{Ptr{Cint}}(p+fieldoffset(DAQPBase.Workspace,10)),pointer(sense))
end

function isfeasible(p::Ptr{DAQPBase.Workspace}, m=nothing, ms=nothing ;validate=false)::Bool
    # Update A and bupper/blower dimensions 
    !isnothing(m)  && unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQPBase.Workspace,3)),m)
    !isnothing(ms) && unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQPBase.Workspace,4)),ms)

    exitflag =ccall((:daqp_ldp,DAQPBase.libdaqp), Int32, (Ptr{Cvoid},),p);

    if(validate && exitflag == -1)
        daqp_ws = unsafe_load(Ptr{DAQPBase.Workspace}(p))
        m,ms,n = daqp_ws.m, daqp_ws.ms, daqp_ws.n
        A = unsafe_wrap(Matrix{Cdouble}, daqp_ws.M, (n,m), own=false)
        b= unsafe_wrap(Vector{Cdouble}, daqp_ws.dupper, m+ms, own=false)
        AS = copy(unsafe_wrap(Vector{Cint}, daqp_ws.WS, daqp_ws.n_active, own=false))
        AS .+= 1 # Offset for C -> Julia
        lam_star = unsafe_wrap(Vector{Cdouble}, daqp_ws.lam_star, daqp_ws.n_active, own=false)
        err  = dot(b[AS],lam_star)+norm(A[:,AS]*lam_star)
        if(err>0)
            @warn "Couldn't validate infeas. with Frakas (err:$(err), fval=$(daqp_ws.fval))"
        end
    end
    # Make sure workspace is clean for next solve
    reset(p)
    return exitflag == 1
end

