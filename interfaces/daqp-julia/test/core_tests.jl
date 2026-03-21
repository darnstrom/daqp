using LinearAlgebra
using Random
using DAQPBase
using DAQP_jll
include(joinpath(dirname(@__FILE__), "utils.jl"))
global templib

# Use local libdaqp if available
_libdaqp = joinpath(pkgdir(DAQPBase),"libdaqp."*Libc.Libdl.dlext)
if isfile(_libdaqp)
    local_lib = true
    @info "Using local libdaqp"
    DAQP_jll.libdaqp = _libdaqp
else
    local_lib = false
end

# API Tests
nQPs,nLPs = 100,100;
n = 100; m = 500; ms = 50;
nAct = 80
kappa = 1e2
tol = 1e-4
@testset "Quadprog (C)" begin
    for nQP in 1:nQPs
        xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa);
        x,fval,exitflag,info = quadprog(H,f,A,bupper,blower,sense);
        @test norm(xref-x) < tol;
    end
    # Test quadprog interface by passing settings 
    xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa);
    s = settings(DAQPBase.Model(),Dict(:iter_limit=>1))
    x,fval,exitflag,info = quadprog(H,f,A,bupper,blower,sense;settings=s)
    @test exitflag == -4
end

@testset "Quadprog (one-sided)" begin
    for nQP in 1:10
        _,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa);
        blower = fill(-1e30,length(bupper))
        xref,fval,exitflag,info = quadprog(H,f,A,bupper,blower,sense);
        x,fval,exitflag,info = quadprog(H,f,A,bupper);
        @test norm(xref-x) < tol;
    end
end
@testset "Quadprog (C)" begin
    for nQP in 1:10
        xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa);
        x,fval,exitflag,info = DAQPBase.quadprog(DAQPBase.QPj(H,f,A,bupper,blower,sense));
        @test norm(xref-x) < tol;
    end
end

@testset "Linprog (C)" begin
    for nQP in 1:nQPs
        xref,f,A,bupper,blower,sense = generate_test_LP(n,m,ms);
        x,fval,exitflag,info = linprog(f,A,bupper,blower,sense);
        @test abs(f'*(xref-x)) < tol;
    end
end

@testset "Linprog (one-sided)" begin
    for nQP in 1:10
        _,f,A,bupper,blower,sense = generate_test_LP(n,m,ms);
        blower = fill(-1e30,length(bupper))
        xref,fval,exitflag,info = linprog(f,A,bupper,blower,sense);
        x,fval,exitflag,info = linprog(f,A,bupper);
        @test abs(f'*(xref-x)) < tol;
    end
end

@testset "Quadprog (JL)" begin
    qpj = DAQPBase.QPj()
    for selection_rule in [DAQPBase.DANTZIG, DAQPBase.BLAND]
        for nQP in 1:10
            xref,H,f,A,bupper,blower,sense = generate_test_QP(20,100,0,16,1e2);
            x,lam,AS,J,iter= DAQPBase.daqp_jl(H,f,[A;-A],[bupper;-blower],[sense;sense],Int64[];selection_rule);
            @test norm(xref-x) < tol;
        end
    end
    # Test warm start
    xref,H,f,A,bupper,blower,sense = generate_test_QP(20,100,0,16,1e2);
    x,lam,AS,J,iter= DAQPBase.daqp_jl(H,f,[A;-A],[bupper;-blower],[sense;sense],Int64[1]);
    @test norm(xref-x) < tol;
    # Test infeasible problem
    x,lam,AS,J,iter= DAQPBase.daqp_jl([1.0 0; 0 1],zeros(2),[1.0 0;-1 0],[1;-2],zeros(Cint,2),Int64[]);
    @test isinf(J)
    # Test unconstrained problem
    x,lam,AS,J,iter= DAQPBase.daqp_jl((H=[1.0 0; 0 1], f=zeros(2),
                                       A=[1.0 0;0 1], b=[1;1],
                                       senses=zeros(Cint,2)),Int64[]);
    @test isempty(AS)
end

@testset "BnB" begin
    nb = 10
    ϵb= 1e-5
    nMIQPs = 5
    for nMIQP = 1:nMIQPs
        M = randn(n,n)
        H = M'*M
        A = randn(m-ms,n);
        bupper = 20*rand(m); blower = -20*rand(m) # ensure that origin is feasible
        f = 100*randn(n); f[1:nb].=-abs.(f[1:nb]) # make it lucrative to avoid origin
        # Make first nb variables binary
        bupper[1:nb].=1
        blower[1:nb].=0
        sense = zeros(Cint,m)
        sense[1:nb].=DAQPBase.BINARY
        x,fval,exitflag,info = DAQPBase.quadprog(H,f,A,bupper,blower,sense);
        @test exitflag == 1 # was able to solve problem
        @test all((abs.(x[1:nb].-1.0).<ϵb) .| (abs.(x[1:nb]).<ϵb)) # is binary feasible
    end

    H = [1 0.5 0; 0.5 1 0.5; 0 0.5 1]
    f = [1.0;0;0]
    A = [1.0 2 3;  1 1 0]
    bu = [1.0;1;1;1e30;1e30]
    bl = [0.0;0;0;4;1]
    sense = Cint.([DAQPBase.BINARY;DAQPBase.BINARY;DAQPBase.BINARY;0;0])
    x,_,_,info = DAQPBase.quadprog(H,f,A,bu,bl,sense)
    @test norm(x-[0;1;1]) < tol 

end

@testset "Model interface" begin
    # Setup model and solve problem
    d = DAQPBase.Model()
    xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)
    setup(d,H,f,A,bupper,blower,sense)
    x,fval,exitflag,info = solve(d)
    @test norm(xref-x) < tol
    @test norm(H*x+[I(n)[1:ms,:];A]'*info.λ+f) < tol

    # Test access settings
    s = settings(d)
    @test s.primal_tol==1e-6
    settings(d,Dict(:primal_tol=>1e-5))
    s = settings(d)
    @test s.primal_tol==1e-5

    # Update existing model with new problem
    xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)
    update(d,H,f,A,bupper,blower,sense)
    x,fval,exitflag,info = solve(d)
    @test norm(xref-x) < tol
    DAQPBase.delete!(d); # Test manually running destructor
end

@testset "C LDP interface" begin
    # Setup model and solve problem
    n = 10; m = 50; ms = 5; nAct =0;
    xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)
    p=DAQPBase.setup_c_workspace(n)
    A = A'[:,:] # since row major...
    DAQPBase.init_c_workspace_ldp(p,A,bupper,blower,sense;max_radius=1e30) 
    @test isfeasible(p,m,ms)
    buold = bupper[1]
    bupper[1] = -1e30 #Make trivially infeasible
    @test !isfeasible(p,m,ms;validate=true)
    bupper[1] = buold;
    settings(p,Dict(:iter_limit => 1))
    @test !isfeasible(p,m,ms;validate=false)
    settings(p,Dict(:iter_limit => 1e4))
    @test isfeasible(p,m,ms;validate=false)

    work = unsafe_load(Ptr{DAQPBase.Workspace}(p));
    @test work.n == n
    DAQPBase.free_c_workspace(p)
end

@testset "Code generation" begin
    n = 5; m = 5; ms = 5; nAct =2;
    xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)
    sense[1] = DAQPBase.BINARY
    d = DAQPBase.Model()
    DAQPBase.setup(d,H,f,A,bupper,blower,sense)
    srcdir = tempname();
    DAQPBase.codegen(d,dir=srcdir,src=!local_lib)
    local_lib && get_local_sources(srcdir)
    src = [f for f in readdir(srcdir) if last(f,1) == "c"]
    if(!isnothing(Sys.which("gcc")))
        testlib = "daqptestlib."* Base.Libc.Libdl.dlext
        run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))
        @test isfile(joinpath(srcdir,testlib))
    end
    rm(srcdir,recursive=true)

    # Try to also get global source...
    try
        DAQPBase.codegen(d,dir=srcdir,src=true)
        @test isfile(joinpath(srcdir,"daqp.c"))
        rm(srcdir,recursive=true)
    catch
    end

    # Test special case with diagonal H 
    d = DAQPBase.Model()
    DAQPBase.setup(d,diagm(5*ones(n)),f,A,bupper,blower,sense)
    srcdir = tempname();
    DAQPBase.codegen(d,dir=srcdir,src=!local_lib)
    local_lib && get_local_sources(srcdir)
    src = [f for f in readdir(srcdir) if last(f,1) == "c"]
    if(!isnothing(Sys.which("gcc")))
        testlib = "daqptestlib."* Base.Libc.Libdl.dlext
        run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))
        @test isfile(joinpath(srcdir,testlib))
    end
    rm(srcdir,recursive=true)
end

@testset "Hierarchical QP" begin
    A = [1.0 1 1; 1 -1 0; 3 1 -1]
    bu = [ones(3);1;0.5;20]
    bl = [-ones(3);-1e30;0.5;10]
    sense = zeros(Cint,6)
    xref = [1; 0.5; -1]
    d = DAQPBase.Model()
    DAQPBase.setup(d,zeros(0,0),zeros(0),A,bu,bl,sense;break_points = [3;4;5;6])
    x,fval,exitflag,info = solve(d)
    @test norm(xref-x) < tol

    # Test codegen
    srcdir = tempname();
    DAQPBase.codegen(d,dir=srcdir,src=!local_lib)
    local_lib && get_local_sources(srcdir)
    src = [f for f in readdir(srcdir) if last(f,1) == "c"]

    if(!isnothing(Sys.which("gcc")))
        testlib = "daqptestlib."* Base.Libc.Libdl.dlext
        run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))
        @test isfile(joinpath(srcdir,testlib))
        global templib = joinpath(srcdir,testlib)
        daqp_work_ptr = cglobal((:daqp_work,templib),DAQPBase.Workspace)
        ws = unsafe_load(daqp_work_ptr)

        # Setup correct LDP 
        dupper = unsafe_wrap(Vector{Cdouble}, ws.dupper, ws.m, own=false)
        dlower = unsafe_wrap(Vector{Cdouble}, ws.dlower, ws.m, own=false)
        sense = unsafe_wrap(Vector{Cint}, ws.sense, ws.m, own=false)
        Anorm = [ones(3);[norm(A[i,:],2) for i in 1:size(A,2)]]
        dupper .= bu ./ Anorm
        dlower .= bl ./ Anorm

        # Solve
        lambda = zeros(ws.m) 
        exitflag = ccall((:daqp_hiqp,templib),Cint,(Ptr{DAQPBase.Workspace},Ptr{Cdouble}),daqp_work_ptr,lambda);
        ccall((:ldp2qp_solution,templib),Cvoid,(Ptr{DAQPBase.Workspace},),daqp_work_ptr);
        xgen = copy(unsafe_wrap(Vector{Cdouble}, ws.x, ws.n, own=false))
        @test norm(xref-xgen) < tol
    end

    # Degenerate
    H = [10.5 4.0 2.0; 4.0 5.5 0.5; 2.0 0.5 2.0]
    f = [-53.0; -30; -11.5]
    A = [1.0 0 0; 1 1 0; 0 0 0; 1 0 0];
    bu = [3*ones(3);7.5;7.5;5.0;10.0]
    bl = [-3*ones(3);4.5;4.5;2;7]
    sense = zeros(Cint,7)
    sense[6] = 4
    d = DAQPBase.Model()
    DAQPBase.setup(d,zeros(0,0),zeros(0),A,bu,bl,sense;break_points = Cint.([3;5;7]))
    x,fval,exitflag,info = solve(d)
    @test exitflag > 0

    # Degenerate #2
    A = [1.0 0; 1 0; 0 1]
    bu = [4.0;8;1]
    bl = [4.0;8;1]
    sense = zeros(Cint,3)
    d = DAQPBase.Model()
    DAQPBase.setup(d,zeros(0,0),zeros(0),A,bu,bl,sense;break_points = [0,2,3])
    x,fval,exitflag,info = solve(d)
    @test norm([6.0;1]-x) < tol

end

@testset "Trivial infeasible" begin
    H = [6.837677669279314 1.3993262799977795 1.9781574256330445 0.7988389688453156;
         1.3993262799977795 4.91607513347457 0.8347008717503388 0.964319980996552;
         1.9781574256330445 0.8347008717503388 5.6186371819867755 0.23421356059787485;
         0.7988389688453156 0.964319980996552 0.23421356059787485 5.512564828518534]
    f = [-0.9018168388545096, 1.3888380439021342, -3.2050167583822065, 6.2604158413126205]

    A = [0.0 0.0 0.0 0.0;
         1.0 0.0 1.0 0.0;
         0.0 0.0 0.0 0.0;
         -1.0 0.0 -1.0 0.0;
         1.0 0.0 1.0 0.0;
         0.5 1.0 0.5 1.0;
         -1.0 0.0 -1.0 0.0;
         -0.5 -1.0 -0.5 -1.0;
         1.0 0.0 0.0 0.0;
         -1.0 0.0 0.0 0.0;
         0.0 1.0 0.0 0.0;
         0.0 -1.0 0.0 0.0;
         0.0 0.0 1.0 0.0;
         0.0 0.0 -1.0 0.0;
         0.0 0.0 0.0 1.0;
         0.0 0.0 0.0 -1.0]
    b = [2.2693082025353517, 1.3445735938597536, -0.2693082025353519, 0.6554264061402464, 2.2330893356345, 1.172286796929877, -0.23308933563449985, 0.8277132030701232, 1.0, 0.5, 1.0, 0.5, 0.5, 2.0, 0.5, 2.0]

    d = DAQPBase.Model()

    exitflag,_ = DAQPBase.setup(d,H,f,A,b,Float64[],zeros(Cint, length(b)))
    @test exitflag == -1
    _, _, exitflag, _ = DAQPBase.solve(d)
    @test exitflag < 0
end

@testset "Affine variational inequality" begin
    for _ in 1:10
        n = 100; m = 500
        xref,H,f,A,b = generate_test_avi(n,m);
        # C api
        x,λ,info = DAQPBase.avi(H,f,A,b);
        @test norm(xref-x) < tol;

        # Pure Julia implemenetation 
        x,λ,info = DAQPBase.solve_avi_jl(H,f,A,b);
        @test norm(xref-x) < tol;

        # With setup
        d = DAQPBase.Model()
        DAQPBase.setup(d,H,f,A,b;is_avi=true)
        x,fval,exitflag,info = solve(d)
        @test norm(xref-x) < tol;
        x,fval,exitflag,info = solve(d)
        @test info.iterations == 5 
        @test norm(xref-x) < tol;
    end

    # Test that update does not cause a segfault for AVIs
    n = 10; m = 50
    xref,H,f,A,b = generate_test_avi(n,m);
    d = DAQPBase.Model()
    DAQPBase.setup(d,H,f,A,b;is_avi=true)
    x,fval,exitflag,info = solve(d)
    @test exitflag > 0
    update(d,nothing,-f,nothing,nothing,nothing)
    GC.gc(); GC.gc(); GC.gc() # Force GC to expose stale pointer bugs
    x2,fval2,exitflag2,info2 = solve(d)
    @test exitflag2 > 0
end

@testset "Prefactorized Hessian" begin
    xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa);
    C = cholesky(H)
    # Quadprog API
    x,fval,exitflag,info = quadprog(C,f,A,bupper,blower,sense)
    @test norm(xref-x) < tol;

    # Test Model API 
    d = DAQPBase.Model()
    setup(d,C,f,A,bupper,blower,sense)
    x,fval,exitflag,info = solve(d)
    @test norm(xref-x) < tol;

    # Test diagonal
    H = diagm(rand(n))
    C = cholesky(H)
    xref,fval,exitflag,info = quadprog(H,f,A,bupper,blower,sense)
    x,fval,exitflag,info = quadprog(C,f,A,bupper,blower,sense)
    @test norm(xref-x) < tol;

end

@testset "Setting Warm Start" begin
    # Test primal start 
    d = DAQPBase.Model()
    xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)
    setup(d,H,f,A,bupper,blower,sense;primal_start=xref)
    x,fval,exitflag,info = solve(d)
    @test norm(xref-x) < tol
    @test info.iterations==1

    x,fval,exitflag,info = quadprog(H,f,A,bupper,blower,sense;primal_start=xref)
    @test norm(xref-x) < tol
    @test info.iterations==1

    # Test dual start
    λstar = info.λ
    d = DAQPBase.Model()
    setup(d,H,f,A,bupper,blower,sense;dual_start=λstar)
    x,fval,exitflag,info = solve(d)
    @test norm(xref-x) < tol
    @test info.iterations==1

    x,fval,exitflag,info = quadprog(H,f,A,bupper,blower,sense;dual_start=λstar)
    @test norm(xref-x) < tol
    @test info.iterations==1

    # Ensure solver recovers from degenerate starting point
    d = DAQPBase.Model()
    H,f = [1.0 0; 0 1.0], zeros(2);
    b = [1.0;1.0;2]
    A = ones(1,2);
    setup(d,H,f,A,b;primal_start=[1.0;1.0])
    x,fval,exitflag,info = solve(d)
    @test norm(x-zeros(2)) < tol;
    @test info.iterations > 1;

    # Test warm start for LPs (prox iters...)
    xref,f,A,bupper,blower,sense = generate_test_LP(n,m,ms);
    xcold,fval,exitflag,info_cold = linprog(f,A,bupper,blower,sense)
    xwarm,fval,exitflag,info_warm = linprog(f,A,bupper,blower,sense;primal_start = 0.95*xref)
    @test info_cold.iterations > info_warm.iterations 
    @test norm(xwarm-xref) < tol

    # Test warm start for AVIs
    xref,H,f,A,b = generate_test_avi(n,m);
    xcold,fval,exitflag,info_cold = DAQPBase.avi(H,f,A,b)
    xwarm,fval,exitflag,info_warm = DAQPBase.avi(H,f,A,b;primal_start = 0.95*xref)
    @test norm(xwarm-xref) < tol

end

@testset "Time limit" begin
    # Use a large problem so the solver takes multiple iterations
    xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)

    # Setting a tiny time limit (1 nanosecond) should trigger TIMELIMIT exit
    s = settings(DAQPBase.Model(), Dict(:time_limit => 1e-9))
    x,fval,exitflag,info = quadprog(H,f,A,bupper,blower,sense; settings=s)
    @test exitflag == DAQPBase.TIMELIMIT

    # Setting a generous time limit should allow the solver to find the optimum
    s = settings(DAQPBase.Model(), Dict(:time_limit => 100.0))
    x,fval,exitflag,info = quadprog(H,f,A,bupper,blower,sense; settings=s)
    @test exitflag == DAQPBase.OPTIMAL
    @test norm(xref-x) < tol
end

@testset "Semi-proximal method" begin
    # --- PD Hessian: n_prox must be 0, result matches eps_prox=0 solve ---
    n2 = 5; m2 = 20; ms2 = 5; nAct2 = 4
    xref,H,f,A,bupper,blower,sense = generate_test_QP(n2,m2,ms2,nAct2,1e2)

    # Reference (no proximal)
    xref2,fval_ref,ef_ref,_ = quadprog(H,f,A,bupper,blower,sense)
    @test ef_ref == DAQPBase.OPTIMAL

    # With eps_prox > 0 on PD H: solver should still reach the optimum.
    s_prox = settings(DAQPBase.Model(), Dict(:eps_prox => 1e-4))
    x_prox,_,ef_prox,_ = quadprog(H,f,A,bupper,blower,sense; settings=s_prox)
    @test ef_prox == DAQPBase.OPTIMAL
    @test norm(xref2 - x_prox) < tol

    # Check that the Workspace struct layout is correct by reading n_prox via
    # unsafe_load -- it must be 0 for a PD Hessian (no regularisation needed).
    d = DAQPBase.Model()
    DAQPBase.settings(d, Dict(:eps_prox => 1e-4))
    DAQPBase.setup(d, H, f, A, bupper, blower, sense)
    p = d.work  # raw workspace pointer (Ptr{Cvoid})
    ws = unsafe_load(Ptr{DAQPBase.Workspace}(p))
    @test ws.n_prox == 0

    # --- Rank-1 Hessian: x2 direction is singular -> n_prox == 1 ---
    H_sing = [1.0 0.0; 0.0 0.0]
    f_sing = [1.0; 1.0]
    A_sing = zeros(0, 2)
    bu_sing = [2.0; 2.0]
    bl_sing = [-2.0; -2.0]
    sense_sing = zeros(Cint, 2)

    d2 = DAQPBase.Model()
    DAQPBase.settings(d2, Dict(:eps_prox => 1e-3))
    DAQPBase.setup(d2, H_sing, f_sing, A_sing, bu_sing, bl_sing, sense_sing)
    ws2 = unsafe_load(Ptr{DAQPBase.Workspace}(d2.work))
    @test ws2.n_prox == 1   # only x2 direction needed regularisation

    x2,_,ef2,_ = DAQPBase.solve(d2)
    @test ef2 == DAQPBase.OPTIMAL
    @test abs(x2[1] - (-1.0)) < tol   # x1* = -1
    @test abs(x2[2] - (-2.0)) < tol   # x2* = -2 (lower bound)

    # --- Zero Hessian: all directions singular -> n_prox == n ---
    n3 = 3
    H_zero = zeros(n3, n3)
    f_zero = ones(n3)
    A_zero = zeros(0, n3)
    bu_zero =  5.0*ones(n3)
    bl_zero = -5.0*ones(n3)
    sense_zero = zeros(Cint, n3)

    d3 = DAQPBase.Model()
    DAQPBase.settings(d3, Dict(:eps_prox => 1e-2))
    DAQPBase.setup(d3, H_zero, f_zero, A_zero, bu_zero, bl_zero, sense_zero)
    ws3 = unsafe_load(Ptr{DAQPBase.Workspace}(d3.work))
    @test ws3.n_prox == n3

    x3,_,ef3,_ = DAQPBase.solve(d3)
    @test ef3 > 0
    @test norm(x3 .- (-5.0)) < tol  # all at lower bound
end

