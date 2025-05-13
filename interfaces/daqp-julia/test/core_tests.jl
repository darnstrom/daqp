using LinearAlgebra
using Random
using DAQPBase
using DAQP_jll
include(joinpath(dirname(@__FILE__), "utils.jl"))

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
    if(local_lib)
        # Get local source
        daqp_dir = joinpath(dirname(@__FILE__), "..","..","..","..")
        cfiles = ["daqp.c","auxiliary.c","factorization.c", "bnb.c", "hierarchical.c"]
        hfiles = ["daqp.h","auxiliary.h","factorization.h", "bnb.h", "hierarchical.h","constants.h", "types.h"]
        for cf in cfiles
            cp(joinpath(daqp_dir,"src",cf), joinpath(srcdir,cf))
        end
        for hf in hfiles
            cp(joinpath(daqp_dir,"include",hf), joinpath(srcdir,hf))
        end
    end
    src = [f for f in readdir(srcdir) if last(f,1) == "c"]
    if(!isnothing(Sys.which("gcc")))
        testlib = "daqptestlib."* Base.Libc.Libdl.dlext
        run(Cmd(`gcc -lm -fPIC -O3 -msse3 -xc -shared -o $testlib $src`; dir=srcdir))
        @test isfile(joinpath(srcdir,testlib))
    end
    rm(srcdir,recursive=true)

    # Try to also get global source...
    DAQPBase.codegen(d,dir=srcdir,src=true)
    @test isfile(joinpath(srcdir,"daqp.c"))
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
end
