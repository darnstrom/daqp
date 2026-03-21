"""
    benchmark.jl

Performance benchmark script for DAQP solver.
Measures solve time for representative QP and LP problems at different scales.
Results are saved to CSV for easy comparison and regression detection.

Also benchmarks the semi-proximal approach vs the old full-proximal (eps*I)
approach for QPs with rank-deficient Hessians.

Usage:
    julia benchmark.jl                          # Run benchmarks with default settings
    julia benchmark.jl --output results.csv    # Specify output file
    julia benchmark.jl --suite small            # Run only small problems
    julia benchmark.jl --prox                   # Run semi-proximal vs full-proximal comparison
"""

using LinearAlgebra
using Random
using Statistics
using Printf
using DAQPBase
using DAQP_jll
using Dates

include("utils.jl")

# Set random seed for reproducibility
Random.seed!(42)

# Problem sizes: (n_variables, n_constraints, n_simple_bounds, n_active, condition_number)
PROBLEM_SIZES = Dict(
    "small" => [(10, 50, 5, 8, 1e2)],
    "medium" => [(50, 250, 25, 40, 1e2), (100, 500, 50, 80, 1e2)],
    "large" => [(200, 1000, 100, 160, 1e2), (500, 2500, 250, 400, 1e2)],
)

# Tolerance for solution correctness
CORRECTNESS_TOL = 1e-4

function benchmark_qp(n, m, ms, nAct, kappa; num_runs=5)
    """
    Benchmark quadratic programming with given problem size.
    Tracks setup_time and solve_time separately.
    Returns statistics for both times and iterations.
    """
    setup_times = Float64[]
    solve_times = Float64[]
    iterations = Int[]
    
    for run in 1:num_runs
        # Generate test problem with known solution
        xref, H, f, A, bupper, blower, sense = generate_test_QP(n, m, ms, nAct, kappa)
        
        # Solve and extract timing from info struct
        x, fval, exitflag, info = quadprog(H, f, A, bupper, blower, sense)
        
        # Verify correctness
        if norm(xref - x) > CORRECTNESS_TOL
            @warn "Solution quality issue: norm(xref - x) = $(norm(xref - x))"
        end
        
        push!(setup_times, info.setup_time)
        push!(solve_times, info.solve_time)
        push!(iterations, info.iterations)
    end
    
    return (
        setup_mean=mean(setup_times),
        setup_std=std(setup_times),
        setup_min=minimum(setup_times),
        setup_max=maximum(setup_times),
        solve_mean=mean(solve_times),
        solve_std=std(solve_times),
        solve_min=minimum(solve_times),
        solve_max=maximum(solve_times),
        total_mean=mean(setup_times .+ solve_times),
        total_std=std(setup_times .+ solve_times),
        iter_mean=mean(iterations),
        iter_std=std(iterations),
        iter_min=minimum(iterations),
        iter_max=maximum(iterations),
        iterations=iterations
    )
end

function benchmark_lp(n, m, ms; num_runs=5)
    """
    Benchmark linear programming with given problem size.
    Tracks setup_time and solve_time separately.
    Returns statistics for both times and iterations.
    """
    setup_times = Float64[]
    solve_times = Float64[]
    iterations = Int[]

    # Build an LP settings struct with eps_prox=1.  We use the old-style
    # quadprog(QPj; settings=...) API so that eps_prox is passed directly to
    # the C function (daqp_quadprog) rather than through the workspace struct.
    # This makes the benchmark work correctly with both the current C library
    # (which dispatches LPs via n_prox) and older libraries (which dispatch via
    # eps_prox != 0), avoiding a struct-layout mismatch when swapping .so files.
    _default_s = DAQPBase.DAQPSettings()
    _lp_settings = DAQPBase.DAQPSettings(
        [f == :eps_prox ? 1.0 : getfield(_default_s, f)
         for f in fieldnames(DAQPBase.DAQPSettings)]...)

    for run in 1:num_runs
        # Generate test problem with known solution
        xref, f, A, bupper, blower, sense = generate_test_LP(n, m, ms)

        # Solve and extract timing from info struct
        qpj = DAQPBase.QPj(zeros(0,0), f, A, bupper, blower, sense)
        x, fval, exitflag, info = DAQPBase.quadprog(qpj; settings=_lp_settings)
        
        # Verify correctness
        if abs(f' * (xref - x)) > CORRECTNESS_TOL
            @warn "Solution quality issue: |f'*(xref - x)| = $(abs(f' * (xref - x)))"
        end
        
        push!(setup_times, info.setup_time)
        push!(solve_times, info.solve_time)
        push!(iterations, info.iterations)
    end
    
    return (
        setup_mean=mean(setup_times),
        setup_std=std(setup_times),
        setup_min=minimum(setup_times),
        setup_max=maximum(setup_times),
        solve_mean=mean(solve_times),
        solve_std=std(solve_times),
        solve_min=minimum(solve_times),
        solve_max=maximum(solve_times),
        total_mean=mean(setup_times .+ solve_times),
        total_std=std(setup_times .+ solve_times),
        iter_mean=mean(iterations),
        iter_std=std(iterations),
        iter_min=minimum(iterations),
        iter_max=maximum(iterations),
        iterations=iterations
    )
end

"""
    generate_rank_deficient_QP(n, m, ms, rank, nAct, kappa)

Generate a QP whose Hessian has exactly `rank` positive eigenvalues; the
remaining (n - rank) eigenvalues are zero.  The linear term f and constraints
are retained from the PD reference problem so that the constrained minimum is
well-defined (the constraints clip the null-space directions).

Returns (H, f, A, bupper, blower, sense) — no `xref` since the true solution
of the rank-deficient problem differs from the PD one.
"""
function generate_rank_deficient_QP(n, m, ms, rank, nAct, kappa)
    _, H_pd, f, A, bupper, blower, sense = generate_test_QP(n, m, ms, nAct, kappa)
    F = eigen(Symmetric(H_pd))
    vals = max.(F.values, 0.0)      # remove floating-point noise
    vals[rank+1:end] .= 0.0         # zero out the (n - rank) smallest eigenvalues
    H = Symmetric(F.vectors * Diagonal(vals) * F.vectors')
    return Matrix(H), f, A, bupper, blower, sense
end

"""
    benchmark_prox_comparison(n, m, ms, rank, nAct, kappa; num_runs=10)

Compare the semi-proximal method (new: eps applied only to singular directions,
selected by `prox_mask`) against the old full-proximal method (eps·I applied to
ALL variables regardless of whether they are singular).

The old behaviour is simulated by adding a small `η·I` to H before factorisation
so that all Cholesky diagonals are strictly positive and `prox_mask` is all-zero
— meaning n_prox = 0 and daqp_ldp would run once without any outer proximal
loop.  That single-shot solve corresponds to what the old code produced when the
user had to supply a pre-regularised H.

The difference between the methods therefore becomes visible through the *number
of outer proximal iterations* the solver performs: the semi-proximal method
perturbs only the truly singular directions, keeping the inner QP closer to the
original problem, which typically requires fewer outer iterations to converge.
"""
function benchmark_prox_comparison(n, m, ms, rank, nAct, kappa; num_runs=10)
    semi_times  = Float64[]
    old_times   = Float64[]
    semi_iters  = Int[]
    old_iters   = Int[]
    n_prox_vals = Int[]

    for _ in 1:num_runs
        H, f, A, bupper, blower, sense = generate_rank_deficient_QP(n, m, ms, rank, nAct, kappa)

        # --- Semi-proximal (new): n_prox identifies and perturbs only singular dirs ---
        x_semi, _, ef_semi, info_semi = quadprog(H, f, A, bupper, blower, sense)
        push!(semi_times, info_semi.solve_time)
        push!(semi_iters, info_semi.iterations)

        # Check n_prox for this problem (all runs should be identical rank)
        if isempty(n_prox_vals)
            d_tmp = DAQPBase.Model()
            setup(d_tmp, H, f, A, bupper, blower, sense)
            ws_tmp = unsafe_load(Ptr{DAQPBase.Workspace}(d_tmp.work))
            push!(n_prox_vals, ws_tmp.n_prox)
        end

        # --- Old full-proximal (simulated): pre-regularise H with a small η·I
        # so all directions look PD to the Cholesky.  This mimics the old code
        # which applied eps to every direction unconditionally.  The resulting
        # problem is then solved with eps_prox = 0 (direct daqp_ldp), giving the
        # "one-shot" result the old code would have produced. ---
        η = 1.0  # regularisation amount added to H (matches the eps used in the old full-proximal approach)
        H_reg = H + η * I
        s_old = settings(DAQPBase.Model(), Dict(:eps_prox => 0.0))
        x_old, _, ef_old, info_old = quadprog(H_reg, f, A, bupper, blower, sense; settings=s_old)
        push!(old_times, info_old.solve_time)
        push!(old_iters, info_old.iterations)
    end

    return (
        semi_time_mean = mean(semi_times),
        semi_time_std  = std(semi_times),
        semi_iter_mean = mean(semi_iters),
        semi_iter_std  = std(semi_iters),
        old_time_mean  = mean(old_times),
        old_time_std   = std(old_times),
        old_iter_mean  = mean(old_iters),
        old_iter_std   = std(old_iters),
        n_prox         = isempty(n_prox_vals) ? -1 : n_prox_vals[1],
        num_runs       = num_runs,
    )
end

"""
    run_prox_benchmark(; output_file="prox_comparison.csv", use_local=false)

Run semi-proximal vs old full-proximal comparison and print a summary table.

Columns:
  n          — problem dimension
  rank       — number of positive Hessian eigenvalues
  n_sing     — number of singular directions (n - rank)
  n_prox     — directions actually regularised by semi-proximal
  semi_µs    — semi-proximal solve time (µs)
  old_µs     — old full-proximal solve time (µs)
  semi_iters — outer iterations for semi-proximal
  old_iters  — outer iterations for old approach (always 1 since it's a
               direct single solve of the pre-regularised problem)
"""
function run_prox_benchmark(; output_file="prox_comparison.csv", use_local=false)
    _libdaqp = joinpath(pkgdir(DAQPBase),"libdaqp."*Libc.Libdl.dlext)
    if isfile(_libdaqp) && use_local
        @info "Using local libdaqp"
        DAQP_jll.libdaqp = _libdaqp
    end

    # (n, m, ms, rank, nAct, kappa)
    cases = [
        (20,  100,  10,  15,  16,  1e2),   # small,  1 singular direction
        (20,  100,  10,  10,  16,  1e2),   # small,  half-rank
        (20,  100,  10,   5,  16,  1e2),   # small,  quarter-rank
        (50,  250,  25,  40,  40,  1e2),   # medium, 10 singular directions
        (50,  250,  25,  25,  40,  1e2),   # medium, half-rank
        (100, 500,  50,  80,  80,  1e2),   # large,  20 singular directions
    ]

    println("\n=== Semi-proximal vs Old Full-proximal Benchmark ===")
    println("  semi: new code — eps only on singular directions (n_prox ≤ n)")
    println("  old:  old code — pre-regularised H with eps·I on all directions")
    println()
    header = @sprintf("%-6s %-5s %-6s %-7s  %10s %10s  %11s %11s",
                      "n", "rank", "n_sing", "n_prox",
                      "semi_µs", "old_µs",
                      "semi_iters", "old_iters")
    println(header)
    println("-"^length(header))

    csv_rows = String["n,rank,n_sing,n_prox,semi_time_mean_us,semi_time_std_us,old_time_mean_us,old_time_std_us,semi_iter_mean,old_iter_mean,num_runs"]

    for (n, m, ms, rank, nAct, kappa) in cases
        r = benchmark_prox_comparison(n, m, ms, rank, nAct, kappa; num_runs=10)
        n_sing = n - rank
        println(@sprintf("%-6d %-5d %-6d %-7d  %10.1f %10.1f  %11.1f %11.1f",
                         n, rank, n_sing, r.n_prox,
                         r.semi_time_mean*1e6, r.old_time_mean*1e6,
                         r.semi_iter_mean, r.old_iter_mean))
        push!(csv_rows, join([n, rank, n_sing, r.n_prox,
                              r.semi_time_mean*1e6, r.semi_time_std*1e6,
                              r.old_time_mean*1e6,  r.old_time_std*1e6,
                              r.semi_iter_mean, r.old_iter_mean, r.num_runs], ","))
    end

    open(output_file, "w") do io
        foreach(l -> println(io, l), csv_rows)
    end
    println("\n✓ Results saved to: $output_file")
end

function run_benchmarks(; suite="all", output_file="daqp_benchmark_results.csv", use_local=false)
    """
    Run all benchmarks and save results to CSV.
    suite: "small", "medium", "large", or "all"
    """

    # Use local libdaqp if available
    _libdaqp = joinpath(pkgdir(DAQPBase),"libdaqp."*Libc.Libdl.dlext)
    if isfile(_libdaqp) && use_local
        local_lib = true
        @info "Using local libdaqp"
        println("local")
        DAQP_jll.libdaqp = _libdaqp
    else
        @info "Using DAQP_jll libdaqp"
        println("not local")
        local_lib = false
    end
    
    # Determine which sizes to run
    sizes_to_run = if suite == "all"
        vcat(values(PROBLEM_SIZES)...)
    elseif haskey(PROBLEM_SIZES, suite)
        PROBLEM_SIZES[suite]
    else
        @error "Unknown suite: $suite. Must be one of: all, small, medium, large"
        return
    end
    
    # Prepare CSV header and data
    csv_lines = [
        "timestamp,daqp_version,problem_type,problem_id,n_variables,n_constraints,n_simple_bounds,condition_number,setup_time_mean_s,setup_time_std_s,solve_time_mean_s,solve_time_std_s,total_time_mean_s,total_time_std_s,iter_mean,iter_std,iter_min,iter_max,num_runs"
    ]
    
    timestamp = string(now())
    version = string(pkgversion(DAQPBase))
    
    @info "Starting DAQP performance benchmarks..."
    @info "Running suite: $suite"
    
    # Run QP benchmarks
    println("\n=== Quadratic Programming Benchmarks ===")
    for (i, (n, m, ms, nAct, kappa)) in enumerate(sizes_to_run)
        problem_id = "qp_$(n)_$(m)_$(ms)_$(nAct)_$(Int(log10(kappa)))"
        println("  [$i/$(length(sizes_to_run))] QP: n=$n, m=$m, ms=$ms, nAct=$nAct, κ=$kappa")
        
        stats = benchmark_qp(n, m, ms, nAct, kappa)
        
        csv_row = join([
            timestamp,
            version,
            "QP",
            problem_id,
            n,
            m,
            ms,
            kappa,
            stats.setup_mean,
            stats.setup_std,
            stats.solve_mean,
            stats.solve_std,
            stats.total_mean,
            stats.total_std,
            stats.iter_mean,
            stats.iter_std,
            stats.iter_min,
            stats.iter_max,
            5  # num_runs
        ], ",")
        push!(csv_lines, csv_row)
        
        println("    Setup: $(round(stats.setup_mean*1e6; digits=2))±$(round(stats.setup_std*1e6; digits=2))µs | Solve: $(round(stats.solve_mean*1e6; digits=2))±$(round(stats.solve_std*1e6; digits=2))µs | Iters: $(round(stats.iter_mean; digits=1))±$(round(stats.iter_std; digits=1))")
    end
    
    # Run LP benchmarks
    println("\n=== Linear Programming Benchmarks ===")
    for (i, (n, m, ms, _, _)) in enumerate(sizes_to_run)
        problem_id = "lp_$(n)_$(m)_$(ms)"
        println("  [$i/$(length(sizes_to_run))] LP: n=$n, m=$m, ms=$ms")
        
        stats = benchmark_lp(n, m, ms)
        
        csv_row = join([
            timestamp,
            version,
            "LP",
            problem_id,
            n,
            m,
            ms,
            "",  # condition_number (N/A for LP)
            stats.setup_mean,
            stats.setup_std,
            stats.solve_mean,
            stats.solve_std,
            stats.total_mean,
            stats.total_std,
            stats.iter_mean,
            stats.iter_std,
            stats.iter_min,
            stats.iter_max,
            5  # num_runs
        ], ",")
        push!(csv_lines, csv_row)
        
        println("    Setup: $(round(stats.setup_mean*1e6; digits=2))±$(round(stats.setup_std*1e6; digits=2))µs | Solve: $(round(stats.solve_mean*1e6; digits=2))±$(round(stats.solve_std*1e6; digits=2))µs | Iters: $(round(stats.iter_mean; digits=1))±$(round(stats.iter_std; digits=1))")
    end
    
    # Save to CSV
    open(output_file, "w") do f
        for line in csv_lines
            println(f, line)
        end
    end
    println("\n✓ Results saved to: $output_file")
    
    return csv_lines
end
if abspath(PROGRAM_FILE) == @__FILE__
    # Parse command line arguments
    local suite = "all"
    local output_file = "daqp_benchmark_results.csv"
    local use_local = false
    local run_prox = false
    local i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--output" && i < length(ARGS)
            output_file = ARGS[i+1]
            i += 2
        elseif ARGS[i] == "--suite" && i < length(ARGS)
            suite = ARGS[i+1]
            i += 2
        elseif ARGS[i] == "--local"
            use_local = true 
            i += 1 
        elseif ARGS[i] == "--prox"
            run_prox = true
            i += 1
        else
            i += 1
        end
    end
    
    if run_prox
        run_prox_benchmark(; output_file="prox_comparison.csv", use_local=use_local)
    else
        run_benchmarks(; suite=suite, output_file=output_file, use_local=use_local)
    end
end
