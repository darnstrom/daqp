"""
    benchmark.jl

Performance benchmark script for DAQP solver.
Measures solve time for representative QP and LP problems at different scales.
Results are saved to CSV for easy comparison and regression detection.

Usage:
    julia benchmark.jl                          # Run benchmarks with default settings
    julia benchmark.jl --output results.csv    # Specify output file
    julia benchmark.jl --suite small            # Run only small problems
"""

using LinearAlgebra
using Random
using Statistics
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

function benchmark_qp(n, m, ms, nAct, kappa; num_runs=11)
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
    
    total_times = setup_times .+ solve_times
    return (
        setup_median=median(setup_times),
        setup_min=minimum(setup_times),
        setup_max=maximum(setup_times),
        solve_median=median(solve_times),
        solve_min=minimum(solve_times),
        solve_max=maximum(solve_times),
        total_median=median(total_times),
        iter_median=median(iterations),
        iter_min=minimum(iterations),
        iter_max=maximum(iterations),
        iterations=iterations,
        num_runs=num_runs
    )
end

function benchmark_lp(n, m, ms; num_runs=11)
    """
    Benchmark linear programming with given problem size.
    Tracks setup_time and solve_time separately.
    Returns statistics for both times and iterations.
    """
    setup_times = Float64[]
    solve_times = Float64[]
    iterations = Int[]
    
    for run in 1:num_runs
        # Generate test problem with known solution
        xref, f, A, bupper, blower, sense = generate_test_LP(n, m, ms)
        
        # Solve and extract timing from info struct
        x, fval, exitflag, info = linprog(f, A, bupper, blower, sense)
        
        # Verify correctness
        if abs(f' * (xref - x)) > CORRECTNESS_TOL
            @warn "Solution quality issue: |f'*(xref - x)| = $(abs(f' * (xref - x)))"
        end
        
        push!(setup_times, info.setup_time)
        push!(solve_times, info.solve_time)
        push!(iterations, info.iterations)
    end
    
    total_times = setup_times .+ solve_times
    return (
        setup_median=median(setup_times),
        setup_min=minimum(setup_times),
        setup_max=maximum(setup_times),
        solve_median=median(solve_times),
        solve_min=minimum(solve_times),
        solve_max=maximum(solve_times),
        total_median=median(total_times),
        iter_median=median(iterations),
        iter_min=minimum(iterations),
        iter_max=maximum(iterations),
        iterations=iterations,
        num_runs=num_runs
    )
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
        "timestamp,daqp_version,problem_type,problem_id,n_variables,n_constraints,n_simple_bounds,condition_number,setup_time_median_s,solve_time_median_s,total_time_median_s,iter_median,iter_min,iter_max,num_runs"
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
            stats.setup_median,
            stats.solve_median,
            stats.total_median,
            stats.iter_median,
            stats.iter_min,
            stats.iter_max,
            stats.num_runs
        ], ",")
        push!(csv_lines, csv_row)
        
        println("    Setup: $(round(stats.setup_median*1e6; digits=2))µs (median) | Solve: $(round(stats.solve_median*1e6; digits=2))µs (median) | Iters: $(round(stats.iter_median; digits=1)) (median)")
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
            stats.setup_median,
            stats.solve_median,
            stats.total_median,
            stats.iter_median,
            stats.iter_min,
            stats.iter_max,
            stats.num_runs
        ], ",")
        push!(csv_lines, csv_row)
        
        println("    Setup: $(round(stats.setup_median*1e6; digits=2))µs (median) | Solve: $(round(stats.solve_median*1e6; digits=2))µs (median) | Iters: $(round(stats.iter_median; digits=1)) (median)")
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
        else
            i += 1
        end
    end
    
    run_benchmarks(; suite=suite, output_file=output_file, use_local=use_local)
end
