"""
    compare_benchmarks.jl

Compare benchmark results to detect performance regressions.
Loads two CSV benchmark files and reports differences.

Usage:
    julia compare_benchmarks.jl baseline.csv current.csv          # Compare two files
    julia compare_benchmarks.jl --baseline baseline.csv current.csv --threshold 5
"""

function load_benchmarks(filename::String)
    """Load benchmark results from CSV file."""
    if !isfile(filename)
        error("File not found: $filename")
    end

    lines = readlines(filename)
    if isempty(lines)
        error("Empty benchmark file: $filename")
    end

    # Parse header
    header = split(lines[1], ",")

    # Parse data rows
    benchmarks = []
    for line in lines[2:end]
        if isempty(strip(line))
            continue
        end
        values = split(line, ",")
        row = Dict(header[i] => values[i] for i in 1:length(header))
        push!(benchmarks, row)
    end

    return benchmarks
end

function parse_float(s::Union{String, SubString})
    """Safely parse a string to float, returning nothing if empty or invalid."""
    s_str = string(s)  # Convert SubString to String if needed
    if isempty(strip(s_str))
        return nothing
    end
    try
        return parse(Float64, s_str)
    catch
        return nothing
    end
end

function compare_benchmarks(baseline_file::String, current_file::String;
                           regression_threshold::Float64=5.0)
    """
    Compare two benchmark result files and report regressions.
    regression_threshold: percentage slowdown to flag as regression (default 5%)
    """

    baseline = load_benchmarks(baseline_file)
    current = load_benchmarks(current_file)

    # Index benchmarks by problem_id for easy lookup
    baseline_dict = Dict(b["problem_id"] => b for b in baseline)
    current_dict = Dict(b["problem_id"] => b for b in current)

    println("\n" * "="^70)
    println("DAQP Performance Comparison")
    println("="^70)
    println("Baseline: $baseline_file")
    println("Current:  $current_file")
    println("Regression threshold: $regression_threshold%\n")

    # Track statistics
    regressions = []
    improvements = []
    unchanged = []
    missing_current = []
    new_in_current = []

    # Compare common benchmarks
    for problem_id in sort(collect(keys(baseline_dict)))
        if !haskey(current_dict, problem_id)
            push!(missing_current, problem_id)
            continue
        end

        base_row = baseline_dict[problem_id]
        curr_row = current_dict[problem_id]

        # Parse values
        base_setup = parse_float(base_row["setup_time_median_s"])
        curr_setup = parse_float(curr_row["setup_time_median_s"])
        base_solve = parse_float(base_row["solve_time_median_s"])
        curr_solve = parse_float(curr_row["solve_time_median_s"])
        base_total = parse_float(base_row["total_time_median_s"])
        curr_total = parse_float(curr_row["total_time_median_s"])
        base_iters = parse_float(base_row["iter_median"])
        curr_iters = parse_float(curr_row["iter_median"])

        # Calculate percentage change (positive = slower)
        if base_total === nothing || curr_total === nothing
            continue
        end

        pct_change = ((curr_total - base_total) / base_total) * 100

        problem_desc = "$(base_row["problem_type"]): n=$(base_row["n_variables"]), m=$(base_row["n_constraints"])"

        if pct_change > regression_threshold
            push!(regressions, (
                problem_id=problem_id,
                desc=problem_desc,
                base_setup=base_setup,
                curr_setup=curr_setup,
                base_solve=base_solve,
                curr_solve=curr_solve,
                base_total=base_total,
                curr_total=curr_total,
                pct_change=pct_change,
                base_iters=base_iters,
                curr_iters=curr_iters
            ))
        elseif pct_change < -regression_threshold
            push!(improvements, (
                problem_id=problem_id,
                desc=problem_desc,
                base_setup=base_setup,
                curr_setup=curr_setup,
                base_solve=base_solve,
                curr_solve=curr_solve,
                base_total=base_total,
                curr_total=curr_total,
                pct_change=pct_change,
                base_iters=base_iters,
                curr_iters=curr_iters
            ))
        else
            push!(unchanged, (
                problem_id=problem_id,
                desc=problem_desc,
                base_setup=base_setup,
                curr_setup=curr_setup,
                base_solve=base_solve,
                curr_solve=curr_solve,
                base_total=base_total,
                curr_total=curr_total,
                pct_change=pct_change,
                base_iters=base_iters,
                curr_iters=curr_iters
            ))
        end
    end

    # Find new benchmarks in current
    for problem_id in sort(collect(keys(current_dict)))
        if !haskey(baseline_dict, problem_id)
            push!(new_in_current, problem_id)
        end
    end

    # Print results
    if !isempty(regressions)
        println("⚠️  PERFORMANCE REGRESSIONS (>$(regression_threshold)% slower):")
        println(repeat("-", 70))
        for reg in sort(regressions, by=x -> -x.pct_change)
            println("  $(reg.desc)")
            baseline_total_ms = reg.base_total * 1000
            current_total_ms = reg.curr_total * 1000
            println("    Total: $(round(baseline_total_ms; digits=3)) ms → $(round(current_total_ms; digits=3)) ms ($(round(reg.pct_change; digits=1))%)")

            if reg.base_setup !== nothing && reg.curr_setup !== nothing
                setup_change = ((reg.curr_setup - reg.base_setup) / reg.base_setup) * 100
                println("    Setup: $(round(reg.base_setup*1e6; digits=1))µs → $(round(reg.curr_setup*1e6; digits=1))µs ($(round(setup_change; digits=1))%)")
            end
            if reg.base_solve !== nothing && reg.curr_solve !== nothing
                solve_change = ((reg.curr_solve - reg.base_solve) / reg.base_solve) * 100
                println("    Solve: $(round(reg.base_solve*1e6; digits=1))µs → $(round(reg.curr_solve*1e6; digits=1))µs ($(round(solve_change; digits=1))%)")
            end
            if reg.base_iters !== nothing && reg.curr_iters !== nothing
                iter_change = ((reg.curr_iters - reg.base_iters) / reg.base_iters) * 100
                println("    Iters: $(round(reg.base_iters; digits=1)) → $(round(reg.curr_iters; digits=1)) ($(round(iter_change; digits=1))%)")
            end
        end
        println()
    end

    if !isempty(improvements)
        println("✓ PERFORMANCE IMPROVEMENTS (>$(regression_threshold)% faster):")
        println(repeat("-", 70))
        for imp in sort(improvements, by=x -> x.pct_change)
            println("  $(imp.desc)")
            baseline_total_ms = imp.base_total * 1000
            current_total_ms = imp.curr_total * 1000
            println("    Total: $(round(baseline_total_ms; digits=3)) ms → $(round(current_total_ms; digits=3)) ms ($(round(imp.pct_change; digits=1))%)")

            if imp.base_setup !== nothing && imp.curr_setup !== nothing
                setup_change = ((imp.curr_setup - imp.base_setup) / imp.base_setup) * 100
                println("    Setup: $(round(imp.base_setup*1e6; digits=1))µs → $(round(imp.curr_setup*1e6; digits=1))µs ($(round(setup_change; digits=1))%)")
            end
            if imp.base_solve !== nothing && imp.curr_solve !== nothing
                solve_change = ((imp.curr_solve - imp.base_solve) / imp.base_solve) * 100
                println("    Solve: $(round(imp.base_solve*1e6; digits=1))µs → $(round(imp.curr_solve*1e6; digits=1))µs ($(round(solve_change; digits=1))%)")
            end
            if imp.base_iters !== nothing && imp.curr_iters !== nothing
                iter_change = ((imp.curr_iters - imp.base_iters) / imp.base_iters) * 100
                println("    Iters: $(round(imp.base_iters; digits=1)) → $(round(imp.curr_iters; digits=1)) ($(round(iter_change; digits=1))%)")
            end
        end
        println()
    end

    if !isempty(unchanged)
        println("≈ UNCHANGED (within ±$(regression_threshold)%):")
        println(repeat("-", 70))
        for unch in sort(unchanged, by=x -> abs(x.pct_change), rev=true)
            time_str = "$(round(unch.pct_change; digits=2))%"
            details = []

            if unch.base_setup !== nothing && unch.curr_setup !== nothing
                setup_change = ((unch.curr_setup - unch.base_setup) / unch.base_setup) * 100
                push!(details, "Setup: $(round(setup_change; digits=1))%")
            end
            if unch.base_solve !== nothing && unch.curr_solve !== nothing
                solve_change = ((unch.curr_solve - unch.base_solve) / unch.base_solve) * 100
                push!(details, "Solve: $(round(solve_change; digits=1))%")
            end
            if unch.base_iters !== nothing && unch.curr_iters !== nothing
                iter_change = ((unch.curr_iters - unch.base_iters) / unch.base_iters) * 100
                push!(details, "Iters: $(round(iter_change; digits=1))%")
            end

            detail_str = isempty(details) ? "" : " | " * join(details, ", ")
            println("  $(unch.desc): $(time_str)$(detail_str)")
        end
        println()
    end

    if !isempty(missing_current)
        println("⚠️  MISSING IN CURRENT RESULTS:")
        println(repeat("-", 70))
        for missing in missing_current
            println("  $(missing)")
        end
        println()
    end

    if !isempty(new_in_current)
        println("★ NEW BENCHMARKS IN CURRENT:")
        println(repeat("-", 70))
        for new in new_in_current
            println("  $(new)")
        end
        println()
    end

    # Summary
    println("="^70)
    println("SUMMARY:")
    println("  Total compared: $(length(unchanged) + length(regressions) + length(improvements))")
    println("  Regressions:    $(length(regressions))")
    println("  Improvements:   $(length(improvements))")
    println("  Unchanged:      $(length(unchanged))")
    println("  Missing:        $(length(missing_current))")
    println("  New:            $(length(new_in_current))")

    if !isempty(regressions)
        println("\n⚠️  WARNING: $(length(regressions)) performance regression(s) detected!")
        return false
    else
        println("\n✓ No regressions detected.")
        return true
    end
end

function print_benchmark_info(filename::String)
    """Print summary information about a benchmark file."""
    if !isfile(filename)
        println("File not found: $filename")
        return
    end

    data = load_benchmarks(filename)

    if isempty(data)
        println("Benchmark file: $filename (empty)")
        return
    end

    # Get version and timestamp from first row
    first_row = data[1]
    version = get(first_row, "daqp_version", "unknown")
    timestamp = get(first_row, "timestamp", "unknown")

    println("\nBenchmark file: $filename")
    println("  Version:   $version")
    println("  Timestamp: $timestamp")
    println("  Benchmarks: $(length(data))")

    # Summary stats
    times = [parse_float(b["total_time_median_s"]) for b in data]
    times = filter(x -> x !== nothing, times)
    if !isempty(times)
        println("  Time range: $(round(minimum(times)*1000; digits=3)) - $(round(maximum(times)*1000; digits=3)) ms")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Parse command line arguments
    if length(ARGS) < 1
        println("""
        Usage:
            julia compare_benchmarks.jl baseline.csv current.csv
            julia compare_benchmarks.jl --baseline baseline.csv current.csv --threshold 5
        """)
        exit(1)
    end

    local baseline_file = ""
    local current_file = ""
    local threshold = 5.0

    local i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--baseline" && i < length(ARGS)
            baseline_file = ARGS[i+1]
            i += 2
        elseif ARGS[i] == "--threshold" && i < length(ARGS)
            threshold = parse(Float64, ARGS[i+1])
            i += 2
        elseif baseline_file == ""
            baseline_file = ARGS[i]
            i += 1
        elseif current_file == ""
            current_file = ARGS[i]
            i += 1
        else
            i += 1
        end
    end

    # Validate arguments
    if baseline_file == "" || current_file == ""
        @error "Must provide both baseline and current benchmark files"
        exit(1)
    end

    if !isfile(baseline_file)
        @error "Baseline file not found: $baseline_file"
        exit(1)
    end

    if !isfile(current_file)
        @error "Current file not found: $current_file"
        exit(1)
    end

    print_benchmark_info(baseline_file)
    print_benchmark_info(current_file)

    success = compare_benchmarks(baseline_file, current_file; regression_threshold=threshold)
    exit(success ? 0 : 1)
end
