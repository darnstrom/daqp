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

"""
    pct_diff(base, curr)

Percentage change from `base` to `curr` (positive = larger), or `nothing` if
either value is missing or the baseline is zero.
"""
function pct_diff(base, curr)
    (base === nothing || curr === nothing || base == 0) && return nothing
    return ((curr - base) / base) * 100
end

"""
    has_branching(entry)

Whether the entry comes from a problem that actually branched. Only MIQPs visit
more than one node, so this keeps the node counts out of the other reports.
"""
function has_branching(entry)
    entry.base_nodes === nothing && return false
    entry.curr_nodes === nothing && return false
    return max(entry.base_nodes, entry.curr_nodes) > 1
end

"""
    print_entry_details(e)

Print the multi-line before/after breakdown used for regressions and improvements.
"""
function print_entry_details(e)
    println("  $(e.desc)")
    println("    Total: $(round(e.base_total*1000; digits=3)) ms → $(round(e.curr_total*1000; digits=3)) ms ($(round(e.pct_change; digits=1))%)")

    setup_change = pct_diff(e.base_setup, e.curr_setup)
    if setup_change !== nothing
        println("    Setup: $(round(e.base_setup*1e6; digits=1))µs → $(round(e.curr_setup*1e6; digits=1))µs ($(round(setup_change; digits=1))%)")
    end
    solve_change = pct_diff(e.base_solve, e.curr_solve)
    if solve_change !== nothing
        println("    Solve: $(round(e.base_solve*1e6; digits=1))µs → $(round(e.curr_solve*1e6; digits=1))µs ($(round(solve_change; digits=1))%)")
    end
    iter_change = pct_diff(e.base_iters, e.curr_iters)
    if iter_change !== nothing
        println("    Iters: $(round(e.base_iters; digits=1)) → $(round(e.curr_iters; digits=1)) ($(round(iter_change; digits=1))%)")
    end
    node_change = pct_diff(e.base_nodes, e.curr_nodes)
    if node_change !== nothing && has_branching(e)
        println("    Nodes: $(round(e.base_nodes; digits=1)) → $(round(e.curr_nodes; digits=1)) ($(round(node_change; digits=1))%)")
    end
end

function compare_benchmarks(baseline_file::String, current_file::String;
                           regression_threshold::Float64=5.0,
                           work_threshold::Float64=5.0)
    """
    Compare two benchmark result files and report regressions.

    regression_threshold: percentage slowdown in wall time to flag as a
                          regression (default 5%). Timings on shared machines
                          are noisy, so this is the coarse criterion.
    work_threshold:       percentage increase in iterations or branch and bound
                          nodes to flag as a regression (default 5%). Both runs
                          solve identical problems from a seeded generator, so
                          these counts carry no measurement noise at all.
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
    println("Time regression threshold: $regression_threshold%")
    println("Work (iterations/nodes) regression threshold: $work_threshold%\n")

    # Track statistics
    regressions = []
    work_regressions = []
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
        # nodes_median is absent from CSVs produced before MIQPs were benchmarked
        base_nodes = parse_float(get(base_row, "nodes_median", ""))
        curr_nodes = parse_float(get(curr_row, "nodes_median", ""))

        # Calculate percentage change (positive = slower)
        if base_total === nothing || curr_total === nothing
            continue
        end

        pct_change = ((curr_total - base_total) / base_total) * 100

        problem_desc = "$(base_row["problem_type"]): n=$(base_row["n_variables"]), m=$(base_row["n_constraints"])"

        entry = (
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
            curr_iters=curr_iters,
            base_nodes=base_nodes,
            curr_nodes=curr_nodes
        )

        iter_change = pct_diff(base_iters, curr_iters)
        node_change = pct_diff(base_nodes, curr_nodes)
        if (iter_change !== nothing && iter_change > work_threshold) ||
           (has_branching(entry) && node_change !== nothing && node_change > work_threshold)
            push!(work_regressions, entry)
        end

        if pct_change > regression_threshold
            push!(regressions, entry)
        elseif pct_change < -regression_threshold
            push!(improvements, entry)
        else
            push!(unchanged, entry)
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
            print_entry_details(reg)
        end
        println()
    end

    if !isempty(improvements)
        println("✓ PERFORMANCE IMPROVEMENTS (>$(regression_threshold)% faster):")
        println(repeat("-", 70))
        for imp in sort(improvements, by=x -> x.pct_change)
            print_entry_details(imp)
        end
        println()
    end

    if !isempty(unchanged)
        println("≈ UNCHANGED (within ±$(regression_threshold)%):")
        println(repeat("-", 70))
        for unch in sort(unchanged, by=x -> abs(x.pct_change), rev=true)
            time_str = "$(round(unch.pct_change; digits=2))%"
            details = String[]

            setup_change = pct_diff(unch.base_setup, unch.curr_setup)
            if setup_change !== nothing
                push!(details, "Setup: $(round(setup_change; digits=1))%")
            end
            solve_change = pct_diff(unch.base_solve, unch.curr_solve)
            if solve_change !== nothing
                push!(details, "Solve: $(round(solve_change; digits=1))%")
            end
            iter_change = pct_diff(unch.base_iters, unch.curr_iters)
            if iter_change !== nothing
                push!(details, "Iters: $(round(iter_change; digits=1))%")
            end
            node_change = pct_diff(unch.base_nodes, unch.curr_nodes)
            if node_change !== nothing && has_branching(unch)
                push!(details, "Nodes: $(round(node_change; digits=1))%")
            end

            detail_str = isempty(details) ? "" : " | " * join(details, ", ")
            println("  $(unch.desc): $(time_str)$(detail_str)")
        end
        println()
    end

    if !isempty(work_regressions)
        println("⚠️  MORE WORK PER SOLVE (>$(work_threshold)% more iterations/nodes):")
        println(repeat("-", 70))
        println("  These counts are exact -- both runs solve identical problems.")
        for wr in sort(work_regressions, by=x -> -something(pct_diff(x.base_iters, x.curr_iters), 0.0))
            println("  $(wr.desc)")
            iter_change = pct_diff(wr.base_iters, wr.curr_iters)
            if iter_change !== nothing
                println("    Iters: $(round(wr.base_iters; digits=1)) → $(round(wr.curr_iters; digits=1)) ($(round(iter_change; digits=1))%)")
            end
            node_change = pct_diff(wr.base_nodes, wr.curr_nodes)
            if node_change !== nothing && has_branching(wr)
                println("    Nodes: $(round(wr.base_nodes; digits=1)) → $(round(wr.curr_nodes; digits=1)) ($(round(node_change; digits=1))%)")
            end
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
    println("  Total compared:  $(length(unchanged) + length(regressions) + length(improvements))")
    println("  Time regressions: $(length(regressions))")
    println("  Work regressions: $(length(work_regressions))")
    println("  Improvements:     $(length(improvements))")
    println("  Unchanged:        $(length(unchanged))")
    println("  Missing:          $(length(missing_current))")
    println("  New:              $(length(new_in_current))")

    if !isempty(regressions) || !isempty(work_regressions)
        msgs = String[]
        isempty(regressions)      || push!(msgs, "$(length(regressions)) time regression(s)")
        isempty(work_regressions) || push!(msgs, "$(length(work_regressions)) work regression(s)")
        println("\n⚠️  WARNING: " * join(msgs, " and ") * " detected!")
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
            julia compare_benchmarks.jl baseline.csv current.csv --threshold 25 --work-threshold 5

        --threshold       percentage slowdown in wall time to flag (noisy measure)
        --work-threshold  percentage increase in iterations/nodes to flag (exact measure)
        """)
        exit(1)
    end

    local baseline_file = ""
    local current_file = ""
    local threshold = 5.0
    local work_threshold = 5.0

    local i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--baseline" && i < length(ARGS)
            baseline_file = ARGS[i+1]
            i += 2
        elseif ARGS[i] == "--threshold" && i < length(ARGS)
            threshold = parse(Float64, ARGS[i+1])
            i += 2
        elseif ARGS[i] == "--work-threshold" && i < length(ARGS)
            work_threshold = parse(Float64, ARGS[i+1])
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

    success = compare_benchmarks(baseline_file, current_file;
                                 regression_threshold=threshold,
                                 work_threshold=work_threshold)
    exit(success ? 0 : 1)
end
