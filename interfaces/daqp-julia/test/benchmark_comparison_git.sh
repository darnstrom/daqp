#!/bin/bash
# benchmark_comparison_git.sh
# Helper script for performance benchmark comparison using git versions.
# Builds libdaqp from a baseline git ref (tag, branch or SHA) and compares it
# against the currently built libdaqp.
#
# Only the shared library is swapped -- the Julia project, benchmark script and
# problem generators always come from the working tree, so both runs solve the
# exact same problems and the CSV problem_ids line up.
#
# Usage:
#   benchmark_comparison_git.sh <julia_project> [ref] [suite] [threshold] [outdir] [fail_on_regression] [work_threshold]

set -e

JULIA_PROJECT="$1"
GIT_REF="${2:-v0.7.1}"
BENCHMARK_SUITE="${3:-small}"
REGRESSION_THRESHOLD="${4:-5}"
OUTPUT_DIR="${5:-.}"
FAIL_ON_REGRESSION="${6:-1}"
WORK_THRESHOLD="${7:-5}"

# Convert to absolute paths
JULIA_PROJECT="$(cd "$JULIA_PROJECT" && pwd)"
OUTPUT_DIR="$(mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR" && pwd)"
BENCHMARK_SCRIPT="$JULIA_PROJECT/test/benchmark.jl"
COMPARE_SCRIPT="$JULIA_PROJECT/test/compare_benchmarks.jl"

# Get the current git repo root
REPO_ROOT="$(git rev-parse --show-toplevel)"

# Safe-for-filenames version of the ref (refs may contain '/')
GIT_REF_SLUG="$(echo "$GIT_REF" | tr '/' '_')"

echo "=========================================="
echo "DAQP Git-Based Benchmark Comparison"
echo "=========================================="
echo "Julia project: $JULIA_PROJECT"
echo "Baseline ref: $GIT_REF"
echo "Benchmark suite: $BENCHMARK_SUITE"
echo "Time regression threshold: $REGRESSION_THRESHOLD%"
echo "Work regression threshold: $WORK_THRESHOLD%"
echo "Output directory: $OUTPUT_DIR"
echo ""

# The benchmark script looks for libdaqp.so in pkgdir(DAQPBase), which is the Julia project dir
CURRENT_LIBDAQP="$JULIA_PROJECT/libdaqp.so"
CURRENT_LIBDAQP_BAK="$OUTPUT_DIR/libdaqp_current.so"
OLD_LIBDAQP_BAK="$OUTPUT_DIR/libdaqp_$GIT_REF_SLUG.so"

# Step 1: Benchmark current version
echo "Step 1: Benchmarking current version..."
if [ ! -f "$CURRENT_LIBDAQP" ]; then
    echo "ERROR: Current libdaqp.so not found at $CURRENT_LIBDAQP"
    echo "Make sure the current version is built with CMake (make) first"
    exit 1
fi

julia --project="$JULIA_PROJECT" "$BENCHMARK_SCRIPT" \
    --suite "$BENCHMARK_SUITE" \
    --output "$OUTPUT_DIR/current_dev.csv" \
    --local \
    2>&1 | grep -E "(Starting|Mean|Results|Using)" || true

if [ ! -f "$OUTPUT_DIR/current_dev.csv" ]; then
    echo "ERROR: Failed to generate current version benchmark"
    exit 1
fi

# Backup current libdaqp.so
cp "$CURRENT_LIBDAQP" "$CURRENT_LIBDAQP_BAK"

# Step 2: Clone and build the baseline version
echo ""
echo "Step 2: Cloning and building baseline $GIT_REF..."
TEMP_REPO=$(mktemp -d)
trap "rm -rf $TEMP_REPO; cp '$CURRENT_LIBDAQP_BAK' '$CURRENT_LIBDAQP' 2>/dev/null || true" EXIT

cd "$TEMP_REPO"
# Clone the local repo (cheap, shares objects) so that tags, branches and bare
# SHAs can all be checked out through the same code path.
git clone --no-checkout "$REPO_ROOT" daqp_old > /dev/null 2>&1 || {
    echo "ERROR: Failed to clone DAQP from $REPO_ROOT"
    exit 1
}

cd daqp_old
if ! git checkout --detach "$GIT_REF" > /dev/null 2>&1; then
    # The ref may not be advertised (e.g. a bare SHA only reachable from a
    # remote-tracking branch) -- fetch it explicitly before giving up.
    if git fetch --depth=1 origin "$GIT_REF" > /dev/null 2>&1 &&
       git checkout --detach FETCH_HEAD > /dev/null 2>&1; then
        :
    else
        echo "ERROR: Failed to check out ref $GIT_REF"
        echo "       (in CI, make sure the checkout has enough history, e.g. fetch-depth: 0)"
        exit 1
    fi
fi
echo "  Baseline commit: $(git rev-parse --short HEAD) $(git log -1 --pretty=%s)"

# Build libdaqp with CMake, using the same build type as the current build
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1 || {
    echo "ERROR: CMake configuration failed for baseline $GIT_REF"
    exit 1
}
cmake --build . -- -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)" > /dev/null 2>&1 || {
    echo "ERROR: Build failed for baseline $GIT_REF"
    exit 1
}

# Find the built libdaqp - check both build dir and interfaces subdir
if [ -f "./libdaqp.so" ]; then
    OLD_LIBDAQP_BUILT="./libdaqp.so"
elif [ -f "./interfaces/daqp-julia/libdaqp.so" ]; then
    OLD_LIBDAQP_BUILT="./interfaces/daqp-julia/libdaqp.so"
else
    OLD_LIBDAQP_BUILT=$(find . -name "libdaqp.so" -type f 2>/dev/null | head -1)
fi

if [ -z "$OLD_LIBDAQP_BUILT" ] || [ ! -f "$OLD_LIBDAQP_BUILT" ]; then
    echo "ERROR: Could not find built libdaqp.so for baseline $GIT_REF"
    exit 1
fi

# Copy it to output and then to the current project location
cp "$OLD_LIBDAQP_BUILT" "$OLD_LIBDAQP_BAK"
cp "$OLD_LIBDAQP_BAK" "$CURRENT_LIBDAQP"

# Go back to original repo
cd "$REPO_ROOT"

# Step 3: Benchmark baseline version
echo "Step 3: Running benchmarks with baseline $GIT_REF..."
julia --project="$JULIA_PROJECT" "$BENCHMARK_SCRIPT" \
    --suite "$BENCHMARK_SUITE" \
    --output "$OUTPUT_DIR/old_version.csv" \
    --local \
    2>&1 | grep -E "(Starting|Mean|Results|Using)" || true

if [ ! -f "$OUTPUT_DIR/old_version.csv" ]; then
    echo "ERROR: Failed to generate baseline benchmark"
    exit 1
fi

# Restore current libdaqp.so
cp "$CURRENT_LIBDAQP_BAK" "$CURRENT_LIBDAQP"

# Step 4: Compare results
echo ""
echo "Step 4: Comparing results (current vs $GIT_REF)..."
echo "=================================================================="
# compare_benchmarks.jl exits non-zero on regression; capture that without
# tripping `set -e` so that the summary below is always printed.
COMPARE_RESULT=0
julia --project="$JULIA_PROJECT" "$COMPARE_SCRIPT" \
    "$OUTPUT_DIR/old_version.csv" \
    "$OUTPUT_DIR/current_dev.csv" \
    --threshold "$REGRESSION_THRESHOLD" \
    --work-threshold "$WORK_THRESHOLD" \
    > "$OUTPUT_DIR/comparison.txt" 2>&1 || COMPARE_RESULT=$?
cat "$OUTPUT_DIR/comparison.txt"

echo ""
echo "Results saved in: $OUTPUT_DIR"
echo "  - current_dev.csv"
echo "  - old_version.csv"
echo "  - comparison.txt"

if [ "$COMPARE_RESULT" -eq 0 ]; then
    echo ""
    echo "✓ Test PASSED: No regressions detected"
    exit 0
elif [ "$FAIL_ON_REGRESSION" != "1" ]; then
    echo ""
    echo "⚠️  Performance regressions detected (not failing: fail_on_regression=$FAIL_ON_REGRESSION)"
    exit 0
else
    echo ""
    echo "⚠️  Test FAILED: Performance regressions detected"
    exit 1
fi
