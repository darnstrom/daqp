#!/bin/bash
# benchmark_comparison_git.sh
# CTest helper script for performance benchmark comparison using git versions
# Builds and compares local libdaqp from different git tags

set -e

JULIA_PROJECT="$1"
GIT_TAG="${2:-v0.7.1}"
BENCHMARK_SUITE="${3:-small}"
REGRESSION_THRESHOLD="${4:-5}"
OUTPUT_DIR="${5:-.}"

# Convert to absolute paths
JULIA_PROJECT="$(cd "$JULIA_PROJECT" && pwd)"
OUTPUT_DIR="$(mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR" && pwd)"
BENCHMARK_SCRIPT="$JULIA_PROJECT/test/benchmark.jl"
COMPARE_SCRIPT="$JULIA_PROJECT/test/compare_benchmarks.jl"

# Get the current git repo root
REPO_ROOT="$(git rev-parse --show-toplevel)"

echo "=========================================="
echo "DAQP Git-Based Benchmark Comparison"
echo "=========================================="
echo "Julia project: $JULIA_PROJECT"
echo "Comparison tag: $GIT_TAG"
echo "Benchmark suite: $BENCHMARK_SUITE"
echo "Regression threshold: $REGRESSION_THRESHOLD%"
echo "Output directory: $OUTPUT_DIR"
echo ""

# The benchmark script looks for libdaqp.so in pkgdir(DAQPBase), which is the Julia project dir
CURRENT_LIBDAQP="$JULIA_PROJECT/libdaqp.so"
CURRENT_LIBDAQP_BAK="$OUTPUT_DIR/libdaqp_current.so"
OLD_LIBDAQP_BAK="$OUTPUT_DIR/libdaqp_$GIT_TAG.so"

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

# Step 2: Clone and build the old version
echo ""
echo "Step 2: Cloning and building version $GIT_TAG..."
TEMP_REPO=$(mktemp -d)
trap "rm -rf $TEMP_REPO; cp '$CURRENT_LIBDAQP_BAK' '$CURRENT_LIBDAQP' 2>/dev/null || true" EXIT

cd "$TEMP_REPO"
git clone --depth 1 --branch "$GIT_TAG" "$REPO_ROOT" daqp_old > /dev/null 2>&1 || {
    echo "ERROR: Failed to clone DAQP at tag $GIT_TAG"
    exit 1
}

cd daqp_old
# Build libdaqp with CMake
mkdir -p build
cd build
cmake .. > /dev/null 2>&1 || {
    echo "ERROR: CMake configuration failed for version $GIT_TAG"
    exit 1
}
make > /dev/null 2>&1 || {
    echo "ERROR: Build failed for version $GIT_TAG"
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
    echo "ERROR: Could not find built libdaqp.so in version $GIT_TAG"
    exit 1
fi

# Copy it to output and then to the current project location
cp "$OLD_LIBDAQP_BUILT" "$OLD_LIBDAQP_BAK"
cp "$OLD_LIBDAQP_BAK" "$CURRENT_LIBDAQP"

# Go back to original repo
cd "$REPO_ROOT"

# Step 3: Benchmark old version
echo "Step 3: Running benchmarks with version $GIT_TAG..."
julia --project="$JULIA_PROJECT" "$BENCHMARK_SCRIPT" \
    --suite "$BENCHMARK_SUITE" \
    --output "$OUTPUT_DIR/old_version.csv" \
    --local \
    2>&1 | grep -E "(Starting|Mean|Results|Using)" || true

if [ ! -f "$OUTPUT_DIR/old_version.csv" ]; then
    echo "ERROR: Failed to generate old version benchmark"
    exit 1
fi

# Restore current libdaqp.so
cp "$CURRENT_LIBDAQP_BAK" "$CURRENT_LIBDAQP"

# Step 4: Compare results
echo ""
echo "Step 4: Comparing results (current vs $GIT_TAG)..."
echo "=================================================================="
julia --project="$JULIA_PROJECT" "$COMPARE_SCRIPT" \
    "$OUTPUT_DIR/old_version.csv" \
    "$OUTPUT_DIR/current_dev.csv" \
    --threshold "$REGRESSION_THRESHOLD"

COMPARE_RESULT=$?

echo ""
echo "Results saved in: $OUTPUT_DIR"
echo "  - current_dev.csv"
echo "  - old_version.csv"

if [ $COMPARE_RESULT -eq 0 ]; then
    echo ""
    echo "✓ Test PASSED: No regressions detected"
    exit 0
else
    echo ""
    echo "⚠️  Test FAILED: Performance regressions detected"
    exit 1
fi
