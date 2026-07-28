#!/bin/bash
# Compare the current DAQP library against a git ref on Maros-Meszaros dense QPs.

set -euo pipefail

if [ "$#" -lt 4 ]; then
    echo "usage: $0 <current-build> <base-ref> <maros-data-dir> <output-dir> [threshold] [repeats]"
    exit 2
fi

CURRENT_BUILD="$(cd "$1" && pwd)"
BASE_REF="$2"
MAROS_DATA="$(cd "$3" && pwd)"
OUTPUT_DIR="$(mkdir -p "$4" && cd "$4" && pwd)"
THRESHOLD="${5:-25}"
REPEATS="${6:-3}"
REPO_ROOT="$(git rev-parse --show-toplevel)"
SCRIPT_DIR="$REPO_ROOT/.github/benchmarks"
EXPORTED_DIR="$OUTPUT_DIR/exported"
CURRENT_LIB="$CURRENT_BUILD/interfaces/daqp-julia/libdaqp.so"

if [ ! -f "$CURRENT_LIB" ]; then
    echo "Current libdaqp.so not found at $CURRENT_LIB"
    exit 1
fi

python3 "$SCRIPT_DIR/export_maros_meszaros.py" "$MAROS_DATA" "$EXPORTED_DIR"

cc -O3 -I"$REPO_ROOT/include" "$SCRIPT_DIR/maros_meszaros_runner.c" \
    -L"$(dirname "$CURRENT_LIB")" -ldaqp -lm -o "$OUTPUT_DIR/maros_meszaros_runner"

echo "Benchmarking current library"
LD_LIBRARY_PATH="$(dirname "$CURRENT_LIB")${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
    "$OUTPUT_DIR/maros_meszaros_runner" "$EXPORTED_DIR" \
    "$OUTPUT_DIR/current.csv" "$REPEATS"

TEMP_REPO="$(mktemp -d)"
trap 'rm -rf "$TEMP_REPO"' EXIT
git clone --no-checkout "$REPO_ROOT" "$TEMP_REPO/daqp_base" >/dev/null 2>&1
git -C "$TEMP_REPO/daqp_base" checkout --detach "$BASE_REF" >/dev/null
cmake -S "$TEMP_REPO/daqp_base" -B "$TEMP_REPO/daqp_base/build" \
    -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build "$TEMP_REPO/daqp_base/build" -- -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)" >/dev/null

BASE_LIB="$(find "$TEMP_REPO/daqp_base/build" -name libdaqp.so -type f -print -quit)"
if [ -z "$BASE_LIB" ]; then
    echo "Baseline build did not produce libdaqp.so"
    exit 1
fi

echo "Benchmarking baseline library $BASE_REF"
LD_LIBRARY_PATH="$(dirname "$BASE_LIB")${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
    "$OUTPUT_DIR/maros_meszaros_runner" "$EXPORTED_DIR" \
    "$OUTPUT_DIR/master.csv" "$REPEATS"

COMPARE_RESULT=0
python3 "$SCRIPT_DIR/compare_maros_meszaros.py" \
    "$OUTPUT_DIR/master.csv" "$OUTPUT_DIR/current.csv" \
    "$OUTPUT_DIR/comparison.md" --threshold "$THRESHOLD" || COMPARE_RESULT=$?
exit "$COMPARE_RESULT"
