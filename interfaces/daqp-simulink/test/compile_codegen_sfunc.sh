#!/bin/sh
# Compile DAQP_sfunc.c as in generated code (i.e., without MATLAB_MEX_FILE)
# against the Simulink headers of a MATLAB installation. This only checks
# that the code compiles (no Simulink Coder license is needed).
#
# Usage: compile_codegen_sfunc.sh <matlabroot>
set -e
ROOT="$1"
HERE="$(cd "$(dirname "$0")" && pwd)"
SFUNC_DIR="$HERE/.."
DAQP_DIR="$HERE/../../.."
BUILD="$(mktemp -d)"

# rtwtypes.h is normally generated with the model => use the MATLAB types
RTWTYPES="$(find "$ROOT/rtw" "$ROOT/simulink/include" -name rtwtypes.h 2>/dev/null | head -n 1)"
if [ -z "$RTWTYPES" ]; then
    printf '#include "tmwtypes.h"\n' > "$BUILD/rtwtypes.h"
    echo "Using a stub rtwtypes.h"
else
    echo "Using $RTWTYPES"
fi

set -x
for src in DAQP_sfunc.c daqp_sfunc_core.c; do
    gcc -c -std=c99 -Wall -Werror=implicit-function-declaration \
        -DRT -DUSE_RTMODEL -DMODEL=daqp_codegen_check -DNUMST=1 -DNCSTATES=0 \
        -DHAVESTDIO -DMT=0 -DTID01EQ=0 -DCLASSIC_INTERFACE=0 -DALLOCATIONFCN=0 \
        -DMULTI_INSTANCE_CODE=0 -DINTEGER_CODE=0 \
        -I"$BUILD" -I"$ROOT/extern/include" -I"$ROOT/simulink/include" \
        -I"$ROOT/rtw/c/src" -I"$DAQP_DIR/include" -I"$SFUNC_DIR" \
        "$SFUNC_DIR/$src" -o "$BUILD/${src%.c}.o"
done
