# DAQP C Benchmark

A self-contained C benchmark for measuring DAQP solve performance.

## Build

From the repository root, after building DAQP with CMake:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFASTER_MATH=ON
cmake --build build

gcc -O3 -fassociative-math -fno-signed-zeros -fno-trapping-math \
    -I include -I codegen \
    benchmark/bench_qp.c build/libdaqpstat.a -lm -o bench_qp
./bench_qp
```

## What it measures

Random dense QPs with a diagonally-dominant Hessian, at five problem sizes
(n = 10, 20, 40, 60, 80 decision variables).  Each size runs 3000 solves per
pass and 7 passes; the best-of-7 is reported to reduce OS scheduling noise.

## Comparing two builds

Build each library separately and link the benchmark against each:

```bash
# baseline (e.g. v0.8.0)
git stash && git checkout v0.8.0
cmake -S . -B build_v080 -DCMAKE_BUILD_TYPE=Release -DFASTER_MATH=ON
cmake --build build_v080
gcc -O3 -fassociative-math -fno-signed-zeros -fno-trapping-math \
    -I include -I codegen \
    benchmark/bench_qp.c build_v080/libdaqpstat.a -lm -o bench_v080
git checkout -

# current branch
cmake -S . -B build_head -DCMAKE_BUILD_TYPE=Release -DFASTER_MATH=ON
cmake --build build_head
gcc -O3 -fassociative-math -fno-signed-zeros -fno-trapping-math \
    -I include -I codegen \
    benchmark/bench_qp.c build_head/libdaqpstat.a -lm -o bench_head

./bench_v080
./bench_head
```
