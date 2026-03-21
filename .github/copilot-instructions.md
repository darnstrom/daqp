# DAQP – Copilot Instructions

DAQP is a library-free C implementation of a dual active-set solver for convex quadratic programs (QP), mixed-integer QPs (via branch-and-bound), hierarchical QPs, and AVI problems. It exposes interfaces to Python, Julia, MATLAB, and C++ (Eigen).

## Build

```bash
# Core library only
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

# With Eigen C++ interface and tests
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DEIGEN=ON -DTEST=ON
cmake --build build

# Other interface flags: -DJULIA=ON, -DMATLAB=ON, -DPYTHON=ON
# Other option flags:    -DPROFILING=ON (default), -DFASTER_MATH=ON (default), -DSOFT_WEIGHTS=OFF (default)
```

Produces `daqp` (shared) and `daqpstat` (static) libraries.

## Tests

```bash
# All tests (after building with -DTEST=ON)
cd build && ctest --output-on-failure

# Single Eigen/C++ test
cd build && ctest -R test_00_basic_qp -V

# Python (no CMake required)
cd interfaces/daqp-python && pip install .
python -m pytest test/example_test.py -v
# Single Python test
python -m pytest test/example_test.py::Testing::test_python_demo -v

# Julia (after building with -DJULIA=ON)
julia --color=yes --project=build/interfaces/daqp-julia -e "using Pkg; Pkg.test()"
```

## Architecture

The entry point for all solvers is `src/api.c`: `daqp_solve()` selects the algorithm based on workspace configuration, then calls into:

| Module | File | Role |
|---|---|---|
| LDP solver | `src/daqp.c` | Core dual active-set iteration |
| Factorization | `src/factorization.c` | Recursive LDL^T add/remove updates |
| Auxiliary | `src/auxiliary.c` | CSP computation, constraint add/remove |
| Branch & Bound | `src/bnb.c` | MIQP via BnB over binary constraints |
| Hierarchical | `src/hierarchical.c` | Multi-priority HiQP |
| AVI | `src/avi.c` | Absolute value inequality solver |
| Proximal | `src/daqp_prox.c` | Proximal-point algorithm variant |
| Codegen | `codegen/codegen.c` | Emit standalone C for embedded use |

Key data structures are in `include/types.h`: `DAQPProblem` (problem data), `DAQPWorkspace` (solver state, LDL factors, working set), `DAQPSettings` (tolerances), `DAQPResult` (output).

The typical call sequence is:
1. `setup_daqp(qp, work, ...)` — allocate workspace and factor the Hessian
2. `daqp_solve(res, work)` — run the solver
3. Optionally call `setup_daqp_bnb` / `setup_daqp_hiqp` / `setup_daqp_avi` before step 2 for specialized variants

## Conventions

**Naming**
- Functions: `daqp_` prefix, `snake_case` — e.g., `daqp_ldp`, `daqp_solve`, `daqp_bnb`
- Macros/constants: `DAQP_` prefix, `UPPER_CASE` — e.g., `DAQP_EXIT_OPTIMAL`, `DAQP_INF`
- Types: `DAQP` prefix, `PascalCase` — e.g., `DAQPWorkspace`, `DAQPResult`
- Floating-point type is always `c_float` (defaults to `double`; define `DAQP_SINGLE_PRECISION` for `float`)

**Constraint state** is encoded as bit flags in the `sense` array (one `int` per constraint). Use the provided macros from `include/constants.h` rather than manipulating bits directly:
```c
DAQP_IS_ACTIVE(i), DAQP_SET_ACTIVE(i), DAQP_SET_INACTIVE(i)
DAQP_IS_LOWER(i),  DAQP_SET_LOWER(i),  DAQP_SET_UPPER(i)
DAQP_IS_SOFT(i),   DAQP_SET_SOFT(i)
DAQP_IS_BINARY(i)
DAQP_IS_IMMUTABLE(i), DAQP_SET_IMMUTABLE(i)
```

**Exit flags** (returned in `DAQPResult.exitflag` and as function return values):
- `DAQP_EXIT_OPTIMAL (1)`, `DAQP_EXIT_SOFT_OPTIMAL (2)` — success
- `DAQP_EXIT_INFEASIBLE (-1)`, `DAQP_EXIT_CYCLE (-2)`, `DAQP_EXIT_UNBOUNDED (-3)`, `DAQP_EXIT_ITERLIMIT (-4)`, `DAQP_EXIT_TIMELIMIT (-7)` — failure

**Memory layout**: constraints are stored as `[simple bounds (ms); general constraints (m-ms)]`. `blower`/`bupper` are flat arrays of length `m` covering both. `A` has `m - ms` rows.

**Profiling**: timing code is guarded by `#ifdef PROFILING`. The `PROFILING` flag is on by default; it is how `solve_time` and `setup_time` are populated in `DAQPResult`.

**Interfaces** live in `interfaces/daqp-{python,julia,matlab,eigen,simulink}/` and are self-contained (own build scripts, tests, examples). The Python interface is a Cython wrapper (`daqp.pyx`). The Julia interface uses `ccall` into the shared library.
