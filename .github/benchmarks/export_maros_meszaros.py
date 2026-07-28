#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2022 Stéphane Caron
# Copyright 2023-2024 Inria
# Adapted from qpsolvers/maros_meszaros_qpbenchmark.
"""Export the Maros-Meszaros dense subset to DAQP's native dense layout."""

import argparse
import os
import struct

import numpy as np
import scipy.io as spio
import scipy.sparse as spa
from qpbenchmark.utils import is_posdef


DAQP_INF = 1e30
EQ_TOL = 1e-10


def load_mat(path):
    mat = spio.loadmat(path)
    P = mat["P"].astype(float).tocsc()
    q = mat["q"].T.flatten().astype(float)
    A = mat["A"].astype(float).tocsc()
    lower = mat["l"].T.flatten().astype(float)
    upper = mat["u"].T.flatten().astype(float)
    n = int(mat["n"].T.flatten()[0])
    m = int(mat["m"].T.flatten()[0])
    assert A.shape == (m, n)

    A = A.copy()
    A.data[A.data > 9e19] = np.inf
    A.data[A.data < -9e19] = -np.inf
    lower[lower > 9e19] = np.inf
    lower[lower < -9e19] = -np.inf
    upper[upper > 9e19] = np.inf
    upper[upper < -9e19] = -np.inf

    # The MAT files store A = vstack([C, eye(n)]).
    return P, q, A[:-n], lower[:-n], upper[:-n], lower[-n:], upper[-n:]


def converted_constraint_count(C, lower, upper, box_lower):
    """Match MarosMeszaros.count_constraints after its format conversion."""
    equal = upper - lower < EQ_TOL
    inequality_rows = np.asarray(np.logical_not(equal)).nonzero()
    G = spa.vstack([C[inequality_rows], -C[inequality_rows]], format="csc")
    h = np.hstack([upper[inequality_rows], -lower[inequality_rows]])
    finite = h < np.inf
    inequalities = G[finite].shape[0] if G.size > 0 else 0
    equalities = C[np.asarray(equal).nonzero()].shape[0]
    return inequalities + equalities + box_lower.shape[0]


def export_problem(output_dir, name, P, q, C, lower, upper, box_lower, box_upper):
    H = np.asarray(P.todense(), dtype=np.float64)
    H = 0.5 * (H + H.T)
    A = np.asarray(C.todense(), dtype=np.float64)
    blower = np.concatenate([box_lower, lower]).astype(np.float64)
    bupper = np.concatenate([box_upper, upper]).astype(np.float64)
    n = H.shape[0]
    m = blower.size
    sense = np.zeros(m, dtype=np.int32)

    equal = bupper - blower < EQ_TOL
    sense[equal] = 5  # DAQP_ACTIVE | DAQP_IMMUTABLE
    blower[equal] = bupper[equal]
    blower = np.clip(blower, -DAQP_INF, DAQP_INF)
    bupper = np.clip(bupper, -DAQP_INF, DAQP_INF)

    with open(os.path.join(output_dir, f"{name}.bin"), "wb") as output:
        output.write(struct.pack("<iii", n, m, n))
        H.astype("<f8").tofile(output)
        q.astype("<f8").tofile(output)
        A.astype("<f8").tofile(output)
        bupper.astype("<f8").tofile(output)
        blower.astype("<f8").tofile(output)
        sense.astype("<i4").tofile(output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    index = []
    for filename in sorted(os.listdir(args.data_dir)):
        if not filename.endswith(".mat"):
            continue
        name = filename[:-4]
        problem = load_mat(os.path.join(args.data_dir, filename))
        P, q, C, lower, upper, box_lower, box_upper = problem
        n = P.shape[0]
        converted_m = converted_constraint_count(C, lower, upper, box_lower)
        if n > 1000 or converted_m > 1000:
            continue

        export_problem(args.output_dir, name, *problem)
        # MarosMeszarosDensePosdef applies this check after `to_dense()`.
        index.append(
            (
                name,
                n,
                box_lower.size + C.shape[0],
                bool(is_posdef(np.asarray(P.todense()))),
            )
        )

    with open(os.path.join(args.output_dir, "index.txt"), "w", encoding="utf-8") as output:
        for name, n, m, posdef in index:
            output.write(f"{name} {n} {m} {int(posdef)}\n")

    posdef_count = sum(row[3] for row in index)
    print(f"Exported {len(index)} dense problems ({posdef_count} positive definite)")


if __name__ == "__main__":
    main()
