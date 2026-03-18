"""
Test the python interface to daqp.

Runs in Github CI
"""
import unittest
from ctypes import c_double, c_int

import daqp

import numpy as np


class Testing(unittest.TestCase):
    """Testing class for daqp python interface."""

    def test_python_demo(self):
        """Python demo code test."""
        H = np.array([[1, 0], [0, 1]], dtype=c_double)
        f = np.array([1, 1], dtype=c_double)
        A = np.array([[1, 1], [1, -1]], dtype=c_double)
        bupper = np.array([1, 2, 3, 4], dtype=c_double)
        blower = np.array([-1, -2, -3, -4], dtype=c_double)
        sense = np.array([0, 0, 0, 0], dtype=c_int)
        xstar, fval, exitflag, lam = daqp.solve(H, f, A, bupper, blower, sense)
        self.assertEqual(exitflag, 1)

    def test_warm_start_dual(self):
        """Dual warm start produces the same optimal solution as a cold start."""
        H = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=c_double)
        f = np.array([1.0, 1.0], dtype=c_double)
        A = np.array([[1.0, 1.0]], dtype=c_double)
        bupper = np.array([2.0, 2.0, 10.0], dtype=c_double)
        blower = np.array([-2.0, -2.0, -1.0], dtype=c_double)
        sense = np.array([0, 0, 0], dtype=c_int)

        # Cold start – reference solution
        x_cold, fval_cold, ef_cold, info_cold = daqp.solve(
            H, f, A, bupper, blower, sense)
        self.assertEqual(ef_cold, 1)

        sense_before = sense.copy()

        # Dual warm start using the multipliers from the cold solve
        lam_cold = info_cold['lam']
        x_warm, fval_warm, ef_warm, info_warm = daqp.solve(
            H, f, A, bupper, blower, sense, dual_start=lam_cold)

        self.assertEqual(ef_warm, 1)
        np.testing.assert_allclose(x_warm, x_cold, atol=1e-8)
        np.testing.assert_allclose(fval_warm, fval_cold, atol=1e-8)
        # Warm start should not require more iterations than cold start
        self.assertLessEqual(info_warm['iterations'], info_cold['iterations'])
        # sense array must not be mutated
        np.testing.assert_array_equal(sense, sense_before)

    def test_warm_start_primal(self):
        """Primal warm start produces the same optimal solution as a cold start."""
        H = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=c_double)
        f = np.array([1.0, 1.0], dtype=c_double)
        A = np.array([[1.0, 1.0]], dtype=c_double)
        bupper = np.array([2.0, 2.0, 10.0], dtype=c_double)
        blower = np.array([-2.0, -2.0, -1.0], dtype=c_double)
        sense = np.array([0, 0, 0], dtype=c_int)

        # Cold start – reference solution
        x_cold, fval_cold, ef_cold, info_cold = daqp.solve(
            H, f, A, bupper, blower, sense)
        self.assertEqual(ef_cold, 1)

        sense_before = sense.copy()

        # Primal warm start using the solution from the cold solve
        x_warm, fval_warm, ef_warm, info_warm = daqp.solve(
            H, f, A, bupper, blower, sense, primal_start=x_cold)

        self.assertEqual(ef_warm, 1)
        np.testing.assert_allclose(x_warm, x_cold, atol=1e-8)
        np.testing.assert_allclose(fval_warm, fval_cold, atol=1e-8)
        # Warm start should not require more iterations than cold start
        self.assertLessEqual(info_warm['iterations'], info_cold['iterations'])
        # sense array must not be mutated
        np.testing.assert_array_equal(sense, sense_before)

    def test_warm_start_does_not_modify_sense(self):
        """Neither primal nor dual warm start should modify the user's sense array."""
        H = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=c_double)
        f = np.array([1.0, 1.0], dtype=c_double)
        A = np.array([[1.0, 1.0]], dtype=c_double)
        bupper = np.array([2.0, 2.0, 10.0], dtype=c_double)
        blower = np.array([-2.0, -2.0, -1.0], dtype=c_double)
        sense = np.array([0, 0, 0], dtype=c_int)

        x_ref, _, _, info_ref = daqp.solve(H, f, A, bupper, blower, sense)
        lam_ref = info_ref['lam']
        sense_original = sense.copy()

        daqp.solve(H, f, A, bupper, blower, sense, primal_start=x_ref)
        np.testing.assert_array_equal(sense, sense_original,
                                      err_msg="primal_start must not mutate sense")

        daqp.solve(H, f, A, bupper, blower, sense, dual_start=lam_ref)
        np.testing.assert_array_equal(sense, sense_original,
                                      err_msg="dual_start must not mutate sense")


if __name__ == '__main__':
    unittest.main()
