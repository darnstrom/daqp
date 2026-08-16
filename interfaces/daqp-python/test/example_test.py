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

    def test_symmetric_avi_uses_qp_fast_path(self):
        """An ill-conditioned symmetric AVI is solved as its equivalent QP."""
        Q = np.array([[1.0, 1.0], [-1.0, 1.0]], dtype=c_double) / np.sqrt(2.0)
        H = Q @ np.diag([1.0e-6, 1.0e6]) @ Q.T
        f = np.array([1.0, 2.0], dtype=c_double)
        A = np.eye(2, dtype=c_double)
        bupper = np.ones(2, dtype=c_double)
        blower = -np.ones(2, dtype=c_double)
        sense = np.zeros(2, dtype=c_int)

        x_qp, _, flag_qp, _ = daqp.solve(
            H, f, A, bupper, blower, sense, iter_limit=10)
        x_avi, _, flag_avi, _ = daqp.solve(
            H, f, A, bupper, blower, sense, is_avi=True, iter_limit=10)

        self.assertEqual(flag_qp, 1)
        self.assertEqual(flag_avi, 1)
        np.testing.assert_allclose(x_avi, x_qp, atol=1.0e-8)

    def test_avi_model_dispatches_on_symmetry(self):
        """Persistent AVI models select the QP or DR path during setup."""
        H_symmetric = np.array([[2.0, 0.5], [0.5, 1.0]], dtype=c_double)
        H_asymmetric = np.array([[2.0, 1.0], [0.0, 1.0]], dtype=c_double)
        f = np.array([1.0, 2.0], dtype=c_double)
        A = np.eye(2, dtype=c_double)
        bupper = np.ones(2, dtype=c_double)
        blower = -np.ones(2, dtype=c_double)
        sense = np.zeros(2, dtype=c_int)

        for H in (H_symmetric, H_asymmetric, H_symmetric):
            model = daqp.Model()
            setup_flag, _ = model.setup(
                H, f, A, bupper, blower, sense, is_avi=True)
            self.assertEqual(setup_flag, 1)
            x_model, _, flag_model, _ = model.solve()
            x_ref, _, flag_ref, _ = daqp.solve(
                H, f, A, bupper, blower, sense, is_avi=True)
            self.assertEqual(flag_model, 1)
            self.assertEqual(flag_ref, 1)
            np.testing.assert_allclose(x_model, x_ref, atol=1.0e-8)

    def test_ill_conditioned_asymmetric_avi(self):
        """An asymmetric AVI remains accurate with a hidden weak direction."""
        Q = np.array([[1.0, 1.0], [-1.0, 1.0]], dtype=c_double) / np.sqrt(2.0)
        H = Q @ np.diag([1.0e-6, 1.0]) @ Q.T
        H += np.array([[0.0, 1.0e-7], [-1.0e-7, 0.0]], dtype=c_double)
        f = np.array([0.75, -0.25], dtype=c_double)
        A = np.eye(2, dtype=c_double)
        bupper = np.ones(2, dtype=c_double)
        blower = -np.ones(2, dtype=c_double)
        sense = np.zeros(2, dtype=c_int)

        x, _, flag, _ = daqp.solve(
            H, f, A, bupper, blower, sense, is_avi=True)

        self.assertEqual(flag, 1)
        residual = H @ x + f
        for xi, ri in zip(x, residual):
            if xi <= -1.0 + 1.0e-7:
                self.assertGreaterEqual(ri, -1.0e-7)
            elif xi >= 1.0 - 1.0e-7:
                self.assertLessEqual(ri, 1.0e-7)
            else:
                self.assertAlmostEqual(ri, 0.0, delta=1.0e-7)
                
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


class TestModel(unittest.TestCase):
    """Tests for the daqp.Model workspace class."""

    def _make_qp(self):
        """Return a simple QP for testing."""
        H = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=c_double)
        f = np.array([2.0, 2.0], dtype=c_double)
        A = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=c_double)
        bupper = np.array([1.0, 1.0], dtype=c_double)
        blower = np.array([-1.0, -1.0], dtype=c_double)
        sense = np.array([0, 0], dtype=c_int)
        return H, f, A, bupper, blower, sense

    def test_model_setup_and_solve(self):
        """Model setup followed by solve returns an optimal solution."""
        H, f, A, bupper, blower, sense = self._make_qp()
        d = daqp.Model()
        exitflag, _ = d.setup(H, f, A, bupper, blower, sense)
        self.assertGreaterEqual(exitflag, 0)
        x, fval, ef, info = d.solve()
        self.assertEqual(ef, 1)
        np.testing.assert_allclose(x, [-1.0, -1.0], atol=1e-6)

    def test_model_solve_before_setup_raises(self):
        """Calling solve before setup raises RuntimeError."""
        d = daqp.Model()
        with self.assertRaises(RuntimeError):
            d.solve()

    def test_model_update_before_setup_raises(self):
        """Calling update before setup raises RuntimeError."""
        d = daqp.Model()
        with self.assertRaises(RuntimeError):
            d.update(f=np.array([1.0, 1.0], dtype=c_double))

    def test_model_update_cost_vector(self):
        """Updating f via update() changes the optimal solution."""
        H, f, A, bupper, blower, sense = self._make_qp()
        d = daqp.Model()
        d.setup(H, f, A, bupper, blower, sense)

        x1, _, ef1, _ = d.solve()
        self.assertEqual(ef1, 1)
        np.testing.assert_allclose(x1, [-1.0, -1.0], atol=1e-6)

        # Flip the sign of f — optimal should flip to [1, 1]
        f_new = np.array([-2.0, -2.0], dtype=c_double)
        exitflag = d.update(f=f_new)
        self.assertEqual(exitflag, 0)

        x2, _, ef2, _ = d.solve()
        self.assertEqual(ef2, 1)
        np.testing.assert_allclose(x2, [1.0, 1.0], atol=1e-6)

    def test_model_update_bounds(self):
        """Updating bounds via update() changes the feasible set."""
        H, f, A, bupper, blower, sense = self._make_qp()
        d = daqp.Model()
        d.setup(H, f, A, bupper, blower, sense)

        bu_new = np.array([0.5, 0.5], dtype=c_double)
        bl_new = np.array([-0.5, -0.5], dtype=c_double)
        d.update(bupper=bu_new, blower=bl_new)
        x, _, ef, _ = d.solve()
        self.assertEqual(ef, 1)
        np.testing.assert_allclose(x, [-0.5, -0.5], atol=1e-6)

    def test_model_keeps_equalities_by_default(self):
        """Persistent models do not eliminate equalities by default."""
        n = 10
        neq = 6
        H = np.eye(n, dtype=c_double)
        f = np.zeros(n, dtype=c_double)
        A = np.vstack((np.eye(n, dtype=c_double)[:neq], np.ones((1, n))))
        bupper = np.concatenate((
            np.full(n, 10.0), np.zeros(neq), np.array([10.0])))
        blower = np.concatenate((
            np.full(n, -10.0), np.zeros(neq), np.array([-10.0])))
        sense = np.concatenate((
            np.zeros(n, dtype=c_int), np.full(neq, 5, dtype=c_int),
            np.zeros(1, dtype=c_int)))

        d = daqp.Model()
        setup_flag, _ = d.setup(H, f, A, bupper, blower, sense)
        self.assertGreater(setup_flag, 0)

        x, _, exitflag, _ = d.solve()
        self.assertEqual(exitflag, 1)
        np.testing.assert_allclose(x, np.zeros(n), atol=1e-8)

        bupper[-1] = 9.0
        blower[-1] = -9.0
        update_flag = d.update(bupper=bupper, blower=blower)
        self.assertEqual(update_flag, 0)

    def test_model_solve_returns_copies(self):
        """Each solve() call returns independent copies of x and lam."""
        H, f, A, bupper, blower, sense = self._make_qp()
        d = daqp.Model()
        d.setup(H, f, A, bupper, blower, sense)

        x1, _, _, info1 = d.solve()
        x1[0] = 999.0  # mutate the returned array
        x2, _, _, info2 = d.solve()
        self.assertNotEqual(x2[0], 999.0)

    def test_model_settings_get(self):
        """settings getter returns a dict with all expected keys."""
        d = daqp.Model()
        s = d.settings
        self.assertIn('iter_limit', s)
        self.assertIn('primal_tol', s)
        self.assertIn('eps_prox', s)
        self.assertIn('time_limit', s)
        self.assertEqual(s['time_limit'], 0.0)

    def test_model_settings_set(self):
        """settings setter updates only the specified keys."""
        d = daqp.Model()
        original_primal_tol = d.settings['primal_tol']
        d.settings = {'iter_limit': 42, 'time_limit': 5.0}
        self.assertEqual(d.settings['iter_limit'], 42)
        self.assertAlmostEqual(d.settings['time_limit'], 5.0)
        # Other settings should be unchanged
        self.assertAlmostEqual(d.settings['primal_tol'], original_primal_tol)

    def test_automatic_prox_tolerance_setting(self):
        """Negative eta remains automatic as dual tolerance changes."""
        d = daqp.Model()
        self.assertAlmostEqual(d.settings['eta_prox'], -1.0)
        d.settings = {'dual_tol': 1e-9}
        self.assertAlmostEqual(d.settings['eta_prox'], -1.0)

        d.settings = {'dual_tol': 1e-10, 'eta_prox': 5e-8}
        self.assertAlmostEqual(d.settings['eta_prox'], 5e-8)

    def test_model_settings_preserved_across_setup(self):
        """Settings modified before setup are preserved after setup."""
        H, f, A, bupper, blower, sense = self._make_qp()
        d = daqp.Model()
        d.settings = {'iter_limit': 123}
        d.setup(H, f, A, bupper, blower, sense)
        self.assertEqual(d.settings['iter_limit'], 123)

    def test_model_settings_preserved_across_reuse(self):
        """Settings are preserved when setup is called a second time."""
        H, f, A, bupper, blower, sense = self._make_qp()
        d = daqp.Model()
        d.setup(H, f, A, bupper, blower, sense)
        d.settings = {'iter_limit': 77}
        # Re-setup with slightly different cost
        f2 = np.array([-2.0, -2.0], dtype=c_double)
        d.setup(H, f2, A, bupper, blower, sense)
        self.assertEqual(d.settings['iter_limit'], 77)
        x, _, ef, _ = d.solve()
        self.assertEqual(ef, 1)

    def test_model_matches_quadprog(self):
        """Model produces the same result as daqp.solve for the same problem."""
        H, f, A, bupper, blower, sense = self._make_qp()

        x_ref, fval_ref, ef_ref, _ = daqp.solve(H, f, A, bupper, blower, sense)
        self.assertEqual(ef_ref, 1)

        d = daqp.Model()
        d.setup(H, f, A, bupper, blower, sense)
        x_ws, fval_ws, ef_ws, _ = d.solve()

        self.assertEqual(ef_ws, 1)
        np.testing.assert_allclose(x_ws, x_ref, atol=1e-8)
        np.testing.assert_allclose(fval_ws, fval_ref, atol=1e-8)

    def test_model_warm_start_dual(self):
        """Dual warm start via setup does not require more iterations than cold start."""
        H, f, A, bupper, blower, sense = self._make_qp()
        d = daqp.Model()
        d.setup(H, f, A, bupper, blower, sense)
        x_cold, _, ef_cold, info_cold = d.solve()
        self.assertEqual(ef_cold, 1)

        lam_cold = info_cold['lam']
        d2 = daqp.Model()
        d2.setup(H, f, A, bupper, blower, sense, dual_start=lam_cold)
        x_warm, _, ef_warm, info_warm = d2.solve()
        self.assertEqual(ef_warm, 1)
        np.testing.assert_allclose(x_warm, x_cold, atol=1e-8)
        self.assertLessEqual(info_warm['iterations'], info_cold['iterations'])

    def test_model_setup_does_not_mutate_sense(self):
        """setup() with dual/primal_start does not modify the caller's sense array."""
        H, f, A, bupper, blower, sense = self._make_qp()
        sense_orig = sense.copy()

        d = daqp.Model()
        d.setup(H, f, A, bupper, blower, sense)
        x_cold, _, _, info_cold = d.solve()

        # dual warm start
        d2 = daqp.Model()
        d2.setup(H, f, A, bupper, blower, sense,
                 dual_start=info_cold['lam'])
        np.testing.assert_array_equal(sense, sense_orig,
                                      err_msg="setup dual_start must not mutate sense")

        # primal warm start
        d3 = daqp.Model()
        d3.setup(H, f, A, bupper, blower, sense, primal_start=x_cold)
        np.testing.assert_array_equal(sense, sense_orig,
                                      err_msg="setup primal_start must not mutate sense")


class TestSemiProximal(unittest.TestCase):
    """Tests for automatic and forced proximal modes."""

    def test_pd_hessian_n_prox_zero(self):
        """Automatic mode leaves a positive-definite Hessian unshifted."""
        H = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=c_double)
        f = np.array([1.0, 1.0], dtype=c_double)
        A = np.zeros((0, 2), dtype=c_double)
        bupper = np.array([1.0, 1.0], dtype=c_double)
        blower = np.array([-1.0, -1.0], dtype=c_double)
        sense = np.array([0, 0], dtype=c_int)

        d = daqp.Model()
        d.settings = {'eps_prox': -1e-4}
        d.setup(H, f, A, bupper, blower, sense)
        x, fval, ef, info = d.solve()

        self.assertEqual(ef, 1)
        np.testing.assert_allclose(x, [-1.0, -1.0], atol=1e-4)

    def test_singular_hessian_semi_proximal(self):
        """Rank-1 Hessian: only the singular direction gets regularised."""
        # H = diag(1, 0): x2 direction is singular
        H = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=c_double)
        f = np.array([1.0, 1.0], dtype=c_double)
        # Use simple-bound columns as A rows so A is non-empty (ms=0, mA=2)
        A = np.eye(2, dtype=c_double)
        bupper = np.array([2.0, 2.0], dtype=c_double)
        blower = np.array([-2.0, -2.0], dtype=c_double)
        sense = np.array([0, 0], dtype=c_int)

        d = daqp.Model()
        d.settings = {'eps_prox': -1e-3}
        d.setup(H, f, A, bupper, blower, sense)
        x, fval, ef, info = d.solve()

        self.assertEqual(ef, 1)
        # x1 minimises 0.5*x1^2 + x1 unconstrained -> x1* = -1
        # x2 is singular; regularisation + f[1]=1 pushes to lower bound
        self.assertAlmostEqual(x[0], -1.0, places=3)
        self.assertAlmostEqual(x[1], -2.0, places=3)

    def test_zero_hessian_all_singular(self):
        """Zero Hessian: all directions are singular, full proximal applied."""
        H = np.array([[0.0, 0.0], [0.0, 0.0]], dtype=c_double)
        f = np.array([1.0, 2.0], dtype=c_double)
        # Use simple-bound columns as A rows so A is non-empty
        A = np.eye(2, dtype=c_double)
        bupper = np.array([3.0, 3.0], dtype=c_double)
        blower = np.array([-3.0, -3.0], dtype=c_double)
        sense = np.array([0, 0], dtype=c_int)

        d = daqp.Model()
        d.settings = {'eps_prox': -1e-2}
        d.setup(H, f, A, bupper, blower, sense)
        x, fval, ef, info = d.solve()

        self.assertGreater(ef, 0)
        # Regularised problem: min eps/2*||x||^2 + f'*x  -> x* = -f/eps (clipped to bounds)
        np.testing.assert_allclose(x, [-3.0, -3.0], atol=1e-3)

    def test_zero_eps_disables_proximal_regularization(self):
        """A singular Hessian is rejected when eps_prox is exactly zero."""
        H = np.diag([1.0, 0.0]).astype(c_double)
        f = np.ones(2, dtype=c_double)
        A = np.zeros((0, 2), dtype=c_double)
        bupper = np.ones(2, dtype=c_double)
        blower = -bupper
        sense = np.zeros(2, dtype=c_int)

        _, _, exitflag, _ = daqp.solve(
            H, f, A, bupper, blower, sense, eps_prox=0.0)

        self.assertEqual(exitflag, -5)

    def test_pd_hessian_matches_no_prox(self):
        """With a PD Hessian, eps_prox > 0 still gives the correct solution."""
        H = np.array([[4.0, 1.0], [1.0, 3.0]], dtype=c_double)
        f = np.array([1.0, 2.0], dtype=c_double)
        A = np.zeros((0, 2), dtype=c_double)
        bupper = np.array([5.0, 5.0], dtype=c_double)
        blower = np.array([-5.0, -5.0], dtype=c_double)
        sense = np.array([0, 0], dtype=c_int)

        # Reference: solve without proximal
        x_ref, _, ef_ref, _ = daqp.solve(H, f, A, bupper, blower, sense)
        self.assertEqual(ef_ref, 1)

        # Solve with proximal on PD H -- should converge to same solution
        d = daqp.Model()
        d.settings = {'eps_prox': 1e-4}
        d.setup(H, f, A, bupper, blower, sense)
        x_prox, _, ef_prox, _ = d.solve()
        self.assertEqual(ef_prox, 1)
        np.testing.assert_allclose(x_prox, x_ref, atol=1e-3)

    def test_positive_eps_forces_full_proximal(self):
        """Positive eps_prox forces the full shift and preserves the QP solution."""
        H = np.array([[4.0, 1.0], [1.0, 3.0]], dtype=c_double)
        f = np.array([1.0, 2.0], dtype=c_double)
        A = np.zeros((0, 2), dtype=c_double)
        bupper = np.array([5.0, 5.0], dtype=c_double)
        blower = np.array([-5.0, -5.0], dtype=c_double)
        sense = np.array([0, 0], dtype=c_int)

        x_ref, _, flag_ref, _ = daqp.solve(
            H, f, A, bupper, blower, sense)
        x_force, _, flag_force, info_force = daqp.solve(
            H, f, A, bupper, blower, sense, eps_prox=1e-2)

        self.assertEqual(flag_ref, 1)
        self.assertEqual(flag_force, 1)
        self.assertGreater(info_force['iterations'], 1)
        np.testing.assert_allclose(x_force, x_ref, atol=1e-6)

    def test_consistent_redundant_equalities_are_ignored(self):
        """Dependent consistent equalities do not make setup fail."""
        H = np.eye(2, dtype=c_double)
        f = np.zeros(2, dtype=c_double)
        A = np.array([[1.0, 0.0], [2.0, 0.0]], dtype=c_double)
        bounds = np.array([1.0, 2.0], dtype=c_double)
        sense = np.array([5, 5], dtype=c_int)

        x, _, exitflag, _ = daqp.solve(
            H, f, A, bounds, bounds, sense
        )

        self.assertEqual(exitflag, 1)
        self.assertAlmostEqual(x[0], 1.0)

    def test_inconsistent_redundant_equalities_still_fail(self):
        """A dependent equality with an inconsistent RHS is not discarded."""
        H = np.eye(2, dtype=c_double)
        f = np.zeros(2, dtype=c_double)
        A = np.array([[1.0, 0.0], [2.0, 0.0]], dtype=c_double)
        bounds = np.array([1.0, 3.0], dtype=c_double)
        sense = np.array([5, 5], dtype=c_int)

        _, _, exitflag, _ = daqp.solve(
            H, f, A, bounds, bounds, sense
        )

        self.assertLess(exitflag, 0)


if __name__ == '__main__':
    unittest.main()
