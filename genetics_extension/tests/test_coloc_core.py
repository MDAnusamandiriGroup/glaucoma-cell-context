import unittest
import numpy as np
from scipy.stats import t
from coloc_core import derived_standard_error, log_abf, combine_log_abf


class ScientificContracts(unittest.TestCase):
    def test_standard_error_recovery_from_known_regression(self):
        beta = np.array([0.2, -0.4, 0.001, 0.6])
        se = np.array([0.06, 0.1, 0.08, 0.06])
        df = 404.0
        p = 2 * t.sf(np.abs(beta / se), df)
        np.testing.assert_allclose(derived_standard_error(beta, p, df), se, rtol=2e-8)

    def test_direct_enumeration_of_distinct_variant_hypotheses(self):
        a = np.array([2., 4., 1.])
        b = np.array([5., 1., 3.])
        p1, p2, p12 = 1e-4, 1e-4, 1e-5
        weights = np.array([1., p1*a.sum(), p2*b.sum(),
            p1*p2*sum(a[i]*b[j] for i in range(3) for j in range(3) if i != j),
            p12*np.dot(a, b)])
        actual, conditional = combine_log_abf(np.log(a), np.log(b), p1, p2, p12)
        np.testing.assert_allclose(actual, weights/weights.sum(), rtol=1e-13)
        np.testing.assert_allclose(conditional, a*b/np.dot(a,b), rtol=1e-13)

    def test_allele_flip_and_row_order_invariance(self):
        beta=np.array([0.2,-0.1,0.08]); se=np.array([0.05,0.04,0.06])
        l1=log_abf(beta,se,0.2); l2=log_abf(beta*2,se,0.15)
        np.testing.assert_array_equal(l1,log_abf(-beta,se,0.2))
        p,_=combine_log_abf(l1,l2)
        q,_=combine_log_abf(l1[::-1],l2[::-1])
        np.testing.assert_allclose(p,q,rtol=1e-13)

    def test_extreme_shared_and_distinct_signals_remain_finite(self):
        for a,b in [(np.array([1000.,0.,0.]),np.array([1000.,0.,0.])),
                    (np.array([1000.,0.,0.]),np.array([0.,1000.,0.]))]:
            p,c=combine_log_abf(a,b)
            self.assertTrue(np.isfinite(p).all())
            self.assertAlmostEqual(float(p.sum()),1.)
            self.assertAlmostEqual(float(c.sum()),1.)
        self.assertGreater(combine_log_abf([1000.,0.,0.],[1000.,0.,0.])[0][4],0.99)
        self.assertGreater(combine_log_abf([1000.,0.,0.],[0.,1000.,0.])[0][3],0.99)

    def test_invalid_input_is_rejected(self):
        for p in [0.,1.,np.nan]:
            with self.assertRaises(ValueError):derived_standard_error([0.2],[p],404)
        with self.assertRaises(ValueError):log_abf([0.2],[0.],0.2)
        with self.assertRaises(ValueError):combine_log_abf([0.],[0.])


if __name__ == '__main__':
    unittest.main()
