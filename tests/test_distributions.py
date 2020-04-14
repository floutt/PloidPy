import unittest
import PloidPy.nbinom as nb
import PloidPy.ploidy_model as pm
import numpy as np

class test_prediction(unittest.TestCase):
    def test_pred(self):
        SIZE = 60000
        reads = np.random.poisson(np.random.randint(200, 300), size=SIZE)
        ploidy = np.random.randint(2, 9)
        mac = np.zeros(SIZE)
        lh = np.zeros((7, SIZE))
        for i in range(2, ploidy//2):
            unit = SIZE // (ploidy // 2)
            interval = np.arange((i - 1) * unit, i * unit, dtype=np.int64)
            mac[interval] = np.random.binomial(reads[interval], 1/i)
        x = np.array([mac, reads]).T
        r, p_nb = nb.fit_nbinom(x[:,1])
        p_err = 0
        llh, aic, w = pm.get_Log_Likelihood_AIC(x, np.array(range(2, 9)), r, p_nb, p_err)
        self.assertEqual(np.argmin(aic), ploidy - 2)
