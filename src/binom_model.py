import numpy as np
from scipy.stats import binom

# since we are using calculations based off of the minor allele frequency, we
# have to use a truncated binomial instead of the traditional one. In order to
# adjust for this we normalize the data based off of this. The maximum possible
# value in our case will be (0.5 * x) and the minimum possible value is 1
def truncated_binom_pmf(x, n, p):
    if n < 1:
        return 0
    else:
        return binom.pmf(x, n, p) / (binom.cdf(n/2, n, p) - binom.pmf(0, n, p))


# calculates the likelihood of each value in x based off of
def binom_mix_pmf(x, p, dens, n_val):
    lh = np.zeros_like(x)
    for i in range(len(n_val)):
        lh = lh + (dens[i] * truncated_binom_pmf(x, n_val[i], p))
    return lh


# calculates a matrix of the binom_mix for a vector of p values
def get_Likelihood(x, p, dens, n_val):
    # likelihood of p
    lh = np.ones((len(p), len(x)))
    for i in range(len(p)):
        lh[i] = binom_mix_pmf(x, p[i], dens, n_val)
    return lh


# uses expectation maximization to get the weights of each subpopulation
# model given a set of fixed distributions. Calculates the weights from
# likelihood data
def get_Weights(lh):
    return np.sum(lh, axis=1) / np.sum(lh)
