import numpy as np
import scipy.stats as stats

# pmf function for the thinned poisson distribution. p is the probability of
# success for each random Poisson variable, with each success being distributed
# according to Bern(p). Lambda is the mean of the underlining poisson
# distribution. This essentially results in a Pois(p * lambda) distribution
def tPois_pmf(x, p, lam):
    return stats.poisson.pmf(x, p*lam)


# calculates the pmf of a certain variable given MIXTURE of thinned
# poisson distribribution. p and lam are now matrices with weights
# corresponding to the value in w
def mixed_tPois_pmf(x, p, lam, w):
    lh = None
    if isinstance(x, int):
        lh = np.zeros(len(p))
    else:
        lh = np.zeros((len(p), len(x)))

    for i in range(len(p)):
        lh[i] = tPois_pmf(x, p[i], lam)
    return np.dot(lh.T, w.T)


def get_Likelihood(x, p, lam):
    # we assume a fixed mean across each data. This is very relevant to our
    # ploidy model, which tries to reduce bias coming from different coverages
    # for each heterozygous model. lambda here is the mean of the underlying
    # poisson distribution BEFORE doing the Bernoulli processes on it

    lh = np.ones((len(p), len(x)))
    for i in range(len(p)):
        lh[i] = tPois_pmf(x, p[i], lam)
    return lh


# uses expectation maximization to get the weights of each thinned Poisson
# model given a set of fixed distributions. Calculates the weights from
# likelihood data
def get_Weights(lh):
    return np.sum(lh, axis=1) / np.sum(lh)
