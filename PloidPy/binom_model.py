import numpy as np
from scipy.stats import binom, nbinom, randint

EPS = np.finfo(np.float64).tiny


# since we are using calculations based off of the minor allele frequency, we
# have to use a truncated binomial instead of the traditional one. In order to
# adjust for this we normalize the data based off of this. The maximum possible
# value in our case will be (0.5 * x) and the minimum possible value is 1
def truncated_binom_pmf(x, n, p):
    return binom.pmf(x, n, p) / (binom.cdf(n/2, n, p) - binom.pmf(0, n, p))


# calculates the likelihood of each value in the joint distribution x given an
# underlying negative binomial distribution distribution for reads and a
# conditional truncated binomial distribution for the minor allele
def compound_nb_binom_pmf(x, p, r, p_nb):
    lh_nb = nbinom.pmf(x[:, 1], r, p_nb)
    lh_b = truncated_binom_pmf(x[:, 0], x[:, 1], np.ones_like(x[:, 1]) * p)
    return lh_nb * lh_b


# probability mass function for the uniform noise distribution of the data.
# truncated at 0.5 * x
def uniform_pmf(x, r, p_nb):
    lh_nb = nbinom.pmf(x[:, 1], r, p_nb)
    # truncated uniform component
    lh_unfm = randint.pmf(x[:, 0], 1, np.floor(0.5 * x[:, 1]))
    return lh_nb * lh_unfm


# calculates a matrix of the binom_mix for a vector of p values
def get_Likelihood(x, p, r, p_nb, p_err, uniform_com=True):
    extra_param = 0
    if uniform_com:
        extra_param += 1
    if p_err > 0:
        extra_param += 1
    # likelihood of p
    lh = np.ones((len(p) + extra_param, len(x)))
    for i in range(len(p)):
        lh[i] = compound_nb_binom_pmf(x, p[i], r, p_nb)
    if extra_param == 1:
        if p_err > 0:
            lh[-1] = compound_nb_binom_pmf(x, p_err, r, p_nb)
        else:
            lh[-1] = uniform_pmf(x, r, p_nb)
        lh[-1][np.isnan(lh[-1])] = EPS
    elif extra_param == 2:
        lh[-1] = uniform_pmf(x, r, p_nb)
        lh[-1][np.isnan(lh[-1])] = EPS
        lh[-2] = compound_nb_binom_pmf(x, p_err, r, p_nb)
        lh[-2][np.isnan(lh[-2])] = EPS
    return lh


# uses expectation maximization to get the weights of each subpopulation
# model given a set of fixed distributions. Calculates the weights from
# likelihood data
def get_Weights(lh, p_err, uniform=False, uniform_com=False):
    size = len(lh)
    w = np.ones(size)/size
    if uniform and p_err == 0 and not uniform_com:
        return w

    def calc_p():
        p0 = np.multiply(lh, w[:, np.newaxis])
        p0[p0 == 0] = EPS
        return p0

    old_p = calc_p()
    w = np.nanmean(old_p / np.sum(old_p, axis=0), axis=1)
    p = calc_p()
    if not np.sum(np.log(p)) < np.sum(np.log(old_p)):
        # stop when the function reaches a maximum
        while np.sum(np.log(p)) > np.sum(np.log(old_p)):
            w = np.nanmean(p / np.sum(p, axis=0), axis=1)
            old_p = p
            p = calc_p()

    if not uniform:
        return w
    elif p_err and uniform_com:
        w[:-2] = (np.sum(w) - np.sum(w[-2:])) / (len(w) - 2)
        return w
    else:
        w[:-1] = (np.sum(w) - np.sum(w[-1:])) / (len(w) - 1)
        return w
