import numpy
import thinned_poisson_mixture_model as tp


# calculates the likelihood values for each heterozygous state for the x matrix
# for the n ploidy model. This model does NOT calculate the WEIGHTED
# likelihood, just the likelihood of each value for each model.
def ploidy_Likelihood(x, n, lam):
    het_p = np.arange(1, np.floor(n/2) + 1) / n
    return tp.get_Likelihood(x, het_p, lam)


def weighted_Ploidy_Log_Likelihood(lh):
    w = tp.get_Weights(lh)
    return np.sum(np.log(np.sum(np.multiply(lh, w[:, np.newaxis]), axis = 0)))
