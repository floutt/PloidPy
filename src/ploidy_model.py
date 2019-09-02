import numpy as np
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


# Calculates the  Akaike Information Criterion (AIC) value of x when given a
# list of ploidy models. Returns a tuple with the log likelihood values and the
# AIC values
def get_Log_Likelihood_AIC(x, models, lam):
    w_lh = np.zeros(np.shape(models))
    k = np.floor(models / 2) + 2
    for i in range(models):
        w_lh[i] = weighted_Ploidy_Log_Likelihood(
            ploidy_Likelihood(x, models[i], lam))
    return (w_lh, (2 * k) - (2 * w_lh))


# selects best model based on the AIC value
def select_model(aic, model):
    return model[np.argmin(aic)]
