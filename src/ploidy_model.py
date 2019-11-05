import numpy as np
import binom_model as bm


EPS = np.finfo(np.float64).tiny


# calculates the likelihood values for each heterozygous state for the x matrix
# for the n ploidy model. This model does NOT calculate the WEIGHTED
# likelihood, just the likelihood of each value for each model.
def ploidy_Likelihood(x, n, r, p_nb, p_err):
    het_p = np.arange(1, np.floor(n/2) + 1) / n
    return bm.get_Likelihood(x, het_p, r, p_nb, p_err)


def weighted_Ploidy_Log_Likelihood(lh):
    w = bm.get_Weights(lh)
    a = np.multiply(lh, w[:, np.newaxis])
    a[a == 0] = EPS
    return np.sum(np.log(np.sum(a, axis = 0))), w


# Calculates the  Akaike Information Criterion (AIC) value of x when given a
# list of ploidy models. Returns a tuple with the log likelihood values and the
# AIC values
def get_Log_Likelihood_AIC(x, models, r, p_nb, p_err):
    w_lh = np.zeros(np.shape(models))
    w = []
    k = (np.floor(models / 2) * 2) + 4
    for i in range(len(models)):
        w_lh[i], w0 = weighted_Ploidy_Log_Likelihood(
            ploidy_Likelihood(x, models[i], r, p_nb, p_err))
        w += [w0]
    return (w_lh, (2 * k) - (2 * w_lh), w)


# selects best model based on the AIC value
def select_model(aic, model):
    return model[np.argmin(aic)]
