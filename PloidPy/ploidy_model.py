import numpy as np
import PloidPy.binom_model as bm

EPS = np.finfo(np.float64).tiny


# calculates the likelihood values for each heterozygous state for the x matrix
# for the n ploidy model. This model does NOT calculate the WEIGHTED
# likelihood, just the likelihood of each value for each model.
def ploidy_Likelihood(x, n, r, p_nb, p_err, uniform_com=True):
    het_p = np.arange(1, np.floor(n/2) + 1) / n
    return bm.get_Likelihood(x, het_p, r, p_nb, p_err, uniform_com=uniform_com)


def weighted_Ploidy_Log_Likelihood(lh, p_err, uniform_weights=False,
                                   uniform_com=False):
    w = bm.get_Weights(lh, p_err, uniform=uniform_weights,
                       uniform_com=uniform_com)
    a = np.multiply(lh, w[:, np.newaxis])
    a[a == 0] = EPS
    return np.sum(np.log(np.sum(a, axis=0))), w


# Calculates the  Akaike Information Criterion (AIC) value of x when given a
# list of ploidy models. Returns a tuple with the log likelihood values and the
# AIC values
def get_Log_Likelihood_AIC(x, models, r, p_nb, p_err, uniform_weights=False,
                           uniform_com=True):
    add_param = 2
    if p_err == 0:
        add_param += 1
    if uniform_com:
        add_param += 1
    w_lh = np.zeros(np.shape(models))
    w = []
    k = (np.floor(models / 2) * 2) + add_param
    for i in range(len(models)):
        w_lh[i], w0 = weighted_Ploidy_Log_Likelihood(
            ploidy_Likelihood(x, models[i], r, p_nb, p_err,
                              uniform_com=uniform_com), p_err,
            uniform_weights=uniform_weights, uniform_com=uniform_com)
        w += [w0]
    return (w_lh, (2 * k) - (2 * w_lh), w)


# selects best model based on the AIC value
def select_model(aic, model):
    return model[np.argmin(aic)]
