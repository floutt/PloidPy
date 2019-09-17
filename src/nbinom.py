import numpy as np
import scipy.stats as stats
import statsmodels.discrete.discrete_model as dm


# fits data onto a negative binomial distribution by approximately maximizing
# the parameters
def fit_nbinom(x):
    params = dm.NegativeBinomial(x, np.ones_like(x)).fit().params
    mu = np.exp(params[0])
    alpha = params[1]
    r = alpha ** -1
    p = r / (r + mu)
    return r, p
