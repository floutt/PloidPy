import numpy as np
import scipy.stats as stats
import statsmodels.discrete.discrete_model as dm


# fits data onto a negative binomial distribution by approximately maximizing
# the parameters
def fit_nbinom(x, cdf_cutoff = 0.95):
    # filter out the outliers by filtering out data after a certain point
    uniq = np.unique(x, return_counts = True)
    cdf = np.cumsum(uniq[1] / np.sum(uniq[1]))
    max_val = uniq[0][np.where(cdf > cdf_cutoff)[0][0]]
    x0 = x[x <= max_val]
    params = dm.NegativeBinomial(x0, np.ones_like(x0)).fit(maxiter=200000).params
    mu = np.exp(params[0])
    alpha = params[1]
    r = alpha ** -1
    p = r / (r + mu)
    return r, p
