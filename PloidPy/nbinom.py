import numpy as np
import statsmodels.discrete.discrete_model as dm


# fits data onto a negative binomial distribution by approximately maximizing
# the parameters
def fit_nbinom(x):
    print("Fitting total read count to negative binomial distribution...")
    params = dm.NegativeBinomial(x, np.ones_like(x)).fit(maxiter=200000,
                                                         disp=0).params
    mu = np.exp(params[0])
    alpha = params[1]
    r = alpha ** -1
    p = r / (r + mu)
    return r, p
