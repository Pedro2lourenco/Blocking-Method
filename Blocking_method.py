import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import chi2


def blocking_method(X):
    """
    Blocking Method for error estimation in correlated time series.

    This method reduces autocorrelation effects by recursively averaging
    neighboring data points (blocking transformation) and applies a
    chi-square test to determine the optimal blocking level.

    Parameters
    ----------
    X : array-like
        Time series data. Ideally the length should be a power of two (2^p).

    Returns
    -------
    var_new : float
        Estimated variance corrected for autocorrelation.
    """

    # Convert input to numpy array
    X = np.asarray(X)

    # Number of measurements
    n = len(X)

    # Mean value of the series
    mu = np.mean(X)

    # Degrees of freedom used for the chi-square test
    ddof = np.arange(1, 48)

    # Critical chi-square values for 95% confidence level
    q = chi2.ppf(0.95, ddof)

    # Maximum number of blocking transformations possible
    d = int(np.floor(np.log2(n)))

    # Arrays to store variance and autocovariance estimates
    s = np.zeros(d)
    gamma = np.zeros(d)

    i = 0

    # Perform blocking transformations
    while n > 1:

        # Centered series
        x = X - mu

        # Estimate autocovariance
        gamma[i] = np.sum(x[1:n] * x[0:n-1]) / n

        # Estimate variance
        s[i] = np.sum(x * x) / n

        # Create new blocked series
        y = np.zeros(n // 2)

        # Average neighboring data points
        for j in range(n // 2):
            y[j] = 0.5 * (X[2*j] + X[2*j+1])

        # Replace series with blocked version
        X = y

        # Reduce number of points
        n //= 2

        i += 1

    # Compute the test statistic M_k
    M = np.zeros(d)

    for k in range(d):

        n = 2**(d - k)

        M[k] = n * ((((n - 1) * s[k] / n**2) + gamma[k])**2) / (s[k]**2)

    # Cumulative sum (reverse order)
    M = np.cumsum(M[::-1])[::-1]

    # Determine optimal blocking level
    for p in range(d):

        if M[p] < q[p]:

            optimal = p
            break

    # Final variance estimate corrected for autocorrelation
    var_new = s[optimal] / (2**(d - optimal))

    return var_new

'''
############ Aplication #################

# New Variance without autocorrelation
VarNew = blocking_method(X) 

# Autocorrelation time
tau = 0.5 * (((VarE * len(X)) / np.var(X)) - 1)

'''