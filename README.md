# Blocking Method for Correlated Time Series

This repository implements the **Blocking Method** for estimating statistical errors in **correlated time series**.

Many real-world datasets exhibit temporal correlations, meaning that successive observations are not statistically independent. In such cases, the naive variance estimator can significantly underestimate the true uncertainty of the mean.

The blocking method provides a systematic approach to estimate the variance of correlated data by progressively reducing autocorrelation through a sequence of data transformations.

---

## Applications

The blocking method can be applied to a wide range of time series data, including:

- Monte Carlo simulations
- Markov Chain Monte Carlo (MCMC)
- Physical simulations
- Financial time series
- Climate and environmental data
- Experimental measurements
- Any temporally correlated dataset

---

## Method

For correlated data, the naive variance estimator often underestimates the real statistical uncertainty.

The **blocking transformation** reduces correlations by iteratively averaging neighboring data points:

x'_i = (x_{2i} + x_{2i+1}) / 2


At each blocking level, statistical quantities such as variance and autocovariance are computed. A chi-square test is then used to determine the optimal blocking level at which the remaining data can be considered effectively uncorrelated.

The final result is a corrected estimate of the variance of the sample mean.

---

## Requirements
numpy
scipy
matplotlib

## Example

```python
import numpy as np
from blocking import blocking_method

data = np.loadtxt("timeseries.dat")

var_corrected = blocking_method(data)

print("Corrected variance:", var_corrected)

```
## Reference

Jonsson, M. (2018).  
Standard error estimation by an automated blocking method.  
Physical Review E, 98(4), 043304.  
https://doi.org/10.1103/PhysRevE.98.043304
