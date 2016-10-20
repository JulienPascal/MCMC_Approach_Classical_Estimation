# MCMC_Approach_Classical_Estimation
This repository illustrates how to use LTEs (Laplace type estimators) to identify the coeffients of a median regression using the LAD estimator (Powell, 1984). LTEs were developped by V. Chernozhukov and H. Hong in "An MCMC approach to classical estimation",Journal of Econometrics, Vol. 115, N°2, pp. 293-346, 2003

This code is still in its preliminary stage. One main problem is that the rejection rate is too high. A multivariate transition kernel is used, and the selection of the "jumping" factor is done manually so far. An adaptive MCMC scheme will be implemented in the future.

## The Problem at hand
The data is generated according to: y = Xb + u 
Yet, only z = max(0,y) is observable.

## Example
This scripts generate y and z for a sample of 1000 individuals such that the truncation rate is appoximately 40%.

[insert graph here]


## Results
b can estimated using the MCMC algorithm described by V. Chernozhukov and H. Hong (2003)

## Presentation
