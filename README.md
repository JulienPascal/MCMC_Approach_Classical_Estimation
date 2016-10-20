# MCMC_Approach_Classical_Estimation
This repository illustrates how to use LTEs (Laplace type estimators) to identify the coeffients of a median regression using the LAD estimator (Powell, 1984). LTEs were developped by V. Chernozhukov and H. Hong in "An MCMC approach to classical estimation", Journal of Econometrics, Vol. 115, N°2, pp. 293-346, 2003.

This code is still in its preliminary stage. One main problem is that the rejection rate is too high. A multivariate transition kernel is used, and the selection of the "jumping" factor is done manually so far. An adaptive MCMC scheme will be implemented in the future.

## The Problem at hand
The data is generated according to: y_star = Xθ + u. Yet, only y = max(0,y) is observable to the econometrician. 

## Example
This scripts generate y_star and y for a sample of 1000 individuals such that the truncation rate is appoximately 40%.

![GitHub Logo](https://github.com/JulienPascal/MCMC_Approach_Classical_Estimation/blob/master/output/censoring.png)


## Results
The parameter b can estimated using the MCMC algorithm described by V. Chernozhukov and H. Hong (2003). In this example, the true value for θ is (1 3 3 3). The median of the quasi-posterior distribution can be used as an estimator of θ. Quantiles of the quasi-posterior distribution can also be used to build confidence intervals. 

![GitHub Logo](https://github.com/JulienPascal/MCMC_Approach_Classical_Estimation/blob/master/output/draws_censored_regression.png)

![GitHub Logo](https://github.com/JulienPascal/MCMC_Approach_Classical_Estimation/blob/master/output/histograms_censored_regression.png)

|           |                     |                     |                     |                    | 
|-----------|---------------------|---------------------|---------------------|--------------------| 
| variable  |  θ0                 |  θ1                 |  θ2                 |  θ3                | 
| P10       | 0.9729222463036142  | 2.823425986798605   | 2.9694823610212064  | 2.9817850591413575 | 
| Median    | 1.0000193577408492  | 2.886062386934549   | 2.9972225679108764  | 3.0101184734177098 | 
| P90       | 1.0291370017480144  | 2.9474441003672367  | 3.02585751230858    | 3.041841506192087  | 
Note: The true value for  θ = (θ1  θ2  θ3  θ4) is (1 3 3 3). This table reports the 10th, the median and te 90th percentiles of the quasi-posterior distribution.

## Presentation 
More on MCMC and LTEs:
[Read my presentation at Sciences Po](https://github.com/JulienPascal/MCMC_Approach_Classical_Estimation/blob/master/main.pdf)

## Requirements
Julia Version 0.4.5


