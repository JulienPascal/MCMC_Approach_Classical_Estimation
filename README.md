# An MCMC Approach to Classical Estimation

This repository illustrates how to use **Laplace type estimators** (LTEs) in the
context of a **median regression**. See the [slides.](https://github.com/JulienPascal/MCMC_Approach_Classical_Estimation/blob/julia.1.6/main.pdf) for more details.

---

## The Problem at hand
The data is generated according to:

<img src="https://render.githubusercontent.com/render/math?math=y*=X\beta%2Bu">

But we only observe:

<img src="https://render.githubusercontent.com/render/math?math=y=\max(0,y*)">

Can we estimate <img src="https://render.githubusercontent.com/render/math?math=\beta">?

---

## Example
See the script `main.jl`. Approximately 40% of observations are truncated.

![alt text](https://github.com/JulienPascal/MCMC_Approach_Classical_Estimation/blob/julia.1.6/output/censoring.png)

---

## Results

![alt text](https://github.com/JulienPascal/MCMC_Approach_Classical_Estimation/blob/julia.1.6/output/chains.png)

![alt text](https://github.com/JulienPascal/MCMC_Approach_Classical_Estimation/blob/julia.1.6/output/histograms.png)

| variable   | beta0    | beta1    | beta2    | beta3    |
|------------|----------|----------|----------|----------|
| P10        | 0.972387 | 1.920641 | 2.966820 | 3.975281 |
| Median     | 1.002964 | 1.974698 | 2.994887 | 4.003844 |
| P90        | 1.034758 | 2.031147 | 3.024863 | 4.034115 |
| True value | 1        | 2        | 3        | 4        |


Note: This table reports the 10th, the median and te 90th percentiles of the quasi-posterior distribution.

---

## References
* J. L. Powell (1984), "Least Absolute Deviations estimation for the censored regression model", Journal of Econometrics 25, 303-325
* V. Chernozhukov and H. Hong in "An MCMC approach to classical estimation", Journal of Econometrics, Vol. 115, N°2, pp. 293-346, 2003.
