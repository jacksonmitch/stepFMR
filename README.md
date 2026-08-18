
# stepFMR

<!-- badges: start -->

[![R-CMD-check](https://github.com/jacksonmitch/stepFMR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jacksonmitch/stepFMR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Finite Mixture Regression models allows relationships between the
response ($\mathbf{Y}$) and the predictors
$(\mathbf{X}_1,\dots, \mathbf{X}_p)$ to differ between a finite number
of subpopulations. An R package, stepFMR implements stepwise procedures
that select predictors having an effect on the response and determines
whether each of the selected predictors has homogeneous or heterogeneous
effects. The stepwise procedures are based on hypothesis testing where a
weighted significance test (WEST) is used. In mixture regression, the
unknown number of the subpopulations $G$ (referred as order of the
mixture), poses challenges in testing the association between a
predictor and the response. With WEST, the significance of a predictor’s
effects on $y$ is the weighted sum of the significance level of all
candidate $G$’s. The package supports the Gaussian, Poisson, and
Binomial distribution families, and returns the effect types of each
predictor, estimate of their effects, as well as the best-fitting $G$.

## Installation

You can install the development version of stepFMR from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("jacksonmitch/stepFMR")
```

## Example

These are basic examples of the function `stepFMR` being used with
different distributions, showcasing the different s3 methods to be used
with the result:

``` r
library(stepFMR)

sim <- fmr_gaussian_example

fit <- stepFMR(sim$formula, sim$data, G_max = 3)
print(fit)
#> stepFMR | family: gaussian | G: 2,3 | procedure: complete
#> 
#> Variable selection
#>   Direction: forward  |  Alpha: 0.05
#>   Selected:      x1, x2
#>   Not selected:  x3
#>   Final formula: y ~ x1 + x2
#> 
#> Effect-type determination
#>   Direction: forward  |  Alpha: 0.05
#>   Heterogeneous: x1, x2
#>   Homogeneous:   none
#>   Final formula: y ~ (x1 + x2 | group)
#> 
#> Best fit (G: 2)
#> 
#> Mixing proportions:
#>     g1     g2 
#> 0.3739 0.6261 
#> 
#> Component standard deviations:
#>     g1     g2 
#> 0.5452 0.5134 
#> 
#> Heterogeneous coefficients:
#>    (Intercept)      x1     x2
#> g1     -2.0123 -2.9878 0.9595
#> g2      1.9927  3.0204 1.0594
# Align fitted components to the known true groups
aligned <- match_groups(fit, sim$true_group)
aligned$best_fit$parameter_values[["beta_g"]]
#>    (Intercept)       x1        x2
#> g1   -2.012291 -2.98783 0.9595001
#> g2    1.992688  3.02040 1.0593727
```

``` r
# Poisson
sim_pois <- fmr_poisson_example
fit_pois <- stepFMR(sim_pois$formula, sim_pois$data, G_max = 3, family = "poisson")
summary(fit_pois)
#> stepFMR summary | family: poisson | G: 2,3 | procedure: complete
#> 
#> Variable selection
#>   Direction: forward  |  Alpha: 0.05
#>   Selected:      x1, x2
#>   Not selected:  x3
#>   Final formula: y ~ x1 + x2
#> 
#> Step 1  (none included)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x1                    6.97e-40   *
#>   x2                    6.69e-07    
#>   x3                    4.75e-03    
#>   Chosen: x1
#> 
#> Step 2  (included: x1  |  excluded: x2, x3)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x2                    4.97e-20   *
#>   x3                    2.19e-01    
#>   Chosen: x2
#> 
#> Step 3  (included: x1, x2  |  excluded: x3)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x3                    6.68e-02    
#>   No eligible candidate -- stopping.
#> 
#> Effect-type determination
#>   Direction: forward  |  Alpha: 0.05
#>   Heterogeneous: x1
#>   Homogeneous:   x2
#>   Final formula: y ~ x2 + (x1 | group)
#> 
#> Step 1  (all homogeneous)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x1                    6.29e-24   *
#>   x2                    1.03e-02    
#>   Chosen: x1
#> 
#> Step 2  (heterogeneous: x1  |  homogeneous: x2)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x2                    9.29e-01    
#>   No eligible candidate -- stopping.
#> 
#> Best fit (G: 2)
#> 
#> Mixing proportions:
#>     g1     g2 
#> 0.3812 0.6188 
#> 
#> Heterogeneous coefficients:
#>    (Intercept)      x1
#> g1     -0.3487 -0.9570
#> g2      0.5004  0.9483
#> 
#> Homogeneous coefficients:
#>    x2 
#> 0.365
```

``` r
# Binomial (25 trials per observation)
sim_bin <- fmr_binomial_example
fit_bin <- stepFMR(sim_bin$formula, sim_bin$data, G_max = 3, family = "binomial")
print_steps(fit_bin$variable_selection)
#> Step 1  (none included)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x1                    1.94e-114  *
#>   x2                    3.31e-04    
#>   x3                    3.73e-04    
#>   Chosen: x1
#> 
#> Step 2  (included: x1  |  excluded: x2, x3)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x2                    3.12e-72   *
#>   x3                    5.36e-03    
#>   Chosen: x2
#> 
#> Step 3  (included: x1, x2  |  excluded: x3)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x3                    2.19e-03   *
#>   Chosen: x3
print_steps(fit_bin$effect_determination)
#> Step 1  (all homogeneous)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x1                    1.77e-182  *
#>   x2                    3.90e-01    
#>   x3                    6.30e-03    
#>   Chosen: x1
#> 
#> Step 2  (heterogeneous: x1  |  homogeneous: x3, x2)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x2                    1.47e-01    
#>   x3                    2.09e-03   *
#>   Chosen: x3
#> 
#> Step 3  (heterogeneous: x1, x3  |  homogeneous: x2)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x2                    1.12e-01    
#>   No eligible candidate -- stopping.
print(fit_bin)
#> stepFMR | family: binomial | G: 2,3 | procedure: complete
#> 
#> Variable selection
#>   Direction: forward  |  Alpha: 0.05
#>   Selected:      x1, x2, x3
#>   Not selected:  none
#>   Final formula: y ~ x1 + x2 + x3
#> 
#> Effect-type determination
#>   Direction: forward  |  Alpha: 0.05
#>   Heterogeneous: x1, x3
#>   Homogeneous:   x2
#>   Final formula: y ~ x2 + (x1 + x3 | group)
#> 
#> Best fit (G: 2)
#> 
#> Mixing proportions:
#>     g1     g2 
#> 0.5949 0.4051 
#> 
#> Heterogeneous coefficients:
#>    (Intercept)      x1      x3
#> g1      1.0408  1.5209 -0.0222
#> g2     -0.9977 -1.5368  0.1357
#> 
#> Homogeneous coefficients:
#>    x2 
#> 0.553
```
