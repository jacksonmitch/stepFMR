
# stepFMR

<!-- badges: start -->

[![R-CMD-check](https://github.com/jacksonmitch/stepFMR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jacksonmitch/stepFMR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Finite Mixture Regression models allows relationships between the
response ($\mathbf{Y}$) and the predictors
$(\mathbf{X}_1,\dots, \mathbf{X}_p)$ to differ between a finite number
of subpopulations. The R package, stepFMR implements stepwise procedures
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

This example simulates data from a two-component mixture where the true
answer is known, runs the stepwise WEST procedure, and then fits a model
with a specified structure. The vignette (`vignette("stepFMR")`)
explains each step in more detail.

### 1. Data generation

`simulate_fmr()` creates data from a mixture regression with
coefficients that you choose, so the fitted results can be checked
against the truth. Each row of `betas` holds the coefficients of one
component, with columns for the intercept and each predictor. Here the
first covariate, `x1`, is heterogeneous, `x2` is homogeneous, and `x3`
has no effect. `pi` is the mixing proportions of each component and
`sigma` is the group specific standard deviations.

``` r
suppressPackageStartupMessages(library(stepFMR))

betas <- rbind(
  c(-2, -3, 1, 0),
  c(2, 3, 1, 0)
)
colnames(betas) <- c("(Intercept)", "x1", "x2", "x3")

sim <- simulate_fmr(
  n = 1000, betas = betas, pi = c(0.3, 0.7), sigma = c(0.5, 0.25),
  family = "gaussian", seed = 1
)
```

The result is a list. `sim$data` is a data frame with the response `y`
and the predictors `x1`, `x2` and `x3`. `sim$formula` is the model
formula `y ~ x1 + x2 + x3`, which names the response on the left and the
predictors on the right, just as in `lm()`. `sim$true_group` records
which component each observation came from. With your own data, pass
your own formula and data frame instead, for example
`stepFMR(y ~ x1 + x2 + x3, data = my_data)`.

### 2. WEST procedure

The number of components $G$ is unknown, and whether a predictor matters
can depend on it. `stepFMR()` therefore tests each predictor under every
number of components from 2 to `G_max` and combines the results into one
weighted p-value (WEST). It runs in two stages. Variable selection
decides which predictors have an effect. Effect determination then
decides whether each selected predictor is heterogeneous or homogeneous.
`G_max = 3` means that models with 2 and 3 components are compared, and
the one with the lowest BIC is reported.

``` r
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
#>   Heterogeneous: x1
#>   Homogeneous:   x2
#>   Final formula: y ~ x2 + (x1 | group)
#> 
#> Best fit (G: 2)
#> 
#> Mixing proportions:
#>     g1     g2 
#> 0.3047 0.6953 
#> 
#> Component standard deviations:
#>     g1     g2 
#> 0.5141 0.2440 
#> 
#> Heterogeneous coefficients:
#>    (Intercept)      x1
#> g1     -2.0488 -2.9971
#> g2      2.0010  3.0164
#> 
#> Homogeneous coefficients:
#>     x2 
#> 0.9977
```

The printed result lists the selected predictors, how each effect was
classified, and the estimates at the best number of components. To see
the WEST procedure, `print_steps()` shows the weighted p-value of every
candidate at each step. A candidate is added when its p-value is below
`alpha` (0.05 by default).

``` r
print_steps(fit$variable_selection)
#> Step 1  (none included)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x1                    0.00e+00   *
#>   x2                    2.34e-18    
#>   x3                    5.41e-01    
#>   Chosen: x1
#> 
#> Step 2  (included: x1  |  excluded: x2, x3)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x2                    0.00e+00   *
#>   x3                    9.33e-01    
#>   Chosen: x2
#> 
#> Step 3  (included: x1, x2  |  excluded: x3)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x3                    1.34e-01    
#>   No eligible candidate -- stopping.
print_steps(fit$effect_determination)
#> Step 1  (all homogeneous)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x1                    2.61e-311  *
#>   x2                    2.85e-01    
#>   Chosen: x1
#> 
#> Step 2  (heterogeneous: x1  |  homogeneous: x2)
#>   Candidate             Weighted p-value
#>   ---------             -------     
#>   x2                    7.67e-02    
#>   No eligible candidate -- stopping.
```

### 3. Fitting

If you already know which predictors have a common effect, `fitFMR()`
fits that structure directly, without the stepwise search.
`common = "x2"` makes `x2` homogeneous, and `G = 2:3` compares models
with 2 and 3 components by BIC.

``` r
fit_fixed <- fitFMR(sim$formula, sim$data, G = 2:3, common = "x2")
fit_fixed
#> fitFMR | family: gaussian | G: 2,3
#>   Formula:       y ~ x1 + x2 + x3
#>   Heterogeneous: x1, x3
#>   Homogeneous:   x2
#> 
#>  G  loglik num_parameters     bic  
#>  2 -804.01             10 1677.10 *
#>  3 -794.38             15 1692.38  
#> 
#> Best fit (G: 2)
#> 
#> Mixing proportions:
#>     g1     g2 
#> 0.3049 0.6951 
#> 
#> Component standard deviations:
#>     g1     g2 
#> 0.5134 0.2432 
#> 
#> Heterogeneous coefficients:
#>    (Intercept)      x1     x3
#> g1     -2.0484 -2.9993 0.0235
#> g2      1.9998  3.0158 0.0168
#> 
#> Homogeneous coefficients:
#>     x2 
#> 0.9977
```

Mixture components have no natural order, so the fitted components can
be numbered differently from the simulated ones. Because the true groups
are known here, `match_groups()` relabels the fit to agree with them.
The rows of the matrix below can then be compared directly with the rows
of `betas`.

``` r
aligned <- match_groups(fit, sim$true_group)
aligned$best_fit$parameter_values[["beta_g"]]
#>    (Intercept)        x1
#> g1   -2.048839 -2.997093
#> g2    2.000963  3.016419
```

### Other families

The package supports `family = "poisson"` and `family = "binomial"`. For
binomial data the response is `cbind(successes, failures)`, which
`simulate_fmr()` builds for you from the number of trials given by
`size`.

``` r
sim_pois <- simulate_fmr(
  n = 500,
  betas = rbind(c(-0.5, -1, 0.3, 0), c(0.5, 1, 0.3, 0)),
  pi = c(0.4, 0.6), family = "poisson", seed = 2
)
fit_pois <- stepFMR(sim_pois$formula, sim_pois$data, G_max = 3, family = "poisson")
print(fit_pois)
#> stepFMR | family: poisson | G: 2,3 | procedure: complete
#> 
#> Variable selection
#>   Direction: forward  |  Alpha: 0.05
#>   Selected:      x1, x2
#>   Not selected:  x3
#>   Final formula: y ~ x1 + x2
#> 
#> Effect-type determination
#>   Direction: forward  |  Alpha: 0.05
#>   Heterogeneous: x1
#>   Homogeneous:   x2
#>   Final formula: y ~ x2 + (x1 | group)
#> 
#> Best fit (G: 2)
#> 
#> Mixing proportions:
#>     g1     g2 
#> 0.3814 0.6186 
#> 
#> Heterogeneous coefficients:
#>    (Intercept)      x1
#> g1     -0.3451 -0.9554
#> g2      0.4997  0.9490
#> 
#> Homogeneous coefficients:
#>     x2 
#> 0.3647
```

``` r
sim_bin <- simulate_fmr(
  n = 500,
  betas = rbind(c(-1, -3, 0.5, 0), c(1, 3, 0.5, 0)),
  pi = c(0.4, 0.6), family = "binomial", size = 25, seed = 3
)
fit_bin <- stepFMR(sim_bin$formula, sim_bin$data, G_max = 3, family = "binomial")
print(fit_bin)
#> stepFMR | family: binomial | G: 2,3 | procedure: complete
#> 
#> Variable selection
#>   Direction: forward  |  Alpha: 0.05
#>   Selected:      x1, x2
#>   Not selected:  x3
#>   Final formula: y ~ x1 + x2
#> 
#> Effect-type determination
#>   Direction: forward  |  Alpha: 0.05
#>   Heterogeneous: x1
#>   Homogeneous:   x2
#>   Final formula: y ~ x2 + (x1 | group)
#> 
#> Best fit (G: 2)
#> 
#> Mixing proportions:
#>     g1     g2 
#> 0.5896 0.4104 
#> 
#> Heterogeneous coefficients:
#>    (Intercept)      x1
#> g1      1.0062  3.0487
#> g2     -1.0151 -3.0145
#> 
#> Homogeneous coefficients:
#>     x2 
#> 0.5557
```
