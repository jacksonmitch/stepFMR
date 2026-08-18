## Generates the three example datasets shipped with the package
## (fmr_gaussian_example, fmr_poisson_example, fmr_binomial_example),
## used in README/vignette examples.

devtools::load_all()

fmr_gaussian_example <- simulate_fmr(
  n = 500,
  betas = rbind(
    c(-2, -3, 1, 0),
    c(2, 3, 1, 0)
  ),
  pi = c(0.4, 0.6),
  sigma = 0.5,
  family = "gaussian",
  seed = 1
)

fmr_poisson_example <- simulate_fmr(
  n = 500,
  betas = rbind(
    c(-0.5, -1.0, 0.3, 0),
    c(0.5, 1.0, 0.3, 0)
  ),
  pi = c(0.4, 0.6),
  family = "poisson",
  seed = 2
)

fmr_binomial_example <- simulate_fmr(
  n = 500,
  betas = rbind(
    c(-1.0, -1.5, 0.5, 0),
    c(1.0, 1.5, 0.5, 0)
  ),
  pi = c(0.4, 0.6),
  family = "binomial",
  size = 25,
  seed = 3
)

usethis::use_data(fmr_gaussian_example, overwrite = TRUE)
usethis::use_data(fmr_poisson_example, overwrite = TRUE)
usethis::use_data(fmr_binomial_example, overwrite = TRUE)