# Test Gaussian mixture regression

set.seed(10)

n <- 500
G_values <- 2:3

x1 <- rnorm(n)
x2 <- rnorm(n)

true_g <- sample(1:2, n, replace = TRUE, prob = c(0.45, 0.55))

y <- ifelse(
  true_g == 1,
  1 + 2.0 * x1,
  -1 - 2.0 * x1
) + 0.5 * x2 + rnorm(n, sd = 0.5)

dat_gauss <- data.frame(y = y, x1 = x1, x2 = x2)

control_gauss <- build_control(
  n_init = 5,
  n_kmeans_init = 2,
  max_iter = 150,
  verbose = FALSE
)

model_gauss <- FMRModel$new(
  formula = y ~ x1 + x2,
  data = dat_gauss,
  G_values = G_values,
  family = "gaussian",
  control = control_gauss
)

prepared_gauss <- prepare_data(
  model = model_gauss,
  included = model_gauss$predictors,
  common = "x2"
)

fits_gauss <- fit_across_G(
  model = model_gauss,
  prepared_data = prepared_gauss
)

cat("\n================ Gaussian results ================\n")
print(lapply(fits_gauss, function(fit) {
  list(
    G = fit$G,
    loglik = fit$loglik,
    bic = fit$bic,
    converged = fit$converged,
    pi_g = fit$pi_g,
    sigma_g = fit$sigma_g,
    beta_g = fit$beta_g,
    beta = fit$beta
  )
}))


# Test Poisson mixture regression

set.seed(20)

n <- 500
G_values <- 2:3

x1 <- rnorm(n)
x2 <- rnorm(n)

true_g <- sample(1:2, n, replace = TRUE, prob = c(0.45, 0.55))

eta <- ifelse(
  true_g == 1,
  0.4 + 0.8 * x1,
  -0.4 - 0.8 * x1
) + 0.3 * x2

mu <- exp(eta)
y <- rpois(n, lambda = mu)

dat_pois <- data.frame(y = y, x1 = x1, x2 = x2)

control_pois <- build_control(
  n_init = 10,
  n_kmeans_init = 3,
  max_iter = 150,
  irwls_max_iter = 50,
  irwls_tol = 1e-8,
  weight_floor = 1e-10,
  verbose = FALSE
)

model_pois <- FMRModel$new(
  formula = y ~ x1 + x2,
  data = dat_pois,
  G_values = G_values,
  family = "poisson",
  control = control_pois
)

prepared_pois <- prepare_data(
  model = model_pois,
  included = model_pois$predictors,
  common = "x2"
)

fits_pois <- fit_across_G(
  model = model_pois,
  prepared_data = prepared_pois
)

cat("\n================ Poisson results ================\n")
print(lapply(fits_pois, function(fit) {
  list(
    G = fit$G,
    loglik = fit$loglik,
    bic = fit$bic,
    converged = fit$converged,
    pi_g = fit$pi_g,
    sigma_g = fit$sigma_g,
    beta_g = fit$beta_g,
    beta = fit$beta,
    irwls_iterations = fit$irwls_iterations,
    irwls_converged = fit$irwls_converged
  )
}))



# Test Bernoulli-binomial mixture regression

set.seed(30)

n <- 500
G_values <- 2:3

x1 <- rnorm(n)
x2 <- rnorm(n)

binomial_size <- sample(15:30, n, replace = TRUE)

true_g <- sample(
  1:2,
  size = n,
  replace = TRUE,
  prob = c(0.45, 0.55)
)

eta <- ifelse(
  true_g == 1,
  -1.25 + 0.75 * x1,
  1.25 - 0.75 * x1
) + 0.2 * x2

p <- stats::plogis(eta)

y <- stats::rbinom(
  n = n,
  size = binomial_size,
  prob = p
)

dat_bin_grouped <- data.frame(
  y = y,
  .binom_size = binomial_size,
  x1 = x1,
  x2 = x2,
  true_g = true_g
)

control_bin <- build_control(
  n_init = 20,
  n_kmeans_init = 5,
  max_iter = 200,
  irwls_max_iter = 100,
  irwls_tol = 1e-8,
  weight_floor = 1e-10,
  verbose = FALSE
)

model_bin <- FMRModel$new(
  formula = y ~ x1 + x2,
  data = dat_bin_grouped,
  G_values = G_values,
  family = "binomial",
  control = control_bin
)

prepared_bin <- prepare_data(
  model = model_bin,
  included = model_bin$predictors,
  common = "x2"
)

prepared_bin$binomial_size[1:10]
range(prepared_bin$binomial_size)

fits_bin <- fit_across_G(
  model = model_bin,
  prepared_data = prepared_bin
)

print(lapply(fits_bin, function(fit) {
  list(
    G = fit$G,
    loglik = fit$loglik,
    bic = fit$bic,
    converged = fit$converged,
    pi_g = fit$pi_g,
    beta_g = fit$beta_g,
    beta = fit$beta,
    irwls_iterations = fit$irwls_iterations,
    irwls_converged = fit$irwls_converged
  )
}))
