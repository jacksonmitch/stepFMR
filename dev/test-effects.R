rm(list = ls())
devtools::load_all()


test_data <- do.call(simulate_fmr,
                     c(scenarios$two_group_effects,
                       list(n = 400, seed = 123)))$data


model <- FMRModel$new(
  formula = formula(test_data),
  data = test_data,
  G_values = 2:3,
  family = "gaussian",
  control = fmr_control(n_init = 5)
)

forward_fit <- determine_effects(
  model = model,
  direction = "forward"
)
print(forward_fit)
#
# backward_fit <- determine_effects(
#   model = model,
#   direction = "backward"
# )
# print(backward_fit)

bench <- microbenchmark::microbenchmark(
  forward_fit = determine_effects(
    model = model,
    direction = "forward"
  ),
  times = 10
)
print(bench)
