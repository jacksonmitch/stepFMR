rm(list = ls())
devtools::load_all()


four_vars_data <- do.call(simulate_fmr,
                       c(scenarios$three_group_four_variables,
                         list(n = 500, seed = 1)))
two_g_3_vars_data <- do.call(simulate_fmr,
                       c(scenarios$two_group_effects,
                         list(n = 500, seed = 1)))
simulation_study_data_500 <- do.call(simulate_fmr,
                       c(scenarios$four_group_twelve_variables,
                         list(n = 500, seed = 123)))
simulation_study_data_1k <- do.call(simulate_fmr,
                       c(scenarios$four_group_twelve_variables,
                         list(n = 1000, seed = 123)))

test_data <- four_vars_data
result_500 <- stepFMR(test_data$formula, test_data$data, G_max = 4,
    control = fmr_control(parallel = TRUE, direction = "forward"))

result_g3 <- stepFMR(four_vars_data$formula, four_vars_data$data, G_max = 4,
    control = fmr_control(parallel = TRUE, direction = "forward", verbose = TRUE))
result_g2 <- stepFMR(two_g_3_vars_data$formula, two_g_3_vars_data$data, G_max = 4,
    control = fmr_control(parallel = TRUE, direction = "forward"))

bench <- microbenchmark::microbenchmark(
  # sequential = stepFMR(test_data$formula, test_data$data,
  #   control = fmr_control()),
  parallel = stepFMR(test_data$formula, test_data$data,
    control = fmr_control(parallel = TRUE)),
  times = 1
)
print(bench)
