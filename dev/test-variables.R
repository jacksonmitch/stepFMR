rm(list = ls())
devtools::load_all()

test_data <- do.call(simulate_fmr,
                       c(scenarios$three_group_twelve_variables_gaussian,
                         list(n = 500, seed = 1)))$data
# result <- stepFMR(formula(test_data), test_data, control = fmr_control(direction = "backward"))
bench <- microbenchmark::microbenchmark(
  sequential = stepFMR(formula(test_data), test_data,
    control = fmr_control()),
  parallel = stepFMR(formula(test_data), test_data,
    control = fmr_control()),
  times = 1
)
print(bench)
