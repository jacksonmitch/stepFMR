test_data <- do.call(
  simulate_fmr,
  c(scenarios$two_group_effects, list(n = 300, seed = 1))
)$data

model <- FMRModel$new(
  formula = formula(test_data),
  data = test_data,
  G_values = 2:3,
  family = "gaussian",
  control = build_control()
)

test_effects <- determine_effects(
  model = model,
  direction = "forward"
)

test_variables <- select_variables(
  model = model,
  direction = "forward"
)

# Checking first line is an easy way to make sure we're using the right prints
test_that("print.fit_fmr uses the correct method", {
  output <- capture.output(print(test_effects$final_fits[[1]]))
  expect_match(output[[1]], "Best fit")
})
test_that("print.determine_effects uses the correct method", {
  output <- capture.output(print(test_effects))
  expect_match(output[[1]], "Effect-type determination")
})
test_that("select_variables uses the correct method", {
  output <- capture.output(print(test_variables))
  expect_match(output[[1]], "Variable selection")
})
