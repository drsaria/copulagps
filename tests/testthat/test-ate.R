set.seed(7)
sim <- simulate_treatments(n = 150, copula_family = 3, copula_par = 2)
sim$outcome <- 0.4 * sim$y[, 1] - 0.2 * sim$y[, 2] + rnorm(150)
gps <- estimate_gps(sim, margins = c("normal", "normal"), verbose = FALSE)

test_that("estimate_ate returns gps_ate", {
  ate <- estimate_ate(gps, outcome = sim$outcome, verbose = FALSE)
  expect_s3_class(ate, "gps_ate")
  expect_true(is.numeric(ate$ate_estimates[[1]]$ate))
  expect_true(ate$ate_estimates[[1]]$se > 0)
})

test_that("get_augmented_data returns correct structure", {
  d <- get_augmented_data(gps)
  expect_true(all(c("T1", "T2", "w") %in% names(d)))
  expect_equal(nrow(d), 150)
})

test_that("apply_augmented calls user function correctly", {
  result <- apply_augmented(gps, outcome = sim$outcome,
    FUN = function(d) lm(outcome ~ T1 + T2, data = d, weights = d$w))
  expect_s3_class(result, "lm")
  expect_true("T1" %in% names(coef(result)))
})

test_that("apply_augmented errors for non-function FUN", {
  expect_error(apply_augmented(gps, FUN = "not_a_function"), "'FUN' must be a function")
})
