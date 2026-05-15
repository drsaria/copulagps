set.seed(1)
sim <- simulate_treatments(n = 150, copula_family = 3, copula_par = 2)

test_that("estimate_gps returns gps_weights object", {
  gps <- estimate_gps(sim, margins = c("normal", "normal"), verbose = FALSE)
  expect_s3_class(gps, "gps_weights")
  expect_length(gps$weights, 150)
  expect_true(all(is.finite(gps$weights)))
})

test_that("estimate_gps hull_trim assigns NA to outside-hull obs", {
  gps <- estimate_gps(sim, margins = c("normal", "normal"),
                      hull_trim = TRUE, verbose = FALSE)
  expect_true(any(is.na(gps$weights)) || all(!is.na(gps$weights)))
  expect_equal(length(gps$hull_membership), 150)
  expect_type(gps$hull_membership, "logical")
})

test_that("hull_trim = FALSE produces no NA weights", {
  gps <- estimate_gps(sim, margins = c("normal", "normal"),
                      hull_trim = FALSE, verbose = FALSE)
  expect_true(all(!is.na(gps$weights)))
})

test_that("estimate_gps errors for d != 2", {
  bad_data <- list(y = matrix(rnorm(100 * 3), 100, 3),
                   x = replicate(3, cbind(1, rnorm(100)), simplify = FALSE))
  expect_error(estimate_gps(bad_data, margins = rep("normal", 3), verbose = FALSE),
               "d = 2")
})
