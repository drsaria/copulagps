set.seed(42)
sim <- simulate_treatments(n = 120, copula_family = 3, copula_par = 2)
sim$outcome <- 0.5 * sim$y[, 1] - 0.3 * sim$y[, 2] + rnorm(120)

test_that("run_chains returns gps_chains with correct M", {
  chains <- run_chains(sim, margins = c("normal", "normal"), M = 3,
                       seed = 1, verbose = FALSE)
  expect_s3_class(chains, "gps_chains")
  expect_equal(chains$M, 3L)
  expect_equal(chains$mode, "gps")
})

test_that("run_chains in ate mode returns ate objects", {
  chains <- run_chains(sim, margins = c("normal", "normal"), M = 3,
                       seed = 1, outcome = sim$outcome, verbose = FALSE)
  expect_equal(chains$mode, "ate")
  ok <- !vapply(chains$chains, is.null, logical(1))
  expect_true(all(vapply(chains$chains[ok], inherits, logical(1), "gps_ate")))
})

test_that("combine_rubin returns rubin_result", {
  chains <- run_chains(sim, margins = c("normal", "normal"), M = 3,
                       seed = 1, outcome = sim$outcome, verbose = FALSE)
  result <- combine_rubin(chains)
  expect_s3_class(result, "rubin_result")
  expect_true(is.numeric(result$estimate))
  expect_true(result$se > 0)
  expect_true(result$fmi >= 0 && result$fmi <= 1)
})

test_that("combine_rubin requires M >= 2", {
  chains <- run_chains(sim, margins = c("normal", "normal"), M = 3,
                       seed = 1, outcome = sim$outcome, verbose = FALSE)
  chains$chains[[1]] <- NULL
  chains$chains[[2]] <- NULL
  expect_error(combine_rubin(chains), "fewer than 2")
})
