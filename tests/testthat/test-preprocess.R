test_that("prepare_data extracts treatment matrix correctly", {
  df <- data.frame(T1 = rnorm(50), T2 = rnorm(50), age = rnorm(50), Y = rnorm(50))
  d  <- prepare_data(df, treatments = c("T1", "T2"), covariates = "age", outcome = "Y")

  expect_equal(dim(d$y), c(50, 2))
  expect_equal(length(d$outcome), 50)
  expect_length(d$x, 2)
  expect_null(d$lod)
})

test_that("prepare_data respects LOD argument", {
  df <- data.frame(T1 = c(0.1, 0.5, 1.0), T2 = c(0.2, 0.3, 0.8),
                   age = 1:3, Y = rnorm(3))
  d  <- prepare_data(df, treatments = c("T1", "T2"), covariates = "age",
                     lod = c(0.3, 0.25))
  expect_equal(d$lod, c(0.3, 0.25))
})

test_that("prepare_data errors on missing columns", {
  df <- data.frame(T1 = rnorm(10), age = rnorm(10))
  expect_error(prepare_data(df, treatments = c("T1", "T2"), covariates = "age"),
               "Treatment columns not found")
})

test_that("prepare_data errors when treatments overlap covariates", {
  df <- data.frame(T1 = rnorm(10), T2 = rnorm(10))
  expect_error(prepare_data(df, treatments = c("T1", "T2"), covariates = "T1"),
               "appear in both")
})
