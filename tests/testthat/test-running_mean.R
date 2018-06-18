context("Checking a running mean")

test_that("checking mean of runnig mean hasn't changed", {
  set.seed(20180618)
  x <- rnorm(1000)
  output <- RunningMean(x, 100)
  expect_equal(mean(output), -0.00971309884119209)
})

