context("Checking running cov")

test_that("checking errors on running covariance", {
  set.seed(20180618)

  x <- sin(seq(0, 1, length.out = 1000) * 2 * pi * 6)
  y <- x[1:100]
  expect_silent({
    out1 <- RunningCov(x, y, circular = TRUE)
  })
  expect_error(RunningCov(y, x))
})


