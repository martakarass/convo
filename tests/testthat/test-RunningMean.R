context("Checking RunningMean")

test_that("Checking mean of RunningMean hasn't changed  (with circular = TRUE)" , {
  set.seed(20180618)
  x <- rnorm(1000)
  output <- RunningMean(x, 100, circular = TRUE)
  expect_equal(mean(output), -0.00740616601618133)
})

test_that("Checking mean of RunningMean hasn't changed  (with circular = FALSE)", {
  set.seed(20180618)
  x <- rnorm(1000)
  output <- RunningMean(x, 100, circular = FALSE)
  expect_equal(mean(output), -0.00971309884119209)
})



RunningMean.CONV <- function(x, W, circular){

  out <- numeric()
  l_x <- length(x)

  for (i in 1:(l_x - W + 1)){
    segm <- x[i:(i+W-1)]
    out  <- c(out, mean(segm))
  }

  if (circular) {
    for (i in (l_x - W + 2):l_x){
      segm <- c(x[i:l_x], x[1:(W-l_x+i-1)])
      out  <- c(out, mean(segm))
    }
  }

  return(out)
}


test_that("Checking RunningMean via FFT returns what is equivalent with conventional computation (with circular = TRUE)", {
  set.seed(20180618)
  x <- rnorm(1000)
  output      <- RunningMean(x, 100, circular = TRUE)
  output.conv <- RunningMean.CONV(x, 100, circular = TRUE)
  max.diff    <- max(abs(output - output.conv))
  expect_equal(max.diff, 0, tolerance = 1e-8)
})


test_that("Checking RunningMean via FFT returns what is equivalent with conventional computation (with circular = FALSE)", {
  set.seed(20180618)
  x <- rnorm(1000)
  output      <- RunningMean(x, 100, circular = FALSE)
  output.conv <- RunningMean.CONV(x, 100, circular = FALSE)
  max.diff    <- max(abs(output - output.conv))
  expect_equal(max.diff, 0, tolerance = 1e-8)
})

