context("Checking RunningVar")

test_that("Checking mean of RunningVar hasn't changed  (with circular = TRUE)" , {
  set.seed(20180618)
  x <- rnorm(1000)
  output <- RunningVar(x, 100, circular = TRUE)
  expect_equal(mean(output), 0.987951966793899)
})

test_that("Checking mean of RunningVar hasn't changed  (with circular = FALSE)", {
  set.seed(20180618)
  x <- rnorm(1000)
  output <- RunningVar(x, 100, circular = FALSE)
  expect_equal(mean(output), 0.978883021627162)
})



RunningVar.CONV <- function(x, W, circular){

  out <- numeric()
  l_x <- length(x)

  for (i in 1:(l_x - W + 1)){
    segm <- x[i:(i+W-1)]
    out  <- c(out, var(segm))
  }

  if (circular) {
    for (i in (l_x - W + 2):l_x){
      segm <- c(x[i:l_x], x[1:(W-l_x+i-1)])
      out  <- c(out, var(segm))
    }
  }

  return(out)
}


test_that("Checking RunningVar via FFT returns what is equivalent with conventional computation (with circular = TRUE)", {
  set.seed(20180618)
  x <- rnorm(1000)
  output      <- RunningVar(x, 100, circular = TRUE)
  output.conv <- RunningVar.CONV(x, 100, circular = TRUE)
  max.diff    <- max(abs(output - output.conv))
  expect_equal(max.diff, 0, tolerance = 1e-8)
})


test_that("Checking RunningVar via FFT returns what is equivalent with conventional computation (with circular = FALSE)", {
  set.seed(20180618)
  x <- rnorm(1000)
  output      <- RunningVar(x, 100, circular = FALSE)
  output.conv <- RunningVar.CONV(x, 100, circular = FALSE)
  max.diff    <- max(abs(output - output.conv))
  expect_equal(max.diff, 0, tolerance = 1e-8)
})

