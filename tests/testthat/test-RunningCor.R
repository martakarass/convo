context("Checking RunningCor")


RunningCor.CONV <- function(x, y, circular){

  out <- numeric()
  W <- length(y)
  l_x <- length(x)

  out <- sapply(1:(l_x - W + 1), function(i){
    segm <- x[i:(i+W-1)]
    cor(segm, y)
  })

  if (circular) {
    out.tail <- sapply((l_x - W + 2):l_x, function(i){
      segm <- c(x[i:l_x], x[1:(W-l_x+i-1)])
      cor(segm, y)
    })
  } else {
    out.tail <- rep(NA, W-1)
  }

  out <- c(out, out.tail)

  return(out)
}


test_that("Checking RunningCor via FFT returns what is equivalent with conventional computation (with circular = TRUE)", {
  x <- sin(seq(0, 1, length.out = 1000) * 2 * pi * 6)
  y <- x[1:100]
  output      <- RunningCor(x, y, circular = TRUE)
  output.conv <- RunningCor.CONV(x, y, circular = TRUE)
  max.diff    <- max(abs(output - output.conv))
  expect_equal(max.diff, 0, tolerance = 1e-8)
})


test_that("Checking RunningCor via FFT returns what is equivalent with conventional computation (with circular = FALSE)", {

  x <- sin(seq(0, 1, length.out = 1000) * 2 * pi * 6)
  y <- x[1:100]

  output      <- RunningCor(x, y, circular = FALSE)
  output.conv <- RunningCor.CONV(x, y, circular = FALSE)

  max.diff    <- max(abs(output - output.conv), na.rm = TRUE)
  expect_equal(max.diff, 0, tolerance = 1e-8)

  expect_true(all(is.na(output) == is.na(output.conv)))
})

