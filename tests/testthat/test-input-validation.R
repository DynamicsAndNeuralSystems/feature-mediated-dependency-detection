test_that("missing values in source or target give a clear error", {
  set.seed(2)
  x <- rnorm(40); y <- rnorm(40)
  x_na <- x; x_na[17] <- NA
  y_na <- y; y_na[5] <- NaN
  expect_error(suppressMessages(detect_dependency_with_catch22(x_na, y, feature_window = 10, target_hist = 0)),
               "source contains 1 missing or non-finite value.*index 17")
  expect_error(suppressMessages(detect_dependency_with_catch22(x, y_na, feature_window = 10, target_hist = 0)),
               "target contains 1 missing or non-finite value.*index 5")
})

test_that("a flat stretch in the source gives a clear error naming the affected features", {
  set.seed(2)
  x <- rnorm(40); y <- rnorm(40)
  x[11:30] <- 0  # windows of length 10 that fall entirely in this stretch are constant
  expect_error(suppressMessages(detect_dependency_with_catch22(x, y, feature_window = 10, target_hist = 0,
                                                               number_surrogates = 10)),
               "non-finite values.*DN_HistogramMode_5 \\(11 windows\\).*constant stretches")
})

test_that("a constant target gives a clear error", {
  set.seed(2)
  expect_error(suppressMessages(detect_dependency_with_catch22(rnorm(40), rep(1, 40), feature_window = 10,
                                                               target_hist = 0)),
               "target is constant")
})
