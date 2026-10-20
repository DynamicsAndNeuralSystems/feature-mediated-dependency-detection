test_that("circular-shift p-values use each admissible shift once and p = (count + 1)/(n + 1)", {
  set.seed(5)
  n <- 50; tw <- 2
  x <- rnorm(n)
  y <- c(0, x[-n]) + rnorm(n, sd = 0.1)  # y[t + 1] depends strongly on x[t]
  res <- compute_TE(tw, x, y, 0, 1, 1, 1, 1, number_surrogates = 1000)
  # with delay 1 there are N = n - 1 observations, and only N - 2*tw - 1 = 44 distinct circular shifts move
  # every value by more than tw samples; the observed MI beats all of them
  N <- n - 1
  expect_equal(res$nSurrogates, N - 2 * tw - 1)
  expect_equal(res$pValue, 1 / (N - 2 * tw - 1 + 1))

  # independent pair: p is a multiple of 1/(n_surrogates + 1) and never 0
  set.seed(6)
  res0 <- compute_TE(tw, rnorm(n), rnorm(n), 0, 1, 1, 1, 1, number_surrogates = 1000)
  expect_gt(res0$pValue, 0)
  k <- res0$pValue * (res0$nSurrogates + 1)
  expect_equal(k, round(k))

  # with fewer requested surrogates than admissible shifts, the request is honored
  res_small <- compute_TE(tw, x, y, 0, 1, 1, 1, 1, number_surrogates = 20)
  expect_equal(res_small$nSurrogates, 20)
  expect_equal(res_small$pValue, 1 / 21)
})

test_that("random-permutation p-values include the +1 correction", {
  set.seed(5)
  x <- rnorm(60)
  res <- compute_TE(0, x, c(0, x[-60]) + rnorm(60, sd = 0.1), 0, 1, 1, 1, 1, number_surrogates = 99,
                    rotating_surrogates = FALSE)
  expect_equal(res$pValue, 1 / 100)
})

test_that("MI_F does not give false positives under independence at T = 50", {
  # The audit found these 12 independent white-noise pairs gave 5 family-wise false positives with the
  # previous test (circular shifts resampled with replacement and p = count/n)
  n_false_positive <- 0
  for (s in 4001:4012) {
    set.seed(s)
    x <- rnorm(250)[-(1:200)]
    y <- rnorm(250)[-(1:200)]
    expect_warning(res <- suppressMessages(detect_dependency_with_catch22(x, y, feature_window = 10, target_hist = 0)),
                   "no feature can be significant")
    expect_true(all(res$pValue > 0))
    n_false_positive <- n_false_positive + any(res$significant)
  }
  expect_equal(n_false_positive, 0)
})

test_that("MI_F still detects a variance-mediated dependence in a longer series", {
  set.seed(8)
  n <- 620; w <- 20
  x <- rnorm(n)
  z <- sapply(w:(n - 1), function(t) sd(x[(t - w + 1):t]))  # window ending at t drives y[t + 1]
  y <- rnorm(n)
  y[(w + 1):n] <- 0.8 * (z - mean(z)) / sd(z) + 0.2 * rnorm(n - w)
  set.seed(1)
  res <- suppressMessages(detect_dependency_with_catch22(x, y, feature_window = w, target_hist = 0))
  expect_true(res$significant[res$Feature == "DN_Spread_Std"])
})
