test_that("detect_dependency_with_catch22 uses its delay argument", {
  # Target driven by the mean of the source window ending d = 4 samples earlier:
  # y[s + d] = z[s] + noise, where z[s] is the mean of x[(s - w + 1):s].
  set.seed(11)
  n <- 90; w <- 5; d <- 4
  x <- rnorm(n)
  z <- as.numeric(stats::filter(x, rep(1 / w, w), sides = 1))
  y <- rnorm(n, sd = 0.2)
  y[(w + d):n] <- z[w:(n - d)] + rnorm(n - w - d + 1, sd = 0.2)

  run <- function(delay) {
    set.seed(1)
    # (few surrogates: only the MI values are compared, so the warning about p-value resolution is expected)
    suppressWarnings(suppressMessages(detect_dependency_with_catch22(x, y, feature_window = w, target_hist = 0,
                                                                     delay = delay, number_surrogates = 20)))
  }
  mi_delay1 <- run(1)
  mi_delay4 <- run(d)

  mean_mi_1 <- mi_delay1$TE[mi_delay1$Feature == "DN_Mean"]
  mean_mi_4 <- mi_delay4$TE[mi_delay4$Feature == "DN_Mean"]
  # At the true delay the window mean shares a lot of information with the target
  expect_gt(mean_mi_4, 0.5)
  expect_gt(mean_mi_4, mean_mi_1 + 0.3)
})
