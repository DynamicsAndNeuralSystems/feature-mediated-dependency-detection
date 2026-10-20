load_script_functions("analyses/evaluate_methods_performance.R")

test_that("Holm correction uses one family of 24 catch24 MI_f tests per run and capturing timescale", {
  catch24 <- read.csv(file.path(repo_root, "Feature_name_mapping.csv"))$Short_name
  set.seed(3)
  one_run <- function(capturing_timescale, p_features) {
    data.frame(seed = 1, process = "AR3", ts_length = 1000, driving_feature = "autocorrelation_lag1",
               driving_feature_timescale = 100, capturing_feature_timescale = capturing_timescale,
               linear_cor = 0.4,
               capturing_feature = c("MI_s", "null", catch24, "autocorrelation_lag1"),
               inference_type = c("MI_s", "null", rep("MI_f", 25)),
               pValue = c(0.01, 0.03, p_features, 0.02))
  }
  p_a <- c(0.0005, 0.001, runif(22, 0.01, 1))
  p_b <- c(0.0015, runif(23, 0.01, 1))
  res <- rbind(one_run(50, p_a), one_run(75, p_b))

  out <- as.data.frame(adjust_pvalues_holm(res, family_features = catch24))

  for (l in c(50, 75)) {
    sub <- out[out$capturing_feature_timescale == l, ]
    fam <- sub$inference_type == "MI_f" & sub$capturing_feature %in% catch24
    p_raw <- if (l == 50) p_a else p_b
    expect_equal(sub$pValue_holm_corrected[fam], p.adjust(p_raw, method = "holm"))
    # the driving feature added as an extra capturing feature is adjusted as a 25th member of the family
    extra <- sub$capturing_feature == "autocorrelation_lag1"
    expect_equal(sub$pValue_holm_corrected[extra], p.adjust(c(p_raw, 0.02), method = "holm")[25])
    # MI_s and z_null are single tests and are not adjusted
    single <- sub$inference_type %in% c("MI_s", "null")
    expect_equal(sub$pValue_holm_corrected[single], c(0.01, 0.03))
  }
  # A p-value of 0.0015 is significant in a family of 24 (24 * 0.0015 = 0.036), but was not when the
  # family pooled both timescales and the MI_s/null/extra rows (54 tests)
  expect_lt(out$pValue_holm_corrected[out$capturing_feature_timescale == 75 & out$capturing_feature == catch24[1]], 0.05)
})
