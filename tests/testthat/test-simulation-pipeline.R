load_script_functions("analyses/simulate_processes_and_measure_MI.R")

results_template <- data.frame(seed = integer(), process = character(), ts_length = integer(),
                               driving_feature = character(), driving_feature_timescale = integer(),
                               capturing_feature = character(), capturing_feature_timescale = integer(),
                               linear_cor = numeric(), inference_type = character(), TE = numeric(),
                               pValue = numeric())

test_that("random-source simulation writes the results of every driving feature", {
  result_file <- tempfile(fileext = ".csv")
  write.csv(results_template, result_file, row.names = FALSE)
  assign("result_file", result_file, envir = globalenv())

  # Cache the catch22 sliding-window features (recomputed for every driving feature) to keep the test fast
  orig_features <- generate_catch22_features_sliding_window
  cache <- new.env()
  assign("generate_catch22_features_sliding_window", function(ts, timescale, step_size = 1) {
    key <- paste(timescale, length(ts), format(sum(ts), digits = 17))
    if (is.null(cache[[key]])) cache[[key]] <- orig_features(ts, timescale, step_size)
    cache[[key]]
  }, envir = globalenv())
  on.exit(assign("generate_catch22_features_sliding_window", orig_features, envir = globalenv()), add = TRUE)
  # ... and use few surrogates (this does not change which rows are written)
  orig_compute_TE <- compute_TE
  assign("compute_TE", function(..., number_surrogates = 50) orig_compute_TE(..., number_surrogates = 10),
         envir = globalenv())
  on.exit(assign("compute_TE", orig_compute_TE, envir = globalenv()), add = TRUE)

  ts_length <- 40; w <- 10
  invisible(capture.output(
    simulate_process_and_compute_TE(seed = 1, source_type = "random", ts_length = ts_length,
                                    driving_feature_timescale = w, capturing_feature_timescale = w,
                                    source_target_function = "linear", linear_cor = 0.8,
                                    theiler_window_multiplier = 1)))
  written <- read.csv(result_file)

  set.seed(1)
  driving_features <- orig_features(rnorm(ts_length), w)
  varying <- names(driving_features)[sapply(driving_features, sd) != 0]
  expect_gt(length(varying), 1)
  expect_setequal(unique(written$driving_feature), varying)
  # each driving feature: MI_s, z_null and 24 catch24 capturing features
  expect_equal(as.vector(table(written$driving_feature)[varying]), rep(26, length(varying)))
})

test_that("simulation settings follow the paper", {
  # Fig. 4: random source, l = l~ = 10, beta = 0.2
  expect_equal(simulation_settings("random", 10),
               data.frame(driving_feature_timescale = 10, capturing_feature_timescale = 10, linear_cor = 0.8))
  # Figs 5-6: l~ = 100, beta = 0.2-0.8 at l = 100, and l = 50-200 at beta = 0.6
  s100 <- simulation_settings("AR3", 100)
  expect_setequal(s100$linear_cor[s100$capturing_feature_timescale == 100], c(0.2, 0.4, 0.6, 0.8))
  expect_setequal(s100$capturing_feature_timescale[s100$linear_cor == 0.4], seq(50, 200, 25))
  expect_equal(anyDuplicated(s100), 0)
  # Fig. 7: l = l~ = 50, 150, 200 at beta = 0.6
  s_other <- simulation_settings("bimodal_spiking", 50)
  expect_equal(s_other$driving_feature_timescale, c(50, 150, 200))
  expect_equal(s_other$capturing_feature_timescale, c(50, 150, 200))
  expect_equal(unique(s_other$linear_cor), 0.4)
})

test_that("the PBS script calls the existing entry point", {
  pbs <- readLines(file.path(repo_root, "PBS_simulations.pbs"))
  rscript_line <- grep("Rscript", pbs, value = TRUE)
  script <- sub(".*Rscript ([^ ]+).*", "\\1", rscript_line)
  expect_true(file.exists(file.path(repo_root, script)))
})
