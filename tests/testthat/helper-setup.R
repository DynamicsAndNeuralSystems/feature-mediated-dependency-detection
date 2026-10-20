# Minimal setup for the tests: the method code in src.R plus JIDT, without the plotting stack in setup.R.
suppressPackageStartupMessages({
  library(Rcatch22)
  library(tidyr)
  library(dplyr)
  library(rJava)
  library(Rlab)
})

repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))

.jinit()
.jaddClassPath(file.path(repo_root, "infodynamics.jar"))
teCalc <- .jnew("infodynamics/measures/continuous/kraskov/TransferEntropyCalculatorKraskov")
assign("teCalc", teCalc, envir = globalenv())

sys.source(file.path(repo_root, "src.R"), envir = globalenv())

# Load only the function definitions from an analysis script (the scripts have top-level side effects).
load_script_functions <- function(script, envir = globalenv()) {
  for (e in parse(file.path(repo_root, script))) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
        identical(e[[3]][[1]], as.name("function"))) {
      eval(e, envir)
    }
  }
}

# Rcatch22 warns once per session that CO_f1ecac changed type in 0.1.14; trigger it here, silently.
invisible(suppressWarnings(catch22_all(rnorm(50), catch24 = TRUE)))
