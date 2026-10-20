# Run the test suite from the repository root with: Rscript tests/run_tests.R
library(testthat)
test_dir("tests/testthat", stop_on_failure = TRUE)
