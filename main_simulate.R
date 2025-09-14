#-----------------------------------------
#Install required packages and set up
#-----------------------------------------
packages <- readLines("requirements.txt")

# install required packages if they are not installed
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p,repos = "http://cran.us.r-project.org")
  }
}

source("setup.R")

source("src.R")

#-----------------------------------------
# Run simulation studies
#-----------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Find the arguments
for (arg in args) {
  if (grepl("--source_type=", arg)) {
    source_type <- gsub("--source_type=", "", arg)
  }
  if (grepl("--ts_length=", arg)) {
    ts_length <- gsub("--ts_length=", "", arg)
  }
  if (grepl("--driving_feature_timescale=", arg)) {
    driving_feature_timescale <- gsub("--driving_feature_timescale=", "", arg)
  }
  if (grepl("--seed=", arg)) {
    seed <- gsub("--seed=", "", arg)
  }
}

# Call the evaluate-method-performance-simulated-processes.R script with the arguments
if (!is.null(source_type) & ! is.null(ts_length) & ! is.null(driving_feature_timescale)) {
  source("analyses/simulate_processes_and_measure_MI.R")
} else {
  stop("Error: Missing required arguments - source_type and ts_length and driving_feature_timescale ")
}
