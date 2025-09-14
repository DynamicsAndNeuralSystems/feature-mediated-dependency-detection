# Entry point for evaluating simulated processes and generating plots

# Validate arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
output_dir <- args[1]

if (is.null(output_dir)) {
  stop("Usage: Rscript main_evaluate.R <output_dir>")
}

if (!dir.exists(output_dir)) {
  message("Output directory does not exist. Creating: ", output_dir)
  dir.create(output_dir, recursive = TRUE)
}

# Install required packages and set up
#-----------------------------------------
packages <- readLines("requirements.txt")

# install required packages if they are not installed
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p,repos = "http://cran.us.r-project.org")
  }
}

source("setup.R")

# Run evaluation 
#-----------------------------------------
# Pass output_dir into the evaluation script
source("analyses/evaluate_methods_performance.R", local = TRUE)

message("Evaluation complete.")
