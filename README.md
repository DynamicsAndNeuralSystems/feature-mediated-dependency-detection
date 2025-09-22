# feature-mediated-dependency-detection

A method for detecting pairwise interaction based on feature-based adaption of information theoretic metrics

## Table of Contents
- [Background](#background)
- [Installation](#installation)
- [Usage](#usage)
  - [Simulating data](#simulating-data)
  - [Reproducing figures](#reproducing-figures)
  - [Applying the method on a pair of time series](#applying-the-method-on-a-pair-of-time-series)


## Background

This repository contains the code for the paper ["A feature-based information-theoretic approach for detecting interpretable, long-timescale pairwise interactions from time series"](https://arxiv.org/abs/2404.05929). Here we provide the code to reproduce the figures in the paper, as well as the R code to apply this method on a pair of time series data to infer pairwise feature-driven dependency.

## Installation

This repository is written in **R** and requires a small number of dependencies.

### Step 1: Install R packages

All required R packages are listed in `requirements.txt`. They will be installed automatically when you run `main_simulate.R` or `main_evaluate.R` with the parameters required.
If you want to install them manually, run:

``` r
packages <- readLines("requirements.txt")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "http://cran.us.r-project.org")
  }
}
```

### Step 2: [JIDT (Java Information Dynamics Toolkit)](https://github.com/jlizier/jidt)

The code here requires `rJvava` and **JIDT**. - The `infofynamics.jar` file needed for JIDT is already included in this repository.
- `rJava` is listed in `requirements.txt` and will be installed automatically.
- No additional setup is required for reproducing the results in the manuscript.
If you want to use JIDT independently of this repo (for custom experiments), see the [official JIDT documentation](https://github.com/jlizier/jidt/wiki/UseInR).

## Usage

This repository provides everything needed to:
1. Run simulations to generate the dynamical processes studied in the manuscript.
2. Compute mutual information (MI) using both the traditional signal-based approach and the feature-based approach.
3. Evaluate simulation results and reproduce the plots shown in the manuscript.

The repository is structured as follows: 
\* `requirements.txt` --- lists all required R packages. \
\* `setup.R` --- loads libraries and initializes **JIDT**.\
\* `src.R` --- core functions for simulations and analyses.\
\* `main_simulate.R` --- runs simulations, including:\
  \* Generating time-series data for the three dynamical source processes studied.\
  \* Extracting their features.\
  \* Producing target time-series data.\
  \* Computing both signal-based MI and feature-based MI.\
\* `main_evaluate.R` --- evaluates simulation results and generates the manuscript plots.

------------------------------------------------------------------------

### Simulating data

#### Random noise processes

Simulating the random noise process requires intensive computation due to the large number of features and time-series lengths considered.
We provide a job script, `PBS_simulations.pbs`, to run parallelized simulations on a distributed cluster.
You can adapt this script to match your own computing environment.

#### Non-stationary processes

For non-stationary processes, simulation data can be generated directly by running:

``` bash
Rscript main_simulate.R --source_type <AR3|bimodal_spiking> --ts_length <1000> --driving_feature_timescale <50|100|150|200>
```

For example:

``` bash
Rscript main_simulate.R --source_type AR3 --ts_length 1000 --driving_feature_timescale 100
```

### Reproducing figures

## Usage

This code repository includes the code and all setup needed to run the simulations to generate the processes studied in the manuscript, as well as the mutual information measures from traditional approach and from feature-based approach. The code repo also contains the script to evaluate the results from the simulations and generate the plots used in the manuscript. The code is structure as below:

-   All packages required are in requirements.txt

-   Libraries and JIDT initialisation is in setup.R

-   src.R contains functions used in simulations and analysis

-   main_simulate.R is the script to run the simulations (including generating time-series data for the 3 dynamical sources processes studied; their features, the target time series data, signal-based mutual information and feature-based mutual information)

-   main_evaluate.R is the script used to run the analysis on simulation data and generate the plots used in the manuscript

### Simulating data

#### Random noise processes

Simulating the random noise process requires intensive computation due to the large number of features and time-series lengths considered.\
We provide a job script, `PBS_simulations.pbs`, to run parallelized simulations on a distributed cluster that we used. You can adapt this script to match your own computing environment.

#### Non-stationary processes

For non-stationary processes, simulation data can be generated directly by running:

``` bash
Rscript main_simulate.R --source_type <AR3|bimodal_spiking> --ts_length <1000> --driving_feature_timescale <50|100|150|200>
```

If you run the scripts above to generate simulation data, the results will be saved to results/simulation_studies/new_results.csv.

**Note**: The exact simulation dataset used in the manuscript is also provided here: results/simulation_studies/featureBasedDependency_simulation_results.csv

### Reproducing figures

To reproduce the plots from the manuscript, run:

``` bash
Rscript main_evaluate.R --output_dir <path/to/output>
```

All plots will be saved into the specified `output_dir`. If the directory does not exist, it will be created automatically.

#### Mapping of plots to manuscript figures

-   **Figure 4a** → `RandomResultsMatrix.pdf`\
-   **Figures 4b--e** → `RandomResultsLinePlots.pdf`\
-   **Figure 5a** → `AR3_Noise.pdf`\
-   **Figure 5b** → `bimodal_spiking_Noise.pdf`\
-   **Figure 6a** → `AR3_DrivingTimescale.pdf`\
-   **Figure 6b** → `bimodal_DrivingTimescale.pdf`\
-   **Figure 7a** → `AR3_CapturingTimescale.pdf`\
-   **Figure 7b** → `bimodal_CapturingTimescale.pdf`

### Applying the method on a pair of time series

You can use the function `detect_dependency_with_catch22()` from `src.R` to infer feature-driven dependencies between a source time series and a target time series using the feature-based Transfer Entropy approach. The function computes Catch22 features for the source series and evaluates the transfer entropy (TE) from each feature to the target. Note that you will need to run `setup.R` before running this function to call all the libraries required and initilize a JIDT object needed for TE calculation. 

#### Function

``` r
detect_dependency_with_catch22(
  source, 
  target, 
  feature_window, 
  target_hist, 
  delay = 1, 
  number_surrogates = 1000, 
  theiler_window_multiplier = 1, 
  rotating_surrogates = TRUE
)
```

#### Parameters

| Parameter                   | Type           | Description                                                                                           |
|-----------------------------|----------------|-------------------------------------------------------------------------------------------------------|
| `source`                    | numeric vector | The source time series.                                                                               |
| `target`                    | numeric vector | The target time series.                                                                               |
| `feature_window`            | integer        | Window size (in samples) for computing Catch22 features on the source series.                         |
| `target_hist`               | integer        | History length of the target series used in TE calculation. For mutual information, this is set to 0. |
| `delay`                     | integer        | Time delay between source and target for TE calculation. Default is 1.                                |
| `number_surrogates`         | integer        | Number of surrogate time-series to generate for statistical testing of TE. Default is 1000.           |
| `theiler_window_multiplier` | numeric        | Multiplier for the Theiler window for dynamic correlation exclusion. Default is 1.                    |
| `rotating_surrogates`       | logical        | If TRUE, use rotating surrogates for TE significance testing. Default is TRUE.                        |

#### Returns

A data frame (`features_results`) with the following columns:

| Column       | Type    | Description                                                            |
|--------------|---------|------------------------------------------------------------------------|
| Feature      | string  | Name of the feature extracted from the source time series.             |
| TE           | numeric | Transfer entropy value between the feature and the target time series. |
| pValue       | numeric | p-Value from surrogate testing for the TE measurement.                 |
| holm-pValue  | numeric | p-Value adjusted with Holm-Bonferroni method                            |
| significance | boolean | Whether the Holm-adjusted pValue is statistically significant !        |

The pvalue are then adjusted using Holm-Bonferroni method, and if any feature has a Holm-adjusted p-value \< 0.05, the function prints "Statistical dependence detected" along with the list of the features with signficant Holm-adjusted pvalues. Otherwise, it prints "No significant statistical dependence detected."
