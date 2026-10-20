<h1 align="center">feature-mediated-dependency-detection</h1>
<h3 align="center">Detecting pairwise interactions with feature-based mutual information (<em>MI<sub>F</sub></em>)</h3>

<p align="center">
  <a href="https://doi.org/10.1103/qnx2-yp4c"><img src="https://img.shields.io/badge/paper-Phys%20Rev%20Research%20(2025)-informational" alt="Paper"/></a>
  <a href="https://doi.org/10.5281/zenodo.17188607"><img src="https://zenodo.org/badge/1062958385.svg" alt="DOI"/></a>
  <img src="https://img.shields.io/badge/language-R-blue" alt="Language: R"/>
</p>

This repository provides a method for detecting pairwise interactions between time series, based on a feature-based adaptation of information-theoretic measures.
Rather than comparing the raw signals, the method computes interpretable time-series features of the source in a sliding window (the 24 _catch24_ features: the 22 [catch22](https://github.com/DynamicsAndNeuralSystems/catch22) features plus the mean and standard deviation), and tests which feature carries information about the target.
This lets it detect interactions that act over long timescales, and tells you which property of the source (e.g., its variance or autocorrelation) mediates the dependence.
The method is implemented in **R**, using the [Java Information Dynamics Toolkit (JIDT)](https://github.com/jlizier/jidt) for the information-theoretic estimates.

The code accompanies the open-access paper:

> A. Nguyen, O. McMullin, J. T. Lizier, and B. D. Fulcher, "[Feature-based information-theoretic approach for detecting interpretable, long-timescale pairwise interactions from time series](https://doi.org/10.1103/qnx2-yp4c)", _Physical Review Research_ **7**, 043283 (2025) ([arXiv:2404.05929](https://arxiv.org/abs/2404.05929)).

It contains the code to apply the method to a pair of time series, and the code to reproduce the simulations and figures in the paper.

## Contents

- [Installation](#installation)
- [Usage](#usage)
  - [Applying the method to a pair of time series](#applying-the-method-to-a-pair-of-time-series)
  - [Simulating data](#simulating-data)
  - [Reproducing the figures](#reproducing-the-figures)
- [Repository structure](#repository-structure)
- [Citation](#citation)

## Installation

You need [R](https://www.r-project.org) and a Java Development Kit (JDK), which the `rJava` package requires.
Clone the repository, and run all commands below from its top-level folder:

```bash
git clone https://github.com/DynamicsAndNeuralSystems/feature-mediated-dependency-detection.git
cd feature-mediated-dependency-detection
```

### Step 1: R packages

All required R packages are listed in `requirements.txt`.
They are installed automatically when you run `main_simulate.R` or `main_evaluate.R`.
To install them manually, run in R:

```r
packages <- readLines("requirements.txt")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "http://cran.us.r-project.org")
  }
}
```

If R warns that `lib = "/usr/local/lib/R/site-library"` is not writable, install into a folder you can write to:

```bash
mkdir -p <some_writeable_folder>
export R_LIBS_USER="<some_writeable_folder>"
```

### Step 2: JIDT

The JIDT file `infodynamics.jar` is already included in this repository, and `setup.R` loads it.
No additional setup is needed to apply the method or to reproduce the results in the paper.
To use JIDT independently of this repository, see the [JIDT documentation for R](https://github.com/jlizier/jidt/wiki/UseInR).

## Usage

### Applying the method to a pair of time series

Load the libraries and initialize JIDT with `setup.R`, load the functions in `src.R`, and call `detect_dependency_with_catch22()`:

```r
source("setup.R")
source("src.R")

set.seed(1)
source_ts <- rnorm(300)   # the source time series
target_ts <- rnorm(300)   # the target time series

results <- detect_dependency_with_catch22(
  source_ts, target_ts,
  feature_window = 20,
  target_hist = 0,
  number_surrogates = 100
)
head(results)
```

The function computes the catch24 features of the source in a sliding window of `feature_window` samples, and estimates the information between each feature and the target.
For each feature, it assesses significance with a surrogate test and corrects for the number of features using the Holm-Bonferroni method.
It prints a message, either "Statistical dependence detected", followed by the features with a Holm-adjusted p-value below 0.05, or "No significant statistical dependence detected".

```r
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

| Parameter                   | Type           | Description |
|-----------------------------|----------------|-------------|
| `source`                    | numeric vector | The source time series. |
| `target`                    | numeric vector | The target time series. |
| `feature_window`            | integer        | Window size (in samples) for computing the features of the source series. |
| `target_hist`               | integer        | History length of the target series used in the transfer entropy calculation. For mutual information, as used throughout the paper, set this to 0. |
| `delay`                     | integer        | Time delay between source and target. Default is 1. |
| `number_surrogates`         | integer        | Number of surrogate time series generated for the statistical test. Default is 1000. |
| `theiler_window_multiplier` | numeric        | Multiplier for the Theiler window used for dynamic correlation exclusion. Default is 1. |
| `rotating_surrogates`       | logical        | If `TRUE`, use rotating (circular-shift) surrogates for the significance test. Default is `TRUE`. |

#### Returns

A data frame with one row per feature and the following columns:

| Column         | Type    | Description |
|----------------|---------|-------------|
| `Feature`      | string  | Name of the feature computed from the source time series. |
| `TE`           | numeric | Information measure between the feature and the target time series (the mutual information when `target_hist = 0`). |
| `pValue`       | numeric | p-value from surrogate testing. |
| `holm_pValue`  | numeric | p-value adjusted with the Holm-Bonferroni method. |
| `significant`  | logical | Whether the Holm-adjusted p-value is below 0.05. |

### Simulating data

The simulations of the paper generate source processes, extract their features, produce target time series, and compute both the signal-based and the feature-based mutual information.

#### Random noise processes

Simulating the random noise process is computationally intensive, because of the large number of features and time-series lengths considered.
The job script `PBS_simulations.pbs` runs parallelized simulations on a PBS cluster, and you can adapt it to your computing environment.

#### Non-stationary processes

Simulation data for the non-stationary processes can be generated directly with:

```bash
Rscript main_simulate.R --source_type=<AR3|bimodal_spiking> --ts_length=<1000> --driving_feature_timescale=<50|100|150|200>
```

For example:

```bash
Rscript main_simulate.R --source_type=AR3 --ts_length=1000 --driving_feature_timescale=100
```

Each argument is given as `--name=value`.
Each run performs 50 simulations, so it can take a long time.
The results are saved to `results/simulation_studies/new_results.csv`.

The exact simulation dataset used in the paper is also provided, in `results/simulation_studies/featureBasedDependency_simulation_results.csv`.

### Reproducing the figures

To reproduce the plots in the paper, run:

```bash
Rscript main_evaluate.R <path/to/output>
```

All plots are saved into the specified folder, which is created if it does not exist.

| Figure in the paper | Output file |
|---|---|
| Fig. 4a | `RandomResultsMatrix.pdf` |
| Figs. 4b-e | `RandomResultsLinePlots.pdf` |
| Fig. 5a | `AR3_Noise.pdf` |
| Fig. 5b | `bimodal_spiking_Noise.pdf` |
| Fig. 6a | `AR3_DrivingTimescale.pdf` |
| Fig. 6b | `bimodal_spiking_DrivingTimescale.pdf` |
| Fig. 7a | `AR3_CapturingTimescale.pdf` |
| Fig. 7b | `bimodal_spiking_CapturingTimescale.pdf` |

## Repository structure

| File or folder | Contents |
|---|---|
| `requirements.txt` | The R packages required. |
| `setup.R` | Loads the libraries and initializes JIDT. |
| `src.R` | Core functions for the simulations and analyses, including `detect_dependency_with_catch22()`. |
| `main_simulate.R` | Runs the simulations: generates the source processes, extracts their features, produces target time series, and computes the signal-based and feature-based mutual information. |
| `main_evaluate.R` | Evaluates the simulation results and generates the plots of the paper. |
| `analyses/` | The scripts called by `main_simulate.R` and `main_evaluate.R`. |
| `PBS_simulations.pbs` | Job script for running the simulations on a PBS cluster. |
| `results/simulation_studies/` | The simulation results of the paper, and the output of new simulations. |
| `infodynamics.jar` | The JIDT library. |

## Citation

If you use this code, please cite the paper:

```bibtex
@article{Nguyen2025Feature,
  title   = {Feature-based information-theoretic approach for detecting interpretable, long-timescale pairwise interactions from time series},
  author  = {Nguyen, Aria and McMullin, Oscar and Lizier, Joseph T. and Fulcher, Ben D.},
  journal = {Physical Review Research},
  volume  = {7},
  number  = {4},
  pages   = {043283},
  year    = {2025},
  doi     = {10.1103/qnx2-yp4c}
}
```

The code is archived on Zenodo: [10.5281/zenodo.17188607](https://doi.org/10.5281/zenodo.17188607).
