# feature-based-dependency-detection
A method for detecting pairwise interaction based on feature-based adaption of information theoretic metrics 

## Background
This repository contains the code for the paper ["A feature-based information-theoretic approach for detecting interpretable, long-timescale pairwise interactions from time series"](arxiv.org/abs/2404.05929). 
Here we provide the code to reproduce the figures in the paper, as well as the R code to apply this method on a pair of time series data to infer pairwise feature-driven dependency. 

## Installation
### JIDT prerequisite
JIDT ((Java Information Dynamics Toolkit) installation is needed for this.
- Installation instructions for R can be found [here](https://github.com/jlizier/jidt/wiki/UseInR).  
- Make sure JIDT is installed and accessible from your R environment before running the analysis scripts.

## Reproducing figures in the manuscript
### Simulating data
You can generate the simulation data used in the manuscript by running:
```Rscript main_simulate_R --...```
Note: the data is also provided at results/simulation_studies/featureBasedDependency_simulation_results.csv

### Reproducing figures
#### Figure 4--7
Here we provide the code to reproduce figures 4-7 in the manuscript. You can run the following command with your chosen output folder, and the plots will be saved in that folder.
```Rscript main_evaluate.R ...```
Plot names and Where it is in the manuscript
Figure 4a - RandomResultsMatrix.pdf
Figures 4b-e - RandomResultsLinePlots.pdf
Figure 5a - AR3_Noise.pdf
Figure 5b - bimodal_spiking_Noise.pdf
Figure 6a - AR3_DrivingTimescale.pdf
Figure 6b - bimodal_DrivingTimescale.pdf
Figure 7a - AR3_CapturingTimescale.pdf
Figure 7b - bimodal_CapturingTimescale.pdf


## Applying the method on a pair of time series
