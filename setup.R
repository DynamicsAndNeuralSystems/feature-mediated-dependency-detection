#-----------------------------------------
# Load all the packages necessary for the project
#-----------------------------------------
library(Rcatch22)
library(tidyr)
library(hash)
library(ggplot2)
library(dplyr)
library(rJava)
library(ggpubr)
library(gridExtra)
library(stringr)
library(RColorBrewer)
library(cowplot)
library(ggfittext)
library(plotly)
library(Rlab)
library(PearsonDS)
library(e1071)
library(gsl)
library(R.matlab)
library(htmlTable)
library(DT)
library(kableExtra)
library(corrplot)
library(reshape2)
library(here)
library(numDeriv)
library(pracma)
library(ggrepel)
library(viridis)

#-----------------------------------------
# Initialise JIDT 
#-----------------------------------------
.jinit()
.jaddClassPath(here(here("infodynamics.jar")))
teCalc <-.jnew("infodynamics/measures/continuous/kraskov/TransferEntropyCalculatorKraskov")

