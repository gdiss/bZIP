# bZIP
This repository contains all the code used to reproduce the analyses an results presented in "The genetic architecture of the human bZIP interaction network" by Bendel et al., except the code regarding the deep learning models (which can be found at https://github.com/fmi-basel/Deep-learning-bZIPs)

Requirements:
R version 4.6.0 (2026-04-24)  
MoCHI 1.1  

R packages:  
parallel 4.6.0  
ggplot2 4.0.3  
rgl 1.3.36  
stringdist 0.9.17  
mutscan 1.2.0  
Biostrings 2.80.2  
stringr 1.6.0  
data.table 1.18.6.1  
readxl 1.6.0  
ShortRead 1.70.0  
dplyr 1.2.1  
cowplot 1.2.0  
gridExtra 2.3.1  
ggh4x 0.3.1  
gplots 3.3.0  


1) Put all fastq files from the GEO repository in the 000-data folder
2) Put the binder_optimization.csv file from the github repository of the deep learning models in the 000-data folder
3) Put Table S4 to S12 (except S5) from the supplementary method in the 00-data folder
4) install mutscan from within R (devtools::install_github("fmicompbio/mutscan"))
5) install MoCHI (https://github.com/lehner-lab/MoCHI)
6) run 001-bar_var_assoc_bZIP_screen.Rmd
7) run 002-digestFastq.Rmd
8) run 003-binding_scores.Rmd
9) run 004-mutation_effects.Rmd
10) run 005-thermodynamic_models/001-global_model/002-run_mochi.sh (from without the 001-global_model folder)
11) run 005-thermodynamic_models/002-individual_models/ 002-run_mochi_llop.sh (from within the 002-individual_models folder)
12) run 006-synthetic_bZIPs.Rmd
13) run 007-figures_and_stats.R
