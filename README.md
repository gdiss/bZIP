# bZIP
all code regarding the bZIP paper "The genetic architecture of the human bZIP interaction network" by Bendel et al., (submitted as of August 2025)

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
