### ENIGMA analysis
To reproduce part of the simulation and analysis results in our manuscript, here we give several example code for data simulation and analysis.
tips:
1. Before users run following example, user need to source the ENIGMA code at first
```
source("ENIGMA.R")
```
2. Before running the CTS-DE analysis, user need to source related functions
```
source("DEG_analysis_uile_function.R")
```
Both scripts could be downloaded through following links:

[ENIGMA.R](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA_Script/ENIGMA.R)

[DEG_analysis_uile_function.R](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA_Script/DEG_analysis_uile_function.R)

## Code for main text
* [Using tumor scRNA-seq data to simulate pseudo-bulk sample and perform deconvolution](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA_Script/Simulation%20(scRNA-seq).R)

* [Simulate cell type-specific gene expression profile to assess CTS-DEGs detection accuracy](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA_Script/Simulation%20(DEG).R)

  [Benchmark CTS-DEG detection performance](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA_Script/DEG_analysis.R)

* [Identify four latent cell senescence states in Fibroblast](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA_Script/latentCellState.R)

* [Identify cell type specific pseudo-trajectory](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA_Script/ESCO_path.R)

### Real Data analysis
* [Deconvolution analysis for arthritis patients](https://htmlpreview.github.io/?https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/Real_Data_Analysis/pancreas/Beta-cell-type-specific-network-in-pancreas-islet-tissues.html)
## Code for supplementary note
* [Using COVID-19 PBMC single cell datasets to explain the role of parameter alpha](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA_Script/Simulation(mutilPaltforms).R)

* [Attaching noise to reference profile to explain the role of parameter alpha through CTS-DE analysis](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA_Script/ADMM_noise.R)

* [Remove the spurious correlations betweeen inferred CSE profiles and cell type fractions](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA_Script/Normalize_celltype_fractions.R)

* [Post-hoc nonnegative constraint](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA_Script/NegativeValueEffects.R)

* [Gradient Renormalization improve the performance of maximum L2 norm model](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA_Script/Renomarlization_solver_compare_new.R)
  Note: we also write a [document](https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/Why-fixed-renormalized-gradient-norm-size-%3D-200.pdf) to introduce why we use 200 as the new gradient norm size of renormalized gradient.
