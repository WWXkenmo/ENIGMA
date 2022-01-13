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

[ENIGMA.R](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ENIGMA.R)

[DEG_analysis_uile_function.R](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/DEG_analysis_uile_function.R)

## Code for main text
* [Using tumor scRNA-seq data to simulate pseudo-bulk sample and perform deconvolution](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/Simulation%20(scRNA-seq).R)

* [Simulate cell type-specific gene expression profile to assess CTS-DEGs detection accuracy](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/Simulation%20(DEG).R)
   
  [Benchmark CTS-DEG detection performance](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/DEG_analysis.R)

* [Identify four latent cell senescence states in Fibroblast](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/latentCellState.R)

* [Identify cell type specific pseudo-trajectory](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/ESCO_path.R)
## Code for supplementary note
