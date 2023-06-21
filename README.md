## DEconvolution based on Regularized Matrix Completion algorithm (ENIGMA)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5907208.svg)](https://doi.org/10.5281/zenodo.5907208)

**warnings: the package in Zenodo is no longer updated, please install the newest version!**

<img src="https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/Fig.1.png" alt="ENIGMA" width="600" />

## ENIGMA
A method that accurately deconvolute bulk tissue RNA-seq into single cell-type resolution given the knowledge gained from scRNA-seq. ENIGMA applies a matrix completion strategy to minimize the distance between mixture transcriptome and weighted combination of cell-type-specific expression, allowing quantification of cell-type proportions and reconstruction of cell-type-specific transcriptome.

## Notes for installation
our newest version of ENIGMA could be downloaded through following step!
### 1. prepare the required packages of ENIGMA
```
install.packages(c("Matrix","S4Vectors","corpcor","MASS","e1071","ggplot2","cowplot","magrittr","purrr","tibble","nnls","doParallel","tidyr","plyr","vctrs","matrixStats"))
BiocManager::install(c("SingleCellExperiment","scater","Biobase","SummarizedExperiment","sva","preprocessCore"))
```
### 2. install ENIGMA
install the newest version of ENIGMA
```
devtools::install_github("WWXKenmo/ENIGMA_test")
```

## Notes for usage
### When user need to conduct sample-level score calculation (e.g. Gene Set Activity Analysis (GSVA)) and each sample will be treated independently, please used unnormalized CSE profile

### When user need to integrate all samples to perform calculation (e.g. Differential Expression Gene Analysis, Gene Co-expression Network Inference, Sample Clustering Analysis), please used normalized CSE profile
## News
### release v1.6
updated stop criteria
<img src="https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/movie_pca.gif" width="350" height="350"/>

### release v1.5
1. Build FindCSE_DEG function to perform CTS-DEG analysis
```
DEG = FindCSE_DEG(object,y)
# object: an ENIGMA object
# y: a binary phenotype vector represents case(1) and control(0)
```
please refer to the [CTS-DE document](https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/Cell-Type-Specific-Differential-Expression-Gene-Test.pdf) of detailed guidence of CTS-DE analysis with ENIGMA. [link to example datasets](https://github.com/WWXkenmo/ENIGMA/blob/main/ENIGMA_analysis/exampleDatasets.Rdata)

2. Build GeneSigTest function to filter the genes
ENIGMA now provide a function to help user to identify the genes which could be accurately estimated through our algorithm.
```
res = GeneSigTest(object,filtering=TRUE)
head(res$call)
head(res$pval)
egm = res$egm # the filtered ENIGMA object
```
### a simple implementation in python (pyENIGMA)
we have implement the ENIGMA algorithm in python for those people who want to use ENIGMA in python version

[pyENIGMA](https://github.com/WWXkenmo/pyENIGMA/blob/main/pyENIGMA_case.ipynb)

### release v1.3
1. add plotLossCurve to visualize the training
2. set model_tracker parameter to track the trained model
3. add new solvers to trace norm model
4. improve the ENIGMA_class function

### release v1.1
1. Fixed the bugs in batch_correct
2. Add new functions to re-normalized inferred CSE
3. Update new tutorial
---------------------------------------------------------
## Usage
Please refer to the [document](https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/A-simple-guide-of-ENIGMA.pdf) of ENIGMA for detailed guidence using ENIGMA as a R package.
[link to example datasets](https://github.com/WWXkenmo/ENIGMA/tree/master)
## Tutorial
* [Using ENIGMA to estimate CSE in brain tissue](https://htmlpreview.github.io/?https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/brain_tutorial.html)
* [Apply ENIGMA to resolve latent cell states](https://htmlpreview.github.io/?https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/Identify-cell-type-specific-cell-states.html)

## Note

**Which model users should use and why?**
In summary, both trace norm and maximum L2 norm models show superior performance at different aspects. First, trace norm model poses trace norm regularizer to inferred CSE profiles, and uses low-rank matrix to approximate cell type-specific gene expression, which may help the model to discover better gene variation across samples. Trace norm could also perform better than maximum L2 norm on CTS-DEG identification. Second, maximum L2 norm has assumed that there exist unknown variables (expression of rare cell types or technique variations) in bulk samples, and maximum L2 norm shows better performance on recovering cell type-specific correlation structure even there exists very strong noise in observed bulk expression matrix. So, choosing which model is dependent on what kind of analyses users want to conduct. When users want to define patients/samples subtypes according to cell type-specific gene expression profile (e.g. malignant cell), users could choose the maximum L2 norm model to perform the deconvolution. Besides, when users want to perform cell type-specific analysis of differentially expressed genes, users could choose the trace norm model to perform the deconvolution. Maximum L2 norm is also preferable if users have a large cohort of bulk samples. Finally, the training of maximum L2 norm model is not involved with any inverse matrix calculation or singular value decomposition, so it is very scalable to the large bulk samples. When users want to perform fast deconvolution on the bulk expression dataset with large sample sizes, we suggest to use maximum L2 norm model.

## Contact Author
Author: Weixu Wang, Xiaolan Zhou, Dr. Jun Yao, Prof. Ting Ni

Report bugs by opening a new issue on this Github page

Provide suggestions by sending email to maintainer!

Maintainer: Weixu Wang (ken71198@hotmail.com)

## Citation
Wang W, Yao J, Wang Y, et al. Improved estimation of cell type-specific gene expression through deconvolution of bulk tissues with matrix completion[J]. bioRxiv, 2021.
