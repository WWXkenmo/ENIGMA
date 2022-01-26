## DEconvolution based on Regularized Matrix Completion algorithm (ENIGMA)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5907208.svg)](https://doi.org/10.5281/zenodo.5907208)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5906932.svg)](https://doi.org/10.5281/zenodo.5906932)


![ENIGMA](https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/Fig1.png)

## ENIGMA
ENIGMA has three main steps. First, ENIGMA requires cell type reference expression matrix (signature matrix), which could be derived from either FACS RNA-seq or scRNA-seq datasets through calculating the average expression value of each gene from each cell type. Previous researches have widely used reference matrix curated from different platforms, for instance, Newman et al. used LM22 immune signature matrix which
derived from microarray platform to deconvolute bulk RNA-seq dataset. However, we
have to note that use references from different platforms would introduce unwanted batch
effect between reference and bulk RNA-seq matrix, especially for the reference matrix
derived from low coverage scRNA-seq dataset. To overcome this challenge, we used
previously presented method that is specifically designed for correcting batch effect among
bulk RNA-seq matrix and reference matrix. Second, ENIGMA applied robust
linear regression model to estimate each cell type fractions among samples based on
reference matrix derived from the first step. Third, based on reference matrix and cell type
fraction matrix estimated from step 1 and step 2, ENIGMA applied constrained matrix
completion algorithm to deconvolute bulk RNA-seq matrix into CSE on sample-level. In
order to constraint the model complexity to prevent overfitting, we proposed to use two
different norm penalty functions to regularize resulted CSE. Finally, the returned CSE could
be used to identify cell type-specific DEG, visualize each gene’s expression pattern on the
cell type-specific manifold space (e.g. t-SNE, UMAP), and build the cell type-specific
co-expression network to identify modules that relevant to phenotypes of interest

## Notes for installation
our newest version of ENIGMA could be downloaded through following step!
### 1. prepare the required packages of ENIGMA
```
install.packages(c("Matrix","S4Vectors","corpcor","MASS","e1071","ggplot2","cowplot","magrittr","purrr","tibble","nnls","doParallel","tidyr","plyr"))
BiocManager::install(c("SingleCellExperiment","scater","Biobase","SummarizedExperiment","sva","preprocessCore"))
```
### 2. install ENIGMA
```
devtools::install_github("WWXKenmo/ENIGMA")
```

## News
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
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [Toy example for running ENIGMA](https://htmlpreview.github.io/?https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/ENIGMA_toy2.html)
* [Apply ENIGMA to resolve latent cell states](https://htmlpreview.github.io/?https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/ENIGMA_cell_heterogeneity1.html)

Please refer to the [document](https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/A-simple-guide-of-ENIGMA.pdf) of ENIGMA for detailed guidence using ENIGMA as a R package. 

## Note
Fundamental hypotheses of the two models of ENIGMA

**Less heterogenous hypothesis (trace norm model)** : The inferred cell type-specific gene expression matrix represents the expression profile of a cell type, while the bulk expression represents the mixture of different cell types of heterogeneous tissues or samples, so there exists some variation in bulk expression driven by the latent cell type compositional change but not the gene expression alteration within a cell type. Therefore, it is natural to hypothesize that the CSE profile is less heterogeneous than bulk expression. We also used the bulk expression matrix as the initialization matrix of each cell type, which can ensure our inferred CSE matrices to have lower rank than bulk expression. Second, low-rank model is also widely used in gene expression imputation or prediction algorithms (ref), because it is universally known that in many biological processes, genes do not act in a solitary manner and rather interact with other genes to play the function. Those interactions make the expression levels of genes to be interdependent. The existence of gene collinearity can result in a highly correlated data matrix. So, assuming the gene expression values lie on a low-dimensional linear subspace, the resulting data matrix may be a low-rank matrix. Therefore, we used trace norm regularizer to promoting reduced rank cell type-specific gene expression matrix.

**Hidden variable hypothesis (maximum l_2 norm model)** : Most of cell type deconvolution algorithms, including ours, are reference-based deconvolution. Using reference-based methods could provide a robust and cost-effective in-silico way to understand the heterogeneity of bulk samples. It also assumes the existence of prior knowledge on the types of cells existing in a sample. These methods may fail to perform accurately when the data includes rare or otherwise unknown cell types with no references incorporated in the algorithm. Therefore, our reconstituted bulk expression profile (Xθ^T) may not include the variation from unknown rare cell types. Ideally, the observed matrix O would be more informative in our reconstituted bulk expression profile (Xθ^T). In other word, the observed bulk expression matrix would have higher rank than Xθ^T. Under this hypothesis, we need to reduce the rank of Xθ^T. We have proved mathematically that controlling the trace norm of reconstituted bulk expression profile (Xθ^T) equals to the maximum L2 norm (X) of CSE profile (see loss design section of Supplementary Notes).


## Tutorial Dataset
the datasets could be downloaded from this repository. ([The link to datasets](https://github.com/WWXkenmo/ENIGMA/tree/master))

## Contact Author
Author: Weixu Wang, Xiaolan Zhou, Dr. Jun Yao, Prof. Ting Ni

Report bugs by opening a new issue on this Github page

Provide suggestions by sending email to maintainer!

Maintainer: Weixu Wang (ken71198@hotmail.com)

## Citation
Wang W, Yao J, Wang Y, et al. Improved estimation of cell type-specific gene expression through deconvolution of bulk tissues with matrix completion[J]. bioRxiv, 2021.
