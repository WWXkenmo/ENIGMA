#' @title Perform CTS-DEG analysis
#'
#' @description Analysis the differential expression genes for each cell type
#' 
#' @param object ENIGMA object
#' @param FDR_control
#' if use Wang et al., FDR controlled DEG analysis model. Default: TRUE
#'
#' @param FoldChange
#' if output the FoldChange of each gene, if FALSE, then output the expression difference
#'
#' @param covariate
#' The data.frame object contains the covariate information of each sample
#'
#' @return A list object contain the DEG results of each cell types
#'
#'
#' @examples
#' \dontrun{
#' DEG = FindCSE_DEG(object,y)
#' DEG = FindCSE_DEG(object,y,covariate=covariate)
#' head(DEG$celltype1)
#' }
#'
#' @reference
#' Wang J, Roeder K, Devlin B. Bayesian estimation of cell type–specific gene expression with prior derived from single-cell data[J]. Genome research, 2021, 31(10): 1807-1818.
#'
#' @export
FindCSE_DEG <- function(object,y,FDR_control = TRUE,covariate=NULL,FoldChange = FALSE){
    ####
    #identify the DEG of ENIGMA outputs require normalized profile
    if(is.null(object@result_CSE_normalized)){
        stop('CTS-DEG analysis required normalized profile!')
    }
    #convert sce object into 3-order array
    Exp = sce2array(object)
    cellName = dimnames(Exp)[[3]]
    DEG_list = list()
    if(FDR_control){result <- DEG_test(Exp,y,covariate)}else{result <- DEG_test2(Exp,y,covariate)}
    ES_m <- EffectiveSize(Exp,y,FoldChange)
    for(i in cellName){
        Tab_m <- cbind(ES_m[,i],result$pval[,i],result$qval[,i])
        if(FoldChange){colnames(Tab_m) <- c("FoldChange","pvalue","qvalue")}else{
            colnames(Tab_m) <- c("ExpressionDifference","pvalue","qvalue")
        }
        DEG_list[[i]] <- Tab_m
    }
    DEG_list
}

EffectiveSize <- function(X_array,y,FoldChange=FoldChange){
    if(FoldChange){
        ### convert all profile into pseudo positive
        X_array[X_array<0] <- 0
        ES_m <- NULL
        nCell = dim(X_array)[3]
        cellName = dimnames(X_array)[[3]]
        for(i in 1:nCell){
            G = X_array[,,i]
            FC <- apply(G,1,FCcall,y=y)
            ES_m = cbind(ES_m,FC)
        }
        colnames(ES_m) <- cellName
    }else{
        ES_m <- NULL
        nCell = dim(X_array)[3]
        cellName = dimnames(X_array)[[3]]
        for(i in 1:nCell){
            G = X_array[,,i]
            FC <- apply(G,1,DEcall,y=y)
            ES_m = cbind(ES_m,FC)
        }
        colnames(ES_m) <- cellName
    }
    ES_m
}

FCcall <- function(g,y){
    mean(g[y==1])/mean(g[y==0])
}

DEcall <- function(g,y){
    mean(g[y==1]) - mean(g[y==0])
}


DEG_test <- function(X_array,y,covariate=NULL){
    O = array(0,
                  dim = c( dim(X_array)[1],
                           dim(X_array)[3],
                           dim(X_array)[2]),
                  dimnames = list( dimnames(X_array)[[1]],
                                   dimnames(X_array)[[3]],
                                   dimnames(X_array)[[2]])
    )
	for(i in 1:dim(X_array)[3]){
	  O[,i,] <- X_array[,,i]
	}
	###Using ANOVA+glm model to evaluate prediction performance
	result <- test(O,y,covariate)
    pval <- t(result$pval)
    qval <- t(result$qval)
    colnames(qval) <- colnames(pval) <- dimnames(X_array)[[3]]
    rownames(qval) <- rownames(pval) <- dimnames(X_array)[[1]]
	return(list(qval = qval, pval = pval))
}


get_pval <- function (pval, cell_type, K)
{
    pval0 = rep(NA, K)
    names(pval0) = cell_type
    names = intersect(names(pval), cell_type)
    pval0[names] = pval[names]
    return(pval0)
}

pval2qval <- function (pval, A, y, covariate = NULL)
{
    ng = nrow(A)
    if (is.null(covariate))
        pval1 = sapply(1:ng, function(g) try(summary(manova(t(A[g,
            , ]) ~ y))$stats[1, "Pr(>F)"], silent = T))
    else pval1 = sapply(1:ng, function(g) try(summary(manova(t(A[g,
        , ]) ~ y + covariate))$stats[1, "Pr(>F)"], silent = T))
    pval = pval[, !is.na(as.numeric(pval1))]
    pval1 = na.omit(as.numeric(pval1))
    qval1 = p.adjust(pval1, "fdr")
    qval = pval
    K = ncol(A)
    for (i in 1:ncol(pval)) {
        qval[, i] = 1
        if (min(pval[, i], na.rm = T) < 0.05/K)
            qval[, i][which.min(pval[, i])] = qval1[i]
    }
    return(qval)
}


test <- function (A, y, covariate = NULL)
{
    if (dim(A)[3] != length(y))
        print("CSE estimates and y have different length")
    if (!is.null(covariate))
        if (dim(A)[3] != nrow(covariate))
            print("CSE estimates and covariate have different number of samples/subjects")
        else {
            if (!is.null(rownames(covariate)) & any(rownames(covariate) !=
                dimnames(A)[[3]]))
                covariate = covariate[dimnames(A)[[3]], ]
        }
    K = ncol(A)
    cell_type = colnames(A)
    if (is.null(covariate))
        pval = apply(A, 1, function(x) {
            pval = coef(summary(glm(y ~ ., data = data.frame(t(x)),
                family = "binomial")))[, 4]
            return(get_pval(pval, cell_type, K))
        })
    else pval = apply(A, 1, function(x) {
        pval = coef(summary(glm(y ~ ., data = data.frame(t(x),
            covariate), family = "binomial")))[, 4]
        return(get_pval(pval, cell_type, K))
    })
    qval = pval2qval(pval, A, y, covariate)
    return(list(qval = qval, pval = pval))
}

DEG_test2 <- function(X_array,y,covariate=NULL){
    dims_ct <- dim(X_array)[3]
    cellName <- dimnames(X_array)[[3]]
    if(is.null(covariate)){
    pval <- qval <- NULL
    for(ct in 1:dims_ct){
	mat <- X_array[,,ct]
	p <- NULL
	for(i in 1:nrow(mat)){
	 p <- c(p, summary(lm(mat[i,]~as.numeric(y)))$coefficients[2,4])
	}
	q = p.adjust(p, "fdr")
    pval <- cbind(pval,p)
    qval <- cbind(qval,q)
	}
    colnames(qval) <- colnames(pval) <- cellName
    rownames(qval) <- rownames(pval) <- dimnames(X_array)[[1]]

    }else{
        pval <- qval <- NULL
        cvname <- colnames(covariate)
        dims_ct <- dim(X_array)[3]
            for(ct in 1:dims_ct){
            mat <- X_array[,,ct]
            p <- NULL
            for(i in 1:nrow(mat)){
                dat <- cbind(mat[i,],y,covariate)
                dat <- as.data.frame(dat)
                colnames(dat) <- c("x","y",cvname)
                p <- c(p, summary(lm(x~.,dat=dat))$coefficients[2,4])
            }
            q = p.adjust(p, "fdr")
            pval <- cbind(pval,p)
            qval <- cbind(qval,q)
	    }
    colnames(qval) <- colnames(pval) <- cellName
    rownames(qval) <- rownames(pval) <- dimnames(X_array)[[1]]    
    }

    return(list(qval = qval, pval = pval))
}