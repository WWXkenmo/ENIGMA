#' @title Calculate the proportion of each cell types in bulk samples
#'
#' @param X
#' Bulk gene expression dataset (gene*n) in linear space (non-log).
#'
#' @param ref
#' Cell-type specific gene expression profile (GEP, gene*prop) in linear space (non-log).
#'
#' @return A matrix with proportions of each cell types in bulk samples.
#' @importFrom MASS rlm
#' @export
#'
get_proportion <- function(X, ref) {
    cat( date(), "Calculating cell type proportion of bulk samples... \n" )
    gene_id = intersect( rownames(X), rownames(ref) )
    X_m = X[gene_id,]
    ref_m = ref[gene_id,]
    ref_m <- apply(ref_m,2,scale)
    theta <- NULL
    coefVec <- NULL
    for(i in 1:ncol(X_m)){
        Exp <- as.matrix(X_m[,i])
        rownames(Exp) <- rownames(ref_m)
        colnames(Exp) <- colnames(X_m)[i]
        Exp <- scale(Exp)

        rlm.o <- rlm(Exp ~ as.matrix(ref_m), maxit = 100)
        coef.v <- summary(rlm.o)$coef[2:(ncol(as.matrix(ref_m)) + 1), 1]
        coefVec <- rbind(coefVec,coef.v)
        coef.v[which(coef.v < 0)] <- 0
        total <- sum(coef.v)
        coef.v <- coef.v/total
        theta <- rbind(theta,coef.v)
    }
    colnames(theta) <- colnames(coefVec) <- colnames(ref_m)
    rownames(theta) <- rownames(coefVec) <- colnames(X_m)
    res <- list()
    res$theta <- theta
    res$coef <- coefVec
    return(res)
}
