#' @title Calculate the proportion of each cell types in bulk samples
#'
#' @param object ENIGMA object
#' 
#' @param max.iter
#' the maximum iteration times for robust linear model. Default = 500
#'
#' @param nu.v
#' parameter needed for nu-regression, user could input a vector of nu, and CIBERSORT would automatically select the nu with the best fitness with users' data. 
#'
#' @return the ENIGMA object contains the cell type fractions in object@result_cell_proportion
#'
#' 
#' @examples
#' \dontrun{
#' egm = get_cell_proportion(egm)
#' egm@result_cell_proportion
#' }
#' @export
#'
#'
#' @reference
#' Teschendorff A E, Breeze C E, Zheng S C, et al. A comparison of reference-based algorithms for correcting cell-type heterogeneity in Epigenome-Wide Association Studies[J]. BMC bioinformatics, 2017, 18(1): 1-14.
get_cell_proportion <- function(object,max.iter = 500,method = "RLR",nu.v) {
    if ( !(method %in% c("RLR", "CBS")) | (length(method) != 1) ) {
        stop("Invalid method type. Please input 'RLR','CBS'. ")
    }
    cat( date(), "Calculating cell type proportion of bulk samples... \n" )
	if(method == "RLR"){
	cat( "Using Robust Linear Regression... \n")
    object@result_cell_proportion = get_proportion(X = object@bulk, ref = object@ref, max.iter = max.iter)$theta
	}
	if(method == "CBS"){
	cat( "Using CIBERSORT... \n")
	object@result_cell_proportion = DoCBS(beta.m = object@bulk, ref.m = object@ref, nu.v = nu.v)$estF
	}
    return(object)
}

#' @importFrom MASS rlm
#'
get_proportion <- function(X, ref,max.iter = 500) {
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

        rlm.o <- rlm(Exp ~ as.matrix(ref_m), maxit = max.iter)
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


#' @importFrom e1071 svm
### the CBS code is provided by https://github.com/sjczheng/EpiDISH/blob/master/R/epidish.R
### please cite Teschendorff A E, Breeze C E, Zheng S C, et al. A comparison of reference-based algorithms for correcting cell-type heterogeneity in Epigenome-Wide Association Studies[J]. BMC bioinformatics, 2017, 18(1): 1-14.
DoCBS <- function(beta.m, ref.m, nu.v) {
    map.idx <- match(rownames(ref.m), rownames(beta.m))
    rep.idx <- which(is.na(map.idx) == FALSE)
    
    data2.m <- beta.m[map.idx[rep.idx], ]
    ref2.m <- ref.m[rep.idx, ]
    
    est.lm <- list()
    nui <- 1
    for (nu in nu.v) {
        est.m <- matrix(nrow = ncol(data2.m), ncol = ncol(ref2.m))
        colnames(est.m) <- colnames(ref2.m)
        rownames(est.m) <- colnames(data2.m)
        for (s in seq_len(ncol(data2.m))) {
            svm.o <- svm(x = ref2.m, y = data2.m[, s], scale = TRUE, type = "nu-regression", 
                kernel = "linear", nu = nu)
            coef.v <- t(svm.o$coefs) %*% svm.o$SV
            coef.v[which(coef.v < 0)] <- 0
            total <- sum(coef.v)
            coef.v <- coef.v/total
            est.m[s, ] <- coef.v
        }
        est.lm[[nui]] <- est.m
        nui <- nui + 1
    }
    
    #### select best nu
    rmse.m <- matrix(NA, nrow = ncol(beta.m), ncol = length(nu.v))
    for (nui in seq_along(nu.v)) {
        reconst.m <- ref2.m %*% t(est.lm[[nui]])
        s <- seq_len(ncol(beta.m))
        rmse.m[s, nui] <- sqrt(colMeans((data2.m[, s] - reconst.m[, s])^2))
        message(nui)
    }
    colnames(rmse.m) <- nu.v
    nu.idx <- apply(rmse.m, 1, which.min)
    estF.m <- est.m
    for (s in seq_len(nrow(estF.m))) {
        estF.m[s, ] <- est.lm[[nu.idx[s]]][s, ]
    }
    return(list(estF = estF.m, nu = nu.v[nu.idx], ref = ref2.m, dataREF = data2.m))
}


