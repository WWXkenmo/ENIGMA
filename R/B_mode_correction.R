#' Title
#'
#' @param X
#' @param ref
#' @param fra
#'
#' @return
#' @importFrom sva ComBat
#' @export
#'
#' @examples
B_mode_batch_effect_remove <- function(X,ref,fra){
    cat("Run B-mode to correct batch effect...")
    X_hat <- ref %*% t(fra)
    X_log <- log2(X+1)
    X_hat_log <- log2(X_hat+1)

    #
    cat("\n do Combat...")
    correct <- ComBat(cbind(X_log,X_hat_log),batch=c(rep("1",ncol(X)),rep("2",nrow(fra))))
    X_correct <- 2^correct[,1:ncol(X)]-1

    correct <- 2^correct[,(ncol(X)+1):ncol(correct)]-1

    #
    #cat("/nRestore the gene specific expression...")
    ##Using nnls
    #gene_ref_exp <- NULL
    #for(i in 1:nrow(X_correct)){
    #    gene_exp <- rlm(fra, correct[i,])$coefficients
    #    gene_ref_exp <- rbind(gene_ref_exp,gene_exp)
    #}

    #rownames(gene_ref_exp) <- rownames(ref)
    #colnames(gene_ref_exp) <- colnames(ref)
    #gene_ref_exp[gene_ref_exp<0] <-0

    cat("\n Done")
    #res <- list()
    res <- X_correct
    res
}
