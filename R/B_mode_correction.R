#' @importFrom sva ComBat
#' @importFrom SummarizedExperiment assay
#' @importFrom sva ComBat
remove_batch_effect_B_mode <- function(object){
    X = object@bulk
    ref = object@ref

    frac <- get_proportion(X,ref)
    fra = frac$theta

    cat(date(), "Run B-mode to correct batch effect...\n")
    X_hat <- ref %*% t(fra)
    X_log <- log2(X+1)
    X_hat_log <- log2(X_hat+1)

    #
    cat(date(), "Doing ComBat... \n")
    gene_id = intersect(rownames(X_log), rownames(X_hat_log))
    X_mix_log2 = cbind(X_log[gene_id,], X_hat_log[gene_id,])

    correct <- ComBat(
        X_mix_log2[rowSums(X_mix_log2)>0,],
        batch=c( rep("1",ncol(X)), rep("2",nrow(fra)) )
    )
    X_correct <- 2^correct[,1:ncol(X)] - 1

    correct <- 2^correct[,(ncol(X)+1):ncol(correct)] - 1 

    cat(date(), "Done. ")
    object@ref = ref
    object@bulk = X_correct
    return(object)
}
