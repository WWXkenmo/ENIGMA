#' @title Remove batch effect with bulk and scRNA reference
#'
#' @param bulk_eset
#' An ExpressionSet object of bulk samples.
#'
#' @param ref_eset
#' An ExpressionSet object of scRNA reference.
#'
#' @param varname_main_ct
#' Varname in pData of \code{ref_eset} to tell which column is cell type.
#'
#' @param varnames_sub_ct
#' Default is NULL. If you want to execute sub cell type mode, set this to tell
#' which cell-type(s) is sub.
#'
#' @param n_pseudo_bulk
#' Number of pseudo bulk samples generated from scRNA reference. Default is 1000.
#'
#' @return A list with batch-removed expression matrix of bulk and reference(s)
#' @importFrom Biobase pData
#' @importFrom magrittr %>%
#' @importFrom sva ComBat
#' @export
#'
remove_batch_effect <- function(bulk_eset, ref_eset, varname_main_ct, varnames_sub_ct=NULL, n_pseudo_bulk=1000) {
    cat(date(), "generating pseudo bulk... \n")
    pseudo_bulk_main_ct = generate_pseudo_bulk_from_scRNA( ref_eset, varname_main_ct, n = n_pseudo_bulk )
    colnames(pseudo_bulk_main_ct$mat) = paste(varname_main_ct, seq_len(ncol(pseudo_bulk_main_ct$mat)), sep = "_")

    if(is.null(varnames_sub_ct)) pseudo_bulk_subs = NULL
    else {
        pseudo_bulk_subs = lapply( varnames_sub_ct, function(x) {
            sample_ids = subset( pData(ref_eset), eval(parse(text = paste0("!is.na(", x, ")") )) ) %>% rownames
            tmp = generate_pseudo_bulk_from_scRNA( ref_eset[,sample_ids], x, n = n_pseudo_bulk )
            colnames(tmp$mat) = paste(x, seq_len(ncol(tmp$mat)), sep = "_")
            return(tmp)
        } )
    }

    cat(date(), "do ComBat... \n")
    M_list = c(
        list( exprs(bulk_eset), pseudo_bulk_main_ct$mat ),
        lapply( pseudo_bulk_subs, function(x) x$mat )
    )
    gene_id = M_list %>% lapply(rownames) %>% Reduce(intersect, .)
    M_mix = M_list %>% lapply(function(x) x[gene_id,] ) %>% do.call(cbind, .)
    M_mix_log2 = log2(M_mix+1)

    M_correct = ComBat(
        M_mix_log2[rowSums(M_mix_log2)>0,],
        rep( c("bulk", "pseudo_bulk"), c( ncol(bulk_eset), ncol(M_mix_log2)-ncol(bulk_eset) ) )
    )
    M_correct = 2^M_correct - 1

    cat(date(), "restore reference... \n")
    res = c(
        list(
            M_correct[,colnames(bulk_eset)],
            restore_ref( pseudo_bulk_main_ct$frac, M_correct[,colnames(pseudo_bulk_main_ct$mat)] )
        ),
        lapply( pseudo_bulk_subs, function(x) {
            restore_ref( x$frac, M_correct[,colnames(x$mat)] )
        } )
    )
    names(res) = c( "bulk", varname_main_ct, varnames_sub_ct )
    colnames(res[[2]]) <- colnames(pseudo_bulk_main_ct$frac)
    return(res)
    cat(date(), "done. \n")
}

# generate pseudo bulk from scRNA
#' @importFrom Biobase pData
#' @importFrom Biobase exprs
#' @importFrom magrittr %>%
#' @importFrom purrr map2_df
generate_pseudo_bulk_from_scRNA <- function(ref_eset, ct_varname, n) {
    frac_init = table( pData(ref_eset)[,ct_varname] )/nrow(pData(ref_eset))

    frac = map2_df(frac_init, frac_init*2, rnorm, n=n) %>% sapply(function(x) x)

    # set nagative equal to zero
    frac[frac<0] <- 0
    frac <- frac[rowSums(frac)>0,]
    # normalization
    frac <- t(t(frac) %*% diag(1/rowSums(frac)))

    M_star <- NULL
    for(i in 1:nrow(frac)){
        Ma = lapply( pData(ref_eset)[,ct_varname] %>% unique, function(ct) {
            sample_ids = subset( pData(ref_eset), eval(parse( text = paste0(ct_varname, "==\"", ct, "\"") )) ) %>%
                rownames %>%
                sample( 1000*frac[i,ct], replace = TRUE )
            exprs(ref_eset)[,colnames(ref_eset) %in% sample_ids]
        } ) %>%
            do.call(cbind, .)

        M_star <- cbind(M_star,rowSums(Ma))
    }

    return(list( frac=frac, mat=M_star ))
}

# robust regression calculate adjusted reference
#' @importFrom nnls nnls
#'
restore_ref <- function(frac, pseudo_bulk) {
    ref_exp <- NULL
    for(i in 1:nrow(pseudo_bulk)){
        coef <- nnls(frac, pseudo_bulk[i,])

        coef.v = coef$x

        ref_exp <- rbind(ref_exp,coef.v)
    }
    rownames(ref_exp) <- rownames(pseudo_bulk)
    ref_exp[ref_exp<0] = 0

    return(ref_exp)
}
