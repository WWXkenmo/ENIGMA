#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#' @importFrom sva ComBat
#'
remove_batch_effect_S_mode <- function(object, varname_cell_type, n_pseudo_bulk) {
    bulk_matrix = assay( object@raw_input$bulk, "raw" )
    ref_se = object@raw_input$ref

    cat(date(), "Generating pseudo bulk... \n")
    pseudo_bulk_main_ct = generate_pseudo_bulk_from_scRNA( ref_se, varname_cell_type, n = n_pseudo_bulk )
    colnames(pseudo_bulk_main_ct$mat) = paste(varname_cell_type, seq_len(ncol(pseudo_bulk_main_ct$mat)), sep = "_")

    cat(date(), "Doing ComBat... \n")
    M_list = list( bulk_matrix, pseudo_bulk_main_ct$mat )
    gene_id = M_list %>% lapply(rownames) %>% Reduce(intersect, .)
    M_mix = M_list %>% lapply(function(x) x[gene_id,] ) %>% do.call(cbind, .)
    M_mix_log2 = log2(M_mix+1)

    M_correct = ComBat(
        M_mix_log2[rowSums(M_mix_log2)>0,],
        rep( c("bulk", "pseudo_bulk"), c( ncol(bulk_matrix), ncol(M_mix_log2)-ncol(bulk_matrix) ) )
    )
    M_correct = 2^M_correct

    cat(date(), "Restore reference... \n")
    object@bulk = M_correct[,colnames(bulk_matrix)]
    object@ref = restore_ref( pseudo_bulk_main_ct$frac, M_correct[,colnames(pseudo_bulk_main_ct$mat)])
    colnames(object@ref) = colnames(pseudo_bulk_main_ct$frac)
    return(object)
    cat(date(), "done. \n")
}

# generate pseudo bulk from scRNA
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#' @importFrom purrr map2_df
generate_pseudo_bulk_from_scRNA <- function(ref_se, ct_varname, n) {
    ref_matrix = assay( ref_se, "raw" )
    meta_ref = colData( ref_se )

    frac_init = table( meta_ref[,ct_varname] )/nrow(meta_ref)
    frac = map2_df(frac_init, frac_init*2, rnorm, n=n) %>% sapply(function(x) x)

    # set nagative equal to zero
    frac[frac<0] <- 0
    frac <- frac[rowSums(frac)>0,]
    # normalization
    frac <- t(t(frac) %*% diag(1/rowSums(frac)))

    M_star <- NULL
    for(i in 1:nrow(frac)){
        Ma = lapply( meta_ref[,ct_varname] %>% unique, function(ct) {
            sample_ids = subset( meta_ref, eval(parse( text = paste0(ct_varname, "==\"", ct, "\"") )) ) %>%
                rownames %>%
                sample( 1000*frac[i,ct], replace = TRUE )
            ref_matrix[,colnames(ref_matrix) %in% sample_ids]
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
