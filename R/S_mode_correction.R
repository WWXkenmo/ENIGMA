#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#' @importFrom sva ComBat
#'
remove_batch_effect_S_mode <- function(object, varname_cell_type, n_pseudo_bulk, sample_size = sample_size, ncores = ncores) {
    bulk_matrix = assay( object@raw_input$bulk, "raw" )
    ref_se = object@raw_input$ref

    cat(date(), "Generating pseudo bulk... \n")
    pseudo_bulk_main_ct = generate_pseudo_bulk_from_scRNA( ref_se, varname_cell_type, n = n_pseudo_bulk, sample_size = sample_size, ncores = ncores)
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
    M_correct = 2^M_correct - 1

    cat(date(), "Restore reference... \n")
    object@ref = restore_ref( pseudo_bulk_main_ct$frac, M_correct[,colnames(pseudo_bulk_main_ct$mat)])
    colnames(object@ref) = colnames(pseudo_bulk_main_ct$frac)
	object@bulk = assay( object@raw_input$bulk, "raw" )
    return(object)
    cat(date(), "done. \n")
}



# generate pseudo bulk from scRNA
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#' @importFrom purrr map2_df
generate_pseudo_bulk_from_scRNA <- function(ref_eset, ct_varname, n=1000,sample_size = 1000,ncores=ncores) {
    suppressPackageStartupMessages(require("doParallel"))
	suppressPackageStartupMessages(require("sva"))
	suppressPackageStartupMessages(require("nnls"))
	suppressPackageStartupMessages(require("tibble"))
	suppressPackageStartupMessages(require("purrr"))
	
	ref_matrix = assay( ref_eset, "raw" )
    ref_matrix = as.matrix(assay( ref_eset, "raw" ))
	meta_ref = as.data.frame(colData(ref_eset))
	ref_eset = ExpressionSet(ref_matrix)
	pData(ref_eset) = meta_ref
	
    frac_init = table( pData(ref_eset)[,ct_varname] )/nrow(pData(ref_eset))
    
    frac = map2_df(frac_init, frac_init*2, rnorm, n=n) %>% sapply(function(x) x)
    
    # set nagative equal to zero
    frac[frac<0] <- 0
    frac <- frac[rowSums(frac)>0,]
    # normalization
    frac <- t(t(frac) %*% diag(1/rowSums(frac)))
    
    colnames(frac) <- names(table( pData(ref_eset)[,ct_varname] ))
    
 bulk_per_sample = function(x,ct_varname,ref_eset){
 require("magrittr")
 require("Biobase")
         Ma = lapply( pData(ref_eset)[,ct_varname] %>% unique, function(ct) {
            sample_ids = subset( pData(ref_eset), eval(parse( text = paste0(ct_varname, "==\"", ct, "\"") )) ) %>%
                rownames %>%
                sample( sample_size*x[ct], replace = TRUE )
            exprs(ref_eset)[,colnames(ref_eset) %in% sample_ids]
        } ) %>%
            do.call(cbind, .)
   
   Ma = rowSums(Ma)
   Ma
 }
 ##
 writeLines( paste("using ",ncores," cores...",sep=""))
 cl = makeCluster(ncores)
    registerDoParallel(cl)
    getDoParWorkers()
    #envir %>% appendEnv(parent_envir)
    ns = nrow(frac)
 result = foreach(j = 1:ns, .errorhandling = 'pass') %dopar% {
 bulk_per_sample(frac[j,],ref_eset = ref_eset,ct_varname = ct_varname)
 }
 result = do.call(cbind, result)
    stopCluster(cl)
    return(list( frac=frac, mat=result ))
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


