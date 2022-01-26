#' @title Remove batch effect with bulk and reference
#' 
#' @description This function would perform cross-platform batch effects correction. if object@ref_type = "single_cell", it would use S-mode batch correction method, if object@ref_type = "sort" or "aggre", it would use B-mode batch correction method. We note that this step is not essential step in deconvolution.
#' 
#' @param object ENIGMA object
#' @param varname_cell_type
#' Varname in colData of reference to tell which column is cell type. Only avaliable when reference is single cell RNA-seq data.
#'
#' @param n_pseudo_bulk
#' Number of pseudo bulk samples generated from scRNA reference. Only avaliable when reference is single cell RNA-seq data. Default is 100.
#' 
#' @param sample_size
#' Number of sampled cell for reconstituting pseudo-bulk samples. Default is 1000.
#'
#' @return For B-mode correction, ENIGMA object where object@bulk contains the corrected bulk gene expression profile. For S-mode correction, ENIGMA object where object@ref contains the corrected reference profile matrix.
#'
#'
#' @examples
#' \dontrun{
#' #S-mode
#' egm = batch_correct(egm, varname_cell_type = "main_celltype", n_pseudo_bulk=1000)
#' #B-mode
#' egm = batch_correct(egm)
#' }
#'
#' @export
#'
batch_correct <- function(object, varname_cell_type, n_pseudo_bulk=100, sample_size = 1000,ncores=6) {
    if (object@ref_type == "single_cell") {
        message(date(), " Reference is from Single Cell RNA-seq, doing batch correction in S mode. ")
        object = remove_batch_effect_S_mode(object, varname_cell_type = varname_cell_type, n_pseudo_bulk=n_pseudo_bulk,sample_size = sample_size,ncores=ncores)
    }
    if (object@ref_type %in% c("sort","aggre")) {
        message(date(), " Reference is from FACS Bulk RNA-seq/microarray, doing batch correction in B mode. ")
        object = remove_batch_effect_B_mode(object)
    }

    return(object)
}

