#'# Class definitions
#' @importFrom methods setClassUnion
#' @importClassesFrom Matrix dgCMatrix
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))
setClassUnion(name = 'AnyFactor', members = c("factor", "list"))

#' The ENIGMA Class
#'
#' @slot raw_input
#' raw bulk and reference data
#'
#' @slot ref_type
#' the reference type
#'
#' @slot bulk
#' bulk matrix after batch correction
#'
#' @slot ref
#' reference matrix after batch correction
#'
#' @slot result_cell_proportion
#' the proportion of each cell types in bulk samples
#'
#' @slot result_CSE
#' inferred cell type-specific expression (CSE)
#'
#' @exportClass ENIGMA
#' @importFrom methods setClass
#'
ENIGMA = methods::setClass(
    "ENIGMA",
    slots = c(
        raw_input = "list",
        ref_type = "character",
        bulk = "AnyMatrix",
        ref = "AnyMatrix",
        result_cell_proportion = "AnyMatrix",
        result_CSE = "SingleCellExperiment"
    )
)

#' Show method for ENIGMA
#'
setMethod(f = "show", signature = "ENIGMA", definition = function(object) {
    cat("An object of class", class(object))
    invisible(x = NULL)
})

#' Create a new ENIGMA object
#'
#' @param bulk
#' raw bulk data matrix
#'
#' @param ref
#' raw reference data, could be single cell (matrix, SingleCellExperiment object or Seurat object) or FACS bulk
#'
#' @param ref_type
#' Determine the reference type. Should be either "single_cell" or "bulk".
#'
#' @param meta_bulk
#' metadata of bulk matrix
#'
#' @param meta_ref
#' metadata of reference
#'
#' @return an ENIGMA initial object
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#'
create_ENIGMA <- function(bulk, ref, ref_type=c("single_cell", "bulk"), meta_bulk=NULL, meta_ref=NULL, assay=NULL) {
    if ( !inherits(bulk, what = c("matrix", "Matrix", "dgCMatrix")) ) {
        stop("Bulk should be a matrix. ")
    }
    if ( !(ref_type %in% c("single_cell", "bulk")) | (length(ref_type) != 1) ) {
        stop("Invalid reference type. Please input 'single_cell' or 'bulk'. ")
    }

    ## reference is single cell RNA-seq
    if (ref_type == "single_cell") {
        message(date(), " Reference from Single Cell RNA-seq. ")
        ## input is matrix
        if ( inherits(ref, what = c("matrix", "Matrix", "dgCMatrix")) ) {
            message(date(), " Obtain reference from a matrix")
            ref_matrix = ref
            if (!is.null(meta_ref)) meta_ref = DataFrame(meta_ref)
        }

        ## input is SingleCellExperiment object
        if (is(ref, "SingleCellExperiment")) {
            message(date(), " Obtain reference from a SingleCellExperiment object")
            if ("counts" %in% SummarizedExperiment::assayNames(ref)) {
                cat("The `counts` assay is used",'\n')
                ref_matrix = SingleCellExperiment::counts(ref)
            } else {
                stop("SingleCellExperiment object must contain an assay named `counts`")
            }
            if (is.null(meta_ref)) {
                cat("The `colData` assay in the SingleCellExperiment object is used as metadata of reference",'\n')
                meta_ref = DataFrame(SingleCellExperiment::colData(ref))
            }
        }

        ## input is Seurat object
        if (is(ref, "Seurat")) {
            if (!requireNamespace("Seurat", quietly = TRUE)) {
                stop("Seurat installation required for working with Seurat objects")
            }
            message(date(), " Obtain reference from a Seurat object")
            if (is.null(assay)) {
                assay = Seurat::DefaultAssay(ref)
                if (assay == "integrated") {
                    warning("The data in the `integrated` assay is not suitable for ENIGMA analysis. Please use the `RNA` or `SCT` assay. ")
                }
                cat(paste0("The `data` slot in the default assay is used. The default assay is ", assay),'\n')
            }
            ref_matrix = Seurat::GetAssayData(ref, assay = assay, slot = "data") # normalized data matrix
            if (min(ref_matrix) < 0) {
                stop("The reference matrix contains negative values. Please ensure the normalized data matrix is used.")
            }
            if (is.null(meta_ref)) {
                cat("The `meta.data` slot in the Seurat object is used as metadata of reference",'\n')
                meta_ref = ref@meta.data
                meta_ref$ident = Seurat::Idents(ref)
                meta_ref = DataFrame(meta_ref)
            }
        }
    }

    ## reference is FACS Bulk RNA-seq/microarray
    if (ref_type == "bulk") {
        message(date(), " Reference from FACS Bulk RNA-seq/microarray. ")
        ref_matrix = ref
    }

    ## packaging bulk into SE object
    if (max(bulk) < 50) {
        warning("Bulk matrix seems to be in log space. Anti-log all values by 2^x. ")
        bulk = 2^bulk
    }
    bulk_se = SummarizedExperiment(
        assays = list( raw=bulk )
    )
    if (!is.null(meta_bulk)) meta_bulk = DataFrame(meta_bulk)
    if (!is.null(meta_bulk) & is(meta_bulk, "DataFrame")) {
        colData(bulk_se) = meta_bulk
    }

    ## packaging reference into SE object
    if (max(ref_matrix) < 50) {
        warning("Reference matrix seems to be in log space. Anti-log all values by 2^x. ")
        ref_matrix = 2^ref_matrix
    }
    ref_se = SummarizedExperiment(
        assays = list( raw=ref_matrix )
    )
    if (!is.null(meta_ref) & is(meta_ref, "DataFrame")) {
        colData(ref_se) = meta_ref
    }
    if (ref_type == "single_cell" & ncol(meta_ref)==0) {
        warning("Single cell reference has no metadata to determine cell type. ")
    }

    ## packaging raw input
    input_list = list(
        bulk=bulk_se,
        ref=ref_se
    )

    ## packaging all
    object = methods::new(
        Class = "ENIGMA",
        raw_input=input_list,
        ref_type=ref_type
    )
    return(object)
}
