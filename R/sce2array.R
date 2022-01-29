#' @title Convert the output of ENIGMA CSE profiles into array format
#'
#' @description Convert the ENIGMA CSE profile (SingleCellExperiments object) into array format (genes * samples * cell types)
#' 
#' @param object
#' the ENIGMA object
#' 
#' @param model_name
#' name of the model, if null, the object@result_CSE or object@result_CSE_normalized would be used
#'
#'
#' @param norm_output
#' Output the normalized CSE profiles
#'
#' @return an array with three dimensions.
#'
#'
#' @examples
#' \dontrun{
#' egm = ENIGMA_l2max_norm(egm,model_tracker = TRUE, model_name = "Ken")
#' CSE = sce2array(egm, model_name = "Ken")
#' CSE = sce2array(egm)
#' }
#'
#'
#' @export
sce2array <- function(object, model_name = NULL,norm_output = TRUE){
   require("SingleCellExperiment")
   if(!norm_output){
   if(is.null(model_name)){
     Exp = counts(object@result_CSE)
	 CellLabel = colData(object@result_CSE)$cell_type
	}else{
	 Exp = counts(object@model[[model_name]]$result_CSE)
	 CellLabel = colData(object@model[[model_name]]$result_CSE)$cell_type
	}
	}else{
	if(is.null(model_name)){
     Exp = counts(object@result_CSE_normalized)
	 CellLabel = colData(object@result_CSE_normalized)$cell_type
	}else{
	 Exp = counts(object@model[[model_name]]$result_CSE_normalized)
	 CellLabel = colData(object@model[[model_name]]$result_CSE_normalized)$cell_type
	}
	}
	 SampleLabel = colnames(object@bulk)
	 array = matrix(0,nrow = nrow(Exp), ncol= length(SampleLabel))
	 array = array(0,
              dim = c( nrow(Exp),
                       length(SampleLabel),
                       length(table(CellLabel))),
              dimnames = list( rownames(Exp),
                               SampleLabel,
                               names(table(CellLabel)))
     )
	 for(ct in names(table(CellLabel))){
	    array[,,ct] = Exp[,CellLabel %in% ct]
	 }
	 array = array[,,colnames(object@ref)]
	 array
}
   
