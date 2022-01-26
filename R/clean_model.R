#' @title Clean the trained models for ENIGMA objects
#'
#' @description clean the unwanted models in ENIGMA objects
#' 
#' @param object
#' ENIGMA object
#'
#' @param model_list
#' input the names (list) of models
#' 
#' @return ENIGMA object
#'
#'
#' @examples
#' \dontrun{
#' egm = clean_model(egm,model_list = c(1,2)) # preserve the 1st and 2st models
#' egm = clean_model(egm,model_list = c("Wang")) # preserve the model which named as "Wang"
#' egm = clean_model(egm,NULL) # clean all saved model
#' }
#'
#'
#' @export
clean_model <- function(object, model_list){
   if(is.numeric(model_list)){
     object@model_name = object@model_name[model_list,]
	 object@model = object@model[model_list]
   }else{
     object@model_name = subset(object@model_name, `Model Name` %in% model_list)
	 object@model = object@model[model_list]
   }
   gc()
   object   
}


