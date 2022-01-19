#' @title plot multiple loss curves of ENIGMA object
#'
#' @param object ENIGMA object
#'
#' @param alpha
#'
#' @param rlgType
#' plot which type of model. "maximum L2 norm" or "trace norm"
#'
#' @param scale
#' if scale loss value. log or raw
#'
#' @param shape
#' if loss curve dot shape by the models
#'
#' @param name_prefix
#' the title of model comparison
#'
#' @param all
#' if plot the all object function in one figure
#'
#' @param plotTerm
#' plot which object function. The alternative values are 1,2,3,4(1: distance to bulk; 2: distance to reference; 3: regularization term; 4: total loss)
#'
#' @param moi
#' the model of interets to plot, input the name of models
#'
#'
#' @return the ggplot object contains the plot of loss curve
#'
#'
#' @examples
#' \dontrun{
#' plotMultiLossCurve(egm,rlgType = "trace_norm") # plot all saved trace norm models
#' plotMultiLossCurve(egm,rlgType = "trace_norm", scale="log")# plot all saved trace norm models and use log scale
#' }
#'
#' @export
#'
plotLossCurve <- function(object,rlgType = "trace_norm",scale = "log",shape = FALSE,name_prefix="method",all=TRUE,plotTerm = 1,moi=NULL){
	suppressPackageStartupMessages(require("ggplot2"))
    suppressPackageStartupMessages(require("cowplot"))
    if ( !(scale %in% c("log", "raw")) | (length(scale) != 1) ) {
        stop("Invalid reference type. Please input 'raw','scale'. ")
    }
	if ( !(rlgType %in% c("trace_norm", "maximum L2 norm")) | (length(rlgType) != 1) ) {
        stop("Invalid reference type. Please input 'trace_norm','maximum L2 norm'. ")
    }
    if(rlgType == "trace_norm"){
	    modelTab = object@model_name
	    modellist = object@model_name$`Model Name`[object@model_name$`Model Type` == "trace norm model"]
		if(!is.null(moi)){moi_filter = intersect(moi, modellist);
		if(is.null(mod_filter)){stop("ERROR:object didn't contains models belongs to moi list! ")}}else{moi_filter = modellist}
		
		##generate loss_history_list
		loss_history_list = list()
		for(mod in moi_filter){
		  loss_history_list[[mod]] <- object@model[[mod]]$loss_his
		}
		
        df_list <- list()
        for(l in 1:length(loss_history_list)){
            loss_history <- loss_history_list[[l]]
            loss_history_list[[l]] <- loss_history
            df <- data.frame(Iteration = rep(c(1:nrow(loss_history)),(ncol(loss_history)+1)),LossType = 
                                 c(rep("Distance to observed mixture",nrow(loss_history)),
                                   rep("Distance to reference",nrow(loss_history)),
                                   rep("Rank regularization",nrow(loss_history)),
                                   rep("Over All",nrow(loss_history))),Loss = c(as.numeric(as.matrix(loss_history)),rowSums(loss_history)))
            df_list[[l]] <- df
        }
        names(df_list) <- names(loss_history_list)
        name_vec <- NULL
        for(i in 1:length(loss_history_list)) name_vec <- c(name_vec,rep(names(loss_history_list)[i],nrow(loss_history_list[[i]])))
        
        df1 <- df2 <- df3 <- df4 <- NULL
        for(l in 1:length(loss_history_list)){
            df1 <- rbind(df1, subset(df_list[[l]],df_list[[l]]$LossType %in% "Distance to observed mixture"))
            df2 <- rbind(df2, subset(df_list[[l]],df_list[[l]]$LossType %in% "Distance to reference"))
            df3 <- rbind(df3, subset(df_list[[l]],df_list[[l]]$LossType %in% "Rank regularization"))
            df4 <- rbind(df4, subset(df_list[[l]],df_list[[l]]$LossType %in% "Over All"))
        }
        id <- c(colnames(df1),"Method")
        df1 <- as.data.frame(cbind(df1, name_vec))
        df2 <- as.data.frame(cbind(df2, name_vec))
        df3 <- as.data.frame(cbind(df3, name_vec))
        df4 <- as.data.frame(cbind(df4, name_vec))
        colnames(df1) <- colnames(df2) <- colnames(df3) <- colnames(df4) <- id
        
        
        if(scale == "log"){
            df1$Loss <- log(df1$Loss+1)
            df2$Loss <- log(df2$Loss+1)
            df3$Loss <- log(df3$Loss+1)
            df4$Loss <- log(df4$Loss+1)
        }
        
        if(!shape){
            g1 <- ggplot(data=df1, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Distance to observed mixture',colour = name_prefix)
            g2 <- ggplot(data=df2, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Distance to reference',colour = name_prefix)
            g3 <- ggplot(data=df3, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Rank regularization',colour = name_prefix)
            g4 <- ggplot(data=df4, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Overall',colour = name_prefix)
        }else{
            g1 <- ggplot(data=df1, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Distance to observed mixture',colour = name_prefix)
            g2 <- ggplot(data=df2, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Distance to reference',colour = name_prefix)
            g3 <- ggplot(data=df3, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Rank regularization',colour = name_prefix)
            g4 <- ggplot(data=df4, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Overall',colour = name_prefix)
        }
        if(all){
            if(scale == "raw"){p1 <- plot_grid(g1+labs(y = 'Loss'),g2+labs(y = 'Loss'),g3+labs(y = 'Loss'),g4+labs(y = 'Loss'),nrow=2);print(p1)}
            if(scale == "log"){p1 <- plot_grid(g1,g2,g3,g4,nrow=2);print(p1)}
        }else{
            if(plotTerm == 1){
                if(scale == "raw") print(g1+labs(y = 'Loss'))
                if(scale == "log") print(g1)}
            if(plotTerm == 2){
                if(scale == "raw") print(g2+labs(y = 'Loss'))
                if(scale == "log") print(g2)}
            if(plotTerm == 3){
                if(scale == "raw") print(g3+labs(y = 'Loss'))
                if(scale == "log") print(g3)}
            if(plotTerm == 4){
                if(scale == "raw") print(g4+labs(y = 'Loss'))
                if(scale == "log") print(g4)}
        }
    }else if(rlgType == "maximum L2 norm"){
	    modelTab = object@model_name
	    modellist = object@model_name$`Model Name`[object@model_name$`Model Type` == "maximum L2 norm model"]
		if(!is.null(moi)){moi_filter = intersect(moi, modellist);
		if(is.null(mod_filter)){stop("ERROR:object didn't contains models belongs to moi list! ")}}else{moi_filter = modellist}
		
		##generate loss_history_list
		loss_history_list = list()
		for(mod in moi_filter){
		  loss_history_list[[mod]] <- object@model[[mod]]$loss_his
		}
		
        df_list <- list()
        for(l in 1:length(loss_history_list)){
            loss_history <- loss_history_list[[l]]
            loss_history_list[[l]] <- loss_history
            df <- data.frame(Iteration = rep(c(1:nrow(loss_history)),(ncol(loss_history)+1)),LossType = 
                                 c(rep("Distance to observed mixture",nrow(loss_history)),
                                   rep("Distance to reference",nrow(loss_history)),
                                   rep("Maximum L2 norm regularization",nrow(loss_history)),
                                   rep("Over All",nrow(loss_history))),Loss = c(as.numeric(as.matrix(loss_history)),rowSums(loss_history)))
            df_list[[l]] <- df
        }
        names(df_list) <- names(loss_history_list)
        name_vec <- NULL
        for(i in 1:length(loss_history_list)) name_vec <- c(name_vec,rep(names(loss_history_list)[i],nrow(loss_history_list[[i]])))
        
        df1 <- df2 <- df3 <- df4 <- NULL
        for(l in 1:length(loss_history_list)){
            df1 <- rbind(df1, subset(df_list[[l]],df_list[[l]]$LossType %in% "Distance to observed mixture"))
            df2 <- rbind(df2, subset(df_list[[l]],df_list[[l]]$LossType %in% "Distance to reference"))
            df3 <- rbind(df3, subset(df_list[[l]],df_list[[l]]$LossType %in% "Maximum L2 norm regularization"))
            df4 <- rbind(df4, subset(df_list[[l]],df_list[[l]]$LossType %in% "Over All"))
        }
        id <- c(colnames(df1),"Method")
        df1 <- as.data.frame(cbind(df1, name_vec))
        df2 <- as.data.frame(cbind(df2, name_vec))
        df3 <- as.data.frame(cbind(df3, name_vec))
        df4 <- as.data.frame(cbind(df4, name_vec))
        colnames(df1) <- colnames(df2) <- colnames(df3) <- colnames(df4) <- id
        
        
        if(scale == "log"){
            df1$Loss <- log(df1$Loss+1)
            df2$Loss <- log(df2$Loss+1)
            df3$Loss <- log(df3$Loss+1)
            df4$Loss <- log(df4$Loss+1)
        }
        
        if(!shape){
            g1 <- ggplot(data=df1, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Distance to observed mixture',colour = name_prefix)
            g2 <- ggplot(data=df2, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Distance to reference',colour = name_prefix)
            g3 <- ggplot(data=df3, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Maximum L2 norm regularization',colour = name_prefix)
            g4 <- ggplot(data=df4, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Overall',colour = name_prefix)
        }else{
            g1 <- ggplot(data=df1, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Distance to observed mixture',colour = name_prefix)
            g2 <- ggplot(data=df2, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Distance to reference',colour = name_prefix)
            g3 <- ggplot(data=df3, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Maximum L2 norm regularization',colour = name_prefix)
            g4 <- ggplot(data=df4, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Overall',colour = name_prefix)
        }
        if(all){
            if(scale == "raw"){p1 <- plot_grid(g1+labs(y = 'Loss'),g2+labs(y = 'Loss'),g3+labs(y = 'Loss'),g4+labs(y = 'Loss'),nrow=2);print(p1)}
            if(scale == "log"){p1 <- plot_grid(g1,g2,g3,g4,nrow=2);print(p1)}
        }else{
            if(plotTerm == 1){
                if(scale == "raw") print(g1+labs(y = 'Loss'))
                if(scale == "log") print(g1)}
            if(plotTerm == 2){
                if(scale == "raw") print(g2+labs(y = 'Loss'))
                if(scale == "log") print(g2)}
            if(plotTerm == 3){
                if(scale == "raw") print(g3+labs(y = 'Loss'))
                if(scale == "log") print(g3)}
            if(plotTerm == 4){
                if(scale == "raw") print(g4+labs(y = 'Loss'))
                if(scale == "log") print(g4)}
        }
    }
}


