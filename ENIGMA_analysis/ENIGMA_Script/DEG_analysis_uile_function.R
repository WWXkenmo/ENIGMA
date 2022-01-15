#! /usr/bin/env Rscript
library(MIND)
library(EpiDISH)
library(TCA)
library(MASS)


sensitivity <- function(lm.test,list){
  pos <- list[list>0]
  neg <- abs(list[list<0])
  sensitivity <- (length(intersect(which(lm.test==1),pos))+length(intersect(which(lm.test== -1),neg)))/length(list)
  sensitivity
}


specificity <- function(lm.test,list){
  sp <- sum(which(lm.test==0) %in% abs(list) == FALSE)/(500-25)
  sp
}

precision <- function(lm.test,list){
  pos <- list[list>0]
  neg <- abs(list[list<0])
  precision <- (length(intersect(which(lm.test==1),pos))+length(intersect(which(lm.test== -1),neg)))/sum(lm.test!=0)
  precision
}

DEG_test <- function(X_array,y,covariate=NULL){
    O = array(0,
                  dim = c( dim(X_array)[1],
                           dim(X_array)[3],
                           dim(X_array)[2]),
                  dimnames = list( dimnames(X_array)[[1]],
                                   dimnames(X_array)[[3]],
                                   dimnames(X_array)[[2]])
    )
	for(i in 1:dim(X_array)[3]){
	  O[,i,] <- X_array[,,i]
	}
	###Using ANOVA+glm model to evaluate prediction performance
	result <- MIND::test(O,y,covariate)
	result
}

DePRCcalculator <- function(qval,X_array,DEG_list,y,n_sample,reshape=NULL){
   # Plot PR-curve and Evaluate the area under the PR-curve
   require(pracma)
   
   if(!is.null(reshape)){
      if(reshape == "bmind"){
	    O = array(0,
                  dim = c( dim(X_array)[1],
                           dim(X_array)[3],
                           dim(X_array)[2]),
                  dimnames = list( dimnames(X_array)[[1]],
                                   dimnames(X_array)[[3]],
                                   dimnames(X_array)[[2]])
			)
			for(i in 1:dim(X_array)[2]){
			  O[,,i] <- X_array[,i,]
			}
		X_array <- O
	  }
	  if(reshape == "TCA"){
	     O = array(0,
                  dim = c( dim(X_array[[1]])[1],
                           dim(X_array[[1]])[2],
                           length(X_array))
			)
			for(i in 1:length(X_array)){
			  O[,,i] <- X_array[[i]]
		}
		X_array <- O
	  }
   }
   
   dims_ct <- dim(X_array)[3]
   
   AUC = NULL
   recallList <- NULL
   preList <- NULL
   for(ct in 1:dims_ct){
	mat <- X_array[,,ct]
	fc <- NULL
	for(i in 1:nrow(mat)){
	 fc <- c(fc, summary(lm(mat[i,]~as.numeric(y)))$coefficients[2,1])
	}
	qvalue <- qval[ct,]
	
	###select n_sample random points as threshold
	seq <- runif(n_sample)
	seq <- quantile(qvalue,seq)
	seq <- as.numeric(names(table(seq)))
	print(length(seq))
	
	recall <- NULL
	pre <- NULL
	for(thre in seq){
	lm.test <- qvalue
	lm.test[lm.test==1] <- 0.99
	lm.test[lm.test<thre&fc>0] <- 1
	lm.test[lm.test<thre&fc< -0] <- -1
	lm.test[abs(lm.test)!=1] <- 0
	###
	recall <- c(recall, sensitivity(lm.test,DEG_list[[ct]]))
	pre <- c(pre, precision(lm.test,DEG_list[[ct]]))
	}
   AUC = c(AUC, trapz(recall,pre))
   recallList[[ct]] <- recall
   preList[[ct]] <- pre
   }
   
   return(
   list(recall = recallList,
        precision = preList,
		AUC_PR = AUC)
   )
}

DEG_test1 <- function(X_array,y,theta,method = NULL,n_sample,DEG_list,qval=FALSE){
	##calculate AUPRC for each methods
	
	method_list <- c("enigma","bmind","TCA")
	msg <- paste0("should be one of ", paste(method_list, collapse = ", "), 
        ".")
    if (!method %in% method_list) {
        stop("'method' ", msg)
    }
	
	if(method == "bmind"){
		O = array(0,
				  dim = c( dim(X_array)[1],
						   dim(X_array)[3],
						   dim(X_array)[2]),
				  dimnames = list( dimnames(X_array)[[1]],
								   dimnames(X_array)[[3]],
								   dimnames(X_array)[[2]])
			)
			for(i in 1:dim(X_array)[2]){
			  O[,,i] <- X_array[,i,]
			}
		X_array <- O
	  }
	  if(method == "TCA"){
		 O = array(0,
				  dim = c( dim(X_array[[1]])[1],
						   dim(X_array[[1]])[2],
						   length(X_array)),
				  dimnames = list( rownames(X_array[[1]]),
								   colnames(X_array[[1]]),
								   paste0("celltype",1:length(X_array),sep="")
			))
			for(i in 1:length(X_array)){
			  O[,,i] <- X_array[[i]]
		}
		X_array <- O
	 }
	
	
	if(method == "enigma"){
		X_k_norm <- X_k <- X_array
		for(k in 1:dim(X_k)[3]){
		exp <- X_k[,,k]

		PC <- theta[,k]
		exp.norm <- NULL
		for(i in 1:nrow(exp)){
		  lm.model <- lm(exp[i,]~PC)
		  exp.norm <- rbind(exp.norm, (exp[i,] - lm.model$coefficients[2] * PC))
		}
		X_k_norm[,,k] <- exp.norm
		}
		X_array <- X_k_norm
	}
	result <- DEG_test(X_array,y)
	
	###Evaluation
	if(qval){Tab <- result$qval}else{Tab <- result$pval}
	perform <- DePRCcalculator(Tab,X_array,DEG_list,y,n_sample)$AUC_PR

	return(
	mean(perform)
	)
}
#####calculate two cases of ENIGMA
## After ANOVA test (improved precision)
## After gene number control (further improved precision)
## only for compare of ENIGMA and bMIND
## qvalue < 0.05

DEG_test2 <- function(ENIGMA_result,bMIND_result,y,theta,DEG_list){
	O = array(0,
			  dim = c( dim(bMIND_result)[1],
					   dim(bMIND_result)[3],
					   dim(bMIND_result)[2]),
			  dimnames = list( dimnames(bMIND_result)[[1]],
							   dimnames(bMIND_result)[[3]],
							   dimnames(bMIND_result)[[2]])
		)
		for(i in 1:dim(bMIND_result)[2]){
		  O[,,i] <- bMIND_result[,i,]
		}
	bMIND_result <- O;rm(O);gc()
	X_k_norm <- X_k <- ENIGMA_result
	for(k in 1:dim(X_k)[3]){
		exp <- X_k[,,k]

		PC <- theta[,k]
		exp.norm <- NULL
		for(i in 1:nrow(exp)){
		  lm.model <- lm(exp[i,]~PC)
		  exp.norm <- rbind(exp.norm, (exp[i,] - lm.model$coefficients[2] * PC))
		}
		X_k_norm[,,k] <- exp.norm
	}
	ENIGMA_result <- X_k_norm
	result_enigma <- DEG_test(ENIGMA_result,y)
	result_mind <- DEG_test(bMIND_result,y)
	
	##Evaluate ENIGMA and bMIND
	
	##First, evaluate ENIGMA
	se <- sp <- pre <- num_gene <- NULL
	for(ct in 1:dim(bMIND_result)[3]){
	mat <- ENIGMA_result[,,ct]
	fc <- NULL
	for(i in 1:nrow(mat)){
	fc <- c(fc, summary(lm(mat[i,]~as.numeric(y)))$coefficients[2,1])
	}
	lm.test <- result_enigma$qval[ct,]
	lm.test[lm.test==1] <- 0.99
	lm.test[lm.test<0.05&fc>0] <- 1
	lm.test[lm.test<0.05&fc< -0] <- -1
	lm.test[abs(lm.test)!=1] <- 0
	se <- c(se, sensitivity(lm.test,DEG_list[[ct]]))
	sp <- c(sp, sum(which(lm.test==0) %in% abs(DEG_list[[ct]]) == FALSE)/(1000-100) )
	pre <- c(pre, precision(lm.test,DEG_list[[ct]]))
	num_gene <- c(num_gene, sum(lm.test!=0))
	}
	res_enigma <- c(mean(se),mean(sp),mean(pre),mean(num_gene))
	##We could observed that ENIGMA has improved precision
	
	##Next, evaluate bMIND
	se <- sp <- pre <- num_gene <- NULL
	for(ct in 1:dim(bMIND_result)[3]){
	mat <- bMIND_result[,,ct]
	fc <- NULL
	for(i in 1:nrow(mat)){
	fc <- c(fc, summary(lm(mat[i,]~as.numeric(y)))$coefficients[2,1])
	}
	lm.test <- result_mind$qval[ct,]
	lm.test[lm.test==1] <- 0.99
	lm.test[lm.test<0.05&fc>0] <- 1
	lm.test[lm.test<0.05&fc< -0] <- -1
	lm.test[abs(lm.test)!=1] <- 0
	se <- c(se, sensitivity(lm.test,DEG_list[[ct]]))
	sp <- c(sp, sum(which(lm.test==0) %in% abs(DEG_list[[ct]]) == FALSE)/(1000-100) )
	pre <- c(pre, precision(lm.test,DEG_list[[ct]]))
	num_gene <- c(num_gene, sum(lm.test!=0))
	}
	res_bmind <- c(mean(se),mean(sp),mean(pre),mean(num_gene))
	##We could observed that even though using this method bMIND has very high precision, but meanwhile its sensitivity is significantly reduced, even worse that bulk
	
	##We note that ENIGMA has both higher sensitivity and precision than bulk, and also, if we control the significant genes equal to the bMIND, we also could find ENIGMA has pretty high precision values. 
	se <- sp <- pre <- num_gene <- NULL
	for(ct in 1:dim(bMIND_result)[3]){
	mat <- ENIGMA_result[,,ct]
	fc <- NULL
	for(i in 1:nrow(mat)){
	fc <- c(fc, summary(lm(mat[i,]~as.numeric(y)))$coefficients[2,1])
	}
	lm.test <- result_enigma$qval[ct,]
	#
	lm.test.c <- result_mind$qval[ct,]
	vec <- order(lm.test,decreasing=FALSE)
	vec <- vec[(sum(lm.test.c<0.05)+1):length(vec)]
	lm.test[vec] <- 1
	#####

	lm.test[lm.test==1] <- 0.99
	lm.test[lm.test<0.05&fc>0] <- 1
	lm.test[lm.test<0.05&fc< -0] <- -1
	lm.test[abs(lm.test)!=1] <- 0

	###
	se <- c(se, sensitivity(lm.test,DEG_list[[ct]]))
	sp <- c(sp, sum(which(lm.test==0) %in% abs(DEG_list[[ct]]) == FALSE)/(1000-100) )
	pre <- c(pre, precision(lm.test,DEG_list[[ct]]))
	num_gene <- c(num_gene, sum(lm.test!=0))
	}
	res_enigma2 <- c(mean(se),mean(sp),mean(pre),mean(num_gene))
	
	##return the results
	return(
	list(
	ENIGMA_res1 = res_enigma,
	bMIND_res = res_bmind,
	ENIGMA_res2 = res_enigma2
	))
}