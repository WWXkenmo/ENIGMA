library(MIND)
library(EpiDISH)
library(TCA)
library(MASS)
source("ENIGMA.R")

sensitivity <- function(lm.test,list){
  ### sensitivity(recall) is defined as the num(true positive predictions)/num(true) 
  pos <- list[list>0]
  neg <- abs(list[list<0])
  sensitivity <- (length(intersect(which(lm.test==1),pos))+length(intersect(which(lm.test== -1),neg)))/length(list)
  sensitivity
}

precision <- function(lm.test,list){
  ### precision(1-FDR) is defined as the num(true positive predictions)/num(positive prediction) 
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
	##calculate Sensitivity for each methods
	
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


###############################################
#####
#Visualize precision, sensitivity and Sensitivity
source("/mnt/data1/weixu/HiDe/revised/DEG_analysis_uile_function.R")
ES <- c(1.8,2.4,3,3.6,4.2,4.8)
for(es in 1:length(ES)){
##The simulated datasets are provided in https://github.com/WWXkenmo/ENIGMA/tree/main/ENIGMA_analysis/Data/cts_DEG

load(paste("/mnt/data1/weixu/HiDe/resCompare_bi_",ES[es],".Rdata",sep=""))

##Rerunning TCA, bMIND and ENIGMA
##Runing TCA
library(MASS)
AUPRC <- NULL
for(rep in 1:10){
G <- testMatrixG[[rep]]
H1 <- testMatrixH1[[rep]]
Fra_Simulate <- get_proportion(G, t(H1))

###Running TCA
library(TCA)
tca.mdl <- tca(X = G, W = Fra_Simulate$theta, C1 = NULL, C2 = NULL,
                parallel = TRUE,verbose=FALSE)
Z_hat_simulate <- tensor(X = (as.matrix(G)), tca.mdl,verbose=FALSE)

###Running bMIND
library(MIND)
sig <- t(H1)
deconv_simulate = MIND::bMIND(bulk=G, profile = sig, ncore = 3,frac=Fra_Simulate$theta)

###Running ENIGMA
ENIGMA_trace <- cell_deconvolve_trace(O = as.matrix(G),
                                              theta=Fra_Simulate$theta,
                                              R=t(H1),
                                              alpha=0.1,beta=1,solver="admm",
                                              verbose=FALSE,max.iter = 1000,Normalize=FALSE,pos=FALSE)
ENIGMA_l2max <- cell_deconvolve(X=as.matrix(G),
                                        theta=Fra_Simulate$theta,
                                        R=t(H1),
                                        epsilon=0.001,
                                        alpha=0.1,
                                        beta=0.5,tao_k=0.01,max.iter=1000,verbose=FALSE,pos=FALSE)
p <- 100
y <- gl(2, p/2)
res_bmind <- DEG_test1(deconv_simulate$A,y,Fra_Simulate$theta,method = "bmind",10000,DEG_list_all[[rep]],qval=TRUE)
res_tca <- DEG_test1(Z_hat_simulate,y,Fra_Simulate$theta,method = "TCA",10000,DEG_list_all[[rep]],qval=TRUE)
res_enigma <- DEG_test1(ENIGMA_trace$X_k,y,Fra_Simulate$theta,method = "enigma",10000,DEG_list_all[[rep]],qval=TRUE)
res_enigma_l2max <- DEG_test1(ENIGMA_l2max$X_k,y,Fra_Simulate$theta,method = "enigma",10000,DEG_list_all[[rep]],qval=TRUE)
line <- c(res_bmind,res_tca,res_enigma,res_enigma_l2max)
AUPRC <- rbind(AUPRC,c(res_enigma,res_enigma_l2max))
}

save(AUPRC,file=paste("/mnt/data1/weixu/HiDe/revised/DEG/",ES[es],"_AUPRC_enigma.Rdata",sep=""))
}


##########################################################
###plot the AUPRC
library(ggplot2)

AUPRC_dat <- NULL
ES <- c(1.8,2.4,3,3.6,4.2,4.8)
for(es in ES){
file1 <- paste(es,"_AUPRC.Rdata",sep="")
load(file1)
AUPRC_all <- AUPRC;
rownames(AUPRC_all) <- paste("Sample-",1:nrow(AUPRC_all),sep="")
colnames(AUPRC_all) <- c("bMIND","TCA","ENIGMA (trace norm)","ENIGMA") 

mat <- cbind(as.numeric(AUPRC_all),c(rep("bMIND",10),rep("TCA",10),rep("ENIGMA (trace norm)",10),rep("ENIGMA",10)),rep(es,4*10))
AUPRC_dat <- rbind(AUPRC_dat,mat)
}
colnames(AUPRC_dat) <- c("AUPRC","Method","SNR")
AUPRC_dat <- as.data.frame(AUPRC_dat)
AUPRC_dat$AUPRC <- as.numeric(as.matrix(AUPRC_dat$AUPRC))
AUPRC_dat <- AUPRC_dat[is.nan(AUPRC_dat[,1])==FALSE,]
AUPRC_dat$Method <- as.character(AUPRC_dat$Method)
AUPRC_dat$Method <- as.factor(AUPRC_dat$Method)
AUPRC_dat$Method <- factor(AUPRC_dat$Method, levels = c("bMIND","TCA","ENIGMA","ENIGMA (trace norm)"))
p_boxplot_AUPRC_all <- AUPRC_dat %>% 
    mutate(SNR=paste0("SNR=", SNR)) %>% 
    #Rmisc::summarySE(measurevar = "AUPRC", groupvars = c("SNR", "Method")) %>% 
    ggplot(aes(Method, AUPRC, color=Method)) + 
    geom_boxplot() + 
    facet_grid(~SNR) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
    theme(legend.title = element_text(size = 12, color = "black", family = "Arial"), legend.text = element_text(size = 12, color = "black", family = "Arial")) + 
    theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))

png("AUPRC.png",res=300,height=1500,width=3000)
p_boxplot_AUPRC_all
dev.off()

##############################################################
###Look at the FDR of the methods
##Compare the ENIGMA with bMIND and TCA on FDR
for(es in 1:length(ES)){
load(paste("/Path/to/Data/resCompare_bi_",ES[es],".Rdata",sep=""))
library(MASS)
Precision <- Sensitivity <- NULL
for(rep in 1:10){
G <- testMatrixG[[rep]]
H1 <- testMatrixH1[[rep]]
Fra_Simulate <- get_proportion(G, t(H1))

###Running TCA
library(TCA)
tca.mdl <- tca(X = G, W = Fra_Simulate$theta, C1 = NULL, C2 = NULL,
                parallel = TRUE,verbose=FALSE)
Z_hat_simulate <- tensor(X = (as.matrix(G)), tca.mdl,verbose=FALSE)

###Running bMIND
library(MIND)
sig <- t(H1)
deconv_simulate = MIND::bMIND(bulk=G, profile = sig, ncore = 3,frac=Fra_Simulate$theta)

###Running ENIGMA
ENIGMA_trace <- cell_deconvolve_trace(O = as.matrix(G),
                                              theta=Fra_Simulate$theta,
                                              R=t(H1),
                                              alpha=0.1,beta=1,solver="admm",
                                              verbose=FALSE,max.iter = 1000,pos=FALSE,Norm.method = "frac",pre.process="none")
ENIGMA_l2max <- cell_deconvolve(X=as.matrix(G),
                                        theta=Fra_Simulate$theta,
                                        R=t(H1),
                                        epsilon=0.001,
                                        alpha=0.1,
                                        beta=0.5,tao_k=0.01,max.iter=1000,verbose=FALSE,pos=FALSE,Norm.method = "frac")
p <- 100
y <- gl(2, p/2)
Z_Tab <- array(NA,
                  dim = dim(deconv_simulate$A),
				  dimnames = dimnames(deconv_simulate$A))
for(i in 1:dim(Z_Tab)[2]){
 Z_Tab[,i,] <- Z_hat_simulate[[i]]
}
res <- DEG_test2(ENIGMA_trace$X_k,deconv_simulate$A,y,Fra_Simulate$theta,DEG_list_all[[rep]])
res2 <- DEG_test2(ENIGMA_trace$X_k,Z_Tab,y,Fra_Simulate$theta,DEG_list_all[[rep]])

line <- c(res[[2]][3],res2[[2]][3],res2[[1]][3],res2[[3]][3])
Precision <- rbind(Precision,line)
line <- c(res[[2]][1],res2[[2]][1],res2[[1]][1],res2[[3]][1])
Sensitivity <- rbind(Sensitivity,line)

save(Precision,file=paste("/mnt/data1/weixu/HiDe/revised/DEG/",ES[es],"_precision.Rdata",sep=""))
save(Sensitivity,file=paste("/mnt/data1/weixu/HiDe/revised/DEG/",ES[es],"_sensitivity.Rdata",sep=""))
}
}


#####
colnames(Sensitivity) = colnames(Precision) = c("bMIND","TCA","ENIGMA(TCA)","ENIGMA(bMIND)")


library(ggplot2)
Sensitivity_dat <- NULL
ES <- c(2.4,3,3.6,4.2,4.8)
for(es in ES){
file1 <- paste(es,"_sensitivity.Rdata",sep="")
load(file1)
Sensitivity_all <- Sensitivity;
rownames(Sensitivity_all) <- paste("Sample-",1:nrow(Sensitivity_all),sep="")
colnames(Sensitivity_all) <- c("bMIND","TCA","ENIGMA (trace norm)","ENIGMA (trace norm, top DEG)") 

mat <- cbind(as.numeric(Sensitivity_all),c(rep("bMIND",10),rep("TCA",10),rep("ENIGMA (trace norm)",10),rep("ENIGMA (trace norm, top DEG)",10)),rep(es,4*10))
Sensitivity_dat <- rbind(Sensitivity_dat,mat)
}
colnames(Sensitivity_dat) <- c("Sensitivity","Method","SNR")
Sensitivity_dat <- as.data.frame(Sensitivity_dat)
Sensitivity_dat$Sensitivity <- as.numeric(as.matrix(Sensitivity_dat$Sensitivity))
Sensitivity_dat <- Sensitivity_dat[is.nan(Sensitivity_dat[,1])==FALSE,]
Sensitivity_dat$Method <- as.character(Sensitivity_dat$Method)
Sensitivity_dat$Method <- as.factor(Sensitivity_dat$Method)
Sensitivity_dat$Method <- factor(Sensitivity_dat$Method, levels = c("bMIND","TCA","ENIGMA (trace norm)","ENIGMA (trace norm, top DEG)"))
p_boxplot_Sensitivity_all <- Sensitivity_dat %>% 
    mutate(SNR=paste0("SNR=", SNR)) %>% 
    #Rmisc::summarySE(measurevar = "Sensitivity", groupvars = c("SNR", "Method")) %>% 
    ggplot(aes(Method, Sensitivity, color=Method)) + 
    geom_boxplot() + 
    facet_grid(~SNR) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
    theme(legend.title = element_text(size = 12, color = "black", family = "Arial"), legend.text = element_text(size = 12, color = "black", family = "Arial")) + 
    theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))

png("Sensitivity.png",res=300,height=1500,width=3000)
p_boxplot_Sensitivity_all
dev.off()

##########################################################

Precision_dat <- NULL
ES <- c(2.4,3,3.6,4.2,4.8)
for(es in ES){
file1 <- paste(es,"_precision.Rdata",sep="")
load(file1)
Precision_all <- Precision;
rownames(Precision_all) <- paste("Sample-",1:nrow(Precision_all),sep="")
colnames(Precision_all) <- c("bMIND","TCA","ENIGMA (trace norm)","ENIGMA (trace norm, top DEG)") 

mat <- cbind(as.numeric(Precision_all),c(rep("bMIND",10),rep("TCA",10),rep("ENIGMA (trace norm)",10),rep("ENIGMA (trace norm, top DEG)",10)),rep(es,4*10))
Precision_dat <- rbind(Precision_dat,mat)
}
colnames(Precision_dat) <- c("Precision","Method","SNR")
Precision_dat <- as.data.frame(Precision_dat)
Precision_dat$Precision <- as.numeric(as.matrix(Precision_dat$Precision))
Precision_dat$Precision <- 1- Precision_dat$Precision
Precision_dat <- Precision_dat[is.nan(Precision_dat[,1])==FALSE,]
Precision_dat$Method <- as.character(Precision_dat$Method)
Precision_dat$Method <- as.factor(Precision_dat$Method)
Precision_dat$Method <- factor(Precision_dat$Method, levels = c("bMIND","TCA","ENIGMA (trace norm)","ENIGMA (trace norm, top DEG)"))
p_boxplot_Precision_all <- Precision_dat %>% 
    mutate(SNR=paste0("SNR=", SNR)) %>% 
    #Rmisc::summarySE(measurevar = "Precision", groupvars = c("SNR", "Method")) %>% 
    ggplot(aes(Method, Precision, color=Method)) + 
    geom_boxplot() + 
    facet_grid(~SNR) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
    theme(legend.title = element_text(size = 12, color = "black", family = "Arial"), legend.text = element_text(size = 12, color = "black", family = "Arial")) + 
    theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA)) + labs(y = "FDR")

png("FDR.png",res=300,height=1500,width=3000)
p_boxplot_Precision_all
dev.off()