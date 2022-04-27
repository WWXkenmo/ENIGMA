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

SLM <- function(O,y){
	ncelltype <- dim(O)[2]
	celltype_name <- dimnames(O)[[2]]
	gene_name <- dimnames(O)[[1]]
	pval_m <- qval_m <- NULL
	for(ct in 1:ncelltype){
	  P <- NULL
	  mat <- O[,ct,]
	  for(i in 1:nrow(mat)){
	    pval <- summary(lm(mat[i,]~as.numeric(y)))$coefficients[2,4]
	    P <- c(P,pval)
	  }
	  Q <- p.adjust(P,method="fdr")
	  pval_m <- rbind(pval_m, P)
	  qval_m <- rbind(qval_m, Q)
	}
	rownames(pval_m) <- rownames(qval_m) <- celltype_name
	colnames(pval_m) <- colnames(qval_m) <- gene_name
	res <- list()
	res$pval <- pval_m
	res$qval <- qval_m
	res
}

DEG_test <- function(X_array,y,covariate=NULL,FDR_control=TRUE){
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
	if(FDR_control) result <- MIND::test(O,y,covariate)
	if(!FDR_control) result <- SLM(O,y)
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

#To perform the analysis in main figure, simply set FDR_method = FALSE
DEG_test1 <- function(X_array,y,theta,method = NULL,n_sample,DEG_list,qval=FALSE,FDR_method=TRUE){
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
	result <- DEG_test(X_array,y,FDR=FDR_method)
	
	###Evaluation
	if(qval){Tab <- result$qval}else{Tab <- result$pval}
	perform <- DePRCcalculator(Tab,X_array,DEG_list,y,n_sample)$AUC_PR

	return(
	mean(perform)
	)
}
#####calculate two cases of ENIGMA
## After ANOVA test (improved precision)
#To perform the analysis in main figure, simply set FDR_method = FALSE
DEG_test2 <- function(ENIGMA_result,bMIND_result,y,theta,DEG_list,FDR_method = TRUE){
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
	result_enigma <- DEG_test(ENIGMA_result,y,FDR=FDR_method)
	result_mind <- DEG_test(bMIND_result,y,FDR=FDR_method)
	
	##Evaluate ENIGMA and bMIND
	
	##First, evaluate ENIGMA
	se <- sp <- pre <- num_gene <- NULL
	for(ct in 1:dim(bMIND_result)[3]){
	mat <- ENIGMA_result[,,ct]
	fc <- NULL
	for(i in 1:nrow(mat)){
	fc <- c(fc, summary(lm(mat[i,]~as.numeric(y)))$coefficients[2,1]) ## lm model is used for evaluating DE directions (upregulated or downregulated)
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
source("/Path/to/Data/DEG_analysis_uile_function.R")
ES <- c(1.8,2.4,3,3.6,4.2,4.8)
for(es in 1:length(ES)){
##The simulated datasets are provided in https://github.com/WWXkenmo/ENIGMA/tree/main/ENIGMA_analysis/Data/cts_DEG
load(paste("/Path/to/Data/resCompare_bi_",ES[es],".Rdata",sep=""))

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
                                              alpha=0.1,beta=1,solver="proximalpoint",
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
                                              alpha=0.1,beta=1,solver="proximalpoint",
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

# Perform results: bMIND, TCA, and ENIGMA
line <- c(res[[2]][3],res2[[2]][3],res2[[1]][3])
Precision <- rbind(Precision,line)
line <- c(res[[2]][1],res2[[2]][1],res2[[1]][1])
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

#####################################################################################################
###Randomly adapte 6 different SNR genes (num:150), Draw the plots for comparision
k <- 5 # number of cell types
ng <- 1000 # number of genes
p <- 100 # number of samples
ndiff <- 0.1*1000 # number of genes differentially expressed

H1 <- matrix(rnorm(5*ng), ncol=ng)
H2 <- H1
# create differential expression for 3rd cell type
DEG_list_all <- list()
SNR <- c(1.8,2.4,3,3.6,4.2,4.8)
for(snr in SNR){
DEG_list <- list()
for(i in 1:nrow(H2)){
    DEG_id <- sample(1:ncol(H2),6,replace=FALSE)
	add <- sample(c(snr,-1*snr),6,replace=TRUE)
    H2[i,DEG_id] <- H2[i,DEG_id] + add
    DEG_list[[i]] <- DEG_id * sign(add)
    #seq <- seq[seq %in% DEG_id == FALSE]
}
DEG_list_all[[as.character(snr)]] <- DEG_list
}



# cell frequency matrix per sample
cc <- matrix(runif(p*k), ncol=k)
cc <- t(scale(t(cc), center=FALSE, scale=rowSums(cc)))
colnames(cc) <- paste('cellType', 1:ncol(cc), sep="")

# sample classes (2 groups)
y <- gl(2, p/2)

###Add the noise to each cell type and generate a sample specific profile
H1_array <- array(0,
                  dim = c( nrow(H1),
                           ncol(H1),
                           (floor(100/2))))
for(i in 1:(floor(p/2))){
    H1_array[,,i] <- H1 + matrix(rnorm(5*ng), ncol=ng)
}
H2_array <- array(0,
                  dim = c( nrow(H2),
                           ncol(H2),
                           (floor(100/2))))
for(i in 1:(floor(p/2))){
    H2_array[,,i] <- H2 + matrix(rnorm(5*ng), ncol=ng)
}
##evaluate differential expression genes
##calculate 
G <- NULL
for(i in 1:50){
    G <- cbind(G, t(as.matrix(t(as.matrix(cc[i,])) %*% H1_array[,,i])))
}
for(i in 1:50){
    G <- cbind(G, t(as.matrix(t(as.matrix(cc[i+50,])) %*% H2_array[,,i])))
}
G <- G + t(matrix(rnorm(100*ng), ncol=ng))

#########################################################################
#Generate ground truth mean different effects
H1_array_mean <- NULL
for(i in 1:5){
  H1_array_mat <- H1_array[i,,]
  H1_array_mean <- cbind(H1_array_mean, rowMeans(H1_array_mat))
}

H2_array_mean <- NULL
for(i in 1:5){
  H2_array_mat <- H2_array[i,,]
  H2_array_mean <- cbind(H2_array_mean, rowMeans(H2_array_mat))
}

SNR <- c(1.8,2.4,3,3.6,4.2,4.8)
meanDiff <- NULL
for(snr in SNR){
  diff.all <- gene.idx.all <- ct.id.all <- NULL  
  for(ct in 1:5){ 
     diff <- H2_array_mean[abs(DEG_list_all[[as.character(snr)]][[ct]]),ct] - H1_array_mean[abs(DEG_list_all[[as.character(snr)]][[ct]]),ct]
     gene.idx <- abs(DEG_list_all[[as.character(snr)]][[ct]])
	 ct.id <- rep(paste("celltype-",ct,sep=""), length(diff))
	 diff.all <- c(diff.all,diff)
	 gene.idx.all <- c(gene.idx.all,gene.idx)
	 ct.id.all <- c(ct.id.all,ct.id)
  }
  tab <- cbind(diff.all,gene.idx.all,ct.id.all)
  tab <- cbind(tab,rep(snr,nrow(tab)))
  meanDiff <- rbind(meanDiff,tab)
}
meanDiff <- as.data.frame(meanDiff)
colnames(meanDiff) <- c("True expression diff.","gene.idx","celltype","SNR")

meanDiff$gene.idx <- as.integer(meanDiff$gene.idx)
##############
y <- gl(2, p/2)

###estimate fractions
rownames(G) <- colnames(H1) <- colnames(H2) <- paste0("gene",c(1:nrow(G)))
Fra_Simulate <- get_proportion(G, t(H1))

colnames(Fra_Simulate$theta) <- colnames(cc)
rownames(Fra_Simulate$theta) <- colnames(G) <- paste0("Sample",1:ncol(G))
bulk <- G
tca.mdl <- tca(X = bulk, W = Fra_Simulate$theta, C1 = NULL, C2 = NULL,
                parallel = TRUE)
Z_hat_simulate <- tensor(X = (as.matrix(bulk)), tca.mdl)


rownames(G) <- colnames(H1) <- colnames(H2) <- paste0("gene",c(1:nrow(G)))
rownames(H1) <- paste0("celltype",1:nrow(H1))
Fra_Simulate <- get_proportion(G, t(H1))
colnames(Fra_Simulate$theta) <- colnames(cc)
rownames(Fra_Simulate$theta) <- colnames(G) <- paste0("Sample",1:ncol(G))
res_alg_all_simulate <- cell_deconvolve_trace(O = as.matrix(G),
                                              theta=Fra_Simulate$theta,
                                              R=t(H1),
                                              alpha=0.1,beta=1,solver="proximalpoint",
                                              verbose=FALSE,max.iter = 1000,pos=FALSE,Norm.method = "frac",pre.process="none")

sig <- t(H1)
colnames(sig) <- paste0("cellType",1:5)
deconv_simulate = MIND::bMIND(bulk=G, profile = sig, ncore = 3,frac=Fra_Simulate$theta)

##############################
#Analysis
meanDiff <- rbind(subset(meanDiff, celltype %in% "celltype-1"),
                  subset(meanDiff, celltype %in% "celltype-2"),
				  subset(meanDiff, celltype %in% "celltype-3"),
				  subset(meanDiff, celltype %in% "celltype-4"),
				  subset(meanDiff, celltype %in% "celltype-5"))
index <- c(1:5);names(index) <- c("celltype-1","celltype-2","celltype-3","celltype-4","celltype-5")
y <- as.numeric(gl(2, p/2))

meanDiff$`True expression diff.` <- as.numeric(as.matrix(meanDiff$`True expression diff.`))

estDiff <- NULL
for(i in c("celltype-1","celltype-2","celltype-3","celltype-4","celltype-5")){
  geneList <- subset(meanDiff,celltype %in% i)$gene.idx
  estimateDiff <- rowMeans(G[geneList,y==2]) - rowMeans(G[geneList,y==1])
  estDiff <- c(estDiff,estimateDiff)
}

meanDiff$estDiff <- estDiff
print(paste("Correlation between True and Estimation ",cor(meanDiff$`True expression diff.`,meanDiff$estDiff)))



estDiff <- NULL
for(i in c("celltype-1","celltype-2","celltype-3","celltype-4","celltype-5")){
  geneList <- subset(meanDiff,celltype %in% i)$gene.idx
  estimateDiff <- rowMeans(Z_hat_simulate[[index[i]]][geneList,y==2]) - rowMeans(Z_hat_simulate[[index[i]]][geneList,y==1])
  estDiff <- c(estDiff,estimateDiff)
}

meanDiff$estDiff <- estDiff
print(paste("Correlation between True and Estimation ",cor(meanDiff$`True expression diff.`,meanDiff$estDiff)))


estDiff <- NULL
for(i in c("celltype-1","celltype-2","celltype-3","celltype-4","celltype-5")){
  geneList <- subset(meanDiff,celltype %in% i)$gene.idx
  estimateDiff <- rowMeans(deconv_simulate$A[geneList,index[i],y==2]) - rowMeans(deconv_simulate$A[geneList,index[i],y==1])
  estDiff <- c(estDiff,estimateDiff)
}

meanDiff$estDiff <- estDiff
print(paste("Correlation between True and Estimation ",cor(meanDiff$`True expression diff.`,meanDiff$estDiff)))


estDiff <- NULL
for(i in c("celltype-1","celltype-2","celltype-3","celltype-4","celltype-5")){
  geneList <- subset(meanDiff,celltype %in% i)$gene.idx
  estimateDiff <- rowMeans(res_alg_all_simulate$X_k_norm[geneList,y==2,index[i]]) - rowMeans(res_alg_all_simulate$X_k_norm[geneList,y==1,index[i]])
  estDiff <- c(estDiff,estimateDiff)
}

meanDiff$estDiff <- estDiff
print(paste("Correlation between True and Estimation ",cor(meanDiff$`True expression diff.`,meanDiff$estDiff)))

##########################################
mytheme <- readRDS("mytheme.rds")
png("scatter_plot.png",res=300,height=1500,width=1500)
ggplot(meanDiff, aes(x=estDiff, y=`True expression diff.`)) +
    geom_point(aes(color=celltype,shape=SNR))+
    mytheme + labs(y="True expression diff.",x= "Predicted expression diff")
dev.off()



















