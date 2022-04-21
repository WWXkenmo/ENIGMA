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

softmaxP <- function(p){
   expp <- exp(1-p)
   expp <- diag(1/rowSums(expp)) %*% expp
   expp
}

bulkRegTest <- function(Bulk,frac,nMC = 1000,p_threshold = 0.05){
  writeLines("Using Linear Regression To Infer Probability of Expression...")
  p_tab <- NULL
  df <- (nrow(frac)-ncol(frac))
  for(i in 1:nrow(Bulk)){
   exp <- Bulk[i,] 
   rlm.o <- lm(exp ~ as.matrix(frac)-1)
   p <- pt(summary(rlm.o)$coef[,3], df,lower.tail=FALSE)
   p_tab <- rbind(p_tab,p)
  }
  p_tab[is.nan(p_tab)] <- 1
  pvalue.m <- p_tab
  for(i in 1:ncol(pvalue.m)) pvalue.m[,i] <- p.adjust(pvalue.m[,i],method="BH")
  
  ### soft max transformation
  expp <- softmaxP(pvalue.m)
  rownames(expp) <- rownames(pvalue.m) <- rownames(Bulk)
  colnames(expp) <- colnames(pvalue.m) <- colnames(frac)
  
  score <- matrix(NA,nrow = nrow(expp),ncol = ncol(expp))
  rownames(score) <- rownames(expp)
  colnames(score) <- colnames(expp)
  
  writeLines("Bootstrapping...")
  for(i in 1:ncol(score)){
   ss <- expp[,i]/rowSums(expp[,-i])
   ss[is.nan(ss)] <- 0
   
   ### using bootstrapping to calculate pvalue
   tab <- list()
   nf <- nrow(p_tab)
   for(n in 1:nMC){
     ss_pe <- expp[,i]/(rowSums(expp[,-i])[sample(1:nf,nf,replace=FALSE)])
	 tab[[n]] <- ss_pe
   }
   tab <- t(matrix(unlist(tab), nrow=length(expp[,1])))
   
   vec <- NULL
   for(j in 1:nf){
      vec <- c(vec,sum(ss[j]<tab[,j]))
   }
   score[,i] <- vec/nMC
  }
  
  ####### refine the results and only assign the gene to the most significant cell type
  writeLines("Refining...")
  call <- matrix(0,nrow = nrow(expp),ncol = ncol(expp))
  rownames(call) <- rownames(expp)
  colnames(call) <- colnames(expp)
  
  for(ct in 1:ncol(score)){
   gene <- rownames(score)[score[,ct]<p_threshold]
   order <- apply(pvalue.m[gene,],1,which.min)
   gene <- gene[order == ct]
   call[gene,ct] <- 1
  }
  
  return(
   list(
    call = call,
	pval = score,
	pvalue.m = pvalue.m
	)
  )
}  

############################################################################################
source("/mnt/data1/weixu/HiDe/revised/ENIGMA.R")
library(MIND)
library(TCA)

testMatrixG <- testMatrixH1 <- testMatrixDEGList <- list()

for(nMC in 1:10){
SNR <- 4.8 # chose the signal to noise ratio level and generate expression profile for deconvolution
k <- 5 # number of cell types
ng <- 2000 # number of genes
p <- 100 # number of samples
ndiff <- 0.1*2000 # number of genes differentially expressed

H1 <- matrix(rnorm(5*ng), ncol=ng) ## generate the expression profile for control group
H2 <- H1 ## generate the expression profile for case group

# create differential expression for each cell type
DEG_list <- list()
seq <- 1:ncol(H2)
for(i in 1:nrow(H2)){
    DEG_id <- sample(1:ncol(H2),ndiff,replace=FALSE) ## randomly select the DEG
	add <- sample(c(SNR, -SNR),ndiff,replace=TRUE) ## randomly generate up-regulated or down-regulated genes
    H2[i,DEG_id] <- H2[i,DEG_id] + add ## add to the case group
    DEG_list[[i]] <- DEG_id * sign(add) ## generate the DEG list
    #seq <- seq[seq %in% DEG_id == FALSE]
}

# cell frequency matrix per sample
cc <- matrix(runif(p*k), ncol=k)
cc <- t(scale(t(cc), center=FALSE, scale=rowSums(cc)))
colnames(cc) <- paste('cellType', 1:ncol(cc), sep="")

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

##generate bulk GEP
G <- NULL
for(i in 1:50){
    G <- cbind(G, t(as.matrix(t(as.matrix(cc[i,])) %*% H1_array[,,i])))
}
for(i in 1:50){
    G <- cbind(G, t(as.matrix(t(as.matrix(cc[i+50,])) %*% H2_array[,,i])))
}
G <- G + t(matrix(rnorm(100*ng), ncol=ng))
# sample classes (2 groups)
y <- gl(2, p/2)

###estimate fractions
rownames(G) <- colnames(H1) <- colnames(H2) <- paste0("gene",c(1:nrow(G)))
rownames(H1) <- paste0("cellType",1:nrow(H1))
colnames(G) <- paste0("Sample",1:ncol(G))

testMatrixG[[nMC]] <- G
testMatrixH1[[nMC]] <- H1
testMatrixDEGList[[nMC]] <- DEG_list
}


res <- NULL
for(nround in 1:10){
testMatrixG[[nround]] -> G
testMatrixH1[[nround]] -> H1
testMatrixDEGList[[nround]] -> DEG_list

Fra_Simulate <- get_proportion(G, t(H1))
res_alg_all_simulate <- cell_deconvolve_trace(O = as.matrix(G),
                                              theta=Fra_Simulate$theta,
                                              R=t(H1),
                                              alpha=0.1,beta=1,solver="proximalpoint",
                                              verbose=FALSE,max.iter = 1000,pos=FALSE,Norm.method = "frac",pre.process="none")

#########################################################################
####Using CellType=Gene Specificity Test to enhance CTS-DE Analysis

degRes <- DEG_test(res_alg_all_simulate$X_k_norm,as.numeric(y)-1,covariate = NULL)
geneTest <- bulkRegTest(G,Fra_Simulate$theta,p_threshold=0.1)


### calculate FDR
pre <-  NULL
for(ct in 1:5){
mat <- res_alg_all_simulate$X_k_norm[,,ct]
fc <- NULL
for(i in 1:nrow(mat)){
fc <- c(fc, summary(lm(mat[i,]~as.numeric(y)))$coefficients[2,1])
}
lm.test <- degRes$qval[ct,]
lm.test[lm.test==1] <- 0.99
lm.test[lm.test<0.05&fc>0] <- 1
lm.test[lm.test<0.05&fc< -0] <- -1
lm.test[abs(lm.test)!=1] <- 0
pre <- c(pre, precision(lm.test,DEG_list[[ct]]))
}


pre_new <- NULL
for(ct in 1:5){
mat <- res_alg_all_simulate$X_k_norm[,,ct]
fc <- NULL
for(i in 1:nrow(mat)){
fc <- c(fc, summary(lm(mat[i,]~as.numeric(y)))$coefficients[2,1])
}
lm.test <- degRes$qval[ct,]
lm.test[lm.test==1] <- 0.99
lm.test[lm.test<0.05&fc>0&geneTest$call[,ct]==1] <- 1
lm.test[lm.test<0.05&fc< -0&geneTest$call[,ct]==1] <- -1
lm.test[abs(lm.test)!=1] <- 0
pre_new <- c(pre_new, precision(lm.test,DEG_list[[ct]]))
}
pre_new <- pre_new[is.nan(pre_new) == FALSE]
pre_true <- mean(pre_new)

res <- rbind(res,c(mean(pre),mean(pre_true)))
}
FDR_res <- data.frame(class = c(rep("raw",nrow(res)),rep("filtered",nrow(res))),
                      FDR = 1-as.numeric(res))
FDR_res$class <- factor(FDR_res$class,levels = rev(c("filtered","raw")))
png("result_show.png",res=300,height=1500,width=1300)
ggplot(FDR_res, aes(x=class, y=FDR)) + 
  geom_boxplot(size = 1,fill = "gray")+
  labs(x="", y = "FDR")+ theme_minimal() +  geom_hline(yintercept=0.05, linetype="dashed", color = "red",size=1) + mytheme
dev.off()

###########################################################################################
Fra_Simulate <- get_proportion(G, t(H1))
res_alg_all_simulate <- cell_deconvolve_trace(O = as.matrix(G),
                                              theta=Fra_Simulate$theta,
                                              R=t(H1),
                                              alpha=0.1,beta=1,solver="proximalpoint",
                                              verbose=FALSE,max.iter = 1000,pos=FALSE,Norm.method = "frac",pre.process="none")

degRes <- DEG_test(res_alg_all_simulate$X_k_norm,as.numeric(y)-1,covariate = NULL)
geneTest <- bulkRegTest(G,Fra_Simulate$theta,p_threshold=0.1)


### calculate FDR
pre <-  NULL
for(ct in 1:5){
mat <- res_alg_all_simulate$X_k_norm[,,ct]
fc <- NULL
for(i in 1:nrow(mat)){
fc <- c(fc, summary(lm(mat[i,]~as.numeric(y)))$coefficients[2,1])
}
lm.test <- degRes$qval[ct,]
lm.test[lm.test==1] <- 0.99
lm.test[lm.test<0.05&fc>0] <- 1
lm.test[lm.test<0.05&fc< -0] <- -1
lm.test[abs(lm.test)!=1] <- 0
pre <- c(pre, precision(lm.test,DEG_list[[ct]]))
}


pre_new <- NULL
for(ct in 1:5){
mat <- res_alg_all_simulate$X_k_norm[,,ct]
fc <- NULL
for(i in 1:nrow(mat)){
fc <- c(fc, summary(lm(mat[i,]~as.numeric(y)))$coefficients[2,1])
}
lm.test <- degRes$qval[ct,]
lm.test[lm.test==1] <- 0.99
lm.test[lm.test<0.05&fc>0&geneTest$call[,ct]==1] <- 1
lm.test[lm.test<0.05&fc< -0&geneTest$call[,ct]==1] <- -1
lm.test[abs(lm.test)!=1] <- 0
pre_new <- c(pre_new, precision(lm.test,DEG_list[[ct]]))
}
pre_new <- pre_new[is.nan(pre_new) == FALSE]
pre_true <- mean(pre_new)


lm.test.m <- fc.m <- NULL
for(ct in 1:5){
mat <- res_alg_all_simulate$X_k_norm[,,ct]
fc <- NULL
for(i in 1:nrow(mat)){
fc <- c(fc, summary(lm(mat[i,]~as.numeric(y)))$coefficients[2,1])
}
lm.test <- degRes$qval[ct,]
lm.test[lm.test==1] <- 0.99
lm.test.m <- cbind(lm.test.m,lm.test)
fc.m <- cbind(fc.m,fc)
}

pre_d <- NULL
for(nMC in 1:1000){
pre_new <- NULL
for(ct in 1:5){
lm.test <- lm.test.m[,ct]
fc <- fc.m[,ct]
sample_idx <- sample(1:nrow(geneTest$call),sum(geneTest$call[,ct]==1),replace=FALSE)
vec <- rep(0, nrow(geneTest$call))
vec[sample_idx] <- 1
lm.test[lm.test<0.05&fc>0&vec==1] <- 1
lm.test[lm.test<0.05&fc< -0&vec==1] <- -1
lm.test[abs(lm.test)!=1] <- 0
pre_new <- c(pre_new, precision(lm.test,DEG_list[[ct]]))
}
pre_new <- pre_new[is.nan(pre_new) == FALSE]
pre_d <- c(pre_d,mean(pre_new))
print(length(pre_d))
}

df <- data.frame(FDR_random = 1-pre_d)

p<-ggplot(df, aes(x=FDR_random)) +
  geom_density() +
  labs(x="FDR", y = "density")+ theme_classic() +  geom_vline(xintercept=1-pre_true, linetype="dashed", color = "red",size=1)

png("dp_show.png",res=300,height=1500,width=1500)
p + mytheme
dev.off()















