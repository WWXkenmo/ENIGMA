##############################################################################################
##Simulation for latent cell states
###Using ESCO to simulate cell states
library(Seurat)
library(SingleCellExperiment)
library(scater)

source("ENIGMA.R")
# The scRNA-seq datasets could be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119807
PD50 <- read.table("GSM3384107_LowPD50Gy.dge.txt",header=TRUE,row.names=1)
LowPDCtrl <- read.table("GSM3384106_LowPDCtrl.dge.txt",header=TRUE,row.names=1)
HighPDCtrl <- read.table("GSM3384108_HighPDCtrl.dge.txt",header=TRUE,row.names=1)
Senes <- read.table("GSM3384109_senescence.dge.txt",header=TRUE,row.names=1)
	
PD50 <- CreateSeuratObject(counts = as.matrix(PD50), project = "mouse.brain", min.cells = 10)

PD50[["percent.mt"]] <- PercentageFeatureSet(PD50, pattern = "^mt-")
VlnPlot(PD50, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

quantile(PD50@meta.data$nCount_RNA,c(0.025,0.975))
quantile(PD50@meta.data$nFeature_RNA,c(0.025,0.975))
plot(PD50@meta.data$nCount_RNA,PD50@meta.data$nFeature_RNA,pch=16,cex=0.7,bty="n")
abline(h=c(quantile(PD50@meta.data$nFeature_RNA,c(0.025,0.975))[1],quantile(PD50@meta.data$nFeature_RNA,c(0.025,0.975))[2]),
       v=c(quantile(PD50@meta.data$nCount_RNA,c(0.025,0.975))[1],quantile(PD50@meta.data$nCount_RNA,c(0.025,0.975))[2]
	   ),lty=2,lwd=1,col="red")
PD50 <- subset(PD50, subset = nFeature_RNA > quantile(PD50@meta.data$nFeature_RNA,c(0.025,0.975))[1] & nFeature_RNA < quantile(PD50@meta.data$nFeature_RNA,c(0.025,0.975))[2] & percent.mt < 10)

####
LowPDCtrl <- CreateSeuratObject(counts = as.matrix(LowPDCtrl), project = "mouse.brain", min.cells = 10)

LowPDCtrl[["percent.mt"]] <- PercentageFeatureSet(LowPDCtrl, pattern = "^mt-")
VlnPlot(LowPDCtrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

quantile(LowPDCtrl@meta.data$nCount_RNA,c(0.025,0.975))
quantile(LowPDCtrl@meta.data$nFeature_RNA,c(0.025,0.975))
plot(LowPDCtrl@meta.data$nCount_RNA,LowPDCtrl@meta.data$nFeature_RNA,pch=16,cex=0.7,bty="n")
abline(h=c(quantile(LowPDCtrl@meta.data$nFeature_RNA,c(0.025,0.975))[1],quantile(LowPDCtrl@meta.data$nFeature_RNA,c(0.025,0.975))[2]),
       v=c(quantile(LowPDCtrl@meta.data$nCount_RNA,c(0.025,0.975))[1],quantile(LowPDCtrl@meta.data$nCount_RNA,c(0.025,0.975))[2]
	   ),lty=2,lwd=1,col="red")
LowPDCtrl <- subset(LowPDCtrl, subset = nFeature_RNA > quantile(LowPDCtrl@meta.data$nFeature_RNA,c(0.025,0.975))[1] & nFeature_RNA < quantile(LowPDCtrl@meta.data$nFeature_RNA,c(0.025,0.975))[2] & percent.mt < 10)	
	
	
####
HighPDCtrl <- CreateSeuratObject(counts = as.matrix(HighPDCtrl), project = "mouse.brain", min.cells = 10)

HighPDCtrl[["percent.mt"]] <- PercentageFeatureSet(HighPDCtrl, pattern = "^mt-")
VlnPlot(HighPDCtrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

quantile(HighPDCtrl@meta.data$nCount_RNA,c(0.025,0.975))
quantile(HighPDCtrl@meta.data$nFeature_RNA,c(0.025,0.975))
plot(HighPDCtrl@meta.data$nCount_RNA,HighPDCtrl@meta.data$nFeature_RNA,pch=16,cex=0.7,bty="n")
abline(h=c(quantile(HighPDCtrl@meta.data$nFeature_RNA,c(0.025,0.975))[1],quantile(HighPDCtrl@meta.data$nFeature_RNA,c(0.025,0.975))[2]),
       v=c(quantile(HighPDCtrl@meta.data$nCount_RNA,c(0.025,0.975))[1],quantile(HighPDCtrl@meta.data$nCount_RNA,c(0.025,0.975))[2]
	   ),lty=2,lwd=1,col="red")
HighPDCtrl <- subset(HighPDCtrl, subset = nFeature_RNA > quantile(HighPDCtrl@meta.data$nFeature_RNA,c(0.025,0.975))[1] & nFeature_RNA < quantile(HighPDCtrl@meta.data$nFeature_RNA,c(0.025,0.975))[2] & percent.mt < 10)	
	
	
####
Senes <- CreateSeuratObject(counts = as.matrix(Senes), project = "mouse.brain", min.cells = 10)

Senes[["percent.mt"]] <- PercentageFeatureSet(Senes, pattern = "^mt-")
VlnPlot(Senes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

quantile(Senes@meta.data$nCount_RNA,c(0.025,0.975))
quantile(Senes@meta.data$nFeature_RNA,c(0.025,0.975))
plot(Senes@meta.data$nCount_RNA,Senes@meta.data$nFeature_RNA,pch=16,cex=0.7,bty="n")
abline(h=c(quantile(Senes@meta.data$nFeature_RNA,c(0.025,0.975))[1],quantile(Senes@meta.data$nFeature_RNA,c(0.025,0.975))[2]),
       v=c(quantile(Senes@meta.data$nCount_RNA,c(0.025,0.975))[1],quantile(Senes@meta.data$nCount_RNA,c(0.025,0.975))[2]
	   ),lty=2,lwd=1,col="red")
Senes <- subset(Senes, subset = nFeature_RNA > quantile(Senes@meta.data$nFeature_RNA,c(0.025,0.975))[1] & nFeature_RNA < quantile(Senes@meta.data$nFeature_RNA,c(0.025,0.975))[2] & percent.mt < 10)


genes <- table(c(rownames(PD50),rownames(LowPDCtrl),rownames(HighPDCtrl),rownames(Senes)))
genes <- names(genes[genes==4])

PD50 <- GetAssayData(PD50)[genes,]
LowPDCtrl <- GetAssayData(LowPDCtrl)[genes,]
HighPDCtrl <- GetAssayData(HighPDCtrl)[genes,]
Senes <- GetAssayData(Senes)[genes,]

Senes_exp <- cbind(PD50,LowPDCtrl,HighPDCtrl,Senes)
Senes_exp <- as.matrix(Senes_exp)

#########################
library(ESCO)
params <- escoEstimate(Senes_exp)
sim<-escoSimulateGroups(params=params)
sim <- assays(sim)$TrueCounts

###Simulate profile with two groups
sim2 <-escoSimulateGroups(params=params,
                        group.prob = c(0.5, 0.5), deall.prob = 0.3,
                        de.prob = c(0.3, 0.7),
                        de.facLoc = c(1.9, 2.5), withcorr = TRUE, 
                        trials = 1, verbose =TRUE)
group <- sim2$Group
sim2 <- assays(sim2)$TrueCounts

sim3 <-escoSimulateGroups(params=params)
sim3 <- assays(sim3)$TrueCounts

rownames(sim) <- rownames(sim2) <- rownames(sim3) <-rownames(Senes_exp)

#########################
##Randomly seperate the cells into two groups, once for generating CSE and Pseudo-bulk, the other generate the reference
sample.ref.ind = sample(1:ncol(sim2),0.1*ncol(sim2),replace=FALSE)
ref.sim1 = sim[,sample.ref.ind]
ref.sim2 = sim2[,sample.ref.ind];ref2.group = group[sample.ref.ind]
ref.Senes = Senes_exp[,sample.ref.ind];ref.Senes.group = c(rep("PD50",ncol(PD50)),rep("LowPDCtrl",ncol(LowPDCtrl)),rep("HighPDCtrl",ncol(HighPDCtrl)),rep("Senes",ncol(Senes)))[sample.ref.ind]

##
bulk.sim1 = sim[,-sample.ref.ind]
bulk.sim2 = sim2[,-sample.ref.ind];bulk2.group = group[-sample.ref.ind]
bulk.sim3 = sim3[,-sample.ref.ind]
bulk.Senes = Senes_exp[,-sample.ref.ind];bulk.Senes.group = c(rep("PD50",ncol(PD50)),rep("LowPDCtrl",ncol(LowPDCtrl)),rep("HighPDCtrl",ncol(HighPDCtrl)),rep("Senes",ncol(Senes)))[-sample.ref.ind]


##################
#Assign the CSE identity
idCell = matrix(NA,nrow=200,ncol=3)
for(i in 1:nrow(idCell)){
  idCell[i,] = c(sample(c("PD50","LowPDCtrl","HighPDCtrl","Senes"),1),NA,sample(c("Group1","Group2"),1))
}
idCell = as.data.frame(idCell)
colnames(idCell) = c("Fibroblast","CellType2","CellType3")
rownames(idCell) = paste0("Sample-",1:200,sep="")

#################
####Generate CSE profile
H1_array <- array(0,
                  dim = c( 3,
                           8368,
                           200))
for(i in 1:dim(H1_array)[3]){
   ##For each CSE, we fixed sampled 50 cells to generate profile
   bulk.Senes.sub = bulk.Senes[,bulk.Senes.group %in% idCell$Fibroblast[i]]
   bulk.Senes.sub = bulk.Senes.sub[,sample(1:ncol(bulk.Senes.sub),50,replace=TRUE)]
   ####Normalization and calculate means
   bulk.Senes.sub = bulk.Senes.sub %*% diag(10^4/colSums(bulk.Senes.sub))
   
   ##CellType2
   bulk.sim1.sub = bulk.sim1[,sample(1:ncol(bulk.sim1),50,replace=TRUE)]
   ####Normalization and calculate means
   bulk.sim1.sub = bulk.sim1.sub %*% diag(10^4/colSums(bulk.sim1.sub))
   
   ##CellType3
   bulk.sim2.sub = bulk.sim2[,bulk2.group %in% idCell$CellType3[i]]
   bulk.sim2.sub = bulk.sim2.sub[,sample(1:ncol(bulk.sim2.sub),50,replace=TRUE)]
   ####Normalization and calculate means
   bulk.sim2.sub = bulk.sim2.sub %*% diag(10^4/colSums(bulk.sim2.sub))
   
   profile = cbind(rowMeans(bulk.Senes.sub),rowMeans(bulk.sim1.sub),rowMeans(bulk.sim2.sub))
   H1_array[,,i] = t(profile)
}


###simulate the cell type fractions
k <- 3 # number of cell types
ng <- 8368 # number of genes
p <- 200 # number of samples

cc <- matrix(runif(p*k), ncol=k)
cc <- t(scale(t(cc), center=FALSE, scale=rowSums(cc)))
colnames(cc) <- colnames(idCell)
rownames(cc) <- rownames(idCell)

###Generate the pesudo-bulk
PseudoBulk <- NULL
for(i in 1:200){
    PseudoBulk <- cbind(PseudoBulk, t(as.matrix(t(as.matrix(cc[i,])) %*% H1_array[,,i])))
}
rownames(PseudoBulk) = rownames(bulk.Senes.sub)
colnames(PseudoBulk) = rownames(idCell)

###Attach Noise
noise <- t(matrix(rnorm(p*ng)*10, ncol=ng))
noise[noise<0] <- 0
PseudoBulk <- PseudoBulk + noise

#calculate reference
ref = cbind(rowMeans(ref.Senes %*% diag(10^4/colSums(ref.Senes))),rowMeans(ref.sim1 %*% diag(10^4/colSums(ref.sim1))),rowMeans(ref.sim2 %*% diag(10^4/colSums(ref.sim2))))
colnames(ref) = colnames(idCell)

########################################################################################
#Evaluation
##calculate cell fraction matrix
Frac <- get_proportion(PseudoBulk, ref)
### run two TCA model
tca.mdl = tca(X = sqrt(PseudoBulk), W = Frac$theta, C1 = NULL, C2 = NULL,
                parallel = TRUE,num_cores=4,max_iters=20)
TCA_res1 <- tensor(X = (as.matrix(sqrt(PseudoBulk))), tca.mdl)

tca.mdl = tca(X = log2(PseudoBulk+1), W = Frac$theta, C1 = NULL, C2 = NULL,
                parallel = TRUE,num_cores=4,max_iters=20)
TCA_res2 <- tensor(X = (as.matrix(log2(PseudoBulk+1))), tca.mdl)

### run bMIND model
deconv1 = bMIND(bulk=sqrt(PseudoBulk), profile = sqrt(ref), ncore = 6,frac=Frac$theta)
deconv2 = bMIND(bulk=log2(PseudoBulk+1), profile = log2(ref+1), ncore = 6,frac=Frac$theta)

### run ENIGMA model

egm1 <- cell_deconvolve(X=as.matrix(sqrt(PseudoBulk)),
                                    theta=Frac$theta,
                                    R=sqrt(ref),
                                    epsilon=0.001,
                                    alpha=0.1,
                                    beta=0.1,tao_k=0.01,max.iter=50,verbose=TRUE,pre.process = "sqrt",Norm.method = "PC")

egm2 <- cell_deconvolve(X=as.matrix(log2(PseudoBulk+1)),
                                    theta=Frac$theta,
                                    R=log2(ref+1),
                                    epsilon=0.001,
                                    alpha=0.1,
                                    beta=0.1,tao_k=0.01,max.iter=50,verbose=TRUE,pre.process = "log",Norm.method = "PC")
									
### run ENIGMA(trace) model

egm_trace1 <- cell_deconvolve_trace(O = as.matrix(sqrt(PseudoBulk)),
                                                theta=Frac$theta,
                                                R=sqrt(ref),
                                                alpha=0.1,beta=1,solver="admm",gamma = 0.5,
                                                verbose=FALSE,max.iter = 1000,pre.process = "sqrt",Norm.method = "PC")
												
egm_trace2 <- cell_deconvolve_trace(O = as.matrix(log2(PseudoBulk+1)),
                                                theta=Frac$theta,
                                                R=log2(ref+1),
                                                alpha=0.1,beta=1,solver="admm",gamma = 0.5,
                                                verbose=FALSE,max.iter = 1000,pre.process = "log",Norm.method = "PC")

### evaluation												
												
### compare TCA
tca1 <- SingleCellExperiment(assays=list(logcounts =  TCA_res1[[1]]^2))
tca1$Group <- idCell$Fibroblast		 
tca1 <- runPCA(tca1)
tca1 <- runUMAP(tca1,dimred="PCA",n_dimred=10)
tca1 <- runTSNE(tca1,dimred="PCA",n_dimred=10)

ARI_tca1 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(tca1,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(tca1$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_tca1 <- c(ARI_tca1,ARI)
}


tca2 <- SingleCellExperiment(assays=list(logcounts =  2^TCA_res2[[1]]-1))
tca2$Group <- idCell$Fibroblast		 
tca2 <- runPCA(tca2)
tca2 <- runUMAP(tca2,dimred="PCA",n_dimred=10)
tca2 <- runTSNE(tca2,dimred="PCA",n_dimred=10)

ARI_tca2 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(tca2,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(tca2$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_tca2 <- c(ARI_tca2,ARI)
}



bmind1 <- SingleCellExperiment(assays=list(logcounts =  deconv1$A[,1,]^2))
bmind1$Group <- idCell$Fibroblast		 
bmind1 <- runPCA(bmind1)
bmind1 <- runUMAP(bmind1,dimred="PCA",n_dimred=10)
bmind1 <- runTSNE(bmind1,dimred="PCA",n_dimred=10)

ARI_bmind1 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(bmind1,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(bmind1$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_bmind1 <- c(ARI_bmind1,ARI)
}


bmind2 <- SingleCellExperiment(assays=list(logcounts =  2^deconv2$A[,2,]-1)) ## reconstructed bMIND into raw space
bmind2$Group <- idCell$Fibroblast		 
bmind2 <- runPCA(bmind2)
bmind2 <- runUMAP(bmind2,dimred="PCA",n_dimred=10)
bmind2 <- runTSNE(bmind2,dimred="PCA",n_dimred=10)

ARI_bmind2 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(bmind2,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(bmind2$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_bmind2 <- c(ARI_bmind2,ARI)
}

bmind_test <- SingleCellExperiment(assays=list(logcounts =  deconv1$A[,1,])) # using the deconvolved space can not reconstruct the four state
bmind_test$Group <- idCell$Fibroblast		 
bmind_test <- runPCA(bmind_test)
bmind_test <- runUMAP(bmind_test,dimred="PCA",n_dimred=10)
bmind_test <- runTSNE(bmind_test,dimred="PCA",n_dimred=10)

for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(bmind_test,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(bmind_test$Group)),as.numeric(clust$cluster))
print(ARI)
}


enigma1 <- SingleCellExperiment(assays=list(logcounts =  egm1$X_k_norm[,,1]))
enigma1$Group <- idCell$Fibroblast		 
enigma1 <- runPCA(enigma1)
enigma1 <- runUMAP(enigma1,dimred="PCA",n_dimred=10)
enigma1 <- runTSNE(enigma1,dimred="PCA",n_dimred=10)

ARI_enigma1 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(enigma1,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(enigma1$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_enigma1 <- c(ARI_enigma1,ARI)
}


enigma2 <- SingleCellExperiment(assays=list(logcounts =  egm2$X_k_norm[,,1]))
enigma2$Group <- idCell$Fibroblast		 
enigma2 <- runPCA(enigma2)
enigma2 <- runUMAP(enigma2,dimred="PCA",n_dimred=10)
enigma2 <- runTSNE(enigma2,dimred="PCA",n_dimred=10)

ARI_enigma2 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(enigma2,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(enigma2$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_enigma2 <- c(ARI_enigma2,ARI)
}



enigma_trace1 <- SingleCellExperiment(assays=list(logcounts =  egm_trace1$X_k_norm[,,1]))
enigma_trace1$Group <- idCell$Fibroblast		 
enigma_trace1 <- runPCA(enigma_trace1)
enigma_trace1 <- runUMAP(enigma_trace1,dimred="PCA",n_dimred=10)
enigma_trace1 <- runTSNE(enigma_trace1,dimred="PCA",n_dimred=10)

ARI_enigma_trace1 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(enigma_trace1,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(enigma_trace1$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_enigma_trace1 <- c(ARI_enigma_trace1,ARI)
}


enigma_trace2 <- SingleCellExperiment(assays=list(logcounts =  egm_trace2$X_k_norm[,,1]))
enigma_trace2$Group <- idCell$Fibroblast		 
enigma_trace2 <- runPCA(enigma_trace2)
enigma_trace2 <- runUMAP(enigma_trace2,dimred="PCA",n_dimred=10)
enigma_trace2 <- runTSNE(enigma_trace2,dimred="PCA",n_dimred=10)

ARI_enigma_trace2 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(enigma_trace2,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(enigma_trace2$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_enigma_trace2 <- c(ARI_enigma_trace2,ARI)
}

###################################################################
### show the results

color.legendF <- c(LowPDCtrl="orange", PD50="chartreuse4", Senes="magenta", HighPDCtrl="light blue")
colmatF <- col2rgb(color.legendF) 
color.legendF <- setNames(rgb(colmatF[1,], colmatF[2,], colmatF[3,], maxColorValue=255), names(color.legendF))
allcolors <- c(color.legendF[enigma$Group])
Y = reducedDim(enigma1,"UMAP")
png("enigma1(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	   
Y = reducedDim(enigma2,"UMAP")
png("enigma2(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	   


Y = reducedDim(tca1,"UMAP")
png("TCA1(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()
Y = reducedDim(tca2,"UMAP")
png("TCA2(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()


Y = reducedDim(bmind1,"UMAP")
png("bMIND1(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	
Y = reducedDim(bmind2,"UMAP")
png("bMIND2(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	

Y = reducedDim(enigma_trace1,"UMAP")
png("enigma_trace1(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	
Y = reducedDim(enigma_trace2,"UMAP")
png("enigma_trace2(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()


Y = reducedDim(enigma1,"TSNE")
png("enigma1(tsne).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	   
Y = reducedDim(enigma2,"TSNE")
png("enigma2(tsne).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	   


Y = reducedDim(tca1,"TSNE")
png("TCA1(tsne).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()
Y = reducedDim(tca2,"TSNE")
png("TCA2(tsne).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()


Y = reducedDim(bmind1,"TSNE")
png("bMIND1(tsne).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	
Y = reducedDim(bmind2,"TSNE")
png("bMIND2(tsne).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	

Y = reducedDim(enigma_trace1,"TSNE")
png("enigma_trace1(tsne).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	
Y = reducedDim(enigma_trace2,"TSNE")
png("enigma_trace2(tsne).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	

png("celltype_color.png",width=1500,height=1500,res=300)
plot(0,0,type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-1, y=1, legend=names(color.legendF), pch=21, cex=2.5, col=color.legendF, pt.bg=color.legendF, bty="n")
dev.off()

