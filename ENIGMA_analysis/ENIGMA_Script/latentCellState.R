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

##Running TCA
tca.mdl = tca(X = sqrt(PseudoBulk), W = Frac$theta, C1 = NULL, C2 = NULL,
                parallel = TRUE,num_cores=4,max_iters=20)
TCA_res2 <- tensor(X = (as.matrix(sqrt(PseudoBulk))), tca.mdl)

TCA <- SingleCellExperiment(assays=list(logcounts = (TCA_res2[[1]])))
TCA$Group <- idCell$Fibroblast			 
TCA <- runPCA(TCA)
TCA <- runUMAP(TCA,dimred="PCA",n_dimred=10)
png("tca.png",res=300,height=1500,width=1500)
p_TCA_mds
dev.off()


####Running ENIGMA
egm <- cell_deconvolve(X=as.matrix(sqrt(PseudoBulk)),
                                    theta=Frac$theta,
                                    R=sqrt(ref),
                                    epsilon=0.001,
                                    alpha=0.1,
                                    beta=0.5,tao_k=0.01,max.iter=1000,verbose=TRUE,pre.process = "sqrt",Norm.method = "PC")

###Running Trace Norm
egm_trace <- cell_deconvolve_trace(O = as.matrix(sqrt(PseudoBulk)),
                                                theta=Frac$theta,
                                                R=sqrt(ref),
                                                alpha=0.1,beta=1,solver="admm",
                                                verbose=TRUE,max.iter = 1000,pre.process = "sqrt",Norm.method = "PC")

enigma <- SingleCellExperiment(assays=list(logcounts = egm$X_k_norm[,,1]))
enigma$Group <- idCell$Fibroblast				 
enigma <- runPCA(enigma)
enigma <- runUMAP(enigma,dimred="PCA",n_dimred=10)

enigma_trace <- SingleCellExperiment(assays=list(logcounts = egm_trace$X_k_norm[,,1]))
enigma_trace$Group <- idCell$Fibroblast				 
enigma_trace <- runPCA(enigma_trace)
enigma_trace <- runUMAP(enigma_trace,dimred="PCA",n_dimred=10)

#########################
##Using bMIND to evaluate
bmind_res = bMIND(bulk=log2(PseudoBulk+1), profile = log2(ref+1), ncore = 6,frac=Frac$theta)
bmind <- SingleCellExperiment(assays=list(logcounts = (bmind_res$A[,1,])))
bmind$Group <- idCell$Fibroblast		 
bmind <- runPCA(bmind)
bmind <- runUMAP(bmind,dimred="PCA",n_dimred=10)

########################
##Using raw expression profile to compare
bulk_plot <- SingleCellExperiment(assays=list(logcounts = PseudoBulk))
bulk_plot$Group <- idCell$Fibroblast				 
bulk_plot <- runPCA(bulk_plot)
bulk_plot <- runUMAP(bulk_plot,dimred="PCA",n_dimred=10)

########################
##plot the ground truth
groundtruth <- SingleCellExperiment(assays=list(logcounts = H1_array[1,,]))
groundtruth$Group <- idCell$Fibroblast				 
groundtruth <- runPCA(groundtruth)
groundtruth <- runUMAP(groundtruth,dimred="PCA",n_dimred=10)


#######################
#Evaluate through adjusted rand index
ARI_bulk <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(bulk_plot,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(bulk_plot$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_bulk <- c(ARI_bulk,ARI)
}

ARI_tca <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(TCA,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(TCA$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_tca <- c(ARI_tca,ARI)
}

ARI_bmind <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(bmind,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(bmind$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_bmind <- c(ARI_bmind,ARI)
}

ARI_enigma <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(enigma,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(enigma$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_enigma <- c(ARI_enigma,ARI)
}

ARI_enigma <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(enigma_trace,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(enigma_trace$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_enigma <- c(ARI_enigma,ARI)
}

color.legendF <- c(LowPDCtrl="orange", PD50="chartreuse4", Senes="magenta", HighPDCtrl="light blue")
colmatF <- col2rgb(color.legendF) 
color.legendF <- setNames(rgb(colmatF[1,], colmatF[2,], colmatF[3,], maxColorValue=255), names(color.legendF))
allcolors <- c(color.legendF[enigma$Group])
Y = reducedDim(enigma,"UMAP")
png("enigma(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	   

Y = reducedDim(TCA,"UMAP")
png("TCA(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()

Y = reducedDim(bmind,"UMAP")
png("bMIND(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	

Y = reducedDim(enigma_trace,"UMAP")
png("enigma_trace(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	

Y = reducedDim(bulk_plot,"UMAP")
png("bulk(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()

Y = reducedDim(groundtruth,"UMAP")
png("groundtruth(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()

png("celltype_color.png",width=1500,height=1500,res=300)
plot(0,0,type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-1, y=1, legend=names(color.legendF), pch=21, cex=2.5, col=color.legendF, pt.bg=color.legendF, bty="n")
dev.off()

####################################################################################
###plot the bulk to make comparsion
###Using two baseline model to visualize 
###Model1: Each bulk expression profile regress out the respect cell type fractions
###Model2: TCA estimation
###Model3: PseudoBulk

###Model1: Regress out
fra_single <- Frac$theta[,1]
the <- (PseudoBulk %*% as.matrix(fra_single) - length(fra_single) * mean(fra_single) * rowMeans(PseudoBulk)) / (sum(fra_single^2) - length(fra_single)*mean(fra_single)^2)
Exp1 <- PseudoBulk - as.matrix(the) %*% t(as.matrix(fra_single))

###Model2: Sqrt
Exp2 <- sqrt(PseudoBulk)

###Model3: log
Exp3 <- log2(PseudoBulk+1)

###Model4: TCA estimation
sf <- apply(Frac$theta,1,function(x){norm(x,"2")^2})
Exp4 <- PseudoBulk %*% diag(Frac$theta[,1]/sf)

model1 <- SingleCellExperiment(assays=list(logcounts = Exp1))
model1$Group <- idCell$Fibroblast				 
model1 <- runPCA(model1)
model1 <- runUMAP(model1,dimred="PCA",n_dimred=10)
ARI_model1 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(model1,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(model1$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_model1 <- c(ARI_model1,ARI)
}

model2 <- SingleCellExperiment(assays=list(logcounts = Exp2))
model2$Group <- idCell$Fibroblast				 
model2 <- runPCA(model2)
model2 <- runUMAP(model2,dimred="PCA",n_dimred=10)
ARI_model2 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(model2,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(model2$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_model2 <- c(ARI_model2,ARI)
}

model3 <- SingleCellExperiment(assays=list(logcounts = Exp3))
model3$Group <- idCell$Fibroblast				 
model3 <- runPCA(model3)
model3 <- runUMAP(model3,dimred="PCA",n_dimred=10)
ARI_model3 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(model3,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(model3$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_model3 <- c(ARI_model3,ARI)
}

model4 <- SingleCellExperiment(assays=list(logcounts = Exp4))
model4$Group <- idCell$Fibroblast				 
model4 <- runPCA(model4)
model4 <- runUMAP(model4,dimred="PCA",n_dimred=10)
ARI_model4 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(model4,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(model4$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_model4 <- c(ARI_model4,ARI)
}

Y = reducedDim(model1,"UMAP")
png("model1(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	   

Y = reducedDim(model2,"UMAP")
png("model2(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()

Y = reducedDim(model3,"UMAP")
png("model3(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	

Y = reducedDim(model4,"UMAP")
png("model4(umap).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()	
########################################################################
##plot the identity of different cell type
model1$Group <- idCell$CellType3
model2$Group <- idCell$CellType3
model3$Group <- idCell$CellType3
model4$Group <- idCell$CellType3
bulk_plot$Group <- idCell$CellType3

color.legendA <- c(Group1="red", Group2="blue")
colmatA <- col2rgb(color.legendA) 
color.legendA <- setNames(rgb(colmatA[1,], colmatA[2,], colmatA[3,], maxColorValue=255), names(color.legendA))
allcolors <- c(color.legendA[model1$Group])
Y = reducedDim(model1,"UMAP")
png("model1(umap,celltype3).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()


Y = reducedDim(model2,"UMAP")
png("model2(umap,celltype3).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()

Y = reducedDim(model3,"UMAP")
png("model3(umap,celltype3).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()

Y = reducedDim(model4,"UMAP")
png("model4(umap,celltype3).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()

Y = reducedDim(bulk_plot,"UMAP")
png("bulk(umap,celltype3).png",width=1500,height=1500,res=300)
plot(Y[,1], Y[,2], cex=1,
     pch=21, 
     col=allcolors,
     bg=allcolors,xaxt="n",yaxt="n",xlab=NA,ylab=NA) 
dev.off()


png("celltype3_color.png",width=1500,height=1500,res=300)
plot(0,0,type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-1, y=1, legend=names(color.legendA), pch=21, cex=2.5, col=color.legendA, pt.bg=color.legendA, bty="n")
dev.off()
#############################################################################
###Fraction influence
enigma$cell_type_fraction = Frac$theta[,1]
p_enigma_umap <- plotUMAP(enigma, colour_by = "cell_type_fraction",point_size=3)

png("enigma_frac_umap.png",width=1500,height=1500,res=300)
p_enigma_umap
dev.off()

enigma_trace$cell_type_fraction = Frac$theta[,1]
p_enigma_umap <- plotUMAP(enigma_trace, colour_by = "cell_type_fraction",point_size=3)

png("enigma_frac_umap(trace).png",width=1500,height=1500,res=300)
p_enigma_umap
dev.off()
