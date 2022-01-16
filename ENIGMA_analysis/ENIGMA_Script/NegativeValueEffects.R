source("/path/to/Data/ENIGMA.R")
Bulk <- readRDS("/path/to/Data/Bulk.rds")
Reference <- readRDS("/path/to/Data/Reference.rds")
cellLabel <- readRDS("/path/to/Data/CellLabel.rds")
Frac <- get_proportion(Bulk, Reference)


######
#Running ENIGMA (proximal point method solver)
ENIGMA_trace_pps <- cell_deconvolve_trace(O = as.matrix(sqrt(Bulk)),
                                      theta=Frac$theta,
                                      R=sqrt(Reference),
                                      epsilon=0.00001,
                                      alpha=0.1,beta=0.5,solver="proximalpoint",tao_k = 0.5,
                                      verbose=TRUE,max.iter = 500,pre.process="sqrt")
#Testing the first celltype									  
sce <- SingleCellExperiment(assays = list(logcounts = ENIGMA_trace_pps$X_k_norm[,,1]))
sce$cell_type <- cellLabel$c1
sce <- sce[,Frac$theta[,1]>0.05]
sce <- runPCA(sce,50)
sce <- runTSNE(sce,dimred="PCA",n_dimred=5)
label <- sce$cell_type
ARI_ct1 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct1 <- c(ARI_ct1,ARI)
}

#Testing the third celltype
sce2 <- SingleCellExperiment(assays = list(logcounts = ENIGMA_trace_pps$X_k_norm[,,3]))
sce2$cell_type <- cellLabel$c3
sce2 <- sce2[,Frac$theta[,3]>0.05]
sce2 <- runPCA(sce2)
sce2 <- runTSNE(sce2,dimred="PCA",n_dimred=5)
label <- sce2$cell_type
ARI_ct3 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce2,"PCA")[,1:i],2)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct3 <- c(ARI_ct3,ARI)
}

ARI_proximalpoint <- list(ARI_ct1 = ARI_ct1,ARI_ct3 = ARI_ct3)

png("clustering_tsne(final).png",res=300,height=1200,width=2600)
p1 <- plotTSNE(sce, colour_by="cell_type",point_size=3)
p2 <- plotTSNE(sce2, colour_by="cell_type",point_size=3)
cowplot::plot_grid(p1,p2,nrow=1)
dev.off()

##########################################
##test negative value effects
X_k_nonneg <- ENIGMA_trace_pps$X_k_norm[,,1]
X_k_nonneg[X_k_nonneg<0] <- 0
sce <- SingleCellExperiment(assays = list(logcounts = X_k_nonneg))
sce$cell_type <- cellLabel$c1
sce <- sce[,Frac$theta[,1]>0.05]
sce <- runPCA(sce,50)
sce <- runTSNE(sce,dimred="PCA",n_dimred=5)
label <- sce$cell_type
ARI_ct1_nonneg <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct1_nonneg <- c(ARI_ct1_nonneg,ARI)
}

X_k_nonneg <- ENIGMA_trace_pps$X_k_norm[,,3]
X_k_nonneg[X_k_nonneg<0] <- 0
sce2 <- SingleCellExperiment(assays = list(logcounts = X_k_nonneg))
sce2$cell_type <- cellLabel$c3
sce2 <- sce2[,Frac$theta[,3]>0.05]
sce2 <- runPCA(sce2,50)
sce2 <- runTSNE(sce2,dimred="PCA",n_dimred=5)
label <- sce2$cell_type
ARI_ct3_nonneg <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce2,"PCA")[,1:i],2)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct3_nonneg <- c(ARI_ct3_nonneg,ARI)
}
ARI_proximalpoint_nonneg <- list(ARI_ct1 = ARI_ct1_nonneg,ARI_ct3 = ARI_ct3_nonneg)

#######################################################################################
#Running ENIGMA (admm method solver)
ENIGMA_trace <- cell_deconvolve_trace(O = as.matrix(sqrt(Bulk)),
                                      theta=Frac$theta,
                                      R=sqrt(Reference),
                                      alpha=0.1,solver="admm",beta = 0.5, epsilon = 10^-5,gamma=1,
                                      verbose=TRUE,max.iter = 500,pre.process="sqrt")

  
sce <- SingleCellExperiment(assays = list(logcounts = ENIGMA_trace$X_k_norm[,,1]))
sce$cell_type <- cellLabel$c1
sce <- sce[,Frac$theta[,1]>0.05] 
sce <- runPCA(sce,50)
sce <- runTSNE(sce,dimred="PCA",n_dimred=5)
label <- sce$cell_type
ARI_ct1 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct1 <- c(ARI_ct1,ARI)
}
		
sce2 <- SingleCellExperiment(assays = list(logcounts = ENIGMA_trace$X_k_norm[,,3]))
sce2$cell_type <- cellLabel$c3
sce2 <- sce2[,Frac$theta[,3]>0.05]
sce2 <- runPCA(sce2)
sce2 <- runTSNE(sce2,dimred="PCA",n_dimred=5)
label <- sce2$cell_type
ARI_ct3 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce2,"PCA")[,1:i],2)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct3 <- c(ARI_ct3,ARI)
}
ARI_admm <- list(ARI_ct1 = ARI_ct1,ARI_ct3 = ARI_ct3)

##################################################
##test negative value effects
X_k_nonneg <- ENIGMA_trace$X_k_norm[,,1];X_k_nonneg[X_k_nonneg<0] <- 0
sce <- SingleCellExperiment(assays = list(logcounts = X_k_nonneg))
sce$cell_type <- cellLabel$c1
sce <- sce[,Frac$theta[,1]>0.05]
sce <- runPCA(sce,50)
sce <- runTSNE(sce,dimred="PCA",n_dimred=5)
label <- sce$cell_type
ARI_ct1_nonneg <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct1_nonneg <- c(ARI_ct1_nonneg,ARI)
}

X_k_nonneg <- ENIGMA_trace$X_k_norm[,,3];X_k_nonneg[X_k_nonneg<0] <- 0
sce2 <- SingleCellExperiment(assays = list(logcounts = X_k_nonneg))
sce2$cell_type <- cellLabel$c3
sce2 <- sce2[,Frac$theta[,3]>0.05]
sce2 <- runPCA(sce2)
sce2 <- runTSNE(sce2,dimred="PCA",n_dimred=5)
label <- sce2$cell_type
ARI_ct3_nonneg <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce2,"PCA")[,1:i],2)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct3_nonneg <- c(ARI_ct3_nonneg,ARI)
}
ARI_admm_nonneg <- list(ARI_ct1 = ARI_ct1_nonneg,ARI_ct3 = ARI_ct3_nonneg)


png("clustering_admm(final).png",res=300,height=1200,width=2600)
p1 <- plotTSNE(sce, colour_by="cell_type",point_size=3)
p2 <- plotTSNE(sce2, colour_by="cell_type",point_size=3)
cowplot::plot_grid(p1,p2,nrow=1)
dev.off()


res1 <- max(ARI_proximalpoint$ARI_ct1)
res2 <- max(ARI_proximalpoint$ARI_ct3)
res3 <- max(ARI_admm$ARI_ct1)
res4 <- max(ARI_admm$ARI_ct3)
df <- data.frame(ARI = c(res1,res2,res3,res4),Method = c("Proximal Point","Proximal Point","ADMM","ADMM"),CellType = c("CellType1","CellType3","CellType1","CellType3"))
p<-ggplot(df, aes(x=CellType, y=ARI, fill=Method)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=20))+labs(x="")
png("BarplotARI.png",res=300,height=2000,width=1600)
p
dev.off()


res1 <- nrow(ENIGMA_trace$loss_history)
res2 <- nrow(ENIGMA_trace_pps$loss_history)
df <- data.frame(iteraction_num = c(res1,res2),Method = c("ADMM","Proximal Point"))
p<-ggplot(df, aes(x=Method, y=iteraction_num, fill=Method)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=20))+labs(x="",y="Number of Iterations")
png("BarplotARI(iterations_num).png",res=300,height=2000,width=1600)
p
dev.off()
###################################################################		
res1 <- max(ARI_proximalpoint$ARI_ct1)
res2 <- max(ARI_proximalpoint$ARI_ct3)
res3 <- max(ARI_proximalpoint_nonneg$ARI_ct1)
res4 <- max(ARI_proximalpoint_nonneg$ARI_ct3)
df <- data.frame(ARI = c(res1,res2,res3,res4),Method = c("Raw","Raw","Post-hoc","Post-hoc"),CellType = c("CellType1","CellType3","CellType1","CellType3"))
p<-ggplot(df, aes(x=CellType, y=ARI, fill=Method)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=20))+labs(x="")
png("BarplotARI(Proximal Point).png",res=300,height=2000,width=1600)
p
dev.off()

res1 <- max(ARI_admm$ARI_ct1)
res2 <- max(ARI_admm$ARI_ct3)
res3 <- max(ARI_admm_nonneg$ARI_ct1)
res4 <- max(ARI_admm_nonneg$ARI_ct3)
df <- data.frame(ARI = c(res1,res2,res3,res4),Method = c("Raw","Raw","Post-hoc","Post-hoc"),CellType = c("CellType1","CellType3","CellType1","CellType3"))
p<-ggplot(df, aes(x=CellType, y=ARI, fill=Method)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=20))+labs(x="")
png("BarplotARI(ADMM).png",res=300,height=2000,width=1600)
p
dev.off()


##density plot to show the effects of zero entries
neg_entries_1 <- ENIGMA_trace_pps$X_k_norm[,,1][ENIGMA_trace_pps$X_k_norm[,,1]<0]
nonneg_entries_1 <-  ENIGMA_trace_pps$X_k_norm[,,1][ENIGMA_trace_pps$X_k_norm[,,1]>0]
neg_entries_2 <- ENIGMA_trace_pps$X_k_norm[,,2][ENIGMA_trace_pps$X_k_norm[,,2]<0]
nonneg_entries_2 <-  ENIGMA_trace_pps$X_k_norm[,,2][ENIGMA_trace_pps$X_k_norm[,,2]>0]
neg_entries_3 <- ENIGMA_trace_pps$X_k_norm[,,3][ENIGMA_trace_pps$X_k_norm[,,3]<0]
nonneg_entries_3 <-  ENIGMA_trace_pps$X_k_norm[,,3][ENIGMA_trace_pps$X_k_norm[,,3]>0]

df <- data.frame(entries = abs(c(neg_entries_1,neg_entries_2,neg_entries_3,nonneg_entries_1,nonneg_entries_2,nonneg_entries_3)),ind = c(rep("negative entries",length(c(neg_entries_1,neg_entries_2,neg_entries_3))),rep("positive entries",length(c(nonneg_entries_1,nonneg_entries_2,nonneg_entries_3)))),celltype = c(rep("celltype1",length(neg_entries_1)),rep("celltype2",length(neg_entries_2)),rep("celltype3",length(neg_entries_3)),rep("celltype1",length(nonneg_entries_1)),rep("celltype2",length(nonneg_entries_2)),rep("celltype3",length(nonneg_entries_3))))

p<-ggplot(df, aes(x=celltype, y=entries, fill=ind)) +
    geom_boxplot(position=position_dodge(1)) +theme_classic2() + theme(text = element_text(size=20))+labs(x="",y="imputed expression values",fill="")
png("negEffects(boxplot,proximalpoint).png",res=300,height=2000,width=2000)
p
dev.off()

neg_entries_1 <- ENIGMA_trace$X_k_norm[,,1][ENIGMA_trace$X_k_norm[,,1]<0]
nonneg_entries_1 <-  ENIGMA_trace$X_k_norm[,,1][ENIGMA_trace$X_k_norm[,,1]>0]
neg_entries_2 <- ENIGMA_trace$X_k_norm[,,2][ENIGMA_trace$X_k_norm[,,2]<0]
nonneg_entries_2 <-  ENIGMA_trace$X_k_norm[,,2][ENIGMA_trace$X_k_norm[,,2]>0]
neg_entries_3 <- ENIGMA_trace$X_k_norm[,,3][ENIGMA_trace$X_k_norm[,,3]<0]
nonneg_entries_3 <-  ENIGMA_trace$X_k_norm[,,3][ENIGMA_trace$X_k_norm[,,3]>0]

df <- data.frame(entries = abs(c(neg_entries_1,neg_entries_2,neg_entries_3,nonneg_entries_1,nonneg_entries_2,nonneg_entries_3)),ind = c(rep("negative entries",length(c(neg_entries_1,neg_entries_2,neg_entries_3))),rep("positive entries",length(c(nonneg_entries_1,nonneg_entries_2,nonneg_entries_3)))),celltype = c(rep("celltype1",length(neg_entries_1)),rep("celltype2",length(neg_entries_2)),rep("celltype3",length(neg_entries_3)),rep("celltype1",length(nonneg_entries_1)),rep("celltype2",length(nonneg_entries_2)),rep("celltype3",length(nonneg_entries_3))))

p<-ggplot(df, aes(x=celltype, y=entries, fill=ind)) +
    geom_boxplot(position=position_dodge(1)) +theme_classic2() + theme(text = element_text(size=20))+labs(x="",y="imputed expression values",fill="")
png("negEffects(boxplot).png",res=300,height=2000,width=2000)
p
dev.off()

##density plot to show the effects of zero entries (on un-normalized data)
ENIGMA_trace <- cell_deconvolve_trace(O = as.matrix(sqrt(Bulk)),
                                      theta=Frac$theta,
                                      R=sqrt(Reference),
                                      alpha=0.1,solver="admm",beta = 0.5, epsilon = 10^-5,gamma=1,
                                      verbose=TRUE,max.iter = 500,pre.process="sqrt",pos=FALSE)
ENIGMA_trace_pps <- cell_deconvolve_trace(O = as.matrix(sqrt(Bulk)),
                                      theta=Frac$theta,
                                      R=sqrt(Reference),
                                      epsilon=0.00001,
                                      alpha=0.1,beta=0.5,solver="proximalpoint",tao_k = 0.5,
                                      verbose=TRUE,max.iter = 500,pre.process="sqrt",pos=FALSE)
neg_entries_1 <- ENIGMA_trace_pps$X_k[,,1][ENIGMA_trace_pps$X_k[,,1]<0]
nonneg_entries_1 <-  ENIGMA_trace_pps$X_k[,,1][ENIGMA_trace_pps$X_k[,,1]>0]
neg_entries_2 <- ENIGMA_trace_pps$X_k[,,2][ENIGMA_trace_pps$X_k[,,2]<0]
nonneg_entries_2 <-  ENIGMA_trace_pps$X_k[,,2][ENIGMA_trace_pps$X_k[,,2]>0]
neg_entries_3 <- ENIGMA_trace_pps$X_k[,,3][ENIGMA_trace_pps$X_k[,,3]<0]
nonneg_entries_3 <-  ENIGMA_trace_pps$X_k[,,3][ENIGMA_trace_pps$X_k[,,3]>0]

df <- data.frame(entries = abs(c(neg_entries_1,neg_entries_2,neg_entries_3,nonneg_entries_1,nonneg_entries_2,nonneg_entries_3)),ind = c(rep("negative entries",length(c(neg_entries_1,neg_entries_2,neg_entries_3))),rep("positive entries",length(c(nonneg_entries_1,nonneg_entries_2,nonneg_entries_3)))),celltype = c(rep("celltype1",length(neg_entries_1)),rep("celltype2",length(neg_entries_2)),rep("celltype3",length(neg_entries_3)),rep("celltype1",length(nonneg_entries_1)),rep("celltype2",length(nonneg_entries_2)),rep("celltype3",length(nonneg_entries_3))))

p<-ggplot(df, aes(x=celltype, y=entries, fill=ind)) +
    geom_boxplot(position=position_dodge(1)) +theme_classic2() + theme(text = element_text(size=20))+labs(x="",y="imputed expression values",fill="")
png("negEffects(boxplot,proximalpoint,raw).png",res=300,height=2000,width=2000)
p
dev.off()

neg_entries_1 <- ENIGMA_trace$X_k[,,1][ENIGMA_trace$X_k[,,1]<0]
nonneg_entries_1 <-  ENIGMA_trace$X_k[,,1][ENIGMA_trace$X_k[,,1]>0]
neg_entries_2 <- ENIGMA_trace$X_k[,,2][ENIGMA_trace$X_k[,,2]<0]
nonneg_entries_2 <-  ENIGMA_trace$X_k[,,2][ENIGMA_trace$X_k[,,2]>0]
neg_entries_3 <- ENIGMA_trace$X_k[,,3][ENIGMA_trace$X_k[,,3]<0]
nonneg_entries_3 <-  ENIGMA_trace$X_k[,,3][ENIGMA_trace$X_k[,,3]>0]

df <- data.frame(entries = abs(c(neg_entries_1,neg_entries_2,neg_entries_3,nonneg_entries_1,nonneg_entries_2,nonneg_entries_3)),ind = c(rep("negative entries",length(c(neg_entries_1,neg_entries_2,neg_entries_3))),rep("positive entries",length(c(nonneg_entries_1,nonneg_entries_2,nonneg_entries_3)))),celltype = c(rep("celltype1",length(neg_entries_1)),rep("celltype2",length(neg_entries_2)),rep("celltype3",length(neg_entries_3)),rep("celltype1",length(nonneg_entries_1)),rep("celltype2",length(nonneg_entries_2)),rep("celltype3",length(nonneg_entries_3))))

p<-ggplot(df, aes(x=celltype, y=entries, fill=ind)) +
    geom_boxplot(position=position_dodge(1)) +theme_classic2() + theme(text = element_text(size=20))+labs(x="",y="imputed expression values",fill="")
png("negEffects(boxplot,raw).png",res=300,height=2000,width=2000)
p
dev.off()

######################################
###Running the L2max norm
ENIGMA_l2max <- cell_deconvolve(X = as.matrix(sqrt(Bulk)),
                                      theta=Frac$theta,
                                      R=sqrt(Reference),
                                      alpha=0.1,beta = 0.5, epsilon = 0.001,
                                      verbose=TRUE,max.iter = 500,pre.process="sqrt")

sce <- SingleCellExperiment(assays = list(logcounts = ENIGMA_l2max$X_k_norm[,,1]))
sce$cell_type <- cellLabel$c1
sce <- sce[,Frac$theta[,1]>0.05] 
sce <- runPCA(sce,50)
sce <- runTSNE(sce,dimred="PCA",n_dimred=5)
label <- sce$cell_type
ARI_ct1 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct1 <- c(ARI_ct1,ARI)
}
		
sce2 <- SingleCellExperiment(assays = list(logcounts = ENIGMA_l2max$X_k_norm[,,3]))
sce2$cell_type <- cellLabel$c3
sce2 <- sce2[,Frac$theta[,3]>0.05]
sce2 <- runPCA(sce2)
sce2 <- runTSNE(sce2,dimred="PCA",n_dimred=5)
label <- sce2$cell_type
ARI_ct3 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce2,"PCA")[,1:i],2)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct3 <- c(ARI_ct3,ARI)
}
ARI_l2max <- list(ARI_ct1 = ARI_ct1,ARI_ct3 = ARI_ct3)


##test negative value effects
X_k_nonneg <- ENIGMA_l2max$X_k_norm[,,1];X_k_nonneg[X_k_nonneg<0] <- 0
sce <- SingleCellExperiment(assays = list(logcounts = X_k_nonneg))
sce$cell_type <- cellLabel$c1
sce <- sce[,Frac$theta[,1]>0.05]
sce <- runPCA(sce,50)
sce <- runTSNE(sce,dimred="PCA",n_dimred=5)
label <- sce$cell_type
ARI_ct1_nonneg <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct1_nonneg <- c(ARI_ct1_nonneg,ARI)
}

X_k_nonneg <- ENIGMA_l2max$X_k_norm[,,3];X_k_nonneg[X_k_nonneg<0] <- 0
sce2 <- SingleCellExperiment(assays = list(logcounts = X_k_nonneg))
sce2$cell_type <- cellLabel$c3
sce2 <- sce2[,Frac$theta[,3]>0.05]
sce2 <- runPCA(sce2)
sce2 <- runTSNE(sce2,dimred="PCA",n_dimred=5)
label <- sce2$cell_type
ARI_ct3_nonneg <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce2,"PCA")[,1:i],2)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct3_nonneg <- c(ARI_ct3_nonneg,ARI)
}
ARI_l2max_nonneg <- list(ARI_ct1 = ARI_ct1_nonneg,ARI_ct3 = ARI_ct3_nonneg)

res1 <- max(ARI_l2max$ARI_ct1)
res2 <- max(ARI_l2max$ARI_ct3)
res3 <- max(ARI_l2max_nonneg$ARI_ct1)
res4 <- max(ARI_l2max_nonneg$ARI_ct3)
df <- data.frame(ARI = c(res1,res2,res3,res4),Method = c("Raw","Raw","Post-hoc","Post-hoc"),CellType = c("CellType1","CellType3","CellType1","CellType3"))
p<-ggplot(df, aes(x=CellType, y=ARI, fill=Method)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=20))+labs(x="")
png("BarplotARI(l2max).png",res=300,height=2000,width=1600)
p
dev.off()

########
neg_entries_1 <- ENIGMA_l2max$X_k_norm[,,1][ENIGMA_l2max$X_k_norm[,,1]<0]
nonneg_entries_1 <-  ENIGMA_l2max$X_k_norm[,,1][ENIGMA_l2max$X_k_norm[,,1]>0]
neg_entries_2 <- ENIGMA_l2max$X_k_norm[,,2][ENIGMA_l2max$X_k_norm[,,2]<0]
nonneg_entries_2 <-  ENIGMA_l2max$X_k_norm[,,2][ENIGMA_l2max$X_k_norm[,,2]>0]
neg_entries_3 <- ENIGMA_l2max$X_k_norm[,,3][ENIGMA_l2max$X_k_norm[,,3]<0]
nonneg_entries_3 <-  ENIGMA_l2max$X_k_norm[,,3][ENIGMA_l2max$X_k_norm[,,3]>0]

df <- data.frame(entries = abs(c(neg_entries_1,neg_entries_2,neg_entries_3,nonneg_entries_1,nonneg_entries_2,nonneg_entries_3)),ind = c(rep("negative entries",length(c(neg_entries_1,neg_entries_2,neg_entries_3))),rep("positive entries",length(c(nonneg_entries_1,nonneg_entries_2,nonneg_entries_3)))),celltype = c(rep("celltype1",length(neg_entries_1)),rep("celltype2",length(neg_entries_2)),rep("celltype3",length(neg_entries_3)),rep("celltype1",length(nonneg_entries_1)),rep("celltype2",length(nonneg_entries_2)),rep("celltype3",length(nonneg_entries_3))))

p<-ggplot(df, aes(x=celltype, y=entries, fill=ind)) +
    geom_boxplot(position=position_dodge(1)) +theme_classic2() + theme(text = element_text(size=20))+labs(x="",y="imputed expression values",fill="")
png("negEffects(boxplot,l2max).png",res=300,height=2000,width=2000)
p
dev.off()
#########################

ENIGMA_l2max <- cell_deconvolve(X = as.matrix(sqrt(Bulk)),
                                      theta=Frac$theta,
                                      R=sqrt(Reference),
                                      alpha=0.1,beta = 0.5, epsilon = 0.001,
                                      verbose=TRUE,max.iter = 500,pre.process="sqrt",pos=FALSE)
neg_entries_1 <- ENIGMA_l2max$X_k[,,1][ENIGMA_l2max$X_k[,,1]<0]
nonneg_entries_1 <-  ENIGMA_l2max$X_k[,,1][ENIGMA_l2max$X_k[,,1]>0]
neg_entries_2 <- ENIGMA_l2max$X_k[,,2][ENIGMA_l2max$X_k[,,2]<0]
nonneg_entries_2 <-  ENIGMA_l2max$X_k[,,2][ENIGMA_l2max$X_k[,,2]>0]
neg_entries_3 <- ENIGMA_l2max$X_k[,,3][ENIGMA_l2max$X_k[,,3]<0]
nonneg_entries_3 <-  ENIGMA_l2max$X_k[,,3][ENIGMA_l2max$X_k[,,3]>0]

df <- data.frame(entries = abs(c(neg_entries_1,neg_entries_2,neg_entries_3,nonneg_entries_1,nonneg_entries_2,nonneg_entries_3)),ind = c(rep("negative entries",length(c(neg_entries_1,neg_entries_2,neg_entries_3))),rep("positive entries",length(c(nonneg_entries_1,nonneg_entries_2,nonneg_entries_3)))),celltype = c(rep("celltype1",length(neg_entries_1)),rep("celltype2",length(neg_entries_2)),rep("celltype3",length(neg_entries_3)),rep("celltype1",length(nonneg_entries_1)),rep("celltype2",length(nonneg_entries_2)),rep("celltype3",length(nonneg_entries_3))))

p<-ggplot(df, aes(x=celltype, y=entries, fill=ind)) +
    geom_boxplot(position=position_dodge(1)) +theme_classic2() + theme(text = element_text(size=20))+labs(x="",y="imputed expression values",fill="")
png("negEffects(boxplot,l2max,raw).png",res=300,height=2000,width=2000)
p
dev.off()


#########################################################
##Benchmark CTS-DEG detection
load("DEG_test_data_4.8")
source("DEG_analysis_uile_function.R")

#####refixed the negative values of inputted simulated gene expression matrix into positive
value = max(abs(Bulk[Bulk<0]))
Bulk = Bulk + value
Reference = Reference + value
Frac <- get_proportion(Bulk, Reference)
y <- gl(2, 100/2)

ENIGMA_trace <- cell_deconvolve_trace(O = as.matrix(Bulk),
                                              theta=Frac$theta,
                                              R=Reference,
                                              alpha=0.1,beta=1,solver="admm",gamma = 0.1,
                                              verbose=FALSE,max.iter = 1000,pos=FALSE,Norm.method = "frac",pre.process="none")
result <- DEG_test(ENIGMA_trace$X_k_norm,y)
###Evaluation
Tab <- result$pval
perform_frac_admm_neg <- DePRCcalculator(Tab,ENIGMA_trace$X_k_norm,DEG_list,y,10000)$AUC_PR


ENIGMA_trace$X_k_norm[ENIGMA_trace$X_k_norm<0] <- 0
result <- DEG_test(ENIGMA_trace$X_k_norm,y)
###Evaluation
Tab <- result$pval
perform_frac_admm_pos <- DePRCcalculator(Tab,ENIGMA_trace$X_k_norm,DEG_list,y,10000)$AUC_PR



######################################################################
ENIGMA_trace <- cell_deconvolve_trace(O = as.matrix(Bulk),
                                              theta=Frac$theta,
                                              R=Reference,
                                              alpha=0.1,beta=1,solver="proximalpoint",tao_k = 4,
                                              verbose=FALSE,max.iter = 1000,pos=FALSE,Norm.method = "frac",pre.process="none")
result <- DEG_test(ENIGMA_trace$X_k_norm,y)
###Evaluation
Tab <- result$pval
perform_frac_pps_neg <- DePRCcalculator(Tab,ENIGMA_trace$X_k_norm,DEG_list,y,10000)$AUC_PR


ENIGMA_trace$X_k_norm[ENIGMA_trace$X_k_norm<0] <- 0
result <- DEG_test(ENIGMA_trace$X_k_norm,y)
###Evaluation
Tab <- result$pval
perform_frac_pps_pos <- DePRCcalculator(Tab,ENIGMA_trace$X_k_norm,DEG_list,y,10000)$AUC_PR


ENIGMA_l2max <- cell_deconvolve(X=as.matrix(Bulk),
                                        theta=Frac$theta,
                                        R=Reference,
                                        epsilon=0.001,
                                        alpha=0.1,
                                        beta=0.5,tao_k=0.01,max.iter=1000,verbose=FALSE,pos=FALSE,Norm.method = "frac")
										
result <- DEG_test(ENIGMA_l2max$X_k_norm,y)
###Evaluation
Tab <- result$pval
perform_frac_l2norm_neg <- DePRCcalculator(Tab,ENIGMA_l2max$X_k_norm,DEG_list,y,10000)$AUC_PR

ENIGMA_l2max$X_k_norm[ENIGMA_l2max$X_k_norm<0] <- 0
result <- DEG_test(ENIGMA_l2max$X_k_norm,y)
###Evaluation
Tab <- result$pval
perform_frac_l2norm_pos <- DePRCcalculator(Tab,ENIGMA_l2max$X_k_norm,DEG_list,y,10000)$AUC_PR

#####################################

df <- data.frame(AUPRC = c(perform_frac_admm_pos,perform_frac_admm_neg,perform_frac_pps_pos,perform_frac_pps_neg,perform_frac_l2norm_pos,perform_frac_l2norm_neg),Method = rep(c(rep("Raw",5),rep("Post-hoc",5)),3),CellType = rep(rep(paste0("CellType-",1:5,sep=""),2),3),Model = c(rep("Trace(ADMM)",10),rep("Trace(ProximalPoint)",10),rep("L2 norm",10)))
p<-ggplot(df, aes(x=CellType, y=AUPRC, fill=Method)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=20))+labs(x="") +  facet_wrap(~Model, nrow=1)
png("BarplotAUPRC.png",res=300,height=1600,width=4500)
p
dev.off()















