###############################################################
###study the alpha effects and why we conduct normalization to the reference
###load HNSCC and simulate bulk RNA-seq samples of cell state identification

###HNSCC data
source("/Path/to/Data/ENIGMA.R")
library(MASS)
library(scater)
library(SingleCellExperiment)
###calculate cell type fractions matrix
Bulk <- readRDS("/Path/to/Data/Bulk.rds")
Reference <- readRDS("/Path/to/Data/Reference.rds")
CellLabel <- readRDS("/Path/to/Data/CellLabel.rds")
###calculate cell type fractions matrix
Frac <- get_proportion(Bulk,Reference)

cor_gene <- NULL
for(alpha in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)){
res_alg_all <- cell_deconvolve(X=sqrt(Bulk),
                    theta=Frac$theta,
					R=sqrt(Reference),
					epsilon=0.001,
					alpha=alpha,
					beta=0.5,tao_k=0.1,verbose=FALSE,Norm.method = "PC",pre.process = "sqrt")
cor_per_ct <- NULL
for(k in 1:3){
   exp <- res_alg_all$X_k[,,k]
   exp.scale <- t(apply(exp,1,scale))
   PC <- svd(exp.scale)$v[,1:10]
   cor <- apply(PC,2,function(x){cor(x,Frac$theta[,k],method="sp")})
   cor_per_ct <- c(cor_per_ct,max(abs(cor)))
}			
cor_gene <- cbind(cor_gene, cor_per_ct)
}

##################################
#plot PCA
res_alg_all <- cell_deconvolve(X=sqrt(Bulk),
                    theta=Frac$theta,
					R=sqrt(Reference),
					epsilon=0.001,
					alpha=0.3,
					beta=0.5,tao_k=0.01,verbose=TRUE,Norm.method = "PC",pre.process = "sqrt")
					
enigma <- SingleCellExperiment(assays=list(logcounts = res_alg_all$X_k[,,1]^2))
enigma$Group <- CellLabel$c1 
enigma <- runPCA(enigma)
enigma <- runUMAP(enigma,dimred="PCA",n_dimred=10)
PCA_embedding <- reducedDim(enigma, "PCA")
PCA_embedding <- PCA_embedding[,c(1,2)]
colnames(PCA_embedding) <- c("PC1","PC2")
reducedDim(enigma, "PCA") <- PCA_embedding

png("/Path/to/Data/raw_PC1_2.png",res=300,height=1500,width=2000)
plotReducedDim(enigma,colour_by = "Group",dimred = "PCA",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))
dev.off()
png("/Path/to/Data/raw_UMAP.png",res=300,height=1500,width=2000)
plotReducedDim(enigma,colour_by = "Group",dimred = "UMAP",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))
dev.off()

enigma <- runPCA(enigma)
PCA_embedding <- reducedDim(enigma, "PCA")
PCA_embedding <- PCA_embedding[,c(2,3)]
colnames(PCA_embedding) <- c("PC2","PC3")
reducedDim(enigma, "PCA") <- PCA_embedding

png("/Path/to/Data/raw_PC2_3.png",res=300,height=1500,width=2000)
plotReducedDim(enigma,colour_by = "Group",dimred = "PCA",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))+labs(x = "PCA 2", y = "PCA 3")
dev.off()

#############################################
enigma <- SingleCellExperiment(assays=list(logcounts = res_alg_all$X_k_norm[,,1]))
enigma$Group <- CellLabel$c1 
enigma <- runPCA(enigma)
enigma <- runUMAP(enigma,dimred="PCA",n_dimred=10)
PCA_embedding <- reducedDim(enigma, "PCA")
PCA_embedding <- PCA_embedding[,c(1,2)]
colnames(PCA_embedding) <- c("PC1","PC2")
reducedDim(enigma, "PCA") <- PCA_embedding

png("/Path/to/Data/new_PC1_2.png",res=300,height=1500,width=2000)
plotReducedDim(enigma,colour_by = "Group",dimred = "PCA",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))
dev.off()
png("/Path/to/Data/new_UMAP.png",res=300,height=1500,width=2000)
plotReducedDim(enigma,colour_by = "Group",dimred = "UMAP",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))
dev.off()

enigma <- runPCA(enigma)
PCA_embedding <- reducedDim(enigma, "PCA")
PCA_embedding <- PCA_embedding[,c(2,3)]
colnames(PCA_embedding) <- c("PC2","PC3")
reducedDim(enigma, "PCA") <- PCA_embedding

png("/Path/to/Data/new_PC2_3.png",res=300,height=1500,width=2000)
plotReducedDim(enigma,colour_by = "Group",dimred = "PCA",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))+labs(x = "PCA 2", y = "PCA 3")
dev.off()
##############################################


cor_ari_performance1 <- cor_ari_performance2 <- NULL
for(alpha in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)){
res_alg_all <- cell_deconvolve(X=sqrt(Bulk),
                    theta=Frac$theta,
					R=sqrt(Reference),
					epsilon=0.001,
					alpha=alpha,
					beta=0.5,tao_k=0.1,verbose=FALSE,Norm.method = "PC",pre.process = "sqrt")
enigma <- SingleCellExperiment(assays=list(logcounts = res_alg_all$X_k[,,1]))
enigma$Group <- CellLabel$c1 
enigma <- runPCA(enigma)

ARI_enigma <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(enigma,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(enigma$Group)),as.numeric(clust$cluster))
ARI_enigma <- c(ARI_enigma,ARI)
}	
cor_ari_performance1 <- c(cor_ari_performance1,max(ARI_enigma))

enigma <- SingleCellExperiment(assays=list(logcounts = res_alg_all$X_k[,,3]))
enigma$Group <- CellLabel$c3
enigma <- runPCA(enigma)

ARI_enigma <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(enigma,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(enigma$Group)),as.numeric(clust$cluster))
ARI_enigma <- c(ARI_enigma,ARI)
}
cor_ari_performance2 <- c(cor_ari_performance2,max(ARI_enigma))
}


cor_ari_performance1_new <- cor_ari_performance2_new <- NULL
for(alpha in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)){
res_alg_all <- cell_deconvolve(X=sqrt(Bulk),
                    theta=Frac$theta,
					R=sqrt(Reference),
					epsilon=0.001,
					alpha=alpha,
					beta=0.5,tao_k=0.1,verbose=FALSE,Norm.method = "frac",pre.process = "sqrt")
enigma <- SingleCellExperiment(assays=list(logcounts = res_alg_all$X_k_norm[,,1]))
enigma$Group <- CellLabel$c1 
enigma <- runPCA(enigma)

ARI_enigma <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(enigma,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(enigma$Group)),as.numeric(clust$cluster))
ARI_enigma <- c(ARI_enigma,ARI)
}	
cor_ari_performance1_new <- c(cor_ari_performance1_new,max(ARI_enigma))

enigma <- SingleCellExperiment(assays=list(logcounts = res_alg_all$X_k_norm[,,3]))
enigma$Group <- CellLabel$c3
enigma <- runPCA(enigma)

ARI_enigma <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(enigma,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(enigma$Group)),as.numeric(clust$cluster))
ARI_enigma <- c(ARI_enigma,ARI)
}
cor_ari_performance2_new <- c(cor_ari_performance2_new,max(ARI_enigma))
}

corDat = data.frame(cor = as.numeric(t(cor_gene)),alpha = rep(c(1:9)*0.1,3),celltype = c(rep("celltype-1",9),rep("celltype-2",9),rep("celltype-3",9)))

pp=ggplot(corDat,aes(x=alpha,y=cor,fill=celltype,colour=celltype))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        panel.grid=element_blank())+
  ylab("Correlation with cell type fractions")

png("/Path/to/Data/cor_celltype_fraction.png",res=300,height = 1600,width = 2200)
pp
dev.off()


performDat = data.frame(ARI = c(cor_ari_performance1,cor_ari_performance1_new),alpha = rep(c(1:9)*0.1,2),method = c(rep("Before Normalization",9),rep("After Normalization",9)))

pp=ggplot(performDat,aes(x=alpha,y=ARI,fill=method,colour=method))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        panel.grid=element_blank())+
  ylab("Adjusted Rand Index")

png("/Path/to/Data/ARI_ct1(frac).png",res=300,height = 1600,width = 2200)
pp
dev.off()


performDat = data.frame(ARI = c(cor_ari_performance2,cor_ari_performance2_new),alpha = rep(c(1:9)*0.1,2),method = c(rep("Before Normalization",9),rep("After Normalization",9)))

pp=ggplot(performDat,aes(x=alpha,y=ARI,fill=method,colour=method))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        panel.grid=element_blank())+
  ylab("Adjusted Rand Index")

png("/mnt/data1/weixu/HiDe/revised/Why_normalize/ARI_ct2(frac).png",res=300,height = 1600,width = 2200)
pp
dev.off()

#############################################################################

load("/Path/to/Data/DEG_test_data_4.8")
source("/Path/to/Data/DEG_analysis_uile_function.R")
Frac <- get_proportion(Bulk, Reference)
y <- gl(2, 100/2)

ENIGMA_trace <- cell_deconvolve_trace(O = as.matrix(Bulk),
                                              theta=Frac$theta,
                                              R=Reference,
                                              alpha=0.1,beta=1,solver="admm",
                                              verbose=FALSE,max.iter = 1000,pos=FALSE,Norm.method = "frac",pre.process="none")
result <- DEG_test(ENIGMA_trace$X_k_norm,y)
###Evaluation
Tab <- result$qval
perform_frac <- DePRCcalculator(Tab,ENIGMA_trace$X_k_norm,DEG_list,y,10000)$AUC_PR


ENIGMA_trace <- cell_deconvolve_trace(O = as.matrix(Bulk),
                                              theta=Frac$theta,
                                              R=Reference,
                                              alpha=0.1,beta=1,solver="admm",
                                              verbose=FALSE,max.iter = 1000,pos=FALSE,Norm.method = "PC",pre.process="none")
result <- DEG_test(ENIGMA_trace$X_k_norm,y)
###Evaluation
Tab <- result$qval
perform_pc <- DePRCcalculator(Tab,ENIGMA_trace$X_k_norm,DEG_list,y,10000)$AUC_PR

result <- DEG_test(ENIGMA_trace$X_k,y)
###Evaluation
Tab <- result$qval
perform_raw <- DePRCcalculator(Tab,ENIGMA_trace$X_k,DEG_list,y,10000)$AUC_PR

cor_gene <- NULL
auprc <- auprc_norm <- NULL
for(alpha in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)){
ENIGMA_l2max <- cell_deconvolve(X=as.matrix(Bulk),
                                        theta=Frac$theta,
                                        R=Reference,
                                        epsilon=0.001,
                                        alpha=alpha,
                                        beta=0.5,tao_k=0.01,max.iter=1000,verbose=FALSE,pos=FALSE,Norm.method = "frac")
										
cor_per_ct <- NULL
for(k in 1:5){
   exp <- ENIGMA_l2max$X_k[,,k]
   exp.scale <- t(apply(exp,1,scale))
   PC <- as.matrix(svd(exp.scale)$v[,1:10])
   cor <- apply(PC,2,function(x){cor(x,Frac$theta[,k],method="sp")})
   cor_per_ct <- c(cor_per_ct,max(abs(cor)))
}			
cor_gene <- cbind(cor_gene, cor_per_ct)

result <- DEG_test(ENIGMA_l2max$X_k,y)
###Evaluation
Tab <- result$pval
perform_raw <- DePRCcalculator(Tab,ENIGMA_l2max$X_k,DEG_list,y,10000)$AUC_PR
perform_raw[is.nan(perform_raw)] <- 0
auprc <- c(auprc,mean(perform_raw))


result <- DEG_test(ENIGMA_l2max$X_k_norm,y)
###Evaluation
Tab <- result$pval
perform_raw <- DePRCcalculator(Tab,ENIGMA_l2max$X_k_norm,DEG_list,y,10000)$AUC_PR
perform_raw[is.nan(perform_raw)] <- 0
auprc_norm <- c(auprc_norm,mean(perform_raw))
}

auprc_norm2 <- NULL
for(alpha in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)){
ENIGMA_l2max <- cell_deconvolve(X=as.matrix(Bulk),
                                        theta=Frac$theta,
                                        R=Reference,
                                        epsilon=0.001,
                                        alpha=alpha,
                                        beta=0.5,tao_k=0.01,max.iter=1000,verbose=FALSE,pos=FALSE,Norm.method = "PC")
										
result <- DEG_test(ENIGMA_l2max$X_k_norm,y)
###Evaluation
Tab <- result$pval
perform_raw <- DePRCcalculator(Tab,ENIGMA_l2max$X_k_norm,DEG_list,y,10000)$AUC_PR
perform_raw[is.nan(perform_raw)] <- 0
auprc_norm2 <- c(auprc_norm2,mean(perform_raw))
}

corDat = data.frame(cor = as.numeric(t(cor_gene)),alpha = rep(c(1:9)*0.1,5),celltype = c(rep("celltype-1",9),rep("celltype-2",9),rep("celltype-3",9),rep("celltype-4",9),rep("celltype-5",9)))

pp=ggplot(corDat,aes(x=alpha,y=cor,fill=celltype,colour=celltype))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        panel.grid=element_blank())+
  ylab("Correlation with cell type fractions")

png("cor_celltype_fraction_deg_sim.png",res=300,height = 1600,width = 2200)
pp
dev.off()

performDat = data.frame(ARI = c(auprc,auprc_norm,auprc_norm2),alpha = rep(c(1:9)*0.1,3),method = c(rep("Before Normalization",9),rep("After Normalization(Frac)",9),rep("After Normalization(PC)",9)))

pp=ggplot(performDat,aes(x=alpha,y=ARI,fill=method,colour=method))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        panel.grid=element_blank())+
  ylab("AUPRC")

png("AUPRC_l2norm.png",res=300,height = 1600,width = 2200)
pp
dev.off()
##############################

cor_gene_trace <- NULL
auprc_trace <- auprc_norm_trace <- NULL
for(alpha in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)){
ENIGMA_trace <- cell_deconvolve_trace(O = as.matrix(Bulk),
                                              theta=Frac$theta,
                                              R=Reference,
                                              alpha=alpha,beta=1,solver="adaptive_admm",
                                              verbose=FALSE,max.iter = 1000,pos=FALSE,Norm.method = "frac",pre.process="none")
										
cor_per_ct <- NULL
for(k in 1:5){
   exp <- ENIGMA_trace$X_k[,,k]
   exp.scale <- t(apply(exp,1,scale))
   PC <- as.matrix(svd(exp.scale)$v[,1:10])
   cor <- apply(PC,2,function(x){cor(x,Frac$theta[,k],method="sp")})
   cor_per_ct <- c(cor_per_ct,max(abs(cor)))
}			
cor_gene_trace <- cbind(cor_gene_trace, cor_per_ct)

result <- DEG_test(ENIGMA_trace$X_k,y)
###Evaluation
Tab <- result$pval
perform_raw <- DePRCcalculator(Tab,ENIGMA_trace$X_k,DEG_list,y,10000)$AUC_PR
perform_raw[is.nan(perform_raw)] <- 0
auprc_trace <- c(auprc_trace,mean(perform_raw))


result <- DEG_test(ENIGMA_trace$X_k_norm,y)
###Evaluation
Tab <- result$pval
perform_raw <- DePRCcalculator(Tab,ENIGMA_trace$X_k_norm,DEG_list,y,10000)$AUC_PR
perform_raw[is.nan(perform_raw)] <- 0
auprc_norm_trace <- c(auprc_norm_trace,mean(perform_raw))
}

auprc_norm2_trace <- NULL
for(alpha in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)){
ENIGMA_trace <- cell_deconvolve_trace(O = as.matrix(Bulk),
                                              theta=Frac$theta,
                                              R=Reference,
                                              alpha=alpha,beta=1,solver="admm",
                                              verbose=FALSE,max.iter = 1000,pos=FALSE,Norm.method = "PC",pre.process="none")
										
result <- DEG_test(ENIGMA_trace$X_k_norm,y)
###Evaluation
Tab <- result$pval
perform_raw <- DePRCcalculator(Tab,ENIGMA_trace$X_k_norm,DEG_list,y,10000)$AUC_PR
perform_raw[is.nan(perform_raw)] <- 0
auprc_norm2_trace <- c(auprc_norm2_trace,mean(perform_raw))
}

corDat = data.frame(cor = as.numeric(t(cor_gene_trace)),alpha = rep(c(1:9)*0.1,5),celltype = c(rep("celltype-1",9),rep("celltype-2",9),rep("celltype-3",9),rep("celltype-4",9),rep("celltype-5",9)))

pp=ggplot(corDat,aes(x=alpha,y=cor,fill=celltype,colour=celltype))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        panel.grid=element_blank())+
  ylab("Correlation with cell type fractions")

png("/Path/to/Data/cor_celltype_fraction_deg_sim(trace).png",res=300,height = 1600,width = 2200)
pp
dev.off()

performDat = data.frame(ARI = c(auprc_trace,auprc_norm_trace,auprc_norm2_trace),alpha = rep(c(1:9)*0.1,3),method = c(rep("Before Normalization",9),rep("After Normalization(Frac)",9),rep("After Normalization(PC)",9)))

pp=ggplot(performDat,aes(x=alpha,y=ARI,fill=method,colour=method))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        panel.grid=element_blank())+
  ylab("AUPRC")

png("/Path/to/Data/AUPRC_trace.png",res=300,height = 1600,width = 2200)
pp
dev.off()
