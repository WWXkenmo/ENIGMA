#################################################
source("/path/to/save/Data/ENIGMA.R")
######load the cell states infor
Bulk <- readRDS("/path/to/save/Data/Bulk.rds")
Reference <- readRDS("/path/to/save/Data/Reference.rds")
label <- readRDS("/path/to/save/Data/CellLabel.rds")

Frac_Simulate <- get_proportion(Bulk, Reference)

egm <- cell_deconvolve_trace(O = as.matrix(sqrt(Bulk)),
                                                theta=Frac_Simulate$theta,
                                                R=sqrt(Reference),
                                                alpha=0.1,solver ="adaptive_admm",
                                                verbose=TRUE,max.iter = 1000,pre.process = "sqrt",Norm.method = "PC")

enigma_trace2 <- SingleCellExperiment(assays=list(logcounts = egm$X_k_norm[,,1]))
enigma_trace2$Group <- label$c1				 
enigma_trace2 <- enigma_trace2[Frac_Simulate$theta[,1]>0.05]
enigma_trace2 <- runPCA(enigma_trace2)
enigma_trace2 <- runTSNE(enigma_trace2,dimred="PCA",n_dimred=5)
ARI_enigma <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(enigma_trace2,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(enigma_trace2$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_enigma <- c(ARI_enigma,ARI)
}
p1 <- plotTSNE(enigma_trace2, colour_by="Group",point_size=3)


enigma_trace2 <- SingleCellExperiment(assays=list(logcounts = egm$X_k_norm[,,3]))
enigma_trace2$Group <- label$c3				 
enigma_trace2 <- enigma_trace2
enigma_trace2 <- runPCA(enigma_trace2)
enigma_trace2 <- runTSNE(enigma_trace2,dimred="PCA",n_dimred=7)
ARI_enigma <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(enigma_trace2,"PCA")[,1:i],2)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(enigma_trace2$Group)),as.numeric(clust$cluster))
print(ARI)
ARI_enigma <- c(ARI_enigma,ARI)
}
p2 <- plotTSNE(enigma_trace2, colour_by="Group",point_size=3)
png("/path/to/save/Data/CellType1.png",res=300,height=1500,width=1500)
p1
dev.off()
png("/path/to/save/Data/CellType3.png",res=300,height=1500,width=1500)
p2
dev.off()
#######################################################
###Benchmark with other algorithm

ES <- c(1.8,2.4,3,3.6,4.2,4.8)
for(es in 1:length(ES)){
load(paste("/path/to/save/Data/resCompare_bi_",ES[es],".Rdata",sep=""))
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
                parallel = TRUE,verbose=FALSE,num_cores=3)
Z_hat_simulate <- tensor(X = (as.matrix(G)), tca.mdl,verbose=FALSE)

egm <- cell_deconvolve_trace(O = as.matrix(sqrt(Bulk)),
                                                theta=Frac_Simulate$theta,
                                                R=sqrt(Reference),
                                                alpha=0.1,solver ="adaptive_admm",
                                                verbose=TRUE,max.iter = 1000,pre.process = "sqrt",Norm.method = "PC")

p <- 100
y <- gl(2, p/2)
res_tca <- DEG_test1(Z_hat_simulate,y,Fra_Simulate$theta,method = "TCA",10000,DEG_list_all[[rep]],qval=TRUE)
res_enigma <- DEG_test1(ENIGMA_trace$X_k,y,Fra_Simulate$theta,method = "enigma",10000,DEG_list_all[[rep]],qval=TRUE)
line <- c(res_tca,res_enigma)
AUPRC <- rbind(AUPRC,line)
}

save(AUPRC,file=paste("/Path/to/save/",ES[es],"_AUPRC.Rdata",sep=""))
}

#####draw the boxplot
setwd("/Path/to/save/")
file <- list.files(pattern=".Rdata")
Tab <- NULL
for(i in 1:length(file)){
load(file[i])
ES <- gsub("_AUPRC.Rdata","",file[i])
T <- cbind(as.numeric(AUPRC),c(rep("TCA",10),rep("ENIGMA",10)),rep(ES,20))
Tab <- rbind(Tab,T)
}
Tab <- as.data.frame(Tab)
colnames(Tab) <- c("AUPRC","Method","EffectiveSize")
mytheme <- readRDS("/path/to/save/Data/mytheme.rds")
Tab$AUPRC <- as.numeric(as.matrix(Tab$AUPRC))
p_boxplot <- Tab %>% 
    mutate(EffectiveSize=paste0("SNR=", EffectiveSize)) %>% 
    ggplot(aes(Method, AUPRC, color=Method)) + 
    geom_boxplot() + 
    facet_grid(~EffectiveSize) + 
    mytheme + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
    theme(legend.title = element_text(size = 12, color = "black", family = "Arial"), legend.text = element_text(size = 12, color = "black", family = "Arial")) + 
    theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))
png("/path/to/save/Data/boxplot_DEG.png",res=300,height=1200,width=2400)
p_boxplot
dev.off()







