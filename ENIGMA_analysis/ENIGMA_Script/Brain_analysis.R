rm(list=ls())
gc()
#load the required packages
setwd("/path/to/Data/ENIGMA_documents")
load("/path/to/Data/brain_data.Rdata")
source("/path/to/Data/ENIGMA.R")
source("/path/to/Data/mean_fun.R")

library(MIND)
library(Biobase)
library(Seurat)
library(MASS)
library(Biobase)
library(TCA)
library(umap)
library(ggpubr)

################################## Comparison between multiple methods ###########################
##################################
#Using ENIGMA to perform analysis 
Frac <- get_proportion(Bulk,profile)
#ENIGMA-l2
res_alg_all <- cell_deconvolve(X=Bulk,
                               theta=Frac$theta,
                               R=profile,
                               epsilon=0.001,
                               alpha=0.1,pre.process="sqrt",
                               beta=0.5,tao_k=0.01,max.iter=1000,verbose=TRUE,Normalize=TRUE,Norm.method = "frac")

#ENIGMA-trace
res_alg_trace <- cell_deconvolve_trace(O=Bulk,
                                       theta=Frac$theta,
                                       R=profile,
                                       alpha=0.1,pre.process = "sqrt",
                                       verbose=TRUE,max.iter = 100,Norm.method = "frac")

####################################
#Using TCA to perform analysis 
tca.mdl <- tca(X = as.matrix(sqrt(Bulk)), W = Frac$theta, C1 = NULL, C2 = NULL,
               parallel = TRUE,num_cores=3,max_iter = 10)
Z_hat <- tensor(X = as.matrix(sqrt(Bulk)), tca.mdl)

####################################
#Using bMIND to perform analysis
cell_type <- unique(pData(scExp)$cell_type)
K <- length(cell_type)
sc <- exprs(scExp)
profile = matrix(NA, nrow(sc), K)
rownames(profile) = rownames(sc)
colnames(profile) = cell_type
for (i in cell_type) {
  profile[, i] = log2(edgeR::cpm(rowMeans(sc[, pData(scExp)$cell_type ==
                                               i])) + 1)
}
profile <- profile[rownames(Bulk),]
frac = est_frac(sig = profile, bulk = log2(1+Bulk))
deconv3 = bMIND(Bulk,ncore = 7, profile = as.matrix(profile),frac=frac)

#####################
#Benchmark
HiDe <- NULL
HiDe2 <- NULL
bMIND <- NULL
TCA <- NULL
celltype <- NULL
bulk_exp <- NULL
for(ct in colnames(profile)){
  B_cell_pre2 <- deconv3$A[,ct,colnames(Bulk)]
  B_cell <- CSE_array[ct,,colnames(Bulk)]
  B_cell_pre <- res_alg_all$X_k[,colnames(Bulk),ct] 
  B_cell_pre_trace <- res_alg_trace$X_k[,,ct]
  B_cell_pre3 <- Z_hat[[which(colnames(profile) == ct)]][,colnames(Bulk)] 
  cor_pre <- NULL
  cor_pre2 <- NULL
  cor_pre3 <- NULL
  cor_pre4 <- NULL
  cor_bulk <- NULL
  for(i in rownames(B_cell_pre_trace)){
    cor_pre <- c(cor_pre, cor(B_cell[i,],B_cell_pre[i,],method="sp"))
    cor_pre2 <- c(cor_pre2, cor(B_cell[i,],B_cell_pre2[i,],method="sp"))
    cor_pre3 <- c(cor_pre3, cor(B_cell[i,],B_cell_pre3[i,],method="sp"))
    cor_pre4 <- c(cor_pre4, cor(B_cell[i,],B_cell_pre_trace[i,],method="sp"))
    cor_bulk <- c(cor_bulk, cor(B_cell[i,],as.matrix(Bulk_normalize)[i,],method="sp"))
  }
  HiDe <- c(HiDe,cor_pre)
  bMIND <- c(bMIND,cor_pre2)
  TCA <- c(TCA,cor_pre3)
  HiDe2 <- c(HiDe2,cor_pre4)
  bulk_exp <- c(bulk_exp, cor_bulk)
  celltype <- c(celltype,rep(ct,length(cor_pre4)))
}

data_g <- data.frame(method=c(rep("TCA",length(TCA)),
                              rep("bMIND",length(bMIND)),
                              rep("ENIGMA",length(HiDe)),
                              rep("ENIGMA (trace)",length(HiDe2)),
                              rep("bulk",length(bulk_exp))),
                     performance=c(TCA,bMIND,HiDe,HiDe2,bulk_exp),
                     celltype=c(rep(celltype,5)))

data_g %>% group_by(method) %>% summarise(m=mean(performance,na.rm=T))
data_g=data_g %>% filter(method!="bulk")
data_g$method=factor(data_g$method,levels = c("TCA","bMIND","ENIGMA (trace)","ENIGMA"))
ggplot(data_g, aes(x=celltype, y=performance, fill=method)) +
  geom_boxplot()+theme_minimal()+labs(y="Correlation per gene")+facet_grid(~celltype, scales = "free_x") +
  theme(legend.text = element_text(size=12),legend.title = element_text(size=14),legend.position = "top",
        axis.title.y = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size=14),title = element_blank()) + scale_fill_jama(alpha=0.8)+
  theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))+labs(fill="Method")


## per_sample
HiDe <- NULL
HiDe2 <- NULL
bMIND <- NULL
TCA <- NULL
celltype <- NULL
bulk_exp <- NULL
reference <- NULL
for(ct in colnames(profile)){
  B_cell_pre2 <- deconv3$A[,ct,]
  B_cell_pre <- res_alg_all$X_k[,,ct]
  B_cell <- CSE_array[ct,rownames(B_cell_pre),]
  B_cell_pre_trace <- res_alg_trace$X_k[,,ct]
  B_cell_pre3 <- Z_hat[[which(colnames(profile) == ct)]]
  
  cor_pre <- NULL
  cor_pre2 <- NULL
  cor_pre3 <- NULL
  cor_pre4 <- NULL
  cor_bulk <- NULL
  cor_ref <- NULL
  
  for(i in 1:ncol(Bulk)){
    cor_pre <- c(cor_pre, cor(B_cell[,i],B_cell_pre[,i],method="sp"))
    cor_pre2 <- c(cor_pre2, cor(B_cell[,i],B_cell_pre2[,i],method="sp"))
    cor_pre3 <- c(cor_pre3, cor(B_cell[,i],B_cell_pre3[,i],method="sp"))
    cor_pre4 <- c(cor_pre4, cor(B_cell[,i],B_cell_pre_trace[,i],method="sp"))
    cor_bulk <- c(cor_bulk, cor(B_cell[,i],as.matrix(Bulk_normalize)[rownames(B_cell),i],method="sp"))
    cor_ref <- c(cor_ref, cor(B_cell[,i],profile[rownames(B_cell),ct],method="sp"))
  }
  
  HiDe <- c(HiDe,cor_pre)
  bMIND <- c(bMIND,cor_pre2)
  TCA <- c(TCA,cor_pre3)
  HiDe2 <- c(HiDe2,cor_pre4)
  bulk_exp <- c(bulk_exp, cor_bulk)
  reference <- c(reference, cor_ref)
  celltype <- c(celltype,rep(ct,length(cor_pre)))
}

data <- data.frame(method=c(rep("TCA",length(TCA)),
                            rep("bMIND",length(bMIND)),
                            rep("ENIGMA",length(HiDe)),
                            rep("ENIGMA (trace)",length(HiDe2)),
                            rep("bulk",length(bulk_exp))),
                   performance=c(TCA,bMIND,HiDe,HiDe2,bulk_exp),
                   celltype=c(rep(celltype,5)))
data=data %>% filter(method!="bulk")

data$method=factor(data$method,levels = c("TCA","bMIND","ENIGMA (trace)","ENIGMA"))
data %>% group_by(method) %>% summarise(m=mean(performance,na.rm=T))

ggplot(data, aes(x=celltype, y=performance, fill=method)) +
  geom_boxplot()+theme_minimal()+labs(y="Correlation per sample")+facet_grid(~celltype, scales = "free_x") +
  theme(legend.text = element_text(size=12),legend.title = element_text(size=14),legend.position = "top",
        axis.title.y = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size=14),title = element_blank()) + scale_fill_jama(alpha=0.8)+
  theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))+labs(fill="Method")


################################## alpha parameter test ###########################
alpha.v <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
ENIGMA_brain= NULL
for(k in 1:length(alpha.v)){
  alpha<-alpha.v[k]
  print(paste0("alpha=",alpha))
  ENIGMA_trace.v <- cell_deconvolve_trace(O = as.matrix(Bulk),
                                          theta=Frac$theta,
                                          R=profile,pre.process="log",
                                          alpha=alpha,Normalize=T,epsilon=0.0001,Norm.method = "frac",
                                          verbose=T)
  ENIGMA_brain[[k]]=ENIGMA_trace.v
}

enigma.mean=NULL
for (i in 1:length(ENIGMA_brain)) {
  ENIGMA_trace.v = ENIGMA_brain[[i]]
  alpha<-alpha.v[i]
  
  enigma.mean1=mean_fun(cell_decon_result=ENIGMA_trace.v,Fraction=Frac$theta,groundtrue=single_celltype)
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean1)=paste(alpha_name,rownames(enigma.mean1),sep="-")
  enigma.mean=rbind(enigma.mean,enigma.mean1)
}

## plot
cor_alpha=enigma.mean.norm
cor_alpha=cor_alpha %>% rownames_to_column("name")
cor_alpha$celltype=sapply(strsplit(cor_alpha$name,"-"),function(x){x[[2]]})
cor_alpha$alpha=rep(alpha.v[1:length(ENIGMA_brain)],each=5)
#per sample
ggplot(cor_alpha,aes(x=alpha,y=as.numeric(s),color=celltype))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14),
        panel.grid=element_blank())+
  labs(y="Correlation per sample",color="Cell type")+ylim(c(0,max(cor_alpha$s)+0.2))+
  scale_color_npg(alpha=0.9)+scale_x_continuous(limits = c(0.1,0.9),breaks = seq(0,0.9,0.1))

#per gene
ggplot(cor_alpha,aes(x=alpha,y=as.numeric(g),color=celltype))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14),
        panel.grid=element_blank())+
  labs(y="Correlation per gene",color="Cell type")+ylim(c(0,max(cor_alpha$g)+0.2))+
  scale_color_npg(alpha=0.9)+scale_x_continuous(limits = c(0.1,0.9),breaks = seq(0,0.9,0.1))#+theme(legend.position = c(0.85,0.25))





################################## Comparison between different parameters (alpha) ###########################
single_celltype=NULL
for (i in colnames(profile)) {
  single_celltype[[i]]=CSE_array[i,,]
  single_celltype[[i]]=single_celltype[[i]][rownames(Bulk),colnames(Bulk)]
}

Fra_Simulate <- get_proportion(Bulk, profile)

## trace
alpha.v <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
ENIGMA_brain= NULL
for(k in 1:length(alpha.v)){
  alpha<-alpha.v[k]
  print(paste0("alpha=",alpha))
  ENIGMA_trace.v <- cell_deconvolve_trace(O = as.matrix(Bulk),
                                          theta=Fra_Simulate$theta,
                                          R=profile,pre.process="log",
                                          alpha=alpha,Normalize=F,epsilon=0.0001,Norm.method = "frac",
                                          verbose=T)
  ENIGMA_brain[[k]]=ENIGMA_trace.v
}

enigma.mean=NULL
for (i in 1:length(ENIGMA_brain)) {
  ENIGMA_trace.v = ENIGMA_brain[[i]]
  alpha<-alpha.v[i]
  
  enigma.mean1=mean_fun(cell_decon_result=ENIGMA_trace.v,Fraction=Fra_Simulate$theta,groundtrue=single_celltype)
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean1)=paste(alpha_name,rownames(enigma.mean1),sep="-")
  enigma.mean=rbind(enigma.mean,enigma.mean1)
}
saveRDS(enigma.mean,"./enigma.mean.rds")

## L2
ENIGMA_brain_l2= NULL
enigma.mean=NULL
for(k in 1:length(alpha.v)){
  alpha<-alpha.v[k]
  print(paste0("alpha=",alpha))
  ENIGMA_trace.v <- cell_deconvolve(X = Bulk,
                                        theta=Fra_Simulate$theta,
                                        R=profile,pre.process="log",
                                        alpha=alpha,Normalize=F,Norm.method = "frac",
                                        verbose=T)
  ENIGMA_brain_l2[[k]]=ENIGMA_trace.v
  
  enigma.mean1=mean_fun(cell_decon_result=ENIGMA_trace.v,Fraction=Fra_Simulate$theta,groundtrue=single_celltype)
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean1)=paste(alpha_name,rownames(enigma.mean1),sep="-")
  enigma.mean=rbind(enigma.mean,enigma.mean1)
}
saveRDS(enigma.mean,"./enigma.mean_l2.rds")


##### Supplementary Figure S2b
enigma.mean=readRDS("./enigma.mean.rds") # or enigma.mean=readRDS("./enigma.mean_l2.rds")
cor_alpha=enigma.mean
cor_alpha=cor_alpha %>% rownames_to_column("name")
cor_alpha$celltype=sapply(strsplit(cor_alpha$name,"-"),function(x){x[[2]]})
cor_alpha$alpha=rep(alpha.v[1:length(ENIGMA_brain_l2)],each=5)

ggplot(cor_alpha,aes(x=alpha,y=as.numeric(s),color=celltype))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14),
        panel.grid=element_blank())+
  labs(y="Correlation per sample",color="Cell type")+ylim(c(0,max(cor_alpha$s)+0.2))+
  scale_color_npg(alpha=0.9)+scale_x_continuous(limits = c(0.1,0.9),breaks = seq(0,0.9,0.1))#+theme(legend.position = c(0.85,0.25))

ggplot(cor_alpha,aes(x=alpha,y=as.numeric(g),color=celltype))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14),
        panel.grid=element_blank())+
  labs(y="Correlation per gene",color="Cell type")+ylim(c(0,max(cor_alpha$g)+0.2))+
  scale_color_npg(alpha=0.9)+scale_x_continuous(limits = c(0.1,0.9),breaks = seq(0,0.9,0.1))#+theme(legend.position = c(0.85,0.25))


##### UMAP [Supplementary Figure S2a]
Frac <- Fra_Simulate
## L2
res_alg_all <- cell_deconvolve(X=Bulk,
                               theta=Frac$theta,
                               R=profile,
                               epsilon=0.001,
                               alpha=0.1,pre.process="sqrt",
                               beta=0.5,tao_k=0.01,max.iter=1000,verbose=TRUE,Normalize=TRUE,Norm.method = "frac")
umap_tr = umap(apply(res_alg_all$X_k, 1, as.vector))
ggscatter(data.frame(umap = umap_tr$layout, cell = rep(colnames(profile), each = ncol(Bulk))), x = "umap.1", y = "umap.2", color = "cell", size = 1) + ggtitle('Ground truth')

## trace
res_alg_trace <- cell_deconvolve_trace(O=Bulk,
                                       theta=Frac$theta,
                                       R=profile,
                                       epsilon=0.001,
                                       alpha=0.1,pre.process = "sqrt",
                                       verbose=TRUE,max.iter = 100,Norm.method = "frac")
umap_tr = umap(apply(res_alg_trace$X_k, 1, as.vector))
ggscatter(data.frame(umap = umap_tr$layout, cell = rep(colnames(profile), each = ncol(Bulk))), x = "umap.1", y = "umap.2", color = "cell", size = 1) + ggtitle('Ground truth')












