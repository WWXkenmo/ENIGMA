#!/usr/bin/env Rscript
########### 10X platform 
#### 1.Caculate results
#### 2.Caculate Correlation
#### 3.Caculate RMSE
#### 4.Plot

### load data
setwd("/path/to/Data/")
gt=readRDS("single_celltyoe_500s_fivecellty.rds")
Bulk=readRDS("bulk_var_500s_fivecellty.rds")
Bulk=as.matrix(Bulk)
Frac=readRDS("Fra_Simulate_same10x_rmbe_500s_fivecellty.rds")
ref=readRDS("ref_same10x_rmbe_500s.rds")
ref=ref$majorType
inde=c("B","Mono","T_CD8","T_CD4","NK")
ref=ref[,inde]
profile=ref

#### 1.Caculate results
## 1.1 bMIND results
library(MIND)
source("/mnt/data2/zhouxl/ENIGMA2/bMIND.R")##Remove the log preprocessing inside the function
deconv3 = bMIND(sqrt(Bulk),ncore = 7, frac = Frac$theta, profile = as.matrix(sqrt(profile)))
deconv4 = bMIND(log2(Bulk+1),ncore = 7, frac = Frac$theta, profile = as.matrix(log2(profile+1)))
save(deconv3,deconv4,file="./bMIND_same10x_rmlog2.Rdata")

## 1.2 TCA results
library(TCA)
tca.mdl <- tca(X = as.matrix(sqrt(Bulk)), W = Frac$theta, C1 = NULL, C2 = NULL,
               parallel = TRUE,num_cores=3,max_iter = 10)
Z_hat2 <- tensor(X = as.matrix(sqrt(Bulk)), tca.mdl)

tca.mdl <- tca(X = as.matrix(log2(Bulk+1)), W = Frac$theta, C1 = NULL, C2 = NULL,
               parallel = TRUE,num_cores=3,max_iter = 10)
Z_hat3 <- tensor(X = as.matrix(log2(Bulk+1)), tca.mdl)
save(Z_hat2,Z_hat3,file="./TCA_same10x.Rdata")

## 1.3 ENIGMA results
source("/mnt/data2/zhouxl/ENIGMA/ENIGMA2.R")
alpha=0.1
res_alg_all3 <- cell_deconvolve(X=sqrt(Bulk),
                                theta=Frac$theta,
                                R=sqrt(profile),
                                epsilon=0.001,
                                alpha=alpha,pre.process="sqrt",
                                beta=0.5,tao_k=0.1,max.iter=1000,verbose=TRUE,Normalize=TRUE,Norm.method = "frac")

res_alg_all4 <- cell_deconvolve(X=log2(Bulk+1),
                                theta=Frac$theta,
                                R=log2(profile+1),
                                epsilon=0.001,
                                alpha=alpha,pre.process="log",
                                beta=0.5,tao_k=0.1,max.iter=1000,verbose=TRUE,Normalize=TRUE,Norm.method = "frac")
save(res_alg_all3,res_alg_all4,file = "./ENIGMA_l2_same10x.Rdata")

## 1.4 ENIGMA-trace results
alpha=0.1
res_alg_trace3<-cell_deconvolve_trace(O=log2(as.matrix(Bulk)+1),
                                      theta=Frac$theta,
                                      R=log2(profile+1),
                                      epsilon=0.001,
                                      alpha=alpha,beta=1,solver="admm",gamma=1,
                                      verbose=F,max.iter = 500,Normalize=T,pos=T,Norm.method = "frac")

res_alg_trace4<-cell_deconvolve_trace(O=sqrt(as.matrix(Bulk)),
                                      theta=Frac$theta,
                                      R=sqrt(profile),
                                      epsilon=0.001,
                                      alpha=alpha,beta=1,solver="admm",gamma=1,
                                      verbose=F,max.iter = 500,Normalize=T,pos=T,Norm.method = "frac")
save(res_alg_trace3,res_alg_trace4,file = "./ENIGMA_trace_same10x_admm.Rdata")

#### 2.Caculate Correlation
## 2.1 Correlation per gene
# 2.1.1 TCA
load("./TCA_same10x.Rdata")
tca_res2=Z_hat2
tca_res3=Z_hat3

TCA2 <- NULL
TCA3 <- NULL
celltype <- NULL
bulk_exp <- NULL

for(ct in colnames(profile)){
  gt_base <- gt[[ct]][,colnames(Bulk)]
  tca2<- tca_res2[[which(colnames(profile) == ct)]][,colnames(Bulk)]
  tca3<- tca_res3[[which(colnames(profile) == ct)]][,colnames(Bulk)]
  
  cor_tca2 <- NULL
  cor_tca3 <- NULL
  cor_bulk <- NULL
  
  for(i in rownames(gt_base)){
    cor_tca2 <- c(cor_tca2, cor(gt_base[i,],tca2[i,],method="sp"))
    cor_tca3 <- c(cor_tca3, cor(gt_base[i,],tca3[i,],method="sp"))
    
    cor_bulk <- c(cor_bulk, cor(gt_base[i,],as.matrix(Bulk)[i,],method="sp"))
  }
  TCA2 <- c(TCA2,cor_tca2)
  TCA3 <- c(TCA3,cor_tca3)
  bulk_exp <- c(bulk_exp, cor_bulk)
  
  celltype <- c(celltype,rep(ct,length(cor_tca3)))
}

dat <- data.frame(method=c(rep("TCA(sqrt)",length(TCA2)),
                           rep("TCA(log)",length(TCA3)),
                           rep("bulk",length(bulk_exp))),
                  performance=c(TCA2,TCA3,bulk_exp),
                  celltype=c(rep(celltype,3)))
save(dat,file="./cor_tca_gene.Rdata")

# 2.1.2 bMIND
load("./bMIND_same10x_rmlog2.Rdata")
bmind_r0_res=deconv3
bmind_r1_res=deconv4

bMIND_raw0 <-NULL
bMIND_raw1 <-NULL
celltype <- NULL
for(ct in colnames(profile)){
  gt_base <- gt[[ct]][,colnames(Bulk)]##groundtrue
  bmind_r0 <- bmind_r0_res$A[,ct,colnames(Bulk)]
  bmind_r1 <- bmind_r1_res$A[,ct,colnames(Bulk)]
  
  cor_bmind_r0 <- NULL
  cor_bmind_r1 <- NULL
  
  for(i in rownames(gt_base)){
    cor_bmind_r0 <- c(cor_bmind_r0, cor(gt_base[i,],bmind_r0[i,],method="sp"))
    cor_bmind_r1 <- c(cor_bmind_r1, cor(gt_base[i,],bmind_r1[i,],method="sp"))
  }
  bMIND_raw0 <- c(bMIND_raw0,cor_bmind_r0)
  bMIND_raw1 <- c(bMIND_raw1,cor_bmind_r1)
  celltype <- c(celltype,rep(ct,length(cor_bmind_r0)))
}
dat <- data.frame(method=c(rep("bMIND(sqrt)",length(bMIND_raw0)),
                           rep("bMIND(log)",length(bMIND_raw1))),
                  performance=c(bMIND_raw0,bMIND_raw1),
                  celltype=c(rep(celltype,2)))
save(dat,file="./cor_bmind_rmlog2_gene.Rdata")

# 2.1.3 ENIGMA
load("./ENIGMA_l2_same10x.Rdata")
enigma_trace_res3=res_alg_all3
enigma_trace_res4=res_alg_all4

HiDe_t3 <- NULL
HiDe_t4 <- NULL
celltype<-NULL
for(ct in colnames(profile)){
  gt_base <- gt[[ct]][,colnames(Bulk)]##gt
  enigma_t3 <- enigma_trace_res3$X_k[,colnames(Bulk),ct]
  enigma_t4 <- enigma_trace_res4$X_k[,colnames(Bulk),ct]
  
  cor_enigma_t3 <- NULL
  cor_enigma_t4 <- NULL
  
  for(i in rownames(gt_base)){
    cor_enigma_t3 <- c(cor_enigma_t3, cor(gt_base[i,],enigma_t3[i,],method="sp"))
    cor_enigma_t4 <- c(cor_enigma_t4, cor(gt_base[i,],enigma_t4[i,],method="sp"))
  }
  
  HiDe_t3 <- c(HiDe_t3,cor_enigma_t3)
  HiDe_t4 <- c(HiDe_t4,cor_enigma_t4)
  celltype <- c(celltype,rep(ct,length(cor_enigma_t4)))
}

dat <- data.frame(method=c(rep("ENIGMA(sqrt)",length(HiDe_t3)),
                           rep("ENIGMA(log)",length(HiDe_t4))),
                  performance=c(HiDe_t3,
                                HiDe_t4),
                  celltype=c(rep(celltype,2)))
save(dat,file="./cor_ENIGAMA_l2_gene.Rdata")

# 2.1.4 ENIGMA-trace
load("./ENIGMA_trace_same10x_admm.Rdata")
enigma_trace_res3=res_alg_trace3
enigma_trace_res4=res_alg_trace4

HiDe_t3 <- NULL
HiDe_t4 <- NULL
celltype<-NULL
for(ct in colnames(profile)){
  gt_base <- gt[[ct]][,colnames(Bulk)]##gt
  enigma_t3 <- enigma_trace_res3$X_k[,colnames(Bulk),ct]
  enigma_t4 <- enigma_trace_res4$X_k[,colnames(Bulk),ct]
  
  cor_enigma_t3 <- NULL
  cor_enigma_t4 <- NULL
  
  for(i in rownames(gt_base)){
    cor_enigma_t3 <- c(cor_enigma_t3, cor(gt_base[i,],enigma_t3[i,],method="sp"))
    cor_enigma_t4 <- c(cor_enigma_t4, cor(gt_base[i,],enigma_t4[i,],method="sp"))
  }
  
  HiDe_t3 <- c(HiDe_t3,cor_enigma_t3)
  HiDe_t4 <- c(HiDe_t4,cor_enigma_t4)
  celltype <- c(celltype,rep(ct,length(cor_enigma_t4)))
}

dat <- data.frame(method=c(rep("ENIGMA-trace(sqrt)",length(HiDe_t3)),
                           rep("ENIGMA-trace(log)",length(HiDe_t4))),
                  performance=c(HiDe_t3,
                                HiDe_t4),
                  celltype=c(rep(celltype,2)))
save(dat,file="./cor_ENIGAMA_trace_gene_admm.Rdata")

## 2.2 Correlation per sample
# 2.2.1 TCA
load("./TCA_same10x.Rdata")
tca_res2=Z_hat2
tca_res3=Z_hat3

TCA2 <- NULL
TCA3 <- NULL
celltype <- NULL
bulk_exp <- NULL
for(ct in colnames(profile)){
  gt_base <- gt[[ct]][,colnames(Bulk)]##gt
  tca2<- tca_res2[[which(colnames(profile) == ct)]][,colnames(Bulk)]
  tca3<- tca_res3[[which(colnames(profile) == ct)]][,colnames(Bulk)]
  
  cor_tca2 <- NULL
  cor_tca3 <- NULL
  cor_bulk <- NULL
  
  for(i in colnames(gt_base)){
    cor_tca2 <- c(cor_tca2, cor(gt_base[,i],tca2[,i],method="sp"))
    cor_tca3 <- c(cor_tca3, cor(gt_base[,i],tca3[,i],method="sp"))
    
    cor_bulk <- c(cor_bulk, cor(gt_base[,i],as.matrix(Bulk)[,i],method="sp"))
  }
  TCA2 <- c(TCA2,cor_tca2)
  TCA3 <- c(TCA3,cor_tca3)
  bulk_exp <- c(bulk_exp, cor_bulk)
  
  celltype <- c(celltype,rep(ct,length(cor_tca3)))
}

dat <- data.frame(method=c(rep("TCA(sqrt)",length(TCA2)),
                           rep("TCA(log)",length(TCA3)),
                           rep("bulk",length(bulk_exp))),
                  performance=c(TCA2,TCA3,bulk_exp),
                  celltype=c(rep(celltype,3)))
save(dat,file="./cor_tca_sample.Rdata")

# 2.2.2 bMIND
load("./bMIND_same10x_rmlog2.Rdata")
bmind_r0_res=deconv3
bmind_r1_res=deconv4

bMIND_raw0 <-NULL
bMIND_raw1 <-NULL
celltype <- NULL
for(ct in colnames(profile)){
  gt_base <- gt[[ct]][,colnames(Bulk)]##groundtrue
  bmind_r0 <- bmind_r0_res$A[,ct,colnames(Bulk)]
  bmind_r1 <- bmind_r1_res$A[,ct,colnames(Bulk)]
  
  cor_bmind_r0 <- NULL
  cor_bmind_r1 <- NULL
  
  for(i in colnames(gt_base)){
    cor_bmind_r0 <- c(cor_bmind_r0, cor(gt_base[,i],bmind_r0[,i],method="sp"))
    cor_bmind_r1 <- c(cor_bmind_r1, cor(gt_base[,i],bmind_r1[,i],method="sp"))
  }
  bMIND_raw0 <- c(bMIND_raw0,cor_bmind_r0)
  bMIND_raw1 <- c(bMIND_raw1,cor_bmind_r1)
  celltype <- c(celltype,rep(ct,length(cor_bmind_r0)))
}
dat <- data.frame(method=c(rep("bMIND(sqrt)",length(bMIND_raw0)),
                           rep("bMIND(log)",length(bMIND_raw1))),
                  performance=c(bMIND_raw0,bMIND_raw1),
                  celltype=c(rep(celltype,2)))
save(dat,file="./cor_bmind_rmlog2_sample.Rdata")
print("bMIND cor done")

# 2.2.3 ENIGMA
load("./ENIGMA_l2_same10x.Rdata")
enigma_trace_res3=res_alg_all3
enigma_trace_res4=res_alg_all4

HiDe_t3 <- NULL
HiDe_t4 <- NULL

celltype<-NULL
for(ct in colnames(profile)){
  gt_base <- gt[[ct]][,colnames(Bulk)]##gt
  enigma_t3 <- enigma_trace_res3$X_k[,colnames(Bulk),ct]
  enigma_t4 <- enigma_trace_res4$X_k[,colnames(Bulk),ct]
  
  cor_enigma_t3 <- NULL
  cor_enigma_t4 <- NULL
  
  for(i in colnames(gt_base)){
    cor_enigma_t3 <- c(cor_enigma_t3, cor(gt_base[,i],enigma_t3[,i],method="sp"))
    cor_enigma_t4 <- c(cor_enigma_t4, cor(gt_base[,i],enigma_t4[,i],method="sp"))
  }
  HiDe_t3 <- c(HiDe_t3,cor_enigma_t3)
  HiDe_t4 <- c(HiDe_t4,cor_enigma_t4)
  celltype <- c(celltype,rep(ct,length(cor_enigma_t4)))
}

dat <- data.frame(method=c(rep("ENIGMA(sqrt)",length(HiDe_t3)),
                           rep("ENIGMA(log)",length(HiDe_t4))),
                  performance=c(HiDe_t3,
                                HiDe_t4),
                  celltype=c(rep(celltype,2)))
save(dat,file="./cor_ENIGAMA_l2_sample.Rdata")

# 2.2.4 ENIGMA-trace
load("./ENIGMA_trace_same10x_admm.Rdata")
enigma_trace_res3=res_alg_trace3
enigma_trace_res4=res_alg_trace4

HiDe_t3 <- NULL
HiDe_t4 <- NULL
celltype<-NULL
for(ct in colnames(profile)){
  gt_base <- gt[[ct]][,colnames(Bulk)]##gt
  enigma_t3 <- enigma_trace_res3$X_k[,colnames(Bulk),ct]
  enigma_t4 <- enigma_trace_res4$X_k[,colnames(Bulk),ct]
 
  cor_enigma_t3 <- NULL
  cor_enigma_t4 <- NULL
  for(i in colnames(gt_base)){
    cor_enigma_t3 <- c(cor_enigma_t3, cor(gt_base[,i],enigma_t3[,i],method="sp"))
    cor_enigma_t4 <- c(cor_enigma_t4, cor(gt_base[,i],enigma_t4[,i],method="sp"))
  }
  HiDe_t3 <- c(HiDe_t3,cor_enigma_t3)
  HiDe_t4 <- c(HiDe_t4,cor_enigma_t4)
  celltype <- c(celltype,rep(ct,length(cor_enigma_t4)))
}
dat <- data.frame(method=c(rep("ENIGMA-trace-a(log)",length(HiDe_t3)),
                           rep("ENIGMA-trace-a(sqrt)",length(HiDe_t4))),
                  performance=c(HiDe_t3,
                                HiDe_t4),
                  celltype=c(rep(celltype,2)))
save(dat,file="./cor_ENIGAMA_trace_sample_admm.Rdata")

#### 3.Caculate RMSE
source("/mnt/data2/zhouxl/ENIGMA2/RMSE_function.R")
## 3.1 TCA
load("./TCA_same10x.Rdata")
TCA_list <- list(tca_sqrt=Z_hat2,tca_log=Z_hat3)
process <- c("sqrt","log")
res_tca_s2 <- TCA_RMSE1(TCA_list,gt,Bulk,process,level="sample")
res_tca_g2 <- TCA_RMSE1(TCA_list,gt,Bulk,process,level="gene")
save(res_tca_s2,res_tca_g2,file="./TCA_rmlog2_RMSE.Rdata")

## 3.2 bMIND
load("./bMIND_same10x_rmlog2.Rdata")  
bmind_list <- list(bmind_sqrt=deconv3,bmind_log=deconv4)
process <- c("sqrt","log")
res_bmind_s2 <- bMIND_RMSE1(bmind_list,gt,Bulk,process,level="sample")
res_bmind_g2 <- bMIND_RMSE1(bmind_list,gt,Bulk,process,level="gene")
save(res_bmind_s2,res_bmind_g2,file="./bMIND_rmlog2_RMSE.Rdata")

## 3.3 ENIGMA
load("./ENIGMA_l2_same10x.Rdata")
load("./ENIGMA_trace_same10x_admm.Rdata")

ENIGMA_list <- list(l2_sqrt2=res_alg_all3,l2_log2=res_alg_all4,
                    trace_sqrt2=res_alg_trace4,trace_log2=res_alg_trace3)
process <- c("sqrt","log","sqrt","log")
res_enigma_s2 <- ENIGMA_RMSE1(ENIGMA_list,gt,Bulk,process,level="sample")
res_enigma_g2 <- ENIGMA_RMSE1(ENIGMA_list,gt,Bulk,process,level="gene")
save(res_enigma_s2,res_enigma_g2,file="./ENIGMA_rmlog2_RMSE_admm.Rdata")

## 3.4 reference
res_ref_s=ref_RMSE1(ref=profile,Bulk,gt,level="sample")
res_ref_g=ref_RMSE1(ref=profile,Bulk,gt,level="gene") 
save(res_ref_s,res_ref_g,file="./RMSE_samedepth.Rdata")


#### 4.Plot
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggsci)
## 4.1 Correlation
# 4.1.1 Correlation per gene
load("./cor_bmind_rmlog2_gene.Rdata")
data=dat
load(file="./cor_ENIGAMA_l2_gene.Rdata")
data=rbind(data,dat)
load(file="./cor_ENIGAMA_trace_gene_admm.Rdata")
data=rbind(data,dat)
load(file="./cor_tca_gene.Rdata")
data=rbind(data,dat)
data2=data %>% distinct(method,performance,celltype,.keep_all=T)

#best
data3=data2[which(data2$method %in% c("TCA(sqrt)","bMIND(log)","ENIGMA(sqrt)","ENIGMA-trace(sqrt)")),]
colnames(data3)[1]="Method"
data3$Method=factor(data3$Method,level=c("TCA(sqrt)","bMIND(log)","ENIGMA(sqrt)","ENIGMA-trace(sqrt)"))

p=ggplot(data3, aes(x=celltype, y=performance, fill=Method)) +
  geom_boxplot()+theme_minimal()+labs(y="Correlation per gene")+facet_grid(~celltype, scales = "free_x") +
  theme(legend.text = element_text(size=11),legend.title = element_text(size=14),
        axis.title.y = element_text(size=14),axis.text.y = element_text(size=12))+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))

png("./result_figs/fig2_cor_gene_best2.png",res=300,height=1000,width=2300)
p+scale_fill_lancet(alpha=0.9)
dev.off()

# 4.1.2 Correlation per sample
load("./cor_bmind_rmlog2_sample.Rdata")
data=dat
load(file="./cor_ENIGAMA_l2_sample.Rdata")
data=rbind(data,dat)
load(file="./cor_ENIGAMA_trace_sample_admm.Rdata")
data=rbind(data,dat)
load(file="./cor_tca_sample.Rdata")
data=rbind(data,dat)
data2=data %>% distinct(method,performance,celltype,.keep_all=T)

data3=data2[which(data2$method %in% c("TCA(sqrt)","bMIND(log)","ENIGMA(sqrt)","ENIGMA-trace(sqrt)")),]
colnames(data3)[1]="Method"
data3$Method=factor(data3$Method,level=c("TCA(sqrt)","bMIND(log)","ENIGMA(sqrt)","ENIGMA-trace(sqrt)"))

p=ggplot(data3, aes(x=celltype, y=performance, fill=Method)) +
  geom_boxplot()+theme_minimal()+labs(y="Correlation per sample")+facet_grid(~celltype, scales = "free_x") +
  theme(legend.text = element_text(size=11),legend.title = element_text(size=14),
        axis.title.y = element_text(size=14),axis.text.y = element_text(size=12))+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))
png("./result_figs/fig2_cor_sample_ENIGMA09.png",res=300,height=1000,width=2300)
p+scale_fill_lancet(alpha=0.9)
dev.off()

#### 4.2 RMSE
### RMSE=RSE(gene)/sqrt(ncol(Bulk));RMSE=RSE(sample)/sqrt(nrow(Bulk))
# 4.2.1 RMSE per gene
load("./RMSE_samedepth.Rdata")
load("./bMIND_rmlog2_RMSE.Rdata")
load("./TCA_rmlog2_RMSE.Rdata")
load("./ENIGMA_rmlog2_RMSE_admm.Rdata")

dat=reshape2::melt(res_ref_g)
data=data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("Reference",nrow(dat)))

dat=reshape2::melt(res_bmind_g2[[1]])
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("bMIND (sqrt)",nrow(dat))))
dat=reshape2::melt(res_bmind_g2[[2]])
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("bMIND (log)",nrow(dat))))

dat=reshape2::melt(res_tca_g2$tca_sqrt)
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("TCA (sqrt)",nrow(dat))))
dat=reshape2::melt(res_tca_g2$tca_log)
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("TCA (log)",nrow(dat))))

dat=reshape2::melt(res_enigma_g2$l2_sqrt2)
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("ENIGMA (sqrt)",nrow(dat))))
dat=reshape2::melt(res_enigma_g2$l2_log2)
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("ENIGMA (log)",nrow(dat))))

dat=reshape2::melt(res_enigma_g2$trace_sqrt2)
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("ENIGMA-trace (sqrt)",nrow(dat))))
dat=reshape2::melt(res_enigma_g2$trace_log2)
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("ENIGMA-trace (log)",nrow(dat))))

data$performance=data$performance/sqrt(ncol(Bulk))

#best
data2=data[data$Method %in% c("Reference","TCA (sqrt)","bMIND (log)","ENIGMA (sqrt)","ENIGMA-trace (sqrt)"),]
data2$Method=factor(data2$Method,level=c("Reference","TCA (sqrt)","bMIND (log)","ENIGMA (sqrt)","ENIGMA-trace (sqrt)"))
data2$celltype=factor(data2$celltype,level=c("B","Mono","NK","T_CD4","T_CD8"))

p=ggplot(data2, aes(x=celltype, y=performance, fill=Method)) +
  geom_boxplot()+theme_minimal()+labs(y="RSE per gene")+facet_grid(~celltype, scales = "free_x") +
  theme(legend.text = element_text(size=12),legend.title = element_text(size=15),
        axis.title.y = element_text(size=15),axis.text.y = element_text(size=12))+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))

png("./result_figs/fig2_rmse_gene_ENIGMA_trace2.png",res=300,height=1000,width=1900)
p+scale_fill_npg(alpha=0.9)
dev.off()


## 4.2.2 RMSE per sample
dat=reshape2::melt(res_ref_s)
data=data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("Reference",nrow(dat)))

dat=reshape2::melt(res_bmind_s2[[1]])
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("bMIND (sqrt)",nrow(dat))))
dat=reshape2::melt(res_bmind_s2[[2]])
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("bMIND (log)",nrow(dat))))

dat=reshape2::melt(res_tca_s2$tca_sqrt)
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("TCA (sqrt)",nrow(dat))))
dat=reshape2::melt(res_tca_s2$tca_log)
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("TCA (log)",nrow(dat))))

dat=reshape2::melt(res_enigma_s2$l2_sqrt2)
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("ENIGMA (sqrt)",nrow(dat))))
dat=reshape2::melt(res_enigma_s2$l2_log2)
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("ENIGMA (log)",nrow(dat))))

dat=reshape2::melt(res_enigma_s2$trace_sqrt2)
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("ENIGMA-trace (sqrt)",nrow(dat))))
dat=reshape2::melt(res_enigma_s2$trace_log2)
data=rbind(data,data.frame(celltype=dat$Var2,performance=dat$value,Method=rep("ENIGMA-trace (log)",nrow(dat))))

data$performance=data$performance/sqrt(nrow(Bulk))

#best
data2=data[data$Method %in% c("Reference","TCA (sqrt)","bMIND (log)","ENIGMA (sqrt)","ENIGMA-trace (sqrt)"),]
data2$Method=factor(data2$Method,level=c("Reference","TCA (sqrt)","bMIND (log)","ENIGMA (sqrt)","ENIGMA-trace (sqrt)"))
data2$celltype=factor(data2$celltype,level=c("B","Mono","NK","T_CD4","T_CD8"))

p=ggplot(data2, aes(x=celltype, y=performance, fill=Method)) +
  geom_boxplot()+theme_minimal()+labs(y="RMSE per sample")+facet_grid(~celltype, scales = "free_x") +
  theme(legend.text = element_text(size=12),legend.title = element_text(size=15),
        axis.title.y = element_text(size=15),axis.text.y = element_text(size=12))+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))

png("./result_figs/fig2_rmse_sample_best2.png",res=300,height=1000,width=2300)
p+scale_fill_npg(alpha=0.9)
dev.off()

