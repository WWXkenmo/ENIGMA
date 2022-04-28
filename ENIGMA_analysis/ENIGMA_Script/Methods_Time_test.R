####TIME
setwd("/mnt/data2/zhouxl/ENIGMA")
options (warn = -1)
library(dplyr)
source("/mnt/data2/zhouxl/ENIGMA/ENIGMA2.R")

Bulk=readRDS("./pbmc sample/bulk_var_500s_fivecellty.rds")
Bulk=as.matrix(Bulk)

ref_seqwell=readRDS("./data/ref_seqwell_rmbe_500s.rds")
ref=ref_seqwell$cell.type.coarse
inde=c("B","Mono","NK","T_CD4","T_CD8")
ref=ref[,inde]
profile=ref

Frac=readRDS("./pbmc sample/Fra_Simulate_seqwell_rmbe_500s_fivecellty.rds")
Frac$theta=Frac$theta[,inde]
#### 1.bMIND
library(MIND) 
source("/mnt/data2/zhouxl/ENIGMA2/bMIND.R")#Remove the log preprocessing inside the function
### 1.1 sample 100
set.seed(1)
samplegene=sample(rownames(Bulk),100)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t=NULL
t1=system.time({deconv100_2 = bMIND(log2(Bulk1+1),ncore = 2,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)
t1=system.time({deconv100_10 = bMIND(log2(Bulk1+1),ncore = 10,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)
t1=system.time({deconv100_20 = bMIND(log2(Bulk1+1),ncore = 20,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)

### 1.2 sample 500
set.seed(2)
samplegene=sample(rownames(Bulk),500)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t1=system.time({deconv500_2 = bMIND(log2(Bulk1+1),ncore = 2,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)
t1=system.time({deconv500_10 = bMIND(log2(Bulk1+1),ncore = 10,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)
t1=system.time({deconv500_20 = bMIND(log2(Bulk1+1),ncore = 20,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)

### 1.3 sample 1000
set.seed(3)
samplegene=sample(rownames(Bulk),1000)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t1=system.time({deconv1000_2 = bMIND(log2(Bulk1+1),ncore = 2,frac = Frac$theta,profile =log2(profile1+1))})
t=rbind(t,t1)
t1=system.time({deconv1000_10 = bMIND(log2(Bulk1+1),ncore = 10,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)
t1=system.time({deconv1000_20 = bMIND(log2(Bulk1+1),ncore = 20,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)

### 1.4 sample 5000
set.seed(4)
samplegene=sample(rownames(Bulk),5000)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t1=system.time({deconv5000_2 = bMIND(log2(Bulk1+1),ncore = 2,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)
t1=system.time({deconv5000_10 = bMIND(log2(Bulk1+1),ncore = 10,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)
t1=system.time({deconv5000_20 = bMIND(log2(Bulk1+1),ncore = 20,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)

### 1.5 sample 10000
set.seed(5)
samplegene=sample(rownames(Bulk),10000)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t1=system.time({deconv10000_2 = bMIND(log2(Bulk1+1),ncore = 2,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)
t1=system.time({deconv10000_10 = bMIND(log2(Bulk1+1),ncore = 10,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)
t1=system.time({deconv10000_20 = bMIND(log2(Bulk1+1),ncore = 20,frac = Frac$theta,profile = log2(profile1+1))})
t=rbind(t,t1)
t_bMIND=t
rownames(t_bMIND)=paste0(rep(c("sample100","sample500","sample1000","sample5000","sample10000"),each=3),
                         "-",rep(c(2,10,20),5))
save(t_bMIND,file="/mnt/data2/zhouxl/ENIGMA2/Covid19/Seqwell/Time/Time_bMIND.Rdata")

#### 2. TCA
library(TCA)
### 2.1 sample 100
set.seed(1)
samplegene=sample(rownames(Bulk),100)
Bulk1=Bulk[samplegene,] %>% as.matrix()

t=NULL
t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=2,max_iter = 10)
Z_hat2 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=10,max_iter = 10)
Z_hat10 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=20,max_iter = 10)
Z_hat20 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

### 2.2 sample 500
set.seed(2)
samplegene=sample(rownames(Bulk),500)
Bulk1=Bulk[samplegene,] %>% as.matrix()

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=2,max_iter = 10)
Z_hat2 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=10,max_iter = 10)
Z_hat10 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=20,max_iter = 10)
Z_hat20 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

### 2.3 sample 1000
set.seed(3)
samplegene=sample(rownames(Bulk),1000)
Bulk1=Bulk[samplegene,] %>% as.matrix()

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=2,max_iter = 10)
Z_hat2 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=10,max_iter = 10)
Z_hat10 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=20,max_iter = 10)
Z_hat20 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

### 2.4 sample 5000
set.seed(4)
samplegene=sample(rownames(Bulk),5000)
Bulk1=Bulk[samplegene,] %>% as.matrix()

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=2,max_iter = 10)
Z_hat2 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=10,max_iter = 10)
Z_hat10 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=20,max_iter = 10)
Z_hat20 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

### 2.4 sample 10000
set.seed(5)
samplegene=sample(rownames(Bulk),10000)
Bulk1=Bulk[samplegene,] %>% as.matrix()

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=2,max_iter = 10)
Z_hat2 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=10,max_iter = 10)
Z_hat10 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)

t1=system.time({tca.mdl <- tca(X = as.matrix(sqrt(Bulk1)), W = Frac$theta, C1 = NULL, C2 = NULL,
                               parallel = TRUE,num_cores=20,max_iter = 10)
Z_hat20 <- tensor(X = as.matrix(sqrt(Bulk1)), tca.mdl)})
t=rbind(t,t1)
print("Sample 10000 Done")

t_TCA=t
rownames(t_TCA)=paste0(rep(c("sample100","sample500","sample1000","sample5000","sample10000"),each=3),
                       "-",rep(c(2,10,20),5))
save(t_TCA,file="/mnt/data2/zhouxl/ENIGMA2/Covid19/Seqwell/Time/Time_TCA.Rdata")


#### 3. ENIGMA
source("/mnt/data2/zhouxl/ENIGMA/ENIGMA2.R")
### 3.1 sample 100
set.seed(1)
samplegene=sample(rownames(Bulk),100)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

alpha=0.9
t=NULL
t1=system.time({res_alg_all <- cell_deconvolve(X=sqrt(Bulk1),
                                               theta=Frac$theta,
                                               R=sqrt(profile1),
                                               epsilon=0.001,
                                               alpha=alpha,pre.process="sqrt",
                                               beta=0.1,tao_k=0.1,max.iter=1000,verbose=F,Normalize=F)})
t=rbind(t,t1)

### 3.2sample 500
set.seed(2)
samplegene=sample(rownames(Bulk),500)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t1=system.time({res_alg_all <- cell_deconvolve(X=sqrt(Bulk1),
                                               theta=Frac$theta,
                                               R=sqrt(profile1),
                                               epsilon=0.001,
                                               alpha=alpha,pre.process="sqrt",
                                               beta=0.1,tao_k=0.1,max.iter=1000,verbose=F,Normalize=F)})
t=rbind(t,t1)

### 3.3 sample 1000
set.seed(3)
samplegene=sample(rownames(Bulk),1000)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t1=system.time({res_alg_all <- cell_deconvolve(X=sqrt(Bulk1),
                                               theta=Frac$theta,
                                               R=sqrt(profile1),
                                               epsilon=0.001,
                                               alpha=alpha,pre.process="sqrt",
                                               beta=0.1,tao_k=0.1,max.iter=1000,verbose=F,Normalize=F)})
t=rbind(t,t1)

### 3.4 sample 5000
set.seed(4)
samplegene=sample(rownames(Bulk),5000)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t1=system.time({res_alg_all <- cell_deconvolve(X=sqrt(Bulk1),
                                               theta=Frac$theta,
                                               R=sqrt(profile1),
                                               epsilon=0.001,
                                               alpha=alpha,pre.process="sqrt",
                                               beta=0.1,tao_k=0.1,max.iter=1000,verbose=F,Normalize=F)})
t=rbind(t,t1)

### 3.5 sample 10000
set.seed(5)
samplegene=sample(rownames(Bulk),10000)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t1=system.time({res_alg_all <- cell_deconvolve(X=sqrt(Bulk1),
                                               theta=Frac$theta,
                                               R=sqrt(profile1),
                                               epsilon=0.001,
                                               alpha=alpha,pre.process="sqrt",
                                               beta=0.1,tao_k=0.1,max.iter=1000,verbose=F,Normalize=F)})
t=rbind(t,t1)
t_ENIGMA=t
rownames(t_ENIGMA)=c("sample100","sample500","sample1000","sample5000","sample10000")

save(t_ENIGMA,file="/mnt/data2/zhouxl/ENIGMA2/Covid19/Seqwell/Time/Time_ENIGMA_beta0.1.Rdata")

#### 4. ENIGMA-trace
### 4.1 sample 100
set.seed(1)
samplegene=sample(rownames(Bulk),100)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

alpha=0.9
t=NULL
t1=system.time({res_alg_trace<-cell_deconvolve_trace(O=sqrt(as.matrix(Bulk1)),
                                                      theta=Frac$theta,
                                                      R=sqrt(profile1),
                                                      epsilon=0.001,
                                                      alpha=alpha,beta=1,solver="admm",gamma=1,
                                                      verbose=F,max.iter = 500,Normalize=F,pos=T)})
t=rbind(t,t1)

### 4.2 sample 500
set.seed(2)
samplegene=sample(rownames(Bulk),500)
print(head(samplegene))
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t1=system.time({res_alg_trace<-cell_deconvolve_trace(O=sqrt(as.matrix(Bulk1)),
                                                     theta=Frac$theta,
                                                     R=sqrt(profile1),
                                                     epsilon=0.001,
                                                     alpha=alpha,beta=1,solver="admm",gamma=1,
                                                     verbose=F,max.iter = 500,Normalize=F,pos=T)})
t=rbind(t,t1)

### 4.3 sample 1000
set.seed(3)
samplegene=sample(rownames(Bulk),1000)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t1=system.time({res_alg_trace<-cell_deconvolve_trace(O=sqrt(as.matrix(Bulk1)),
                                                     theta=Frac$theta,
                                                     R=sqrt(profile1),
                                                     epsilon=0.001,
                                                     alpha=alpha,beta=1,solver="admm",gamma=1,
                                                     verbose=F,max.iter = 500,Normalize=F,pos=T)})
t=rbind(t,t1)

### 4.4 sample 5000
set.seed(4)
samplegene=sample(rownames(Bulk),5000)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t1=system.time({res_alg_trace<-cell_deconvolve_trace(O=sqrt(as.matrix(Bulk1)),
                                                     theta=Frac$theta,
                                                     R=sqrt(profile1),
                                                     epsilon=0.001,
                                                     alpha=alpha,beta=1,solver="admm",gamma=1,
                                                     verbose=F,max.iter = 500,Normalize=F,pos=T)})
t=rbind(t,t1)

### 4.5 sample 10000
set.seed(5)
samplegene=sample(rownames(Bulk),10000)
Bulk1=Bulk[samplegene,] %>% as.matrix()
profile1=profile[samplegene,]

t1=system.time({res_alg_trace<-cell_deconvolve_trace(O=sqrt(as.matrix(Bulk1)),
                                                     theta=Frac$theta,
                                                     R=sqrt(profile1),
                                                     epsilon=0.001,
                                                     alpha=alpha,beta=1,solver="admm",gamma=1,
                                                     verbose=F,max.iter = 500,Normalize=F,pos=T)})
t=rbind(t,t1)

t_ENIGMA_trace=t
rownames(t_ENIGMA_trace)=c("sample100","sample500","sample1000","sample5000","sample10000")

save(t_ENIGMA_trace,file="/mnt/data2/zhouxl/ENIGMA2/Covid19/Seqwell/Time/Time_ENIGMA_trace.Rdata")

#### 5.Plot
load("/mnt/data2/zhouxl/ENIGMA2/Covid19/Seqwell/Time/Time_ENIGMA_beta0.1.Rdata")
load("/mnt/data2/zhouxl/ENIGMA2/Covid19/Seqwell/Time/Time_ENIGMA_trace.Rdata")
load("/mnt/data2/zhouxl/ENIGMA2/Covid19/Seqwell/Time/Time_TCA.Rdata")
load("/mnt/data2/zhouxl/ENIGMA2/Covid19/Seqwell/Time/Time_bMIND.Rdata")

bMIND=data.frame(Times=t_bMIND[,3],Method=rep("bMIND",nrow(t_bMIND)),
                 Ncores=sapply(strsplit(rownames(t_bMIND),"-"), function(x){x[[2]]}),
                 Gene=rep(c(100,500,1000,5000,10000),each=3))
bMIND$Method=paste0(bMIND$Method," (ncores=",bMIND$Ncores,")")

TCA=data.frame(Times=t_TCA[,3],Method=rep("TCA",nrow(t_TCA)),
               Ncores=sapply(strsplit(rownames(t_TCA),"-"), function(x){x[[2]]}),
               Gene=rep(c(100,500,1000,5000,10000),each=3))
TCA$Method=paste0(TCA$Method," (ncores=",TCA$Ncores,")")

ENIGMA=data.frame(Times=t_ENIGMA_beta01[,3],Method=rep("ENIGMA",nrow(t_ENIGMA_beta01)),
                  Ncores=rep(1,5),
                  Gene=c(100,500,1000,5000,10000))

ENIGMA_t=data.frame(Times=t_ENIGMA_trace[,3],Method=rep("ENIGMA (trace)",nrow(t_ENIGMA_trace)),
                    Ncores=rep(1,5),
                    Gene=c(100,500,1000,5000,10000))

data=rbind(bMIND,TCA,ENIGMA,ENIGMA_t)
data$Method=factor(data$Method,levels = c("bMIND (ncores=2)","bMIND (ncores=10)","bMIND (ncores=20)",
                                          "TCA (ncores=2)","TCA (ncores=10)","TCA (ncores=20)",
                                          "ENIGMA (trace)","ENIGMA"))
data$Gene=factor(data$Gene,levels = c(100,500,1000,5000,10000))

p=ggplot(data,aes(x=Gene,y=Times,fill=Method))+
  geom_col(position = 'dodge')+theme_classic()+
  theme(legend.text = element_text(size=11),legend.title = element_text(size=14),
        axis.title= element_text(size=14),axis.text = element_text(size=12))+
  ylab("Running times (seconds)")+xlab("Number of genes")

png("/mnt/data2/zhouxl/ENIGMA2/Covid19/Seqwell/Time/Time_bar.png",res=300,height=1100,width=2400)
p+scale_fill_nejm()
dev.off()

