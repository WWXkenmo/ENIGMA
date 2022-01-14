#######ES_4.8  noise
setwd("/path/to/Data/ENIGMA_documents")
source("/path/to/Data/ENIGMA.R")
source("/path/to/Data/DEG_analysis_uile_function.R")

library(MASS)
p = 100
y <- gl(2, p/2)

## our simulated datasets are provided in https://github.com/WWXkenmo/ENIGMA/tree/main/ENIGMA_analysis/Data/cts_DEG
load("/path/to/Data/resCompare_bi_4.8.Rdata")
G <- testMatrixG[[1]]
H1 <- testMatrixH1[[1]]
set.seed(123)
noise0.1<-matrix(rnorm(5000)*0.1,nrow = 5,ncol = 1000)
noise0.5<-matrix(rnorm(5000)*0.5,nrow = 5,ncol = 1000)
noise1<-matrix(rnorm(5000)*1,nrow = 5,ncol = 1000)
noise1.5<-matrix(rnorm(5000)*1.5,nrow = 5,ncol = 1000)
noise2<-matrix(rnorm(5000)*2,nrow = 5,ncol = 1000)
noise2.5<-matrix(rnorm(5000)*2.5,nrow = 5,ncol = 1000)
noise3<-matrix(rnorm(5000)*3,nrow = 5,ncol = 1000)

H10.1<-H1+noise0.1
H10.5<-H1+noise0.5
H11<-H1+noise1
H11.5<-H1+noise1.5
H12<-H1+noise2
H12.5<-H1+noise2.5
H13<-H1+noise3

Fra_Simulate <- get_proportion(G, t(H1))
Fra_Simulate0.1 <- get_proportion(G, t(H10.1))
Fra_Simulate0.5 <- get_proportion(G, t(H10.5))
Fra_Simulate1 <- get_proportion(G, t(H11))
Fra_Simulate1.5 <- get_proportion(G, t(H11.5))
Fra_Simulate2 <- get_proportion(G, t(H12))
Fra_Simulate2.5 <- get_proportion(G, t(H12.5))
Fra_Simulate3 <- get_proportion(G, t(H13))


##admm
alpha.v <- c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99)
at=at0.1=at0.5=at1=at1.5=at2=at2.5=at3=NULL
for(alpha in alpha.v){
  ENIGMA_trace.v <- cell_deconvolve_trace(O = as.matrix(G),
                                          theta=Fra_Simulate$theta,
                                          R=t(H1),
                                          epsilon=0.00001,solver ="admm",
                                          alpha=alpha,gamma=0.5,beta=1,Normalize=F,
                                          verbose=FALSE,max.iter = 1000,pos = F)
  res_enigma.v<- DEG_test1(ENIGMA_trace.v$X_k,y,Fra_Simulate$theta,method = "enigma",10000,DEG_list_all[[1]],qval=FALSE)
  at<- c(at, res_enigma.v)
  
  ENIGMA_trace.v0.1 <- cell_deconvolve_trace(O = as.matrix(G),
                                             theta=Fra_Simulate0.1$theta,
                                             R=t(H10.1),
                                             epsilon=0.00001,solver ="admm",
                                             alpha=alpha,gamma=0.5,beta=1,Normalize=F,
                                             verbose=FALSE,max.iter = 1000,pos = F)
  res_enigma.v0.1 <- DEG_test1(ENIGMA_trace.v0.1$X_k,y,Fra_Simulate0.1$theta,method = "enigma",10000,DEG_list_all[[1]],qval=FALSE)
  at0.1 <- c(at0.1, res_enigma.v0.1)
  
  ENIGMA_trace.v0.5 <- cell_deconvolve_trace(O = as.matrix(G),
                                             theta=Fra_Simulate0.5$theta,
                                             R=t(H10.5),
                                             epsilon=0.00001,solver ="admm",
                                             alpha=alpha,gamma=0.5,beta=1,Normalize=F,
                                             verbose=FALSE,max.iter = 1000,pos = F)
  res_enigma.v0.5 <- DEG_test1(ENIGMA_trace.v0.5$X_k,y,Fra_Simulate0.5$theta,method = "enigma",10000,DEG_list_all[[1]],qval=FALSE)
  at0.5 <- c(at0.5, res_enigma.v0.5)
  
  ENIGMA_trace.v1 <- cell_deconvolve_trace(O = as.matrix(G),
                                           theta=Fra_Simulate1$theta,
                                           R=t(H11),
                                           epsilon=0.00001,solver ="admm",
                                           alpha=alpha,gamma=0.5,beta=1,Normalize=F,
                                           verbose=FALSE,max.iter = 1000,pos = F)
  res_enigma.v1 <- DEG_test1(ENIGMA_trace.v1$X_k,y,Fra_Simulate1$theta,method = "enigma",10000,DEG_list_all[[1]],qval=FALSE)
  at1 <- c(at1, res_enigma.v1)
  
  ENIGMA_trace.v1.5 <- cell_deconvolve_trace(O = as.matrix(G),
                                             theta=Fra_Simulate1.5$theta,
                                             R=t(H11.5),
                                             epsilon=0.00001,solver ="admm",
                                             alpha=alpha,gamma=0.5,beta=1,Normalize=F,
                                             verbose=FALSE,max.iter = 1000,pos = F)
  res_enigma.v1.5 <- DEG_test1(ENIGMA_trace.v1.5$X_k,y,Fra_Simulate1.5$theta,method = "enigma",10000,DEG_list_all[[1]],qval=FALSE)
  at1.5 <- c(at1.5, res_enigma.v1.5)
  
  ENIGMA_trace.v2 <- cell_deconvolve_trace(O = as.matrix(G),
                                           theta=Fra_Simulate2$theta,
                                           R=t(H12),
                                           epsilon=0.00001,solver ="admm",
                                           alpha=alpha,gamma=0.5,beta=1,
                                           verbose=FALSE,max.iter = 1000,pos = F)
  res_enigma.v2 <- DEG_test1(ENIGMA_trace.v2$X_k,y,Fra_Simulate2$theta,method = "enigma",10000,DEG_list_all[[1]],qval=FALSE)
  at2 <- c(at2, res_enigma.v2)
  
  ENIGMA_trace.v2.5 <- cell_deconvolve_trace(O = as.matrix(G),
                                             theta=Fra_Simulate2.5$theta,
                                             R=t(H12.5),
                                             epsilon=0.00001,solver ="admm",
                                             alpha=alpha,gamma=0.5,beta=1,Normalize=F,
                                             verbose=FALSE,max.iter = 1000,pos = F)
  res_enigma.v2.5 <- DEG_test1(ENIGMA_trace.v2.5$X_k,y,Fra_Simulate2.5$theta,method = "enigma",10000,DEG_list_all[[1]],qval=FALSE)
  at2.5 <- c(at2.5, res_enigma.v2.5)
  
  ENIGMA_trace.v3 <- cell_deconvolve_trace(O = as.matrix(G),
                                           theta=Fra_Simulate3$theta,
                                           R=t(H13),
                                           epsilon=0.00001,solver ="admm",
                                           alpha=alpha,gamma=0.5,beta=1,Normalize=F,
                                           verbose=FALSE,max.iter = 1000,pos = F)
  res_enigma.v3 <- DEG_test1(ENIGMA_trace.v3$X_k,y,Fra_Simulate3$theta,method = "enigma",10000,DEG_list_all[[1]],qval=FALSE)
  at3 <- c(at3, res_enigma.v3)
}
admm_res=data.frame(alpha=alpha.v,ref_beta1=at,ref_beta1_noise0.1=at0.1,ref_beta1_noise0.5=at0.5,ref_beta1_noise1=at1,
                    ref_beta1_noise1.5=at1.5,ref_beta1_noise2=at2,ref_beta1_noise2.5=at2.5,ref_beta1_noise3=at3)
saveRDS(admm_res,"/mnt/data2/zhouxl/ENIGMA/github/admm_noise/Data/admm_noise_res.rds")


admm_res=readRDS("/Path/to/Data/admm_noise_res.rds")
colnames(admm_res)=c("alpha","without noise","noise_level:0.1","noise_level:0.5","noise_level:1",
                     "noise_level:1.5","noise_level:2","noise_level:2.5","noise_level:3")
admm_data=melt(admm_res,id="alpha")
colnames(admm_data)[2]="Data"
admm_p=ggplot(admm_data,aes(x=alpha,y=value,colour=Data))+
  geom_point()+
  geom_line()+theme_classic()+
  theme(axis.title = element_text(size=17),
        axis.text = element_text(size=13),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12),
        panel.grid=element_blank())+
  ylab("AUPRC")
admm_p+ylim(c(0,0.85))+ scale_x_continuous(limits = c(0.1,1),breaks = seq(0,1,0.1))

