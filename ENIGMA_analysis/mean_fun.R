mean_fun=function(cell_decon_result,Fraction,groundtrue){
  #cell_decon_result=ENIGMA_L2norm_pbmc,Fraction=Fra_Simulate$theta
  ENIGMA_L2norm_PBMC=NULL
  for (i in 1:dim(cell_decon_result$X_k)[3]) {
    ENIGMA_L2norm_PBMC[[i]]=cell_decon_result$X_k[,,i]
  }
  names(ENIGMA_L2norm_PBMC)=dimnames(cell_decon_result$X_k)[3] %>% unlist()
  
  #per sample
  enigmal2norm_cor_sample=cor_sample_fun(groundtrue=groundtrue,CTE=ENIGMA_L2norm_PBMC,Fraction=Fraction)
  enigmal2norm_cor_sample = as.matrix(enigmal2norm_cor_sample)
  mean_sample=apply(enigmal2norm_cor_sample,1,function(x){mean(x,na.rm=T)})
  #mean_sample=mean(enigmal2norm_cor_sample,na.rm=TRUE)
  
  ##per gene
  enigmal2norm_cor_gene=cor_gene_fun(groundtrue=groundtrue,CTE=ENIGMA_L2norm_PBMC,Fraction=Fraction)
  enigmal2norm_cor_gene=as.matrix(enigmal2norm_cor_gene)
  gene_mean=apply(enigmal2norm_cor_gene,2,function(x){mean(x,na.rm=T)})
  #gene_mean=mean(enigmal2norm_cor_gene,na.rm=T)
  
  mean_r=data.frame(s=mean_sample,g=gene_mean)
  return(mean_r)
}

cor_sample_fun=function(groundtrue,CTE,Fraction){
  library(Hmisc)
  corre=NULL
  for (i in names(groundtrue)) {#cell type
    groundtrue[[i]]=as.matrix(groundtrue[[i]])
    corr=NULL
    for (j in colnames(groundtrue[[1]])) {#sample
      if(Fraction[j,i]==0){
        corr=cbind(corr,NA)
      }else{
        a=groundtrue[[i]][,j]
        b=CTE[[i]][,j]
        c=cor(a,b,method = "spearman")
        corr=cbind(corr,c)
      }
    }
    corr=as.data.frame(corr)
    colnames(corr)=colnames(groundtrue[[i]])
    corre[[i]]=corr
  }
  cor_sample=do.call(rbind,corre)
  return(cor_sample)
}

cor_gene_fun=function(groundtrue,CTE,Fraction){
  library(Hmisc)
  corre=NULL
  groundtrue_new=NULL
  CTE_new=NULL
  for (i in names(groundtrue)) {
    corr=NULL
    groundtrue_new[[i]]=groundtrue[[i]][,rownames(Fraction)[which(Fraction[,i]!=0)]] %>% as.matrix()
    CTE_new[[i]]=CTE[[i]][,rownames(Fraction)[which(Fraction[,i]!=0)]] %>% as.matrix()
    
    for (j in rownames(groundtrue_new[[i]])) {
      a=groundtrue_new[[i]][j,]
      b=CTE_new[[i]][j,]
      c=cor(a,b,method = "spearman")
      corr=rbind(corr,c)
    }
    rownames(corr)=rownames(groundtrue_new[[i]])
    corre[[i]]=corr
  }
  
  cor_gene=do.call(cbind,corre)
  colnames(cor_gene)=names(corre)
  return(cor_gene)
}