#######################
#load the required packages
setwd("/path/to/Data/ENIGMA_documents")
source("/path/to/Data/ENIGMA.R")
source("/path/to/Data/mean_fun.R")

library(Seurat)
library(data.table)
library(stringr)

library(ggsci)
library(ggplot2)
library(gridExtra)

################################################
### 1.prepare data
################################################
#Please note that the original data involved in this part takes up a lot of memory and is not uploaded to GitHub, 
#but all the processing results required for subsequent analysis are uploaded. 
#If you need the original data, you can contact us by email.

############## 1.1 simulate bulk and ground true data
#The PBMC single cell RNA-seq data could downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158055
# Because the original data is too large for direct analysis in R, we use separate processing and finally merge.

##step1: Use Python to convert files GSE158055_covid19_counts.mtx.gz,GSE158055_covid19_features.tsv.gz,GSE158055_covid19_barcodes.tsv.gz
##       into an annData object, split it into multiple subsets and store them as h5ad files.
## example coding(python):
# import scanpy as sc
# import pandas as pd
# import anndata as ad
# scdata=sc.read('./GSE158055_covid19_counts.mtx.gz')
# obs=pd.read_csv('./GSE158055_cell_annotation.csv')
# obs=obs.set_index("cellName")
# var=pd.read_csv("./features.tsv")
# var=var.set_index("gene")
# X=scdata.X.T
# adata = ad.AnnData(X, obs=obs, var=var, dtype='int32')
# adata.filename = './Data/PBMC.h5ad'##save

##step2: Using R convert the h5ad file to seurat object.
##example coding(R):
# library(Seurat)
# library(reticulate)
# sc<-import ("scanpy")
# adata <- sc$read_h5ad("./Data/PBMC_sub.h5ad")
# metadata <- adata$obs
# exp <- adata$X
# exp<-Matrix::t(exp)
# colnames(exp) <- rownames(metadata)
# rownames(exp) <- rownames(adata$var)
# Seurat.obj <- CreateSeuratObject(exp)
# Seurat.obj@meta.data <- cbind(Seurat.obj@meta.data, metadata)
# save(Seurat.obj,file="./pbmc_seurat10_1.rds")

##step3: Merge all subsets
#first subset
pbmc=readRDS("./data/pbmc_seurat10_1.rds")
pbmc=NormalizeData(pbmc,normalization.method = "RC")#cpm
metadata=pbmc@meta.data %>% rownames_to_column("cell")
metadata$majorType=str_replace_all(metadata$majorType,c("CD8"="T_CD8","CD4"="T_CD4"))

inde=c("B","Mono","T_CD8","T_CD4","NK")#overlap cell types
metadata=metadata[metadata$majorType %in% inde,]

groundTrue_r=NULL
bulk_r=NULL
samplesid=names(table(metadata$sampleID)[which(table(metadata$sampleID)!=0)])
for(x in samplesid){
  meta=metadata[metadata$sampleID== x,]
  sample_cellid=sample(meta$cell,size=500,replace = FALSE)
  
  ##simulate bulk 
  bulk=apply(GetAssayData(pbmc)[,sample_cellid],1,sum) %>% as.data.frame()
  names(bulk)=x
  bulk_r[[x]]=bulk
  
  ##simulate ground true
  metaid=merge(data.frame(cell=sample_cellid),meta,by="cell")
  groundTrue=do.call(cbind,lapply(unique(metaid$majorType),function(y){
    apply(GetAssayData(pbmc)[,metaid[metaid$majorType==y,]$cell],1,mean)
  }))
  colnames(groundTrue)=unique(metaid$majorType)
  groundTrue_r[[x]]=groundTrue
  writeLines( paste("produced sample id",x,sep=":") )
}

# second subset
pbmc=readRDS("./data/pbmc_seurat10_2.rds")
pbmc=NormalizeData(pbmc,normalization.method = "RC")#cpm
metadata=pbmc@meta.data %>% rownames_to_column("cell")
metadata$majorType=str_replace_all(metadata$majorType,c("CD8"="T_CD8","CD4"="T_CD4"))

inde=c("B","Mono","T_CD8","T_CD4","NK")#overlap celltypes
metadata=metadata[metadata$majorType %in% inde,]

samplesid=names(table(metadata$sampleID)[which(table(metadata$sampleID)!=0)])
for(x in samplesid){
  meta=metadata[metadata$sampleID== x,]
  sample_cellid=sample(meta$cell,size=500,replace = FALSE)
  
  ##simulate bulk
  bulk=apply(GetAssayData(pbmc)[,sample_cellid],1,sum) %>% as.data.frame()
  names(bulk)=x
  bulk_r[[x]]=bulk
  
  ##simulate groundtrue
  metaid=merge(data.frame(cell=sample_cellid),meta,by="cell")
  groundTrue=do.call(cbind,lapply(unique(metaid$majorType),function(y){
    apply(GetAssayData(pbmc)[,metaid[metaid$majorType==y,]$cell],1,mean)
  }))
  colnames(groundTrue)=unique(metaid$majorType)
  groundTrue_r[[x]]=groundTrue
  writeLines( paste("produced sample id",x,sep=":") )
}

## third subset
pbmc=readRDS("./data/pbmc_seurat10_3.rds")
## repeat the code from line89 to line 115...

## fouth subset
pbmc=readRDS("./data/pbmc_seurat10_4.rds")
## repeat the code from line89 to line 115...

bulk=do.call(cbind,bulk_r)#genes x 70sampels
#remove the fresh pbmc sample which is sequencing based 3'-end 10X 
bulk=bulk[,which(colnames(bulk)!="S-S001-2")]
groundTrue_r=groundTrue_r[which(names(groundTrue_r)!="S-S001-2")]

freshpbmc69samples=list(bulk=bulk,groundTrue=groundTrue_r)
saveRDS(freshpbmc69samples,"./freshpbmc69samples_sample500s_fivecellty.rds")

## There are 172 PBMC samples in this data, including 83 fresh PBMC samples. 
## The 10x 3'-end sequencing samples are removed, leaving 82 10x 5'-end sequencing samples, 
## 69 samples for simulating bulk and groundtrue data, and 13 samples for simulating reference data on the same platform(1.2.6 ref six: same 10x 5'-end).

###########################1.2 simulate reference data
#######1.2.1 ref one:seqwell platform
##Data could be downloaded from https://hosted-matrices-prod.s3-us-west-2.amazonaws.com/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection-25/blish_covid.seu.rds
seqwell=readRDS("./seq_well.obj.rds")
DefaultAssay(seqwell)="RNA"
seqwell=NormalizeData(seqwell,normalization.method = "RC")
seqwell=ScaleData(seqwell)

cellid=seqwell@meta.data[seqwell@meta.data$cell.type.coarse %in% c("B","CD4 T","CD8 T","CD14 Monocyte","CD16 Monocyte","NK"),] %>% rownames()
seqwell_sub=subset(seqwell,cells=cellid)
seqwell_sub@meta.data$cell.type.coarse=str_replace_all(seqwell_sub@meta.data$cell.type.coarse,c("CD14 "="","CD16 "="","pDC"="DC",
                                                                                                "CD4 T"="T_CD4","CD8 T"="T_CD8","Monocyte"="Mono"))

Idents(seqwell_sub)=seqwell_sub$cell.type.coarse
pbmc.markers2 <- FindAllMarkers(seqwell_sub, only.pos = TRUE, logfc.threshold = 0.15)
pbmc.markers2 %>% 
  group_by(cluster) %>%  
  top_n(n = 500, wt = avg_log2FC) -> top500_seqwell

seqwell_sub_eset=ExpressionSet(as.matrix(GetAssayData(seqwell_sub)))
pData(seqwell_sub_eset)=metadata
##remove batch effect
ref_seqwell_rmbe=remove_batch_effect(ExpressionSet(bulk_var),seqwell_sub_eset,"cell.type.coarse", n_pseudo_bulk=5000,ncores = 4)
saveRDS(ref_seqwell_rmbe,"ref_seqwell_rmbe_500s.rds")

#######
##Data could be downloaded from https://singlecell.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data?cluster=Harmony%20TSNE&spatialGroups=--&annotation=Experiment--group--study&subsample=all#study-download
pbmc.data=Read10X("./PBMC/")##ln(cpm)
pbmc.data3=exp(pbmc.data)-1

sixplat <- CreateSeuratObject(counts = pbmc.data3, project = "pbmc")
sixplat=NormalizeData(sixplat,normalization.method = "RC")
sixplat=ScaleData(sixplat)

#rename cell type
meta=fread("/path/to/meta.txt") %>% as.data.frame()
rownames(meta)=meta$NAME_TYPE
meta=meta[colnames(pbmc.data),]
sixplat@meta.data=meta
sixplat@meta.data$CellType_group=str_replace_all(sixplat@meta.data$CellType_group,c("B cell"="B","CD14+ monocyte"="Mono","CD16+ monocyte"="Mono","CD4+ T cell"="T_CD4",
                                                                                    "Cytotoxic T cell"="T_CD8","Dendritic cell"="DC",
                                                                                    "Natural killer cell"="NK","Plasmacytoid dendritic cell"="DC"))
cell=c("B","T_CD8","DC","NK","CD4+ T cell","Megakaryocyte","Unassigned")
sixplat@meta.data$CellType_group=ifelse(sixplat@meta.data$CellType_group %in% cell,
                                        sixplat@meta.data$CellType_group,
                                        substring(sixplat@meta.data$CellType_group,7))
cell=c("B","T_CD8","DC","NK","monocyte","Megakaryocyte","Unassigned")
sixplat@meta.data$CellType_group=ifelse(sixplat@meta.data$CellType_group %in% cell,
                                        sixplat@meta.data$CellType_group,
                                        substring(sixplat@meta.data$CellType_group,6))
sixplat@meta.data$CellType_group=gsub('T cell',"T_CD4",sixplat@meta.data$CellType_group)
sixplat@meta.data$CellType_group=gsub("monocyte","Mono",sixplat@meta.data$CellType_group)

inde=c("B","Mono","T_CD8","T_CD4","NK")
cellid=sixplat@meta.data[sixplat@meta.data$CellType_group %in% inde,] %>% rownames()
sixplat=subset(sixplat,cells=cellid)

########1.2.2 ref two:Drop-seq
cellid=sixplat@meta.data[sixplat@meta.data$orig.ident=="pbmc1"&sixplat@meta.data$Method_group=="Drop-seq",]
cellid=cellid[cellid$CellType_group %in% inde,]
sixplat_sub=subset(sixplat,cells=cellid$NAME_TYPE)

Idents(sixplat_sub)=sixplat_sub@meta.data$CellType_group
pbmc.markers2 <- FindAllMarkers(sixplat_sub, only.pos = TRUE, logfc.threshold = 0.05)
pbmc.markers2 %>% 
  group_by(cluster) %>%  
  top_n(n = 500, wt = avg_log2FC) -> top500_drop

Drop=ExpressionSet(as.matrix(GetAssayData(sixplat_sub)))
pData(Drop)=cellid
##remove batch effect
ref_drop_rmbe=remove_batch_effect(ExpressionSet(bulk_var),Drop,"CellType_group", n_pseudo_bulk=5000)
saveRDS(ref_drop_rmbe,"ref_drop_rmbe_500s.rds")

########1.2.3 ref three:10x 3'-end
cellid=sixplat@meta.data[sixplat@meta.data$orig.ident=="pbmc1"&sixplat@meta.data$Method_group=="10x Chromium (v2) B",]
cellid=cellid[cellid$CellType_group %in% inde,]
sixplat_sub=subset(sixplat,cells=cellid$NAME_TYPE)

Idents(sixplat_sub)=sixplat_sub@meta.data$CellType_group
pbmc.markers2 <- FindAllMarkers(sixplat_sub, only.pos = TRUE, logfc.threshold = 0.05)
pbmc.markers2 %>% 
  group_by(cluster) %>%  
  top_n(n = 500, wt = avg_log2FC) -> top500_x10

x10=ExpressionSet(as.matrix(GetAssayData(sixplat_sub)))
pData(x10)=cellid
ref_10x_rmbe=remove_batch_effect(ExpressionSet(bulk_var),x10,"CellType_group", n_pseudo_bulk=5000)
saveRDS(ref_seqwell_rmbe,"ref_seqwell_rmbe_500s.rds")

########1.2.4 ref four:inDrops
cellid=sixplat@meta.data[sixplat@meta.data$orig.ident=="pbmc2"&sixplat@meta.data$Method_group=="inDrops",]
cellid=cellid[cellid$CellType_group %in% inde,]
sixplat_sub=subset(sixplat,cells=cellid$NAME_TYPE)

Idents(sixplat_sub)=sixplat_sub@meta.data$CellType_group
pbmc.markers2 <- FindAllMarkers(sixplat_sub, only.pos = TRUE, logfc.threshold = 0.05)
pbmc.markers2 %>% 
  group_by(cluster) %>%  
  top_n(n = 500, wt = avg_log2FC) -> top500_indrop

inDrops=ExpressionSet(as.matrix(GetAssayData(sixplat_sub)))
pData(inDrops)=cellid
ref_inDrops_rmbe=remove_batch_effect(ExpressionSet(bulk_var),inDrops,"CellType_group", n_pseudo_bulk=5000)
saveRDS(ref_inDrops_rmbe,"ref_indrop_rmbe_500s.rds")

########1.2.5 ref five:Smart-seq2
cellid=sixplat@meta.data[sixplat@meta.data$orig.ident=="pbmc2"&sixplat@meta.data$Method_group=="Smart-seq2",]
cellid=cellid[cellid$CellType_group %in% inde,]
sixplat_sub=subset(sixplat,cells=cellid$NAME_TYPE)

Idents(sixplat_sub)=sixplat_sub@meta.data$CellType_group
pbmc.markers2 <- FindAllMarkers(sixplat_sub, only.pos = TRUE, logfc.threshold = 0.05)
pbmc.markers2 %>% 
  group_by(cluster) %>%  
  top_n(n = 500, wt = avg_log2FC) -> top500_smart

Smart=ExpressionSet(as.matrix(GetAssayData(sixplat_sub)))
pData(Smart)=cellid
ref_Smart_rmbe=remove_batch_effect(ExpressionSet(bulk_var),Smart,"CellType_group", n_pseudo_bulk=5000)
saveRDS(ref_Smart_rmbe,"ref_smart_rmbe_500s.rds")

#######1.2.6 ref six: same 10x 5'-end
##
pbmc13=readRDS("./pbmc_seurat_13_ref.rds")
pbmc13@meta.data$cell=rownames(pbmc13@meta.data)
pbmc13=NormalizeData(pbmc13,normalization.method = "RC")
pbmc13=ScaleData(pbmc13)

pbmc13@meta.data$majorType=str_replace_all(pbmc13@meta.data$majorType,c("CD4"="T_CD4","CD8"="T_CD8"))
metadata=pbmc13@meta.data
metadata=metadata[metadata$majorType %in% inde,]
pbmc13_sub=subset(pbmc13,cells=metadata$cell)

Idents(pbmc13_sub)=pbmc13_sub@meta.data$majorType
pbmc.markers2 <- FindAllMarkers(pbmc13_sub, only.pos = TRUE, logfc.threshold = 0.05)
pbmc.markers2 %>% 
  group_by(cluster) %>%  
  top_n(n = 500, wt = avg_log2FC) -> top500_same10x

pbmc13_sub_eset=ExpressionSet(as.matrix(GetAssayData(pbmc13_sub)))
pData(pbmc13_sub_eset)=metadata
ref_same10x_rmbe=remove_batch_effect(ExpressionSet(bulk_var),pbmc13_sub_eset,"majorType", n_pseudo_bulk=5000,ncores = 4)
saveRDS(ref_same10x_rmbe,"./ref_same10x_rmbe_500s.rds")

top500=list(top500_drop,top500_indrop,top500_seqwell,top500_smart,top500_x10,top500_same10x)
names(top500)=c('top500_drop','top500_indrop','top500_seqwell','top500_smart','top500_x10','top500_same10x')
saveRDS(top500,"./top500_markergenes.rds")

#################################### Pseudo-bulk #########################
###### remove genes which var < 10^(-8)
bulk_var=bulk %>% as.matrix(.)%>% apply(.,1,var)
bulk_var=bulk_var[bulk_var>10^(-8)]
bulk_var=bulk[rownames(bulk) %in% names(bulk_var),]

#overlap genes
seqwell=readRDS("/path/to/Data/seqwell_ref.rds")
sixplat=readRDS("/path/to/Data/pbmc_ref.rds")
geneset=Reduce(intersect,list(r1=rownames(seqwell@assays$RNA),
                              r2=rownames(sixplat@assays$RNA),
                              r3=rownames(bulk_var)))
bulk_var=bulk_var[geneset,]
saveRDS(bulk_var,"./bulk_var_500s_fivecellty.rds")

################################### Ground True ########################
groundTrue_r=freshpbmc69samples$groundTrue
#All cell types covered by sampling
celltype=NULL
for (i in 1:length(groundTrue_r)) {
  celltype=unique(c(celltype,colnames(groundTrue_r[[i]])))
}

#Complete five cell types
for (i in 1:length(groundTrue_r)){
  df=setdiff(celltype,colnames(groundTrue_r[[i]]))
  dfdata=matrix(0,ncol=length(df),nrow=nrow(groundTrue_r[[i]])) %>% as.data.frame()
  colnames(dfdata)=df
  groundTrue_r[[i]]=cbind(groundTrue_r[[i]],dfdata)
}

for (i in names(groundTrue_r)) {
  groundTrue_r[[i]]=groundTrue_r[[i]][geneset,]
  groundTrue_r[[i]]=groundTrue_r[[i]][,inde]
}

single_celltype=NULL
for (j in colnames(groundTrue_r[[1]])) {
  single_celltype[[j]]=do.call(cbind,lapply(groundTrue_r,function(x){
    x[j]
  }))
  colnames(single_celltype[[j]])=names(groundTrue_r)
}
saveRDS(single_celltype,"./single_celltype_500s_fivecellty.rds")

##Please note
#The above code is the reference of the data preparation part (the user needs to download the original data by himself). 
#The user can directly reproduce the results in our manuscript, according to the following code.

#################################################################
##########      2. L2-norm 
#################################################################
##Datasets are provided in /ENIGMA_analysis/Data/mutilPlatforms/
single_celltype=readRDS("/path/to/Data/single_celltype_500s_fivecellty.rds")
bulk_var=readRDS("/path/to/Data/bulk_var_500s_fivecellty.rds")

ref_seqwell=readRDS("/path/to/Data/ref_seqwell_rmbe_500s.rds")
ref_drops=readRDS("/path/to/Data/ref_drop_rmbe_500s.rds")
ref_10x=readRDS("/path/to/Data/ref_10x_rmbe_500s.rds")
ref_indrop=readRDS("/path/to/Data/ref_indrop_rmbe_500s.rds")
ref_smart=readRDS("/path/to/Data/ref_smart_rmbe_500s.rds")
ref_same10x=readRDS("/path/to/Data/ref_same10x_rmbe_500s.rds")

inde=c("B","Mono","T_CD8","T_CD4","NK")
ref_seqwell=ref_seqwell$cell.type.coarse
ref_drops=ref_drops$CellType_group
ref_10x=ref_10x$CellType_group
ref_smart=ref_smart$CellType_group
ref_indrop=ref_indrop$CellType_group
ref_same10x=ref_same10x$majorType

ref_seqwell=ref_seqwell[,inde]
ref_drops=ref_drops[,inde]
ref_10x=ref_10x[,inde]
ref_smart=ref_smart[,inde]
ref_indrop=ref_indrop[,inde]
ref_same10x=ref_same10x[,inde]

Fra_Simulate_10x <- get_proportion(bulk_var, ref_10x)
Fra_Simulate_smart <- get_proportion(bulk_var, ref_smart)
Fra_Simulate_seqwell<- get_proportion(bulk_var, ref_seqwell)
Fra_Simulate_drops <- get_proportion(bulk_var, ref_drops)
Fra_Simulate_indrop <- get_proportion(bulk_var, ref_indrop)
Fra_Simulate_same10x <- get_proportion(bulk_var, ref_same10x)

alpha.v <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
enigma.mean.10x<-enigma.mean.seqwell<-enigma.mean.smart<-enigma.mean.indrop<-enigma.mean.drops<-enigma.mean.same10x<-NULL
ENIGMA_10x<-ENIGMA_seqwell<-ENIGMA_smart<-ENIGMA_indrop<-ENIGMA_drops<-ENIGMA_same10x<-NULL

bulk_var=as.matrix(bulk_var)
for (k in 1:length(alpha.v)) {
  alpha<-alpha.v[k]
  ENIGMA_L2norm_pbmc_10x <- cell_deconvolve(X=log2(bulk_var+1),
                                            theta=Fra_Simulate_10x$theta,
                                            R=log2(ref_10x+1),
                                            epsilon=0.001,
                                            alpha=alpha,
                                            beta=0.3,tao_k=0.4,verbose=TRUE,Normalize=T,pos=F)
  enigma.mean1=mean_fun(cell_decon_result=ENIGMA_L2norm_pbmc_10x,Fraction=Fra_Simulate_10x$theta,groundtrue=single_celltype)
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean1)=paste(alpha_name,rownames(enigma.mean1),sep="-")
  enigma.mean.10x=rbind(enigma.mean.10x,enigma.mean1)
  
  ENIGMA_L2norm_pbmc_seqwell <- cell_deconvolve(X=log2(bulk_var+1),
                                                theta=Fra_Simulate_seqwell$theta,
                                                R=log2(ref_seqwell+1),
                                                epsilon=0.001,
                                                alpha=alpha,
                                                beta=0.3,tao_k=0.4,verbose=TRUE,Normalize=T,pos=F)
  enigma.mean2=mean_fun(cell_decon_result=ENIGMA_L2norm_pbmc_seqwell,Fraction=Fra_Simulate_seqwell$theta,groundtrue=single_celltype)
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean2)=paste(alpha_name,rownames(enigma.mean2),sep="-")
  enigma.mean.seqwell=rbind(enigma.mean.seqwell,enigma.mean2)
  
  ENIGMA_L2norm_pbmc_smart <- cell_deconvolve(X=log2(bulk_var+1),
                                              theta=Fra_Simulate_smart$theta,
                                              R=log2(ref_smart+1),
                                              epsilon=0.001,
                                              alpha=alpha,
                                              beta=0.3,tao_k=0.4,verbose=TRUE,Normalize=T,pos=F)
  enigma.mean3=mean_fun(cell_decon_result=ENIGMA_L2norm_pbmc_smart,Fraction=Fra_Simulate_smart$theta,groundtrue=single_celltype)#
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean3)=paste(alpha_name,rownames(enigma.mean3),sep="-")
  enigma.mean.smart=rbind(enigma.mean.smart,enigma.mean3)
  
  ENIGMA_L2norm_pbmc_indrop<- cell_deconvolve(X=log2(bulk_var+1),
                                              theta=Fra_Simulate_indrop$theta,
                                              R=log2(ref_indrop+1),
                                              epsilon=0.001,
                                              alpha=alpha,
                                              beta=0.3,tao_k=0.4,verbose=TRUE,Normalize=T,pos=F)
  enigma.mean4=mean_fun(cell_decon_result=ENIGMA_L2norm_pbmc_indrop,Fraction=Fra_Simulate_indrop$theta,groundtrue=single_celltype)
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean4)=paste(alpha_name,rownames(enigma.mean4),sep="-")
  enigma.mean.indrop=rbind(enigma.mean.indrop,enigma.mean4)
  
  ENIGMA_L2norm_pbmc_drops <- cell_deconvolve(X=log2(bulk_var+1),
                                              theta=Fra_Simulate_drops$theta,
                                              R=log2(ref_drops+1),
                                              epsilon=0.001,
                                              alpha=alpha,
                                              beta=0.3,tao_k=0.4,verbose=TRUE,Normalize=T,pos=F)
  enigma.mean5=mean_fun(cell_decon_result=ENIGMA_L2norm_pbmc_drops,Fraction=Fra_Simulate_drops$theta,groundtrue=single_celltype)#
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean5)=paste(alpha_name,rownames(enigma.mean5),sep="-")
  enigma.mean.drops=rbind(enigma.mean.drops,enigma.mean5)
  
  ENIGMA_L2norm_pbmc_same10x <- cell_deconvolve(X=log2(bulk_var+1),
                                                theta=Fra_Simulate_same10x$theta,
                                                R=log2(ref_same10x+1),
                                                epsilon=0.001,
                                                alpha=alpha,
                                                beta=0.3,tao_k=0.4,verbose=TRUE,Normalize=T,pos=F)
  enigma.mean6=mean_fun(cell_decon_result=ENIGMA_L2norm_pbmc_same10x,Fraction=Fra_Simulate_same10x$theta,groundtrue=single_celltype)#
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean6)=paste(alpha_name,rownames(enigma.mean6),sep="-")
  enigma.mean.same10x=rbind(enigma.mean.same10x,enigma.mean6)
}
enigma_mean_rmbe=list(enigma.mean.10x=enigma.mean.10x, enigma.mean.indrop=enigma.mean.indrop,
                      enigma.mean.smart=enigma.mean.smart, enigma.mean.seqwell=enigma.mean.seqwell,
                      enigma.mean.drops=enigma.mean.drops,enigma.mean.same10x=enigma.mean.same10x)


################### plot
cor_alpha=do.call(rbind,enigma_mean_rmbe)
cor_alpha=cor_alpha %>% rownames_to_column("name")
cor_alpha$Platforms=sapply(strsplit(cor_alpha$name,"[.]"),function(x){x[[3]]})
cor_alpha$celltype=sapply(strsplit(cor_alpha$name,"-"),function(x){x[[2]]})
cor_alpha$alpha=rep((rep(alpha.v,each=5)),6)
 
pp2=ggplot(cor_alpha,aes(x=alpha,y=s,fill=Platforms,colour=Platforms))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        panel.grid=element_blank())+
  ylab("Correlation per sample")+
  facet_wrap(~celltype, nrow=2)
pp2+scale_color_npg(alpha=0.9)+scale_x_continuous(breaks = seq(0,1,0.1))+theme(legend.position = c(0.85,0.25))


##########################################################
####      3.trace norm
##########################################################
single_celltype=readRDS("single_celltype_500s_fivecellty.rds")
bulk_var=readRDS("bulk_var_500s_fivecellty.rds")
top500=readRDS("top500_markergenes.rds")

ref_seqwell=readRDS("/path/to/Data/ref_seqwell_rmbe_500s.rds")
ref_drops=readRDS("/path/to/Data/ref_drop_rmbe_500s.rds")
ref_10x=readRDS("/path/to/Data/ref_10x_rmbe_500s.rds")
ref_indrop=readRDS("/path/to/Data/ref_indrop_rmbe_500s.rds")
ref_smart=readRDS("/path/to/Data/ref_smart_rmbe_500s.rds")
ref_same10x=readRDS("/path/to/Data/ref_same10x_rmbe_500s.rds")

inde=c("B","Mono","T_CD8","T_CD4","NK")
ref_seqwell=ref_seqwell$cell.type.coarse %>% as.data.frame()
ref_drops=ref_drops$CellType_group %>% as.data.frame()
ref_10x=ref_10x$CellType_group %>% as.data.frame()
ref_smart=ref_smart$CellType_group %>% as.data.frame()
ref_indrop=ref_indrop$CellType_group %>% as.data.frame()
ref_same10x=ref_same10x$majorType %>% as.data.frame()

ref_seqwell=ref_seqwell[intersect(top500$top500_seqwell$gene,rownames(ref_seqwell)),inde] %>% as.matrix()
ref_drops=ref_drops[intersect(top500$top500_drop$gene,rownames(ref_drops)),inde] %>% as.matrix()
ref_10x=ref_10x[intersect(top500$top500_x10$gene,rownames(ref_10x)),inde] %>% as.matrix()
ref_smart=ref_smart[intersect(top500$top500_smart$gene,rownames(ref_smart)),inde] %>% as.matrix()
ref_indrop=ref_indrop[intersect(top500$top500_indrop$gene,rownames(ref_indrop)),inde] %>% as.matrix()
ref_same10x=ref_same10x[intersect(top500$top500_same10x$gene,rownames(ref_same10x)),inde] %>% as.matrix()

single_celltype_10x=single_celltype_seqwell=single_celltype_drops=single_celltype_indrop=single_celltype_smart=single_celltype_same10x=NULL
for (i in names(single_celltype)) {
  single_celltype_10x[[i]]=single_celltype[[i]][rownames(ref_10x),]
  single_celltype_seqwell[[i]]=single_celltype[[i]][rownames(ref_seqwell),]
  single_celltype_drops[[i]]=single_celltype[[i]][rownames(ref_drops),]
  single_celltype_indrop[[i]]=single_celltype[[i]][rownames(ref_indrop),]
  single_celltype_smart[[i]]=single_celltype[[i]][rownames(ref_smart),]
  single_celltype_same10x[[i]]=single_celltype[[i]][rownames(ref_same10x),]
}

Fra_Simulate_10x <- get_proportion(bulk_var[rownames(ref_10x),], ref_10x)
Fra_Simulate_smart <- get_proportion(bulk_var[rownames(ref_smart),], ref_smart)
Fra_Simulate_seqwell<- get_proportion(bulk_var[rownames(ref_seqwell),], ref_seqwell)
Fra_Simulate_drops <- get_proportion(bulk_var[rownames(ref_drops),], ref_drops)
Fra_Simulate_indrop <- get_proportion(bulk_var[rownames(ref_indrop),], ref_indrop)
Fra_Simulate_same10x <- get_proportion(bulk_var[rownames(ref_same10x),], ref_same10x)

alpha.v <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
enigma.mean.10x<-enigma.mean.seqwell<-enigma.mean.smart<-enigma.mean.indrop<-enigma.mean.drops<-enigma.mean.same10x<-NULL
ENIGMA_10x<-ENIGMA_seqwell<-ENIGMA_smart<-ENIGMA_indrop<-ENIGMA_drops<-ENIGMA_same10x<-NULL

bulk_var=as.matrix(bulk_var)
for (k in 1:length(alpha.v)) {
  alpha<-alpha.v[k]
  ENIGMA_L2norm_pbmc_10x <- cell_deconvolve_trace(O=bulk_var[rownames(ref_10x),],
                                                  theta=Fra_Simulate_10x$theta,
                                                  R=ref_10x,Norm.method = "frac",tao_k=1,
                                                  alpha=alpha,verbose=T,Normalize=T)
  ENIGMA_10x[[k]]=ENIGMA_L2norm_pbmc_10x

  enigma.mean1=mean_fun(cell_decon_result=ENIGMA_L2norm_pbmc_10x,Fraction=Fra_Simulate_10x$theta,groundtrue=single_celltype_10x)
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean1)=paste(alpha_name,rownames(enigma.mean1),sep="-")
  enigma.mean.10x=rbind(enigma.mean.10x,enigma.mean1)

  ENIGMA_L2norm_pbmc_seqwell <-cell_deconvolve_trace(O=bulk_var[rownames(ref_seqwell),],
                                                     theta=Fra_Simulate_seqwell$theta,
                                                     R=ref_seqwell,Norm.method = "frac",
                                                     alpha=alpha,tao_k=1,
                                                     verbose=T,Normalize=T)

  ENIGMA_seqwell[[k]]=ENIGMA_L2norm_pbmc_seqwell
  enigma.mean2=mean_fun(cell_decon_result=ENIGMA_L2norm_pbmc_seqwell,Fraction=Fra_Simulate_seqwell$theta,groundtrue=single_celltype_seqwell)
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean2)=paste(alpha_name,rownames(enigma.mean2),sep="-")
  enigma.mean.seqwell=rbind(enigma.mean.seqwell,enigma.mean2)

  ENIGMA_L2norm_pbmc_smart <- cell_deconvolve_trace(O=bulk_var[rownames(ref_smart),],
                                                    theta=Fra_Simulate_smart$theta,
                                                    R=ref_smart,Norm.method = "frac",
                                                    alpha=alpha,tao_k=1,
                                                    verbose=T,Normalize=T)

  ENIGMA_smart[[k]]=ENIGMA_L2norm_pbmc_smart
  enigma.mean3=mean_fun(cell_decon_result=ENIGMA_L2norm_pbmc_smart,Fraction=Fra_Simulate_smart$theta,groundtrue=single_celltype_smart)#
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean3)=paste(alpha_name,rownames(enigma.mean3),sep="-")
  enigma.mean.smart=rbind(enigma.mean.smart,enigma.mean3)

  ENIGMA_L2norm_pbmc_indrop<- cell_deconvolve_trace(O=bulk_var[rownames(ref_indrop),],
                                                    theta=Fra_Simulate_indrop$theta,
                                                    R=ref_indrop,Norm.method = "frac",
                                                    alpha=alpha,tao_k=1,
                                                    verbose=T,Normalize=T)

  ENIGMA_indrop[[k]]=ENIGMA_L2norm_pbmc_indrop
  enigma.mean4=mean_fun(cell_decon_result=ENIGMA_L2norm_pbmc_indrop,Fraction=Fra_Simulate_indrop$theta,groundtrue=single_celltype_indrop)
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean4)=paste(alpha_name,rownames(enigma.mean4),sep="-")
  enigma.mean.indrop=rbind(enigma.mean.indrop,enigma.mean4)

  ENIGMA_L2norm_pbmc_drops <-cell_deconvolve_trace(O=bulk_var[rownames(ref_drops),],
                                                   theta=Fra_Simulate_drops$theta,
                                                   R=ref_drops,Norm.method = "frac",
                                                   alpha=alpha,tao_k=1,
                                                   verbose=T,Normalize=T)

  ENIGMA_drops[[k]]=ENIGMA_L2norm_pbmc_drops
  enigma.mean5=mean_fun(cell_decon_result=ENIGMA_L2norm_pbmc_drops,Fraction=Fra_Simulate_drops$theta,groundtrue=single_celltype_drops)#
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean5)=paste(alpha_name,rownames(enigma.mean5),sep="-")
  enigma.mean.drops=rbind(enigma.mean.drops,enigma.mean5)

  ENIGMA_L2norm_pbmc_same10x <- cell_deconvolve_trace(O=bulk_var[rownames(ref_same10x),], 
                                                      theta=Fra_Simulate_same10x$theta, 
                                                      R=ref_same10x, Norm.method = "frac",
                                                      alpha=alpha,tao_k=0.1,
                                                      verbose=T,Normalize=T) 
  
  ENIGMA_same10x[[k]]=ENIGMA_L2norm_pbmc_same10x
  enigma.mean6=mean_fun(cell_decon_result=ENIGMA_L2norm_pbmc_same10x,Fraction=Fra_Simulate_same10x$theta,groundtrue=single_celltype_same10x)#
  alpha_name=paste("alpha",alpha,sep = "_")
  rownames(enigma.mean6)=paste(alpha_name,rownames(enigma.mean6),sep="-")
  enigma.mean.same10x=rbind(enigma.mean.same10x,enigma.mean6)
}
enigma_trace_mean_rmbe=list(enigma.mean.10x=enigma.mean.10x, enigma.mean.indrop=enigma.mean.indrop,
                            enigma.mean.smart=enigma.mean.smart, enigma.mean.seqwell=enigma.mean.seqwell,
                            enigma.mean.drops=enigma.mean.drops,enigma.mean.same10x=enigma.mean.same10x)

#plot
trace_cor_alpha=do.call(rbind,enigma_trace_mean_rmbe)
trace_cor_alpha=trace_cor_alpha %>% rownames_to_column("name")
trace_cor_alpha$Platforms=sapply(strsplit(trace_cor_alpha$name,"[.]"),function(x){x[[3]]})
trace_cor_alpha$celltype=sapply(strsplit(trace_cor_alpha$name,"-"),function(x){x[[2]]})
trace_cor_alpha$alpha=rep((rep(alpha.v,each=5)),6)

pps=ggplot(trace_cor_alpha,aes(x=alpha,y=s,fill=Platforms,colour=Platforms))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        panel.grid=element_blank())+
  ylab("Correlation per sample")+
  facet_wrap(~celltype, nrow=2)
pps+scale_color_npg(alpha=0.9)+scale_x_continuous(breaks = seq(0,1,0.1))+theme(legend.position = c(0.85,0.25))

ppg=ggplot(trace_cor_alpha,aes(x=alpha,y=g,fill=Platforms,colour=Platforms))+
  geom_point()+
  geom_line()+theme_bw()+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        panel.grid=element_blank())+
  ylab("Correlation per gene")+
  facet_wrap(~celltype, nrow=2)
ppg+scale_color_npg(alpha=0.9)+scale_x_continuous(breaks = seq(0,1,0.1))+theme(legend.position = c(0.85,0.25))














