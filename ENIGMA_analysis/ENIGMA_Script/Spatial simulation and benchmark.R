rm(list = ls())

library(dplyr)
library(MIND)
library(TCA)
library(STdeconvolve)
library(Seurat)
library(ggplot2)
library(reticulate)
library(sceasy)
library(anndata)
########################### Load Simulation dataset ############################
### structed simulation
dataset <- readRDS("dataset_structed.rds")
Bulk = dataset$mixture %>% as.matrix()
### Build deconvolution profile
metadata = dataset$ref_sc_metadata
profile = NULL
for(i in c("Luminal cell","Fibroblast","Myoepithelial")){
  profile <- cbind(profile, rowMeans(dataset$ref_sc[,rownames(subset(metadata,m.CellTypes == i))]))
}
colnames(profile) <- c("LuminalCell","Fibroblast","Myoepithelial")
profile <- profile[rownames(Bulk),]
### Expression matrix truth
gt <- dataset[c("lumi_gep","fibro_gep","myo_gep")]
names(gt) <- c("LuminalCell","Fibroblast","Myoepithelial")
############################### Predict fraction ###############################
source("ENIGMA.R")
source("ENIGMA_revise.R")
source("get_proportion.R")
Frac <- get_proportion(Bulk, profile)
Frac <- Frac$theta
####################################### ENIGMA #################################
### ENIGMA L2
alpha = 0.5
egm.L2 <- cell_deconvolve(X=Bulk,
                            theta=Frac,
                            R=profile,
                            epsilon=0.001,
                            alpha=alpha,
                            pre.process="log",
                            beta=0.5,
                            tao_k=0.1,
                            max.iter=1000,
                            verbose=TRUE,
                            Normalize=TRUE,
                            Norm.method = "frac")
egm.L2.exp <- egm.L2$X_k
### ENIGMA trace
alpha = 0.1
egm <- cell_deconvolve_trace(O=Bulk,
                                theta=Frac$theta,
                                R=profile,
                                alpha=alpha,
                                pre.process="log",
                                Norm.method = "frac",
                                verbose = TRUE,
                                epsilon = 0.001,
                                epsilon_ks = 0.001,
                                tao_k = 1,
                                max_ks = 1)
egm.trace.exp <- egm.trace$X_k
##################################### bMIND ####################################
mind = bMIND(log2(Bulk+1),
             ncore = 7,
             frac = Frac,
             profile = as.matrix(log2(profile+1)))
### predict fraction in bMIND
sc_metadata <- filter(metadata,m.CellTypes %in% c("Luminal cell","Fibroblast","Myoepithelial"))
sc_counts <- dataset$ref_sc[,rownames(sc_metadata)]
dim(sc_counts)
sc_metadata <- sc_metadata[,"m.CellTypes"] %>% data.frame() %>% rename("cell_type" = ".")
rownames(sc_metadata) <- colnames(sc_counts)
sc_metadata$subject <- rownames(sc_metadata)
dim(sc_metadata)
head(sc_metadata)
mind <- bMIND(bulk = log2(Bulk+1),  
              sc_count = as.matrix(sc_counts), 
              sc_meta = sc_metadata,
              frac_method = "Bisque",
              ncore = 40)
bmind.exp<- mind$A
dimnames(bmind.exp)[2]
dimnames(bmind.exp)[2][[1]] <- c("Fibroblast","LuminalCell","Myoepithelial")
################################## TCA #########################################
tca.mdl <- tca(X = as.matrix(log2(Bulk+1)),
               W = Frac,
               C1 = NULL,
               C2 = NULL,
               parallel = TRUE,
               num_cores=20,               
               max_iter = 2)
tca.exp <- tensor(X = as.matrix(log2(Bulk+1)), tca.mdl)
names(tca.exp) <- colnames(profile)
################################ STdeconvolve ##################################
counts <- t(as.matrix(round(dataset$mixture)))
ldas <- fitLDA(counts,perc.rare.thresh = 0.05,plot= F,verbose=TRUE,Ks = seq(2, 8, by = 1),ncores=30)
optLDA <- optimalModel(models = ldas, opt = "min")
results <- getBetaTheta(optLDA,
                        perc.filt = 0.05,
                        betaScale = 1000)
### proportion
deconProp <- results$theta
### profile
deconGexp <- results$beta
head(deconGexp)
### assign cell type
metadata = dataset$ref_sc_metadata
profile <- profile[colnames(deconGexp),]
dim(profile)
mat <- matrix(nrow = nrow(deconGexp),ncol = ncol(profile))
for (i in 1:nrow(deconGexp)) {
  for (j in 1:ncol(profile)){
    mat[i,j] <- cor(deconGexp[i,],profile[,j],method = "sp")
  }
}
rownames(mat) <- rownames(deconGexp)
colnames(mat) <- colnames(profile)
type.use = apply(mat,2,function(x){which.max(x)})
type.use
STd.exp <- deconGexp[type.use,]
rownames(STd.exp) <- names(type.use)
#################################### DestVI ####################################
st_data <- CreateSeuratObject(counts = dataset$mixture,min.cells = 0,min.features = 0,assay = "Spatial")
st_data <- SCTransform(st_data, assay = "Spatial", verbose = FALSE)
st_data <- RunPCA(st_data, assay = "SCT", verbose = FALSE)
st_data <- FindNeighbors(st_data, reduction = "pca", dims = 1:30)
st_data <- FindClusters(st_data, verbose = FALSE)
st_data <- RunUMAP(st_data, reduction = "pca", dims = 1:30)
sc_data <- CreateSeuratObject(counts = dataset$ref_sc,meta.data = dataset$ref_sc_metadata,min.cells = 0,min.features = 0)
sc_data <- subset(sc_data,m.CellTypes %in% c("Luminal cell","Fibroblast","Myoepithelial"))
sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 10000)
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)
intersect <- intersect(rownames(st_data), rownames(sc_data))
sc_data <- sc_data[intersect]
st_data <- st_data[intersect]
use_condaenv(condaenv = "scvi-env",conda = "/home/tools/anaconda3/bin/conda")
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)
sc_adata <- convertFormat(sc_data, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
st_adata <- convertFormat(st_data, from="seurat", to="anndata", assay="Spatial", main_layer="counts", drop_single_values=FALSE)

### trian sc model
scvi$model$CondSCVI$setup_anndata(sc_adata, labels_key="m.CellTypes")
sclvm <- scvi$model$CondSCVI(sc_adata, weight_obs=TRUE)
sclvm$train(max_epochs=as.integer(250),early_stopping = T,early_stopping_monitor = 'train_loss')

### deconvolution with stLVM
scvi$model$DestVI$setup_anndata(st_adata)
stlvm <- scvi$model$DestVI$from_rna_model(st_adata, sclvm)
stlvm$train(max_epochs=as.integer(2500),early_stopping = T,early_stopping_monitor = 'train_loss')

### proportion
st_adata$obsm["proportions"] <- stlvm$get_proportions()
head(py_to_r(st_adata$obsm$get("proportions")))
### profile
destvi.exp <- list()
cts = c("Luminal cell","Fibroblast","Myoepithelial")
for(i in 1: length(cts)) {
  destvi.exp[[i]] <- as.matrix(py_to_r(stlvm$get_scale_for_ct(cts[i])$T))
}
names(destvi.exp) <- c("LuminalCell","Fibroblast","Myoepithelial")
############################# Calculate correlation ############################
Bulk <- Bulk[,rownames(Frac[rowSums(Frac)!=0,])]
dim(Bulk)
## correlation of per gene
cts <- names(gt)
corr_gene <- list()
for (i in 1:length(cts)) {
  t.exp <- tca.exp[[cts[i]]][rownames(Bulk),colnames(Bulk)]
  b.exp <- bmind.exp[,cts[i],][rownames(Bulk),colnames(Bulk)] 
  d.exp <- destvi.exp[[cts[i]]][rownames(Bulk),colnames(Bulk)]
  gt[[cts[i]]] <- gt[[cts[i]]][rownames(Bulk),colnames(Bulk)]
  e.trace.exp <- egm.trace.exp[,,cts[i]][rownames(Bulk),colnames(Bulk)]
  e.L2.exp <- egm.L2.exp[,,cts[i]][rownames(Bulk),colnames(Bulk)]
  t.corr <- c()
  b.corr <- c()
  d.corr <- c()
  e.trace.corr <- c()
  e.L2.corr <- c()
  for (gene in rownames(Bulk)) {
    t.corr <- c(t.corr,cor(gt[[cts[i]]][gene,],t.exp[gene,],method="sp"))
    b.corr <- c(b.corr,cor(gt[[cts[i]]][gene,],b.exp[gene,],method="sp"))
    d.corr <- c(d.corr,cor(gt[[cts[i]]][gene,],d.exp[gene,],method="sp"))
    e.trace.corr <- c(e.trace.corr,cor(gt[[cts[i]]][gene,],e.trace.exp[gene,],method="sp"))
    e.L2.corr <- c(e.L2.corr,cor(gt[[cts[i]]][gene,],e.L2.exp[gene,],method="sp"))
  }
  corr <- cbind(t.corr,b.corr,d.corr,e.trace.corr,e.L2.corr)
  colnames(corr) <- c("TCA","bMIND","DestVI","ENIGMA (trace)","ENIGMA")
  rownames(corr) <- rownames(Bulk)
  corr_gene[[i]] <- corr
}
names(corr_gene) <- cts
## correlation of per pixel
corr_pixel <- list()
for (i in 1:length(cts)) {
  t.exp <- tca.exp[[cts[i]]][rownames(Bulk),colnames(Bulk)]
  b.exp <- bmind.exp[,cts[i],][rownames(Bulk),colnames(Bulk)]
  s.exp <- STd.exp[cts[i],rownames(Bulk)]
  d.exp <- destvi.exp[[cts[i]]][rownames(Bulk),colnames(Bulk)]
  gt[[cts[i]]] <- gt[[cts[i]]][rownames(Bulk),colnames(Bulk)]
  e.trace.exp <- egm.trace.exp[,,cts[i]][rownames(Bulk),colnames(Bulk)]
  e.L2.exp <- egm.L2.exp[,,cts[i]][rownames(Bulk),colnames(Bulk)] 
  t.corr <- c()
  b.corr <- c()
  s.corr <- c()
  d.corr <- c()
  e.trace.corr <- c()
  e.L2.corr <- c()
  for (pixel in colnames(Bulk)) {
    t.corr <- c(t.corr,cor(gt[[cts[i]]][,pixel],t.exp[,pixel],method="sp"))
    b.corr <- c(b.corr,cor(gt[[cts[i]]][,pixel],b.exp[,pixel],method="sp"))
    s.corr <- c(s.corr,cor(gt[[cts[i]]][,pixel],s.exp,method="sp"))
    d.corr <- c(d.corr,cor(gt[[cts[i]]][,pixel],d.exp[,pixel],method="sp"))
    e.trace.corr <- c(e.trace.corr,cor(gt[[cts[i]]][,pixel],e.trace.exp[,pixel],method="sp"))
    e.L2.corr <- c(e.L2.corr,cor(gt[[cts[i]]][,pixel],e.L2.exp[,pixel],method="sp"))
  }
  corr <- cbind(t.corr,b.corr,s.corr,d.corr,e.trace.corr,e.L2.corr)
  colnames(corr) <- c("TCA","bMIND","STdeconvolve","DestVI","ENIGMA (trace)","ENIGMA")
  rownames(corr) <- colnames(Bulk)
  corr_pixel[[i]] <- corr
}
names(corr_pixel) <- cts
############################# plot Figure 2c ############################
library(ggplot2)
library(ggsci)
library(ggpubr)
df_box_gene <- reshape2::melt(corr_gene)
df_box_sample <- reshape2::melt(corr_pixel)
colnames(df_box_gene) <- c("Gene","Method","Value","Celltype")
colnames(df_box_sample) <- c("Gene","Method","Value","Celltype")
df_box_sample$Celltype <- stringi::stri_replace_all_fixed(df_box_sample$Celltype,"LuminalCell","Luminal cell",vectorize_all = T)
df_box_gene$Celltype <- stringi::stri_replace_all_fixed(df_box_gene$Celltype,"LuminalCell","Luminal cell",vectorize_all = T)

p2c.1 <- 
  df_box_sample %>%
  ggplot() +
  geom_boxplot(aes(x = Method,y = Value,fill = Method))+
  theme_minimal()+
  labs(y="Correlation per pixel",fill="Method")+
  theme(legend.text = element_text(size=12),
        legend.key.size = unit(0.6, "cm"),
        legend.title = element_blank(),
        legend.position = "top",
        axis.ticks.y = element_line(),
        axis.line.y  = element_line(),
        legend.direction = "horizontal",
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA),
        strip.text.x = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 2, bycol = TRUE))+
  scale_fill_manual(values = cols)+
  facet_grid(~Celltype, scales = "free_x",switch = "x")+
  scale_y_continuous(limits = c(0,1),breaks = c(0.00,0.25,0.5,0.75,1.00))

p2c.2 <- 
  df_box_gene %>%
  ggplot() +
  geom_boxplot(aes(x = Method,y = Value,fill = Method))+
  theme_minimal()+
  labs(y="Correlation per gene",fill="Method")+
  theme(legend.text = element_text(size=12),
        legend.key.size = unit(0.6, "cm"),
        legend.title = element_blank(),
        legend.position = "top",
        axis.ticks.y = element_line(),
        axis.line.y  = element_line(),
        legend.direction = "horizontal",
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA),
        strip.text.x = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 2, bycol = TRUE))+
  scale_fill_manual(values = cols)+
  facet_grid(~Celltype, scales = "free_x",switch = "x")+
  scale_y_continuous(limits = c(-1,1),breaks = c(-1,-0.50,0.00,0.5,1.00))

fig2c <- ggarrange(p2c.1,p2c.2,ncol = 1,nrow = 2,common.legend = T)+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
