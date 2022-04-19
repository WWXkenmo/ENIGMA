source("/mnt/data1/weixu/HiDe/revised/ENIGMA.R")
library(TCA)
library(MIND)
library(Biobase)
library(SingleCellExperiment)

dataNSCLC <- readRDS("/mnt/data1/weixu/HiDe/dataNSCLC.rds")
ref_sc <- readRDS("/mnt/data1/weixu/HiDe/ref.rds")

#We used the third patients to generate reference
ref_sc_sub <- ref_sc[,ref_sc$PatientID %in% "3" == TRUE]
ref_sc_sub <- ref_sc_sub[,ref_sc_sub$CellFromTumor %in% "1"]
ref_sc_sub <- ref_sc_sub[,ref_sc_sub$main_celltype %in% c("Alveolar","Epi") == FALSE]

Bulk <- as.matrix(dataNSCLC[[5]])
Tumor <- as.matrix(dataNSCLC[[1]])
Immune <- as.matrix(dataNSCLC[[2]])
Endothelial <- as.matrix(dataNSCLC[[3]])
Fibroblast <- as.matrix(dataNSCLC[[4]])
label <- dataNSCLC[[6]]
# The pheno variable contain the label of each samples (LUSC vs LUAD)
names(label) <- colnames(Tumor)

tmp <- readRDS("/mnt/data1/weixu/HiDe/tmp.rds")
colnames(tmp$main_celltype) <- c("B_cell","EC","Fibro","Myeloid","T_cell","Tumor")
Fra <- get_proportion(X=Bulk,ref=tmp$main_celltype)

egm_trace1 <- cell_deconvolve_trace(as.matrix(log2(Bulk+1)[rownames(tmp$main_celltype),]),Fra$theta,as.matrix(log2(tmp$main_celltype+1)),beta=1,verbose=FALSE,gamma=1,alpha=0.9,Normalize=FALSE,max.iter=1000)
egm_l2max1 <- cell_deconvolve(X=as.matrix(log2(Bulk+1)[rownames(tmp$main_celltype),]),
                    theta=Fra$theta,
					R=as.matrix(log2(tmp$main_celltype+1)),
					epsilon = 0.001,
					alpha=0.9,
					beta=0.5,tao_k=0.01,max.iter=1000,verbose=FALSE,Normalize=FALSE)

Bulk_f <- Bulk[rownames(tmp$main_celltype),][apply(Bulk[rownames(tmp$main_celltype),],1,var)>10^-8,]
tca.mdl <- tca(X = log2(Bulk_f+1), W = Fra$theta, C1 = NULL, C2 = NULL,
                parallel = TRUE,num_cores=3);
Z_hat <- tensor(X = log2(Bulk_f+1), tca.mdl)
##########################################################
###Estimation of bMIND
deconv = bMIND(log2(Bulk[rownames(tmp$main_celltype),]+1), frac = Fra$theta, profile=log2(tmp$main_celltype+1),  ncore = 5)

###
# Get the profile matrix through aggregate the profile
metadata = pData(ref_sc_sub)
Exp = exprs(ref_sc_sub)
profileMatrix = NULL
for(ct in names(table(metadata$main_celltype))){
  profileMatrix = cbind(profileMatrix, rowMeans(Exp[,metadata$main_celltype %in% ct]))
}
colnames(profileMatrix) <- names(table(metadata$main_celltype))

Fra2 <- get_proportion(X=Bulk,ref=profileMatrix)
deconv2 = bMIND(log2(Bulk[rownames(tmp$main_celltype),]+1), frac = Fra2$theta, profile=profileMatrix[rownames(tmp$main_celltype),],  ncore = 5)


# Raw pipline of bMIND
# Get the profile matrix through aggregate the profile
metadata <- pData(ref_sc_sub)
cell_type <- names(table(metadata$main_celltype))
K <- length(cell_type)
sc <- exprs(ref_sc_sub)
profile = matrix(NA, nrow(sc), K)
rownames(profile) = rownames(sc)
colnames(profile) = cell_type
for (i in cell_type) {
	profile[, i] = log2(edgeR::cpm(rowMeans(sc[, metadata$main_celltype ==
		i])) + 1)
}
colnames(profile) <- names(table(metadata$main_celltype))
geneID <- intersect(rownames(profile),rownames(Bulk))
colnames(profile) <- c("B_cell","EC","Fibro","Myeloid","T_cell","Tumor")
frac = MIND::est_frac(sig = profile[geneID,], bulk = log2(1+Bulk[geneID,]))
deconv3 = bMIND(log2(Bulk[geneID,]+1),ncore = 7, profile = as.matrix(profile[geneID,]),frac=frac)


##############################################################
##Compare the performance between three bMIND results

tumor_bmind <- NULL
for(i in 1:ncol(Bulk)) tumor_bmind <- c(tumor_bmind, cor(deconv$A[,6,i],Tumor[rownames(deconv$A[,6,]),i],method="sp"))

tumor_bmind2 <- NULL
for(i in 1:ncol(Bulk)) tumor_bmind2 <- c(tumor_bmind2, cor(deconv2$A[,6,i],Tumor[rownames(deconv2$A[,6,]),i],method="sp"))

tumor_bmind3 <- NULL
for(i in 1:ncol(Bulk)) tumor_bmind3 <- c(tumor_bmind3, cor(deconv3$A[,6,i],Tumor[rownames(deconv3$A[,6,]),i],method="sp"))

boxplot(tumor_bmind, tumor_bmind2, tumor_bmind3, col=terrain.colors(3),names = c("Model1","Model2","Model3"))

Tumor.gene <- rownames(Tumor)[apply(Tumor, 1, mean)>2]
Tumor.gene <- intersect(Tumor.gene,rownames(tmp$main_celltype))
gene <- intersect(rownames(deconv$A[,6,]),Tumor.gene)

tumor_bmind <- NULL
for(i in gene) tumor_bmind <- c(tumor_bmind, cor(deconv$A[i,6,],Tumor[i,],method="sp"))
tumor_bmind2 <- NULL
for(i in gene) tumor_bmind2 <- c(tumor_bmind2, cor(deconv2$A[i,6,],Tumor[i,],method="sp"))
tumor_bmind3 <- NULL
for(i in gene) tumor_bmind3 <- c(tumor_bmind3, cor(deconv3$A[i,6,],Tumor[i,],method="sp"))

png("bMIND_compare(gene).png",res=300,height=1500,width=1500)
boxplot(tumor_bmind, tumor_bmind2, tumor_bmind3, col=terrain.colors(3),names = c("Model1","Model2","Model3"))
dev.off()
###############################################################
###Estimate ISOpure
fra <- Fra$theta[,-6]
fra <- diag(1/rowSums(fra)) %*% fra
rownames(fra) <- rownames(Fra$theta)
N1_profile <- tmp$main_celltype[,-ncol(tmp$main_celltype)] %*% t(fra)

ISOpureS1model <- ISOpure.step1.CPE(Bulk[rownames(tmp$main_celltype),], N1_profile, logging.level = 'WARN');
ISOpureS2model <- ISOpure.step2.PPE(
    Bulk[rownames(tmp$main_celltype),],
    N1_profile,
    ISOpureS1model,logging.level = 'WARN'
)

# ISOpureS2model <- readRDS("isopure.rds")
rownames(ISOpureS2model$cc_cancerprofiles) <- rownames(tmp$main_celltype)
#############################################
fra <- Fra$theta[,-6]
rownames(fra) <- rownames(Fra$theta)
N1_profile <- tmp$main_celltype[,-ncol(tmp$main_celltype)] %*% t(fra)
#seq <- rowMeans(Bulk[rownames(tmp$main_celltype),])
#seq2 <- rowMeans(N1_profile[seq,])
sampleID <- colnames(Bulk)
colnames(N1_profile) <- colnames(Bulk) <- paste0("Sample ",1:ncol(Bulk),sep="")

#bulk_se <- SummarizedExperiment(assays=list(counts=log2(Bulk+1)[rownames(tmp$main_celltype),][seq, ][seq2,]))
#N_se <- SummarizedExperiment(assays=list(counts=as.matrix(log2(N1_profile+1)[seq, ][seq2,])))
bulk_se <- SummarizedExperiment(assays=list(counts=Bulk[rownames(tmp$main_celltype),] + 1))
N_se <- SummarizedExperiment(assays=list(counts=as.matrix(N1_profile + 1)))
fra_n <- as.numeric(rowSums(Fra$theta[,-6]))
names(fra_n) <- rownames(Fra$theta)

res.Demixt <- DeMixT_S2(data.Y = bulk_se,
                    data.N1 = N_se,
                    data.N2 = NULL,
                    givenpi = c(fra_n),nthread=5)

# res.Demixt <- readRDS("res.Demixt.rds")					
##################################
####load CIBERSORTx estimation
CBS <- read.table("/mnt/data1/weixu/HiDe/revised/Trace_Norm_test/cbs/CIBERSORTxHiRes_Job38_Tumor_Window8.txt",header=TRUE,row.names=1,sep="\t")
CBS <- na.omit(CBS)
CBS <- as.matrix(CBS)
CBS <- CBS[rowMeans(CBS)!=1,]

####Baseline
sf <- apply(Fra$theta,1,function(x){norm(x,"2")^2})
Bulk_tumor <- Bulk %*% diag(Fra$theta[,6]/sf)
Bulk_fibro <- Bulk %*% diag(Fra$theta[,3]/sf)
Bulk_endo <- Bulk %*% diag(Fra$theta[,2]/sf)

#####################################################################
##Evaluation
##sample
#Reference
tumor_ref <- NULL
for(i in 1:ncol(Bulk)) tumor_ref <- c(tumor_ref, cor(tmp$main_celltype[,6],Tumor[rownames(tmp$main_celltype),i],method="sp"))
endo_ref <- NULL
for(i in 1:ncol(Bulk)) endo_ref <- c(endo_ref, cor(tmp$main_celltype[,2],Endothelial[rownames(tmp$main_celltype),i],method="sp"))
fibro_ref <- NULL
for(i in 1:ncol(Bulk)) fibro_ref <- c(fibro_ref, cor(tmp$main_celltype[,3],Fibroblast[rownames(tmp$main_celltype),i],method="sp"))

#Bulk
tumor_bulk <- NULL
for(i in 1:ncol(Bulk)) tumor_bulk <- c(tumor_bulk, cor(Bulk_tumor[rownames(tmp$main_celltype),i],Tumor[rownames(tmp$main_celltype),i],method="sp"))
endo_bulk <- NULL
for(i in 1:ncol(Bulk)) endo_bulk <- c(endo_bulk, cor(Bulk_endo[rownames(tmp$main_celltype),i],Endothelial[rownames(tmp$main_celltype),i],method="sp"))
fibro_bulk <- NULL
for(i in 1:ncol(Bulk)) fibro_bulk <- c(fibro_bulk, cor(Bulk_fibro[rownames(tmp$main_celltype),i],Fibroblast[rownames(tmp$main_celltype),i],method="sp"))

#ENIGMA (trace and l2max)
tumor_trace1 <- NULL
for(i in 1:ncol(Bulk)) tumor_trace1 <- c(tumor_trace1, cor(egm_trace1$X_k[,i,6],Tumor[rownames(tmp$main_celltype),i],method="sp"))
endo_trace1 <- NULL
for(i in 1:ncol(Bulk)) endo_trace1 <- c(endo_trace1, cor(egm_trace1$X_k[,i,2],Endothelial[rownames(tmp$main_celltype),i],method="sp"))
fibro_trace1 <- NULL
for(i in 1:ncol(Bulk)) fibro_trace1 <- c(fibro_trace1, cor(egm_trace1$X_k[,i,3],Fibroblast[rownames(tmp$main_celltype),i],method="sp"))

tumor_l2max1 <- NULL
for(i in 1:ncol(Bulk)) tumor_l2max1 <- c(tumor_l2max1, cor(egm_l2max1$X_k[,i,6],Tumor[rownames(tmp$main_celltype),i],method="sp"))
endo_l2max1 <- NULL
for(i in 1:ncol(Bulk)) endo_l2max1 <- c(endo_l2max1, cor(egm_l2max1$X_k[,i,2],Endothelial[rownames(tmp$main_celltype),i],method="sp"))
fibro_l2max1 <- NULL
for(i in 1:ncol(Bulk)) fibro_l2max1 <- c(fibro_l2max1, cor(egm_l2max1$X_k[,i,3],Fibroblast[rownames(tmp$main_celltype),i],method="sp"))

#bMIND
tumor_bmind <- NULL
for(i in 1:ncol(Bulk)) tumor_bmind <- c(tumor_bmind, cor(deconv3$A[,6,i],Tumor[rownames(deconv3$A[,6,]),i],method="sp"))
endo_bmind <- NULL
for(i in 1:ncol(Bulk)) endo_bmind <- c(endo_bmind, cor(deconv3$A[,2,i],Endothelial[rownames(deconv3$A[,6,]),i],method="sp"))
fibro_bmind <- NULL
for(i in 1:ncol(Bulk)) fibro_bmind <- c(fibro_bmind, cor(deconv3$A[,3,i],Fibroblast[rownames(deconv3$A[,6,]),i],method="sp"))

#TCA
tumor_tca <- NULL
for(i in 1:ncol(Bulk)) tumor_tca <- c(tumor_tca, cor(Z_hat[[6]][,i],Tumor[rownames(Z_hat[[6]]),i],method="sp"))
endo_tca <- NULL
for(i in 1:ncol(Bulk)) endo_tca <- c(endo_tca, cor(Z_hat[[2]][,i],Endothelial[rownames(Z_hat[[6]]),i],method="sp"))
fibro_tca <- NULL
for(i in 1:ncol(Bulk)) fibro_tca <- c(fibro_tca, cor(Z_hat[[3]][,i],Fibroblast[rownames(Z_hat[[6]]),i],method="sp"))

#CBS
tumor_cbs <- NULL
for(i in 1:ncol(Bulk)) tumor_cbs <- c(tumor_cbs, cor(CBS[,i],Tumor[rownames(CBS),i],method="sp"))

#ISOpure
tumor_isopure <- NULL
for(i in 1:ncol(Bulk)) tumor_isopure <- c(tumor_isopure, cor(ISOpureS2model$cc_cancerprofiles[,i],Tumor[rownames(ISOpureS2model$cc_cancerprofiles),i],method="sp"))

#Demixt
tumor_demixt <- NULL
for(i in 1:ncol(Tumor)) tumor_demixt <- c(tumor_demixt, cor(res.Demixt$decovExprT[,i],Tumor[rownames(res.Demixt$decovExprT),i],method="sp"))

#Collect two table
res <- cbind(tumor_bulk,tumor_trace1,tumor_l2max1,tumor_bmind,tumor_tca,tumor_cbs,tumor_isopure,tumor_demixt,tumor_ref)
colnames(res) <- c("Bulk(Baseline)","ENIGMA(trace norm)","ENIGMA","bMIND","TCA","CIBERSORTx(2-comp)","ISOpure(2-comp)","DeMixT(2-comp)","Reference")
perfor <- as.data.frame(res)

val <- as.numeric(as.matrix(perfor))
method <- c(rep("Bulk(Baseline)",24),rep("ENIGMA(trace norm)",24),rep("ENIGMA",24),rep("bMIND",24),rep("TCA",24),rep("CIBERSORTx(2-comp)",24),rep("ISOpure(2-comp)",24),rep("DeMixT(2-comp)",24),rep("Reference",24))
res <- cbind(method,val)
colnames(res) <- c("method","SCC")
res <- as.data.frame(res)
res$SCC <- as.numeric(as.matrix(res$SCC))
res_all <- res
saveRDS(res_all,file="result_to_plot_tumor.rds")

########################################
##plot all
res <- cbind(fibro_bulk,fibro_trace1,fibro_l2max1,fibro_bmind,fibro_tca,fibro_ref)
colnames(res) <- c("Bulk(Baseline)","ENIGMA(trace norm)","ENIGMA","bMIND","TCA","Reference")
perfor <- as.data.frame(res)

val <- as.numeric(as.matrix(perfor))
method <- c(rep("Bulk(Baseline)",24),rep("ENIGMA(trace norm)",24),rep("ENIGMA",24),rep("bMIND",24),rep("TCA",24),rep("Reference",24))
res <- cbind(method,val)
colnames(res) <- c("method","SCC")
res <- as.data.frame(res)
res$SCC <- as.numeric(as.matrix(res$SCC))
res_all <- res
saveRDS(res_all,file="result_to_plot_fibro.rds")
################################################
res <- cbind(endo_bulk,endo_trace1,endo_l2max1,endo_bmind,endo_tca,endo_ref)
colnames(res) <- c("Bulk(Baseline)","ENIGMA(trace norm)","ENIGMA","bMIND","TCA","Reference")
perfor <- as.data.frame(res)

val <- as.numeric(as.matrix(perfor))
method <- c(rep("Bulk(Baseline)",24),rep("ENIGMA(trace norm)",24),rep("ENIGMA",24),rep("bMIND",24),rep("TCA",24),rep("Reference",24))
res <- cbind(method,val)
colnames(res) <- c("method","SCC")
res <- as.data.frame(res)
res$SCC <- as.numeric(as.matrix(res$SCC))
res_all <- res
saveRDS(res_all,file="result_to_plot_endo.rds")


###########################################
Tumor.gene <- rownames(Tumor)[apply(Tumor, 1, mean)>=0]
Tumor.gene <- intersect(Tumor.gene,rownames(tmp$main_celltype))

Fibro.gene <- rownames(Fibroblast)[apply(Fibroblast, 1, mean)>=0]
Fibro.gene <- intersect(Fibro.gene,rownames(tmp$main_celltype))

Endo.gene <- rownames(Endothelial)[apply(Endothelial, 1, mean)>=0]
Endo.gene <- intersect(Endo.gene,rownames(tmp$main_celltype))

#################################################
##Baseline
tumor_bulk <- NULL
for(i in Tumor.gene) tumor_bulk <- c(tumor_bulk, cor(Bulk[i,],Tumor[i,],method="sp"))
fibro_bulk <- NULL
for(i in Fibro.gene) fibro_bulk <- c(fibro_bulk, cor(Bulk_fibro[i,],Fibroblast[i,],method="sp"))
endo_bulk <- NULL
for(i in Endo.gene) endo_bulk <- c(endo_bulk, cor(Bulk_endo[i,],Endothelial[i,],method="sp"))

####ENIGMA
tumor_trace1 <- NULL
for(i in Tumor.gene) tumor_trace1 <- c(tumor_trace1, cor(egm_trace1$X_k[i,,6],Tumor[i,],method="sp"))
tumor_l2max1 <- NULL
for(i in Tumor.gene) tumor_l2max1 <- c(tumor_l2max1, cor(egm_l2max1$X_k[i,,6],Tumor[i,],method="sp"))

fibro_trace1 <- NULL
for(i in Fibro.gene) fibro_trace1 <- c(fibro_trace1, cor(egm_trace1$X_k[i,,3],Fibroblast[i,],method="sp"))
fibro_l2max1 <- NULL
for(i in Fibro.gene) fibro_l2max1 <- c(fibro_l2max1, cor(egm_l2max1$X_k[i,,3],Fibroblast[i,],method="sp"))

endo_trace1 <- NULL
for(i in Endo.gene) endo_trace1 <- c(endo_trace1, cor(egm_trace1$X_k[i,,2],Endothelial[i,],method="sp"))
endo_l2max1 <- NULL
for(i in Endo.gene) endo_l2max1 <- c(endo_l2max1, cor(egm_l2max1$X_k[i,,2],Endothelial[i,],method="sp"))

###bMIND
gene <- intersect(rownames(deconv3$A[,6,]),Tumor.gene)
tumor_bmind <- NULL
for(i in gene) tumor_bmind <- c(tumor_bmind, cor(deconv3$A[i,6,],Tumor[i,],method="sp"))

gene <- intersect(rownames(deconv3$A[,3,]),Fibro.gene)
fibro_bmind <- NULL
for(i in gene) fibro_bmind <- c(fibro_bmind, cor(deconv3$A[i,3,],Fibroblast[i,],method="sp"))

gene <- intersect(rownames(deconv3$A[,2,]),Endo.gene)
endo_bmind <- NULL
for(i in gene) endo_bmind <- c(endo_bmind, cor(deconv3$A[i,2,],Endothelial[i,],method="sp"))

###tca
gene <- intersect(rownames(Z_hat[[6]]),Tumor.gene)
tumor_tca <- NULL
for(i in gene) tumor_tca <- c(tumor_tca, cor(Z_hat[[6]][i,],Tumor[i,],method="sp"))

gene <- intersect(rownames(Z_hat[[6]]),Fibro.gene)
fibro_tca <- NULL
for(i in gene) fibro_tca <- c(fibro_tca, cor(Z_hat[[3]][i,],Fibroblast[i,],method="sp"))

gene <- intersect(rownames(Z_hat[[6]]),Endo.gene)
endo_tca <- NULL
for(i in gene) endo_tca <- c(endo_tca, cor(Z_hat[[2]][i,],Endothelial[i,],method="sp"))

###CIBERSORTx
gene <- intersect(rownames(CBS),Tumor.gene)
tumor_cbs <- NULL
for(i in gene) tumor_cbs <- c(tumor_cbs, cor(CBS[i,],Tumor[i,],method="sp"))

###ISOpure
gene <- intersect(rownames(ISOpureS2model$cc_cancerprofiles),Tumor.gene)
tumor_isopure <- NULL
for(i in gene) tumor_isopure <- c(tumor_isopure, cor(ISOpureS2model$cc_cancerprofiles[i,],Tumor[i,],method="sp"))

###DeMixT
gene <- intersect(rownames(res.Demixt$decovExprT),Tumor.gene)
tumor_demixt <- NULL
for(i in gene) tumor_demixt <- c(tumor_demixt, cor(res.Demixt$decovExprT[i,],Tumor[i,],method="sp"))

##return
val <- c(tumor_bulk,tumor_trace1,tumor_l2max1,tumor_bmind,tumor_tca,tumor_cbs,tumor_isopure,tumor_demixt)
val[is.na(val)] <- 0
method <- c(rep("Bulk(Baseline)",length(tumor_bulk)),rep("ENIGMA(trace norm)",length(tumor_trace1)),rep("ENIGMA",length(tumor_l2max1)),rep("bMIND",length(tumor_bmind)),rep("TCA",length(tumor_tca)),rep("CIBERSORTx(2-comp)",length(tumor_cbs)),rep("ISOpure(2-comp)",length(tumor_isopure)),rep("DeMixT(2-comp)",length(tumor_demixt)))
res <- cbind(method,val)
colnames(res) <- c("method","SCC")
res <- as.data.frame(res)
res$SCC <- as.numeric(as.matrix(res$SCC))
res_all <- res
saveRDS(res_all,file="result_to_plot_tumor(gene_all).rds")


##############
val <- c(fibro_bulk,fibro_trace1,fibro_l2max1,fibro_bmind,fibro_tca)
val[is.na(val)] <- 0
method <- c(rep("Bulk(Baseline)",length(fibro_bulk)),rep("ENIGMA(trace norm)",length(fibro_trace1)),rep("ENIGMA",length(fibro_l2max1)),rep("bMIND",length(fibro_bmind)),rep("TCA",length(fibro_tca)))
res <- cbind(method,val)
colnames(res) <- c("method","SCC")
res <- as.data.frame(res)
res$SCC <- as.numeric(as.matrix(res$SCC))
res_all <- res
saveRDS(res_all,file="result_to_plot_fibro(gene_all).rds")


##############
val <- c(endo_bulk,endo_trace1,endo_l2max1,endo_bmind,endo_tca)
val[is.na(val)] <- 0
method <- c(rep("Bulk(Baseline)",length(endo_bulk)),rep("ENIGMA(trace norm)",length(endo_trace1)),rep("ENIGMA",length(endo_l2max1)),rep("bMIND",length(endo_bmind)),rep("TCA",length(endo_tca)))
res <- cbind(method,val)
colnames(res) <- c("method","SCC")
res <- as.data.frame(res)
res$SCC <- as.numeric(as.matrix(res$SCC))
res_all <- res
saveRDS(res_all,file="result_to_plot_endo(gene_all).rds")

########################################
B <- readRDS("result_to_plot_tumor.rds")
B$method[B$method == "Bulk(Baseline)"] <- "Baseline"
B$method[B$method == "ENIGMA(trace norm)"] <- "ENIGMA(trace)"
mycols <- c("#a4c639","#9966cc","#F00000","#FE82AB","#ffbf00","#006699","#538b8b","#999999","#F5B184")
method <- names(table(B$method))
names(mycols) <- method
median <- NULL; for(i in method) median <- c(median,median(na.omit(B$SCC[B$method == i])))
names(median) <- method
B$method <- factor(B$method, levels = method[order(median,decreasing=TRUE)])
mycols = mycols[names(table(B$method))]

png("NSCLC_sample_plot.png",res=300,height=1500,width=2500)
ggplot(B, aes(method,SCC,fill = method)) +
  geom_boxplot()  +
  theme_classic() + 
  scale_fill_manual(values= as.character(mycols)) +
  ylab("Spearman correlation against ground truth \n (FACS-purified tumor cells") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + ylim(0,1.1)+
  theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(), axis.text.y = element_text(size = 15,  face = "bold"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 15, face = "bold",vjust = 0.9, hjust = 0.5,angle = 90)) +
  theme(axis.ticks.length = unit(7, "pt"),axis.line = element_line(size = 1), axis.ticks = element_line(colour = "black", size = 1), panel.grid.major = element_line(colour = NA, size = 1), panel.grid.minor = element_line(colour = NA, size = 1)) + 
  theme(legend.text = element_text(size = 16),legend.title = element_text(size = 18), legend.background = element_rect(fill = "white"))
dev.off()


B <- readRDS("result_to_plot_tumor(gene).rds")
B$method[B$method == "Bulk(Baseline)"] <- "Baseline"
B$method[B$method == "ENIGMA(trace norm)"] <- "ENIGMA(trace)"
B$SCC[is.na(B$SCC)] <- 0
method <- names(table(B$method))
median <- NULL; for(i in method) median <- c(median,median(na.omit(B$SCC[B$method == i])))
names(median) <- method 
B$method <- factor(B$method, levels = method[order(median,decreasing=TRUE)])
mycols = mycols[names(table(B$method))]

png("NSCLC_gene_plot(tumor).png",res=300,height=1500,width=2500)
ggplot(B, aes(method,SCC,fill = method)) +
  geom_boxplot()  +
  theme_classic() + 
  scale_fill_manual(values= as.character(mycols)) +
  ylab("Spearman correlation against ground truth \n (FACS-purified fibro cells") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )+
  theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(), axis.text.y = element_text(size = 15,  face = "bold"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 15, face = "bold",vjust = 0.9, hjust = 0.5,angle = 90)) +
  theme(axis.ticks.length = unit(7, "pt"),axis.line = element_line(size = 1), axis.ticks = element_line(colour = "black", size = 1), panel.grid.major = element_line(colour = NA, size = 1), panel.grid.minor = element_line(colour = NA, size = 1)) + 
  theme(legend.text = element_text(size = 16),legend.title = element_text(size = 18), legend.background = element_rect(fill = "white"))
dev.off()
###################################################################

##function
quantile_thr <- function(seq,vec){
    vec <- vec[order(vec,decreasing=TRUE)]
    vec <- vec[1:seq]
    value <- vec[seq]
    value
}
sensitivity <- function(lm.test,list){
    ### sensitivity(recall) is defined as the num(true positive predictions)/num(true) 
    pos <- list[list>0]
    neg <- abs(list[list<0])
    sensitivity <- (length(intersect(which(lm.test==1),pos))+length(intersect(which(lm.test== -1),neg)))/length(list)
    sensitivity
}
specificity <- function(lm.test,list){
    sp <- sum(which(lm.test==0) %in% abs(list) == FALSE)/(1000-100)
    sp
}
precision <- function(lm.test,list){
    ### precision(1-Sensitivity) is defined as the num(true positive predictions)/num(positive prediction) 
    pos <- list[list>0]
    neg <- abs(list[list<0])
    precision <- (length(intersect(which(lm.test==1),pos))+length(intersect(which(lm.test== -1),neg)))/sum(lm.test!=0)
    precision
}
######################################################################
geneID <- rownames(Z_hat[[6]])
Tumor <- as.matrix(Tumor)
Tumor_table <- NULL
for(i in geneID){
   pvalue <- wilcox.test(Tumor[i,label==1],Tumor[i,label==0])$p.value
   Tumor_table <- rbind(Tumor_table,cbind(pvalue, mean(Tumor[i,label==1])/mean(Tumor[i,label==0])))
}
rownames(Tumor_table) <- geneID

colnames(Tumor_table) <- c("pvalue","FoldChange")
Tumor_table <- as.data.frame(Tumor_table)
Tumor_table$p.adjust <- p.adjust(Tumor_table$pvalue,method="BH")

DEG_list <- c(which(Tumor_table$p.adjust<0.01&Tumor_table$FoldChange>1),-which(Tumor_table$p.adjust<0.01&Tumor_table$FoldChange<1))


################################
y <- label
lm.test <- NULL
fc <- NULL
mat <- sqrt(Bulk)[geneID,]
for(i in 1:nrow(mat)){
 if(sum(mat[i,])==0){lm.test <- c(lm.test,1); fc <- c(fc,0)}else{
 pval <- summary(lm(scale(mat[i,])~as.numeric(y)))$coefficients[2,4]
 fc <- c(fc, summary(lm(scale(mat[i,])~as.numeric(y)))$coefficients[2,1])
 lm.test <- c(lm.test,pval)}
}
lm.test <- p.adjust(lm.test,method="fdr")
lm.test[lm.test==1] <- 0.9999
fc_test <- abs(fc)
fc_test[lm.test>0.05] <- 0
qtile <-  apply(as.matrix(c(20,50,100,120,150,min(sum(lm.test<=0.05),length(DEG_list)))),1,quantile_thr,vec=fc_test)
seq_sen <- seq_sep <- seq_pre <- NULL
pval <- lm.test
for(i in 1:length(qtile)){
#lm.test[order(lm.test,decreasing=FALSE)[1:50]] <- 1
lm.test[lm.test<0.05&fc>qtile[i]] <- 1
lm.test[lm.test<0.05&fc< -qtile[i]] <- -1
lm.test[abs(lm.test)!=1] <- 0
seq_sen <- c(seq_sen,sensitivity(lm.test,DEG_list))
seq_sep <- c(seq_sep,sum(which(lm.test==0) %in% abs(DEG_list) == FALSE)/(nrow(Tumor_table)-length(DEG_list)))
seq_pre <- c(seq_pre,precision(lm.test,DEG_list))
lm.test <- pval
}
seq_sen_bulk <- seq_sen
seq_pre_bulk <- seq_pre

####
y <- label
lm.test <- NULL
fc <- NULL
mat <- egm_trace1$X_k[geneID,,6]
for(i in 1:nrow(mat)){
 if(sum(mat[i,])==0){lm.test <- c(lm.test,1); fc <- c(fc,0)}else{
 pval <- summary(lm(scale(mat[i,])~as.numeric(y)))$coefficients[2,4]
 fc <- c(fc, summary(lm(scale(mat[i,])~as.numeric(y)))$coefficients[2,1])
 lm.test <- c(lm.test,pval)}
}
lm.test <- p.adjust(lm.test,method="fdr")
lm.test[lm.test==1] <- 0.9999
fc_test <- abs(fc)
fc_test[lm.test>0.05] <- 0
qtile <-  apply(as.matrix(c(20,50,100,120,150,min(sum(lm.test<=0.05),length(DEG_list)))),1,quantile_thr,vec=fc_test)
seq_sen <- seq_sep <- seq_pre <- NULL
pval <- lm.test
for(i in 1:length(qtile)){
#lm.test[order(lm.test,decreasing=FALSE)[1:50]] <- 1
lm.test[lm.test<0.05&fc>qtile[i]] <- 1
lm.test[lm.test<0.05&fc< -qtile[i]] <- -1
lm.test[abs(lm.test)!=1] <- 0
seq_sen <- c(seq_sen,sensitivity(lm.test,DEG_list))
seq_sep <- c(seq_sep,sum(which(lm.test==0) %in% abs(DEG_list) == FALSE)/(nrow(Tumor_table)-length(DEG_list)))
seq_pre <- c(seq_pre,precision(lm.test,DEG_list))
lm.test <- pval
}
seq_sen_trace <- seq_sen
seq_pre_trace <- seq_pre


####
y <- label
lm.test <- NULL
fc <- NULL
mat <- egm_l2max1$X_k[geneID,,6]
for(i in 1:nrow(mat)){
 if(sum(mat[i,])==0){lm.test <- c(lm.test,1); fc <- c(fc,0)}else{
 pval <- summary(lm(scale(mat[i,])~as.numeric(y)))$coefficients[2,4]
 fc <- c(fc, summary(lm(scale(mat[i,])~as.numeric(y)))$coefficients[2,1])
 lm.test <- c(lm.test,pval)}
}
lm.test <- p.adjust(lm.test,method="fdr")
lm.test[lm.test==1] <- 0.9999
fc_test <- abs(fc)
fc_test[lm.test>0.05] <- 0
qtile <-  apply(as.matrix(c(20,50,100,120,150,min(sum(lm.test<=0.05),length(DEG_list)))),1,quantile_thr,vec=fc_test)
seq_sen <- seq_sep <- seq_pre <- NULL
pval <- lm.test
for(i in 1:length(qtile)){
#lm.test[order(lm.test,decreasing=FALSE)[1:50]] <- 1
lm.test[lm.test<0.05&fc>qtile[i]] <- 1
lm.test[lm.test<0.05&fc< -qtile[i]] <- -1
lm.test[abs(lm.test)!=1] <- 0
seq_sen <- c(seq_sen,sensitivity(lm.test,DEG_list))
seq_sep <- c(seq_sep,sum(which(lm.test==0) %in% abs(DEG_list) == FALSE)/(nrow(Tumor_table)-length(DEG_list)))
seq_pre <- c(seq_pre,precision(lm.test,DEG_list))
lm.test <- pval
}
seq_sen_l2max <- seq_sen
seq_pre_l2max <- seq_pre


####
y <- label
lm.test <- NULL
fc <- NULL
mat <- deconv$A[geneID,6,]
for(i in 1:nrow(mat)){
 if(sum(mat[i,])==0){lm.test <- c(lm.test,1); fc <- c(fc,0)}else{
 pval <- summary(lm(scale(mat[i,])~as.numeric(y)))$coefficients[2,4]
 fc <- c(fc, summary(lm(scale(mat[i,])~as.numeric(y)))$coefficients[2,1])
 lm.test <- c(lm.test,pval)}
}
lm.test <- p.adjust(lm.test,method="fdr")
lm.test[lm.test==1] <- 0.9999
fc_test <- abs(fc)
fc_test[lm.test>0.05] <- 0
qtile <-  apply(as.matrix(c(20,50,100,120,150,min(sum(lm.test<=0.05),length(DEG_list)))),1,quantile_thr,vec=fc_test)
seq_sen <- seq_sep <- seq_pre <- NULL
pval <- lm.test
for(i in 1:length(qtile)){
#lm.test[order(lm.test,decreasing=FALSE)[1:50]] <- 1
lm.test[lm.test<0.05&fc>qtile[i]] <- 1
lm.test[lm.test<0.05&fc< -qtile[i]] <- -1
lm.test[abs(lm.test)!=1] <- 0
seq_sen <- c(seq_sen,sensitivity(lm.test,DEG_list))
seq_sep <- c(seq_sep,sum(which(lm.test==0) %in% abs(DEG_list) == FALSE)/(nrow(Tumor_table)-length(DEG_list)))
seq_pre <- c(seq_pre,precision(lm.test,DEG_list))
lm.test <- pval
}
seq_sen_mind <- seq_sen
seq_pre_mind <- seq_pre

####
y <- label
lm.test <- NULL
fc <- NULL
mat <- Z_hat[[6]][geneID,]
for(i in 1:nrow(mat)){
 if(sum(mat[i,])==0){lm.test <- c(lm.test,1); fc <- c(fc,0)}else{
 pval <- summary(lm(scale(mat[i,])~as.numeric(y)))$coefficients[2,4]
 fc <- c(fc, summary(lm(scale(mat[i,])~as.numeric(y)))$coefficients[2,1])
 lm.test <- c(lm.test,pval)}
}
lm.test <- p.adjust(lm.test,method="fdr")
lm.test[lm.test==1] <- 0.9999
fc_test <- abs(fc)
fc_test[lm.test>0.05] <- 0
qtile <-  apply(as.matrix(c(20,50,100,120,150,min(sum(lm.test<=0.05),length(DEG_list)))),1,quantile_thr,vec=fc_test)
seq_sen <- seq_sep <- seq_pre <- NULL
pval <- lm.test
for(i in 1:length(qtile)){
#lm.test[order(lm.test,decreasing=FALSE)[1:50]] <- 1
lm.test[lm.test<0.05&fc>qtile[i]] <- 1
lm.test[lm.test<0.05&fc< -qtile[i]] <- -1
lm.test[abs(lm.test)!=1] <- 0
seq_sen <- c(seq_sen,sensitivity(lm.test,DEG_list))
seq_sep <- c(seq_sep,sum(which(lm.test==0) %in% abs(DEG_list) == FALSE)/(nrow(Tumor_table)-length(DEG_list)))
seq_pre <- c(seq_pre,precision(lm.test,DEG_list))
lm.test <- pval
}
seq_sen_tca <- seq_sen
seq_pre_tca <- seq_pre
seq_sen_tca[length(seq_sen_tca)] <- seq_sen_tca[length(seq_sen_tca)-1]
######################################################################
###plot the lines
GeneNumber <- c(rep(c("20","50","100","120","150","# DE"),5))
Method <- c(rep("Bulk",6),
            rep("TCA",6),
            rep("bMIND",6),
			rep("ENIGMA",6),
			rep("ENIGMA(trace norm)",6))
Value <- c(seq_sen_bulk,
           seq_sen_tca,
           seq_sen_mind,
		   seq_sen_l2max,
		   seq_sen_trace)
perfor <- data.frame(GeneNumber = GeneNumber,
                     Method = Method,
					 Value = Value)
perfor$GeneNumber <- factor(perfor$GeneNumber,levels=c("20","50","100","120","150","# DE"))
p_sen <-ggplot(data=perfor,aes(x=GeneNumber,y=Value,group=Method)) + geom_line(aes(colour=Method)) + geom_point(size=1,aes(colour=Method)) + xlab("")+ylab("Sensitivity") + theme_classic2() 

tiff("sensitivity.tiff",res=300,height=1500,width=1800)
p_sen
dev.off()














