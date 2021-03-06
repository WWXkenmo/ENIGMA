############################################
source("ENIGMA.R")
library(Seurat)
####Generating single cell gene expression reference profile from Rheumatoid Arthritis Patients
###the single cell datasets, FACS RNA-seq datasets of Arthritis Patients could be downloaded from ImmPort (https://www.immport.org/shared/study/SDY998 and https://www.immport.org/shared/study/SDY999, with study accession codes SDY998 and SDY999
cel_matrix <- read.table("celseq_matrix_ru10_molecules.tsv",header=TRUE,row.names=1,sep="\t")
metadata <- read.table("celseq_meta.tsv",header=TRUE,row.names=1,stringsAsFactors = FALSE,sep="\t")
cel_matrix <- as.matrix(cel_matrix)
cel_matrix[is.na(cel_matrix)] <- 0

sample <- intersect(rownames(metadata),colnames(cel_matrix))
cel_matrix <- cel_matrix[,sample]

cel_RA <- CreateSeuratObject(counts = cel_matrix, project = "pbmc3k", min.cells = 50, min.features = 100)
cel_RA[["percent.mt"]] <- PercentageFeatureSet(cel_RA, pattern = "^MT-")
cel_RA <- subset(cel_RA, subset = percent.mt < 20)

cel_RA <- NormalizeData(cel_RA,normalization.method="RC",scale.factor = 100000)
cel_RA <- as.matrix(GetAssayData(cel_RA))

metadata <- metadata[colnames(cel_RA),]
metadata <- metadata[metadata$disease == "RA",]

cel_RA <- cel_RA[,rownames(metadata)]
cel_RA <- ExpressionSet(cel_RA);pData(cel_RA) <- metadata
saveRDS(cel_RA,file="RA_sc.rds")

###################################################################################
###Generate DEG genes from scRNA-seq
cel_matrix <- read.table("celseq_matrix_ru10_molecules.tsv",header=TRUE,row.names=1,sep="\t")
metadata <- read.table("celseq_meta.tsv",header=TRUE,row.names=1,stringsAsFactors = FALSE,sep="\t")
cel_matrix <- as.matrix(cel_matrix)
cel_matrix[is.na(cel_matrix)] <- 0

sample <- intersect(rownames(metadata),colnames(cel_matrix))
cel_matrix <- cel_matrix[,sample]

cel_RA <- CreateSeuratObject(counts = cel_matrix, project = "pbmc3k", min.cells = 50, min.features = 100)
cel_RA[["percent.mt"]] <- PercentageFeatureSet(cel_RA, pattern = "^MT-")
cel_RA <- subset(cel_RA, subset = percent.mt < 20)

cel_RA <- NormalizeData(cel_RA,normalization.method="RC",scale.factor = 100000)
cel_RA <- as.matrix(GetAssayData(cel_RA))

metadata <- metadata[colnames(cel_RA),]
metadata <- metadata[metadata$disease %in% c("RA","OA"),]

cel_RA <- cel_RA[,rownames(metadata)]
cel_RA <- ExpressionSet(cel_RA);pData(cel_RA) <- metadata
saveRDS(cel_RA,file="RA_OA_scData.rds")
###########################################################################################

RA_scData2 <- readRDS("RA_OA_scData.rds")
disease_state <- RA_scData2$disease

RA <- RA_scData2[,disease_state=="RA"]
OA <- RA_scData2[,disease_state=="OA"]
RA <- RA[rownames(tmp_ra$type),]
OA <- OA[rownames(tmp_ra$type),]
################################
###Generating reference profile matrix from single cell expression data
###The bulk RNA-seq datasets for deconvolution could be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89408
RA_scData <- readRDS("RA_sc.rds")
bulk_ra_new <- read.table("GSE89408_GEO_count_matrix_rename.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)
sample_type <- strsplit(colnames(bulk_ra_new),split="_")
sampleType <- NULL
for(i in 1:length(sample_type)){
  sampleType <- c(sampleType,sample_type[[i]][1])
}
bulk_ra_vs_oa <- bulk_ra_new[,sampleType %in% c("OA","RA")] 
bulk_ra_vs_oa <- ExpressionSet(as.matrix(bulk_ra_vs_oa))

tmp_ra = remove_batch_effect(bulk_ra_vs_oa,RA_scData,"type",n=1000,ncores=4)

####
###Generating reference profile matrix from FACS RNA-seq expression data
library(clusterProfiler)
library(org.Hs.eg.db)
bulk_sort_ra <- read.table("low_input_gene_sample_tpm_matrix.725714.tsv",row.names=1,stringsAsFactors = FALSE)
meta_ra <- read.table("low_input_meta.725715.tsv",header=TRUE,row.names=1,stringsAsFactors = FALSE,sep="\t")

#bulk_sort_ra <- as.matrix(bulk_sort_ra)
meta_ra <- meta_ra[colnames(bulk_sort_ra),]
meta_ra <- subset(meta_ra, Tissue.type %in% c("OA-arthro","RA-arthro"))


genes <- bitr(rownames(bulk_sort_ra),fromType="ENSEMBL",toType="SYMBOL",org.Hs.eg.db)
genes <- genes[!duplicated(genes[,2]),]
bulk_sort_ra <- bulk_sort_ra[genes[,1],rownames(meta_ra)]
rownames(bulk_sort_ra) <- genes[,2]

####to fairly compare, we took the intersect genes shared by scRNA-seq reference and FACS RNA-seq reference
bulk_sort_ra <- bulk_sort_ra[intersect(rownames(tmp_ra$type),rownames(bulk_sort_ra)),]
bulk_sort_ra <- 2^bulk_sort_ra - 1

ref_ra_all <- NULL
for(i in names(table(meta_ra$Cell.type))){
    vec <- bulk_sort_ra[,meta_ra$Cell.type %in% i]
    ref_ra_all <- cbind(ref_ra_all,rowMeans(vec))
}
ref_ra <- NULL
for(i in names(table(meta_ra$Cell.type))){
    vec <- bulk_sort_ra[,meta_ra$Cell.type %in% i & meta_ra$Tissue.type %in% "RA-arthro"]
    ref_ra <- cbind(ref_ra,rowMeans(vec))
}

colnames(ref_ra) <- colnames(ref_ra_all) <- names(table(meta_ra$Cell.type))

###Estimating fraction using both RA and OA expression to make the prediction more accurate
Fra_ra2 <- get_proportion(exprs(bulk_ra_vs_oa),ref_ra_all)
bulk_correct_ra <- B_mode_batch_effect_remove(exprs(bulk_ra_vs_oa)[rownames(ref_ra),],ref_ra,Fra_ra2$theta)
Fra_ra_correct <- get_proportion(bulk_correct_ra,ref_ra)

###Estimating fractions merely based on RA information
Fra_ra3 <- get_proportion(exprs(bulk_ra_vs_oa),ref_ra)

##########maximum L2 norm model
####raw model (FACS reference)
res_alg_all_ra_l2_norm_sort <- cell_deconvolve(X=as.matrix(log2(exprs(bulk_ra_vs_oa)+1)),
                    theta=Fra_ra2$theta,
					R=as.matrix(log2(ref_ra+1)),
					alpha=0.5,
					beta=0.1,tao_k=0.01,max.iter=1000,verbose=TRUE,Normalize=TRUE,Norm.method = "frac")

####B-mode batch effect correction model (FACS reference)
res_alg_all_ra_l2_norm_sort_2 <- cell_deconvolve(X=as.matrix(log2(bulk_correct_ra+1)),
                    theta=Fra_ra_correct$theta,
					R=as.matrix(log2(ref_ra+1)),
					alpha=0.5,
					beta=0.1,tao_k=0.01,max.iter=1000,verbose=TRUE,Normalize=TRUE,Norm.method = "frac")
					
####S-mode batch effect correction model (single cell reference)	
Fra_ra <- get_proportion(exprs(bulk_ra_vs_oa),tmp_ra$type)
res_alg_all_ra_l2_norm <- cell_deconvolve(X=as.matrix(log2(exprs(bulk_ra_vs_oa)[rownames(tmp_ra$type),]+1)),
                    theta=Fra_ra$theta,
					R=as.matrix(log2(tmp_ra$type+1)),
					alpha=0.5,
					beta=0.1,tao_k=0.01,max.iter=1000,verbose=TRUE,Normalize=FALSE)


##########trace norm model
####raw model (FACS reference)
res_alg_all_ra_l2_norm_sort_trace <- cell_deconvolve_trace(O=log2(exprs(bulk_ra_vs_oa)[rownames(ref_ra),]+1), 
                                                  theta=Fra_ra2$theta, 
                                                  R=log2(ref_ra+1), 
                                                  alpha=0.5,beta=1,gamma=1,epsilon=0.0001,
                                                  max.iter=1000,solver = "admm",verbose=FALSE,Normalize=FALSE)


####S-mode batch effect correction model (single cell reference)	
Fra_ra <- get_proportion(exprs(bulk_ra_vs_oa),tmp_ra$type)
res_alg_all_ra_l2_norm_trace <- cell_deconvolve_trace(O=log2(exprs(bulk_ra_vs_oa)[rownames(tmp_ra$type),]+1), 
                                                  theta=Fra_ra$theta, 
                                                  R=log2(tmp_ra$type+1), 
                                                  alpha=0.5,beta=1,gamma=1,epsilon=0.0001,
                                                  max.iter=1000,solver = "admm",verbose=TRUE,Normalize=FALSE)



#########trace norm model
##Using bMIND to find CTE			
sig <- as.matrix(log2(tmp_ra$type+1))
colnames(sig) <- paste0("cellType",1:4)
colnames(Fra_ra$theta) <- colnames(sig)
G <- log2(exprs(bulk_ra_vs_oa)[rownames(tmp_ra$type),]+1)
deconv = bMIND(G, frac = Fra_ra$theta, profile=sig,  ncore = 5, noRE = F)

sig <- as.matrix(ref_ra)
colnames(sig) <- paste0("cellType",1:4)
colnames(Fra_ra2$theta) <- colnames(sig)
G <- log2(exprs(bulk_ra_vs_oa)[rownames(ref_ra),]+1)
deconv2 = bMIND(G, frac = Fra_ra2$theta, profile=sig,  ncore = 5, noRE = F)
				

##Running TCA
G <- log2(exprs(bulk_ra_vs_oa)[rownames(tmp_ra$type),]+1)
tca.mdl <- tca(X = G, W = Fra_ra$theta, C1 = NULL, C2 = NULL,
                parallel = TRUE);
Z_hat <- tensor(X = G, tca.mdl)


G <- log2(exprs(bulk_ra_vs_oa)[rownames(ref_ra),]+1)
tca.mdl <- tca(X = G, W = Fra_ra2$theta, C1 = NULL, C2 = NULL,
                parallel = TRUE);
Z_hat2 <- tensor(X = G, tca.mdl)

			
pheno <- c(rep("OA",22),rep("RA",(174-23+1)))
pheno <- as.factor(pheno)

####################################################################################################
###perform concordance analysis

evaluate_concordance <- function(gene, sort_exp,sort_label,pre_exp,pre_label,n_boot=100){
   seq <- runif(100)
   sort_exp <- sort_exp[gene,]
   pre_exp <- pre_exp[gene,]
   pre_exp[pre_exp<0] <- 0
   table <- NULL
   for(i in gene){
       pvalue <- wilcox.test((pre_exp[i,pre_label=="OA"]),(pre_exp[i,pre_label=="RA"]))$p.value
       table <- rbind(table,cbind(pvalue, mean((pre_exp[i,pre_label=="OA"]))/mean((pre_exp[i,pre_label=="RA"]))))
   }
   rownames(table) <- gene
   predict <- na.omit(table)

   table <- NULL
   for(i in gene){
      pvalue <- wilcox.test(sort_exp[i,sort_label=="OA"],sort_exp[i,sort_label=="RA"])$p.value
      table <- rbind(table,cbind(pvalue, mean(sort_exp[i,sort_label=="OA"])/mean(sort_exp[i,sort_label=="RA"])))
   }
   rownames(table) <- gene

   ES <- abs(log2(predict[predict[,1]<0.01,2]))
   seq <- quantile(ES,seq)
   prob <- NULL
   for(i in 1:length(seq)){   
   predict_filter <- predict[predict[,1]<0.01&abs(log2(predict[,2]))>seq[i],]
   predict_filter <- predict_filter[order(predict_filter[,1],decreasing=FALSE),]

   table_sub <- table[rownames(predict_filter),]
   prob <- c(prob,sum(diag(table(sign(log2(predict_filter[,2])),sign(log2(table_sub[,2])))))/sum(table(sign(log2(predict_filter[,2])),sign(log2(table_sub[,2])))))
   }
   
   prob
}

B <- cbind(exprs(OA)[,OA$type == "B cell"],exprs(RA)[,RA$type == "B cell"])
label_sc <- c(rep("OA",sum(OA$type == "B cell")),rep("RA",sum(RA$type == "B cell")))
label <- as.character(pheno)
genes <- intersect(rownames(B),rownames(ref_ra))

###FACS reference
bulk_B <- evaluate_concordance(genes,B,label_sc,log2(exprs(bulk_ra_vs_oa)[rownames(res_alg_all_ra_l2_norm$X_k),]+1),label)
ENIGMA_B <- evaluate_concordance(genes,B,label_sc,res_alg_all_ra_l2_norm_sort_2$X_k_norm[,,1],label)
ENIGMA_B2 <- evaluate_concordance(genes,B,label_sc,res_alg_all_ra_l2_norm_sort$X_k_norm[,,1],label)
ENIGMA2_B <- evaluate_concordance(genes,B,label_sc,res_alg_all_ra_l2_norm_sort_trace$X_k_norm[,,1],label)
bMIND_B <- evaluate_concordance(genes,B,label_sc,deconv2$A[,1,],label)
TCA_B <- evaluate_concordance(genes,B,label_sc,Z_hat2[[1]],label)

fibro <- cbind(exprs(OA)[,OA$type == "Fibroblast"],exprs(RA)[,RA$type == "Fibroblast"])
label_sc <- c(rep("OA",sum(OA$type == "Fibroblast")),rep("RA",sum(RA$type == "Fibroblast")))
label <- as.character(pheno)
genes <- intersect(rownames(fibro),rownames(ref_ra))

bulk_fibro <- evaluate_concordance(genes,fibro,label_sc,log2(exprs(bulk_ra_vs_oa)[rownames(res_alg_all_ra_l2_norm$X_k),]+1),label)
ENIGMA_fibro <- evaluate_concordance(genes,fibro,label_sc,res_alg_all_ra_l2_norm_sort_2$X_k_norm[,,2],label)
ENIGMA_fibro2 <- evaluate_concordance(genes,fibro,label_sc,res_alg_all_ra_l2_norm_sort$X_k_norm[,,2],label)
ENIGMA2_fibro <- evaluate_concordance(genes,fibro,label_sc,res_alg_all_ra_l2_norm_sort_trace$X_k_norm[,,2],label)
bMIND_fibro <- evaluate_concordance(genes,fibro,label_sc,deconv2$A[,2,],label)
TCA_fibro <- evaluate_concordance(genes,fibro,label_sc,Z_hat2[[2]],label)

mono <- cbind(exprs(OA)[,OA$type == "Monocyte"],exprs(RA)[,RA$type == "Monocyte"])
label_sc <- c(rep("OA",sum(OA$type == "Monocyte")),rep("RA",sum(RA$type == "Monocyte")))
label <- as.character(pheno)
genes <- intersect(rownames(mono),rownames(ref_ra))

bulk_mono <- evaluate_concordance(genes,mono,label_sc,log2(exprs(bulk_ra_vs_oa)[rownames(res_alg_all_ra_l2_norm$X_k),]+1),label)
ENIGMA_mono <- evaluate_concordance(genes,mono,label_sc,res_alg_all_ra_l2_norm_sort_2$X_k_norm[,,3],label)
ENIGMA_mono2 <- evaluate_concordance(genes,mono,label_sc,res_alg_all_ra_l2_norm_sort$X_k_norm[,,3],label)
ENIGMA2_mono <- evaluate_concordance(genes,mono,label_sc,res_alg_all_ra_l2_norm_sort_trace$X_k_norm[,,3],label)
bMIND_mono <- evaluate_concordance(genes,mono,label_sc,deconv2$A[,3,],label)
TCA_mono <- evaluate_concordance(genes,mono,label_sc,Z_hat2[[3]],label)


T <- cbind(exprs(OA)[,OA$type == "T cell"],exprs(RA)[,RA$type == "T cell"])
label_sc <- c(rep("OA",sum(OA$type == "T cell")),rep("RA",sum(RA$type == "T cell")))
label <- as.character(pheno)
genes <- intersect(rownames(T),rownames(ref_ra))


bulk_T <- evaluate_concordance(genes,T,label_sc,log2(exprs(bulk_ra_vs_oa)[rownames(res_alg_all_ra_l2_norm$X_k),]+1),label)
ENIGMA_T <- evaluate_concordance(genes,T,label_sc,res_alg_all_ra_l2_norm_sort_2$X_k_norm[,,4],label)
ENIGMA_T2 <- evaluate_concordance(genes,T,label_sc,res_alg_all_ra_l2_norm_sort$X_k_norm[,,4],label)
ENIGMA2_T <- evaluate_concordance(genes,T,label_sc,res_alg_all_ra_l2_norm_sort_trace$X_k_norm[,,4],label)
bMIND_T <- evaluate_concordance(genes,T,label_sc,deconv2$A[,4,],label)
TCA_T <- evaluate_concordance(genes,T,label_sc,Z_hat2[[4]],label)


##################
##ENIGMA perform better than TCA and bMIND regardless if perform frac-based normalization
boxDat <- data.frame(method=c(rep("Bulk",400),rep("TCA",400),rep("bMIND",400),rep("ENIGMA(trace)",400),rep("ENIGMA",400),rep("ENIGMA(B-mode correction)",400)),
                     Concordance=c(c(bulk_B,bulk_fibro,bulk_mono,bulk_T),
					                  c(TCA_B,TCA_fibro,TCA_mono,TCA_T),
									  c(bMIND_B,bMIND_fibro,bMIND_mono,bMIND_T),
									  c(ENIGMA2_B,ENIGMA2_fibro,ENIGMA2_mono,ENIGMA2_T),
									  c(ENIGMA_B2,ENIGMA_fibro2,ENIGMA_mono2,ENIGMA_T2),
									  c(ENIGMA_B,ENIGMA_fibro,ENIGMA_mono,ENIGMA_T))
					 ,
                     CellType = rep(c(rep("B cell",100),rep("Fibroblast",100),rep("Monocyte",100),rep("T cell",100)),6))

boxDat$method <- as.factor(boxDat$method)
boxDat$method <- factor(boxDat$method,levels = c("Bulk","bMIND","TCA","ENIGMA(trace)","ENIGMA","ENIGMA(B-mode correction)"))
boxDat <- boxDat[is.nan(boxDat$Concordance)==FALSE,]		 
p1 <- ggplot(boxDat, aes(x=CellType, y=Concordance,fill = method)) +
  geom_boxplot()+labs(y="Concordance")+
  theme_classic() +
  facet_grid(~CellType, scales = "free_x",as.table = FALSE) +
  theme(strip.background = element_blank(),strip.text.x = element_blank(),text = element_text(size=20))
png("Normalize_concordance.png",res=300,height=1500,width=3000)
p1
dev.off()

########################################################################
##Measuring the concordance through scRNA-seq
B <- cbind(exprs(OA)[,OA$type == "B cell"],exprs(RA)[,RA$type == "B cell"])
label_sc <- c(rep("OA",sum(OA$type == "B cell")),rep("RA",sum(RA$type == "B cell")))
label <- as.character(pheno)
genes <- intersect(rownames(B),rownames(tmp_ra$type))

###FACS reference
bulk_B <- evaluate_concordance(genes,B,label_sc,log2(exprs(bulk_ra_vs_oa)[rownames(res_alg_all_ra_l2_norm$X_k),]+1),label)
ENIGMA_B <- evaluate_concordance(genes,B,label_sc,res_alg_all_ra_l2_norm$X_k[,,1],label)
ENIGMA2_B <- evaluate_concordance(genes,B,label_sc,res_alg_all_ra_l2_norm_trace$X_k[,,1],label)
bMIND_B <- evaluate_concordance(genes,B,label_sc,deconv$A[,1,],label)
TCA_B <- evaluate_concordance(genes,B,label_sc,Z_hat[[1]],label)

fibro <- cbind(exprs(OA)[,OA$type == "Fibroblast"],exprs(RA)[,RA$type == "Fibroblast"])
label_sc <- c(rep("OA",sum(OA$type == "Fibroblast")),rep("RA",sum(RA$type == "Fibroblast")))
label <- as.character(pheno)
genes <- intersect(rownames(fibro),rownames(tmp_ra$type))

bulk_fibro <- evaluate_concordance(genes,fibro,label_sc,log2(exprs(bulk_ra_vs_oa)[rownames(res_alg_all_ra_l2_norm$X_k),]+1),label)
ENIGMA_fibro <- evaluate_concordance(genes,fibro,label_sc,res_alg_all_ra_l2_norm$X_k[,,2],label)
ENIGMA2_fibro <- evaluate_concordance(genes,fibro,label_sc,res_alg_all_ra_l2_norm_trace$X_k[,,2],label)
bMIND_fibro <- evaluate_concordance(genes,fibro,label_sc,deconv$A[,2,],label)
TCA_fibro <- evaluate_concordance(genes,fibro,label_sc,Z_hat[[2]],label)

mono <- cbind(exprs(OA)[,OA$type == "Monocyte"],exprs(RA)[,RA$type == "Monocyte"])
label_sc <- c(rep("OA",sum(OA$type == "Monocyte")),rep("RA",sum(RA$type == "Monocyte")))
label <- as.character(pheno)
genes <- intersect(rownames(mono),rownames(tmp_ra$type))

bulk_mono <- evaluate_concordance(genes,mono,label_sc,log2(exprs(bulk_ra_vs_oa)[rownames(res_alg_all_ra_l2_norm$X_k),]+1),label)
ENIGMA_mono <- evaluate_concordance(genes,mono,label_sc,res_alg_all_ra_l2_norm$X_k[,,3],label)
ENIGMA2_mono <- evaluate_concordance(genes,mono,label_sc,res_alg_all_ra_l2_norm_trace$X_k[,,3],label)
bMIND_mono <- evaluate_concordance(genes,mono,label_sc,deconv$A[,3,],label)
TCA_mono <- evaluate_concordance(genes,mono,label_sc,Z_hat[[3]],label)


T <- cbind(exprs(OA)[,OA$type == "T cell"],exprs(RA)[,RA$type == "T cell"])
label_sc <- c(rep("OA",sum(OA$type == "T cell")),rep("RA",sum(RA$type == "T cell")))
label <- as.character(pheno)
genes <- intersect(rownames(T),rownames(tmp_ra$type))


bulk_T <- evaluate_concordance(genes,T,label_sc,log2(exprs(bulk_ra_vs_oa)[rownames(res_alg_all_ra_l2_norm$X_k),]+1),label)
ENIGMA_T <- evaluate_concordance(genes,T,label_sc,res_alg_all_ra_l2_norm$X_k[,,4],label)
ENIGMA2_T <- evaluate_concordance(genes,T,label_sc,res_alg_all_ra_l2_norm_trace$X_k[,,4],label)
bMIND_T <- evaluate_concordance(genes,T,label_sc,deconv$A[,4,],label)
TCA_T <- evaluate_concordance(genes,T,label_sc,Z_hat[[4]],label)


##################
##ENIGMA perform better than TCA and bMIND regardless if perform frac-based normalization
boxDat2 <- data.frame(method=c(rep("Bulk",400),rep("TCA",400),rep("bMIND",400),rep("ENIGMA(trace)",400),rep("ENIGMA",400)),
                     Concordance=c(c(bulk_B,bulk_fibro,bulk_mono,bulk_T),
					                  c(TCA_B,TCA_fibro,TCA_mono,TCA_T),
									  c(bMIND_B,bMIND_fibro,bMIND_mono,bMIND_T),
									  c(ENIGMA2_B,ENIGMA2_fibro,ENIGMA2_mono,ENIGMA2_T),
									  c(ENIGMA_B,ENIGMA_fibro,ENIGMA_mono,ENIGMA_T))
					 ,
                     CellType = rep(c(rep("B cell",100),rep("Fibroblast",100),rep("Monocyte",100),rep("T cell",100)),5))

boxDat2$method <- as.factor(boxDat2$method)
boxDat2$method <- factor(boxDat2$method,levels = c("Bulk","bMIND","TCA","ENIGMA(trace)","ENIGMA","ENIGMA(B-mode correction)"))
boxDat2 <- boxDat2[is.nan(boxDat2$Concordance)==FALSE,]		 
p1 <- ggplot(boxDat2, aes(x=CellType, y=Concordance,fill = method)) +
  geom_boxplot()+labs(y="Concordance")+
  theme_classic() +
  facet_grid(~CellType, scales = "free_x",as.table = FALSE) +
  theme(strip.background = element_blank(),strip.text.x = element_blank(),text = element_text(size=20)) 
png("scRNA_ref_concordance_main.png",res=300,height=1500,width=3000)
p1
dev.off()

###############################################################################
###Analysis monocyte trajectory
library(scater)
library(scran)
library(igraph)
library(Seurat)

res_alg_all_ra2 <- cell_deconvolve(X=as.matrix(sqrt(exprs(bulk_ra_vs_oa))),
                    theta=Fra_ra3$theta,
					R=as.matrix(sqrt(ref_ra)),
					epsilon = 0.001,
					alpha=0.5,
					beta=0.1,tao_k=0.01,max.iter=1000,verbose=TRUE,pre.process = "sqrt", Norm.method = "frac")

sce_monocyte2 <- SingleCellExperiment(assays = list(logcounts = res_alg_all_ra2$X_k_norm[,,3]))
sce_monocyte2 <- runPCA(sce_monocyte2,scale=TRUE)
pca <- reducedDims(sce_monocyte2)
plot(attributes(pca$PCA)$varExplained)
#####Running the Seurat to perform analysis

sce_monocyte2 = CreateSeuratObject(res_alg_all_ra2$X_k_norm[,,3])
sce_monocyte2[['pca']] = CreateDimReducObject(embeddings = pca$PCA,key="PC_")
sce_monocyte2 = RunUMAP(sce_monocyte2,dims=1:10)

###save the coordinate of UMAP
sce_monocyte2@reductions$umap@cell.embeddings -> umap
sce_monocyte2@reductions$pca@cell.embeddings -> pca

sce_monocyte2 <- SingleCellExperiment(assays = list(logcounts = res_alg_all_ra2$X_k_norm[,,3]))
reducedDim(sce_monocyte2, "PCA") <- pca
reducedDim(sce_monocyte2, "UMAP") <- umap
g <- buildSNNGraph(sce_monocyte2, k = 10, use.dimred = "PCA")
reducedDim(sce_monocyte2, "SNN") <- as_adjacency_matrix(g, attr = "weight")
sce_monocyte2$cell_clusters <- factor(cluster_louvain(g)$membership)
plotReducedDim(sce_monocyte2, dimred = "UMAP", colour_by = "cell_clusters",point_alpha=1)
assay <- normalize.quantiles(res_alg_all_ra2$X_k[,,3])
rownames(assay) <- rownames(res_alg_all_ra2$X_k[,,3])
colnames(assay) <- colnames(res_alg_all_ra2$X_k[,,3])

################################################
##Using the same pipeline to process the bulk
sce_bulk2 <- SingleCellExperiment(assays = list(logcounts = log2(exprs(bulk_ra_vs_oa)+1)))
sce_bulk2$cell_type <- factor(cluster_louvain(g)$membership)
sce_bulk2 <- runPCA(sce_bulk2,scale=TRUE)
pca <- reducedDims(sce_bulk2)

sce_bulk2 = CreateSeuratObject(log2(exprs(bulk_ra_vs_oa)+1))
sce_bulk2[['pca']] = CreateDimReducObject(embeddings = pca$PCA,key="PC_")
sce_bulk2 = RunUMAP(sce_bulk2,dims=1:10)

###save the coordinate of UMAP
sce_bulk2@reductions$umap@cell.embeddings -> umap
sce_bulk2@reductions$pca@cell.embeddings -> pca

sce_bulk2 <- SingleCellExperiment(assays = list(logcounts = res_alg_all_ra2$X_k_norm[,,3]))
reducedDim(sce_bulk2, "PCA") <- pca
reducedDim(sce_bulk2, "UMAP") <- umap
g <- buildSNNGraph(sce_bulk2, k = 10, use.dimred = "PCA")
reducedDim(sce_bulk2, "SNN") <- as_adjacency_matrix(g, attr = "weight")
sce_bulk2$cell_clusters <- factor(cluster_louvain(g)$membership)
sce_bulk2$cell_type <- sce_monocyte2$cell_clusters
plotReducedDim(sce_bulk2, dimred = "UMAP", colour_by = "cell_type",point_alpha=1)
								
################################################################
###measure cell differentiation potency
library(CytoTRACE)
library(SCENT)
### getPPI_String function could be downloaded from https://github.com/WWXkenmo/ICAnet/blob/master/getPPI_String.R
PPI <- getPPI_String(sce_monocyte2,species=9606)
score <- CompCCAT(res_alg_all_ra2$X_k_norm[,,3],PPI[rownames(PPI)[rownames(PPI)%in%ribosome$V1 == FALSE],rownames(PPI)[rownames(PPI)%in%ribosome$V1 == FALSE]])
rs <- CytoTRACE::CytoTRACE(res_alg_all_ra2$X_k[rownames(res_alg_all_ra2$X_k[,,3])[rownames(res_alg_all_ra2$X_k[,,3])%in%ribosome$V1 == FALSE],,3])
ribosome_exp <- colMeans((res_alg_all_ra2$X_k_norm[intersect(ribosome$V1,rownames(res_alg_all_ra2$X_k_norm[,,3])),,3]))
sce_monocyte2$ribosome <- ribosome_exp
sce_monocyte2$cytotrace <- rs$CytoTRACE
sce_monocyte2$scent <- score


PPI <- getPPI_String(sce_bulk2,species=9606)
score <- CompCCAT(log2(exprs(bulk_ra_vs_oa)+1),PPI[rownames(PPI)[rownames(PPI)%in%ribosome$V1 == FALSE],rownames(PPI)[rownames(PPI)%in%ribosome$V1 == FALSE]])
rs <- CytoTRACE::CytoTRACE(log2(exprs(bulk_ra_vs_oa)+1)[rownames(exprs(bulk_ra_vs_oa))[rownames(res_alg_all_ra2$X_k[,,3])%in%ribosome$V1 == FALSE],])
ribosome_exp <- colMeans((log2(exprs(bulk_ra_vs_oa)+1)[intersect(ribosome$V1,rownames(res_alg_all_ra2$X_k_norm[,,3])),]))
sce_bulk2$ribosome <- ribosome_exp
sce_bulk2$cytotrace <- rs$CytoTRACE
sce_bulk2$scent <- score


png("ribosome.png",res=300,height=1500,width=1800)
plotReducedDim(sce_monocyte2, dimred = "UMAP", colour_by = "ribosome",point_alpha=1,text_by = "cell_clusters") + theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))+ labs(title = "Ribosome average expression")
dev.off()

png("CytoTRACE.png",res=300,height=1500,width=1800)
plotReducedDim(sce_monocyte2, dimred = "UMAP", colour_by = "cytotrace",point_alpha=1,text_by = "cell_clusters") + theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))+ labs(title = "CytoTRACE differentiation potency")
dev.off()

png("SCENT.png",res=300,height=1500,width=1800)
plotReducedDim(sce_monocyte2, dimred = "UMAP", colour_by = "scent",point_alpha=1,text_by = "cell_clusters") + theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15)) + labs(title = "SCENT differentiation potency")
dev.off()

png("ribosome(boxplot).png",res=300,height=1500,width=1800)
plotColData(sce_monocyte2,colour_by = "cell_clusters",x = "cell_clusters", y = c("ribosome"),point_alpha=1)+theme_bw()+theme(text = element_text(size=15))+labs( y = "Ribosome average expression",x="Clusters")
dev.off()

png("CytoTRACE(boxplot).png",res=300,height=1500,width=1800)
plotColData(sce_monocyte2,colour_by = "cell_clusters",x = "cell_clusters", y = c("cytotrace"),point_alpha=1)+theme_bw()+theme(text = element_text(size=15))+labs( y = "CytoTRACE differentiation potency",x="Clusters")
dev.off()

png("SCENT(boxplot).png",res=300,height=1500,width=1800)
plotColData(sce_monocyte2,colour_by = "cell_clusters",x = "cell_clusters", y = c("scent"),point_alpha=1)+theme_bw()+theme(text = element_text(size=15))+labs( y = "SCENT differentiation potency",x="Clusters")
dev.off()

################################################
#fit slingshot
set.seed(1)
rootLabel <- as.character(as.matrix(sce_monocyte2$cell_clusters))
### 26th sample have the highest differentiation score in both SCENT prediction and CytoTRACE prediction
rootLabel[26] <- "root"
sce_monocyte2$rootLabel <- rootLabel

sce_pt <- slingshot(sce_monocyte2, clusterLabels = 'cell_clusters', reducedDim = 'UMAP',start.clus = "3")
plot(reducedDims(sce)$UMAP, col = as.factor(sce$cell_clusters), pch=16, asp = 1)
lines(slingshot::SlingshotDataSet(sce), lwd=2, col = 'black')
lineages <- getLineages(reducedDims(sce_monocyte2)$UMAP, sce_monocyte2$rootLabel, start.clus = "root")
curves <- getCurves(lineages)

##the color
gg_plot <- function(sds, col = NULL,colour=NULL, title = NULL, lineSize = 1, reduction = "UMAP"
                    , titleFsize = 20
                    , line.colors = gray.colors(n = length(sds@curves), start = 0, end = .6, alpha = 1 )
                    , ...) {
  rd <- reducedDim(sds)
  
  if (is.null(col)) {
    cl <- slingClusterLabels(sds)
    if ("matrix" %in% is(cl)) {
      cl <- apply(cl, 1, which.max)
      cl <- as.character(cl)
    }
  } else {
    cl <- col
  }
  
  # Getting the main plot
  df <- data.frame(dim1 = rd[, 1], dim2 = rd[, 2], col = cl)
  p <- ggplot(df, aes(x = dim1, y = dim2, col = col)) +
    geom_point(...) +
    theme_bw() + +theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15)) + scale_colour_manual(values = colour)+
    labs(title = title, col = "") +
    theme(plot.title = element_text(size = titleFsize)) # , face = "bold"
  # Adding the curves
  for (i in seq_along(slingCurves(sds))) {
    curve_i <- slingCurves(sds)[[i]]
    curve_i <- curve_i$s[curve_i$ord, ]
    colnames(curve_i) <- c("dim1", "dim2")
    p <- p + geom_path(data = as.data.frame(curve_i), col = line.colors[i], size = 1) + 
      labs(x = paste(reduction, 1), y = paste(reduction, 1)) 
    
  }
  return(p)
}

col = c(rgb(110,155,205,max=255),rgb(255,154,67,max=255),rgb(97,189,86,max=255),rgb(242,147,140,max=255))
gg_plot(curves,col=sce_monocyte2$label,colour=col,line.colors = rep("black","2"))
##############

topgenes <- rownames(different_end_association)[different_end_association$pvalue<0.1]
pst.ord <- order(sce_pt$slingPseudotime_1, na.last = NA)
##reverse the order
pst.ord <- rev(pst.ord)

heatdata <- assays(sce_pt)$logcounts[topgenes, pst.ord]
heatclus <- sce_pt$cell_clusters[pst.ord]
##reverse the order

pst.ord2 <- order(sce_pt$slingPseudotime_2, na.last = NA)
heatdata2 <- assays(sce_pt)$logcounts[topgenes, pst.ord2]
heatclus2 <- sce_pt$cell_clusters[pst.ord2]

##concatenate
heatdata <- cbind(heatdata,heatdata2)
heatclus <- c(heatclus,heatclus2)
pseudotime <- c(sce_pt$slingPseudotime_1[pst.ord],sce_pt$slingPseudotime_2[pst.ord2])

heatdata_scale <- t(apply(heatdata,1,scale))
rownames(heatdata_scale) <- rownames(heatdata)
colnames(heatdata_scale) <- colnames(heatdata)
heatdata_scale[heatdata_scale>1] <- 1
heatdata_scale[heatdata_scale< -1] <- -1
pheatmap::pheatmap(heatdata_scale, Colv = NA,
                             ColSideColors = RColorBrewer::brewer.pal(4,"Set1")[heatclus],cluster_cols = FALSE,show_rownames = FALSE,show_colnames = FALSE)

bar <- col
names(bar) <-  names(table(heatclus))
color <- circlize::colorRamp2(c(-1,0,1), c("grey","white", "red"), transparency = 0.05, space = "LAB")
top_anno <- ComplexHeatmap::HeatmapAnnotation(cluster = heatclus,pseudotime = pseudotime, col = list(cluster = bar),show_legend = TRUE,show_annotation_name = TRUE)
p <- ComplexHeatmap::Heatmap(heatdata_scale,show_column_names = FALSE,show_row_names = FALSE,col=color,cluster_rows = TRUE,cluster_columns = FALSE,top_annotation = top_anno,heatmap_legend_param = list(title = "Expression"))
p <- draw(p)  #Show the heatmap
rcl.list <- row_order(p)
rcl.list <- rownames(heatdata_scale)[rcl.list]
keyGene_pos <- which(rcl.list %in% c("CD44","PTPRC","S100A8","S100A9","CCL3","CXCL12","CAV1","SEMA3C","DAB2"))
keyGene <- rcl.list[keyGene_pos]
row_anno <- rowAnnotation(keypath = anno_mark(at = keyGene_pos,labels = keyGene))
p <- ComplexHeatmap::Heatmap(heatdata_scale,show_column_names = FALSE,show_row_names = FALSE,col=color,cluster_rows = TRUE,cluster_columns = FALSE,top_annotation = top_anno,right_annotation = row_anno, heatmap_legend_param = list(title = "Expression"))

pdf("heatmap.pdf",height=8,width=12)
p
dev.off()
##################################
###GO analysis
upgene <- rownames(different_end_association)[different_end_association$pvalue<0.1 & different_end_association$logFC1_2>0]
upgene <- bitr(upgene,fromType="SYMBOL",toType="ENTREZID",org.Hs.eg.db)
upgene_go <- enrichGO(upgene[,2], ont="BP",pvalueCutoff = 0.01,qvalueCutoff = 0.05,OrgDb = org.Hs.eg.db,universe = bkgene$ENTREZID,readable=TRUE)


downgene <- rownames(different_end_association)[different_end_association$pvalue<0.1 & different_end_association$logFC1_2<0]
downgene <- bitr(downgene,fromType="SYMBOL",toType="ENTREZID",org.Hs.eg.db)
downgene_go <- enrichGO(downgene[,2], ont="BP",pvalueCutoff = 0.01,qvalueCutoff = 0.05,OrgDb = org.Hs.eg.db,universe = bkgene$ENTREZID,readable=TRUE)

###################################
pdf("upgene_go.pdf",height=6,width=10)
dotplot(upgene_go,showCategory = 7)
dev.off()
pdf("downgene_go.pdf",height=6,width=10)
dotplot(downgene_go,showCategory = 7)
dev.off()

##############################################



