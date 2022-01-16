##########################
#load the required packages
source("/Path/to/Data/ENIGMA.R")
library(MIND)
library(Biobase)
library(Seurat)
library(MASS)

################Using HNSCC to reconsititute profile
##The HNSCC single cell RNA-seq data could downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322
HNSCC <- read.table("/Path/to/Data/GSE103322_HNSCC_all_data.txt",header=TRUE,row.names=1,sep="\t")

###Processed scRNA-seq data
celltype <- as.character(as.matrix(HNSCC[5,]))
celltype[celltype=="0"] <- "Cancer Cell"
celltype[celltype=="-Fibroblast"] <- "Fibroblast"
batch <- as.character(as.matrix(HNSCC[1,]))

HNSCC <- HNSCC[-c(1:5),]
HNSCC <- CreateSeuratObject(counts = HNSCC, project = "HNSCC", min.cells=0.01*ncol(HNSCC))
HNSCC <- GetAssayData(HNSCC)
HNSCC.m <- 2^(HNSCC) - 1
HNSCC.m <- as.matrix(HNSCC.m)
###
patient <- colnames(HNSCC.m)
patient.list <- strsplit(patient, split="[_]")
for(i in 1:length(patient.list)) patient[i] <- patient.list[[i]][1]

###filter the samples which has less 30 cells
p.id <- names(table(patient)[table(patient)>30])
HNSCC.m <- HNSCC.m[,patient %in% p.id]
celltype <- celltype[patient %in% p.id]
patient <- patient[patient %in% p.id]

###only preserve 6 major cell types
HNSCC.m <- HNSCC.m[,celltype %in% c("B cell","Cancer Cell","Endothelial","Fibroblast","Mast","T cell")]
patient <- patient[celltype %in% c("B cell","Cancer Cell","Endothelial","Fibroblast","Mast","T cell")]
celltype <- celltype[celltype %in% c("B cell","Cancer Cell","Endothelial","Fibroblast","Mast","T cell")]
names(patient) <- names(celltype) <- colnames(HNSCC.m)

###randomly select parts of patient sample to build reference
singleCell_set <- HNSCC.m[,patient %in% c("HNSCC","HNSCC16") == TRUE]
bulkSet <- HNSCC.m[,patient %in% c("HNSCC","HNSCC16") == FALSE]
patient.sub <- patient[patient %in% c("HNSCC","HNSCC16") == FALSE]
celltype.sub <- celltype[patient %in% c("HNSCC","HNSCC16") == FALSE]

###using the remain samples to build pesudo-bulk samples
bulkSet_HNSCC <- NULL
HNSCC_Fra <- matrix(0, nrow=length(table(patient.sub)),ncol=length(table(celltype.sub)))
colnames(HNSCC_Fra) <- names(table(celltype.sub))
rownames(HNSCC_Fra) <- names(table(patient.sub))
for(i in names(table(patient.sub))){
  sample <- names(patient.sub)[patient.sub == i]
  sample <- sample(sample,500,replace=TRUE) ### fixed sampling 500 cells and perform aggregation, so that each sample have comparable number of moleculars
  bulkSet_HNSCC <- cbind(bulkSet_HNSCC, rowSums(bulkSet[,sample]))
  HNSCC_Fra[i,names(table(celltype.sub[patient.sub==i]))] <- (table(celltype.sub[patient.sub==i]))/sum(patient.sub==i)
}
colnames(bulkSet_HNSCC) <- rownames(HNSCC_Fra)


###Construct ground truth cell type specific expression profile
celltype.sc <- celltype[patient %in% c("HNSCC","HNSCC16") == TRUE]
P_HNSCC <- array(0, dim=c(nrow(bulkSet),ncol(bulkSet_HNSCC),ncol(HNSCC_Fra)),
                 dimnames =list(rownames(bulkSet),colnames(bulkSet_HNSCC),colnames(HNSCC_Fra)))
names(celltype.sub) <- names(patient.sub) <- colnames(bulkSet)
for(i in names(table(patient.sub))){
  sample.subset <- colnames(bulkSet)[patient.sub %in% i]
  for(j in colnames(HNSCC_Fra)){
  sample.subset.cell <- sample.subset[celltype.sub[sample.subset] %in% j]
  if(length(sample.subset.cell)!=0){
  exp <- rowMeans(as.matrix(bulkSet[,sample.subset.cell]))
  P_HNSCC[,i,j] <- exp
  }
}
}

###build reference profile (profile matrix)
profile_hnscc <- NULL
for(i in names(table(celltype.sc))){
  sig <- rowMeans(as.matrix(singleCell_set[,celltype.sc == i]))
  profile_hnscc <- cbind(profile_hnscc,sig)
}
colnames(profile_hnscc) <- names(table(celltype.sc))

###calculate cell type fractions matrix
fra_HNSCC <- get_proportion(bulkSet_HNSCC,profile_hnscc)
###Running ENIGMA
res_alg_all_HNSCC <- cell_deconvolve(X=log2(bulkSet_HNSCC+1),
                    theta=fra_HNSCC$theta,
					R=log2(profile_hnscc+1),
					epsilon=0.001,
					alpha=0.8,
					beta=0.1,tao_k=0.5,verbose=TRUE,Normalize=FALSE)

###Running ENIGMA (trace norm)					
res_alg_all_HNSCC2 <- cell_deconvolve_trace(O = as.matrix(log2(bulkSet_HNSCC+1)),
                                                  theta=fra_HNSCC$theta,
                                                  R=log2(profile_hnscc[rownames(bulkSet_HNSCC),]+1),
                                                  epsilon=0.001,
                                                  alpha=0.8,beta=1,solver="admm",gamma=1,
                                                  verbose=TRUE,max.iter = 500,Normalize=FALSE)

###Running TCA
tca.mdl <- tca(X = as.matrix(log2(bulkSet_HNSCC+1)), W = fra_HNSCC$theta, C1 = NULL, C2 = NULL,
                parallel = TRUE,num_cores=4)
TCA_HNSCC <- tensor(X = as.matrix(log2(bulkSet_HNSCC+1)), tca.mdl)		

###############################################################################
###Running bMIND
###rename cell names
cellName = colnames(fra_HNSCC$theta)
cellName = gsub(' ', '', cellName)
theta = fra_HNSCC$theta
colnames(theta) = cellName
deconv_HNSCC = bMIND2(log2(bulkSet_HNSCC+1),ncore = 3, frac = theta, profile = log2(profile_hnscc[rownames(bulkSet_HNSCC),]+1), noRE = FALSE)

############################################
##perform evaluation
dimnames(deconv_HNSCC$A)[[2]] <- colnames(profile_hnscc)
HiDe <- NULL
HiDe2 <- NULL
bMIND <- NULL
TCA <- NULL
celltype <- NULL
bulk_exp <- NULL
reference <- NULL
for(ct in colnames(deconv_HNSCC$A[,,1])){
CSE_pre2 <- deconv_HNSCC$A[,ct,]
CSE <- P_HNSCC[rownames(profile_hnscc),,ct]
CSE_pre <- res_alg_all_HNSCC$X_k[,,ct]
CSE_pre_trace <- res_alg_all_HNSCC2$X_k[rownames(CSE_pre2),,ct]
CSE_pre3 <- TCA_HNSCC[[which(colnames(profile_hnscc) == ct)]]
cor_pre <- NULL
cor_pre2 <- NULL
cor_pre3 <- NULL
cor_pre4 <- NULL
cor_bulk <- NULL
cor_ref <- NULL
for(i in which(fra_HNSCC$theta[,ct]!=0)){
 cor_pre <- c(cor_pre, cor(log2(CSE[,i]+1),CSE_pre[,i],method="sp"))
 cor_pre2 <- c(cor_pre2, cor(log2(CSE[,i]+1),CSE_pre2[,i],method="sp"))
 cor_pre3 <- c(cor_pre3, cor(log2(CSE[,i]+1),CSE_pre3[,i],method="sp"))
 cor_pre4 <- c(cor_pre4, cor(log2(CSE[,i]+1),CSE_pre_trace[,i],method="sp"))
 cor_bulk <- c(cor_bulk, cor(log2(CSE[,i]+1),as.matrix(bulkSet_HNSCC)[rownames(CSE),i],method="sp"))
 cor_ref <- c(cor_ref, cor(log2(CSE[,i]+1),profile_hnscc[rownames(CSE),ct],method="sp"))
}
HiDe <- c(HiDe,cor_pre)
bMIND <- c(bMIND,cor_pre2)
TCA <- c(TCA,cor_pre3)
HiDe2 <- c(HiDe2,cor_pre4)
bulk_exp <- c(bulk_exp, cor_bulk)
reference <- c(reference, cor_ref)
celltype <- c(celltype,rep(ct,length(cor_pre)))
}

dat <- data.frame(method=c(rep("ENIGMA",length(HiDe)),
                           rep("ENIGMA(trace)",length(HiDe2)),
                           rep("bMIND",length(bMIND)),
						   rep("TCA",length(TCA)),
						   rep("bulk",length(bulk_exp)),
						   rep("reference",length(reference))),
				  performance=c(HiDe,
				                HiDe2,
				                bMIND,
								TCA,
								bulk_exp,
								reference),
				  celltype=c(rep(celltype,6)))
dat_sample <- dat
dat$method <- factor(dat$method,levels=c("reference","bulk","TCA","bMIND","ENIGMA(trace)","ENIGMA"))
tiff("performance_compare_spearman_HNSCC(sample).tiff",res=300,height=1000,width=2500)
p<-ggplot(dat, aes(x=celltype, y=performance, color=method)) +
    geom_boxplot()+theme_minimal()+labs(y="Correlation per sample")+facet_grid(~celltype, scales = "free_x") + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
    mytheme + theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))+NoLegend()+ylim(c(0,1))
p
dev.off()


HiDe <- NULL
HiDe2 <- NULL
bMIND <- NULL
TCA <- NULL
celltype <- NULL
bulk_exp <- NULL
for(ct in colnames(deconv_HNSCC$A[,,1])){
CSE_pre2 <- deconv_HNSCC$A[rownames(profile_hnscc),ct,]
CSE <- P_HNSCC[rownames(profile_hnscc),,ct]
CSE_pre <- res_alg_all_HNSCC$X_k[rownames(profile_hnscc),,ct]
CSE_pre_trace <- res_alg_all_HNSCC2$X_k[rownames(CSE_pre2),,ct]
CSE_pre3 <- TCA_HNSCC[[which(colnames(profile_hnscc) == ct)]]
cor_pre <- NULL
cor_pre2 <- NULL
cor_pre3 <- NULL
cor_pre4 <- NULL
cor_bulk <- NULL
for(i in rownames(profile_hnscc)){
 cor_pre <- c(cor_pre, cor(log2(CSE[i,which(fra_HNSCC$theta[,ct]!=0)]+1),CSE_pre[i,which(fra_HNSCC$theta[,ct]!=0)],method="sp"))
 cor_pre2 <- c(cor_pre2, cor(log2(CSE[i,which(fra_HNSCC$theta[,ct]!=0)]+1),CSE_pre2[i,which(fra_HNSCC$theta[,ct]!=0)],method="sp"))
 cor_pre3 <- c(cor_pre3, cor(log2(CSE[i,which(fra_HNSCC$theta[,ct]!=0)]+1),CSE_pre3[i,which(fra_HNSCC$theta[,ct]!=0)],method="sp"))
 cor_pre4 <- c(cor_pre4, cor(log2(CSE[i,which(fra_HNSCC$theta[,ct]!=0)]+1),CSE_pre_trace[i,which(fra_HNSCC$theta[,ct]!=0)],method="sp"))
 cor_bulk <- c(cor_bulk, cor(log2(CSE[i,which(fra_HNSCC$theta[,ct]!=0)]+1),log2(as.matrix(bulkSet_HNSCC)[i,which(fra_HNSCC$theta[,ct]!=0)]+1),method="sp"))
}
HiDe <- c(HiDe,cor_pre)
bMIND <- c(bMIND,cor_pre2)
TCA <- c(TCA,cor_pre3)
HiDe2 <- c(HiDe2,cor_pre4)
bulk_exp <- c(bulk_exp, cor_bulk)
celltype <- c(celltype,rep(ct,length(cor_pre)))
}


dat <- data.frame(method=c(rep("ENIGMA",length(HiDe)),
                           rep("ENIGMA(trace)",length(HiDe2)),
                           rep("bMIND",length(bMIND)),
						   rep("TCA",length(TCA)),
						   rep("Bulk",length(bulk_exp))),
				  performance=c(HiDe,
				                HiDe2,
				                bMIND,
								TCA,
								bulk_exp),
				  celltype=c(rep(celltype,5)))
dat$method <- factor(dat$method,levels=c("Bulk","TCA","bMIND","ENIGMA(trace)","ENIGMA"))
dat_genes <- dat
tiff("performance_compare_spearman_HNSCC(gene).tiff",res=300,height=1000,width=2500)
p<-ggplot(dat, aes(x=celltype, y=performance, color=method)) +
    geom_boxplot()+theme_minimal()+labs(y="Correlation per sample")+facet_grid(~celltype, scales = "free_x") + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
    mytheme + theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))+NoLegend()
p
dev.off()
############################################################
#Simulation code of melanoma
######################################################################
##Running the melanoma
# The melanoma datasets could be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056
melanoma <- read.table("/Path/to/Data/GSE72056_melanoma_single_cell_revised_v2.txt",sep="\t",header=TRUE)
sample <- read.csv("/Path/to/Data/sample.csv",row.names=1,stringsAsFactors=FALSE)
patients.id <- sample[,1]
patients.id <- strsplit(patients.id,split="[_]")
patients <- rep("empty",length(patients.id))


for(i in 1:length(patients.id)){
  if(length(patients.id[[i]])>1){
    patients[i] <- patients.id[[i]][1]  
  }
}
patients <- gsub("cy","Cy",patients)
patients <- gsub("CY","Cy",patients)

###replace monika,SS2
patients[patients=="monika"] <- "Cy75"
patients[patients=="SS2"] <- "Cy75"

patients[patients=="Cy88CD45"] <- "Cy88"
patients[patients=="Cy88CD45POS"] <- "Cy88"
patients[patients=="Cy89A"] <- "Cy89"
patients[patients=="Cy89CORE11"] <- "Cy89"
patients[patients=="Cy89COREQ1"] <- "Cy89"
patients[patients=="Cy89FNA"] <- "Cy89"
patients[patients=="Cy89FNAQ2"] <- "Cy89"
patients[patients=="Cy89NEG"] <- "Cy89"
patients[patients=="Cy94CD45POS"] <- "Cy94"

###
metadata <- t(melanoma[c(2,3),-1])
metadata <- metadata[metadata[,1]!=0,]
metadata.sub <- metadata[metadata[,2]==0,]
metadata.sub <- metadata.sub[metadata.sub[,1]==2,]
#######
metadata.sub2 <- metadata[metadata[,2]!=0,]
metadata.sub2 <- metadata.sub2[metadata.sub2[,1]!=2,]
metadata <- rbind(metadata.sub,metadata.sub2)
celltype <- c(rep("Tumor",nrow(metadata.sub)),metadata.sub2[,2])
celltype[celltype=="1"] <- "T cell"
celltype[celltype=="2"] <- "B cell"
celltype[celltype=="3"] <- "Macro. cell"
celltype[celltype=="4"] <- "Endo. cell"
celltype[celltype=="5"] <- "CAF cell"
celltype[celltype=="6"] <- "NK cell"
#####
names(patients) <- colnames(melanoma)[-1]
patients <- patients[c(rownames(metadata.sub),rownames(metadata.sub2))]
melanoma <- melanoma[!duplicated(melanoma[,1]),]
label <- melanoma[,1]
melanoma.sub <- melanoma[-c(1:3),c(rownames(metadata.sub),rownames(metadata.sub2))]
rownames(melanoma.sub) <- label[-c(1:3)]

#####
#melanoma.sub <- 2^(melanoma.sub)-1
singleCell_set <- melanoma.sub[,patients %in% c("Cy94","Cy80") == TRUE]
bulkSet <- melanoma.sub[,patients %in% c("Cy94","Cy80") == FALSE]
patients.sub <- patients[patients %in% c("Cy94","Cy80") == FALSE]
celltype.sub <- celltype[patients %in% c("Cy94","Cy80") == FALSE]
names(celltype.sub) <- names(patients.sub)

###names
names(patients.sub) <- names(celltype.sub) <- colnames(bulkSet)

###generate bulk sample
bulkSet_melanoma <- NULL
melanoma_Fra <- matrix(0, nrow=length(table(patients.sub)),ncol=length(table(celltype.sub)))
colnames(melanoma_Fra) <- names(table(celltype.sub))
rownames(melanoma_Fra) <- names(table(patients.sub))
for(i in names(table(patients.sub))){
  sample <- names(patients.sub)[patients.sub == i] ## extract the samples belong to the specific patients
  sample <- sample(sample,500,replace=TRUE) ## fix sample 500 cells so that all bulk sample have similar "sequencing depth"
  bulkSet_melanoma <- cbind(bulkSet_melanoma, rowSums(bulkSet[,sample]))
  melanoma_Fra[i,names(table(celltype.sub[patients.sub==i]))] <- (table(celltype.sub[patients.sub==i]))/sum(patients.sub==i)
}
colnames(bulkSet_melanoma) <- rownames(melanoma_Fra)

P_melanoma <- array(0, dim=c(nrow(bulkSet),ncol(bulkSet_melanoma),ncol(melanoma_Fra)),
                 dimnames =list(rownames(bulkSet),colnames(bulkSet_melanoma),colnames(melanoma_Fra)))
for(i in colnames(bulkSet_melanoma)){
  sample.subset <- names(patients.sub)[patients.sub %in% i]
  for(j in colnames(melanoma_Fra)){
  sample.subset.cell <- sample.subset[celltype.sub[sample.subset] %in% j]
  if(length(sample.subset.cell)!=0){
  exp <- rowMeans(as.matrix(bulkSet[,sample.subset.cell]))
  P_melanoma[,i,j] <- exp
  }
}
}

celltype.sc <- celltype[patients %in% c("Cy94","Cy80") == TRUE]
profile_melanoma <- NULL
for(i in names(table(celltype.sc))){
  sig <- rowMeans(as.matrix(singleCell_set[,celltype.sc == i]))
  profile_melanoma <- cbind(profile_melanoma,sig)
}
colnames(profile_melanoma) <- names(table(celltype.sc))
###filter out endo cells
#profile_melanoma <- profile_melanoma[,-3]
fra_melanoma <- get_proportion(bulkSet_melanoma,profile_melanoma)

##filter bulk
genes <- rownames(bulkSet_melanoma)[apply(bulkSet_melanoma,1,var)>10^-8]
profileMat <- profile_melanoma[genes,]
bulk <- log2(bulkSet_melanoma[genes,]+1)

####################################################
###Running bMIND
###rename cell names
cellName = colnames(fra_melanoma$theta)
cellName = gsub(' ', '', cellName)
theta = fra_melanoma$theta
colnames(theta) = cellName
deconv_melanoma = bMIND2(bulk,ncore = 3, frac = theta, profile = profileMat, noRE = FALSE)

########################################
tca.mdl <- tca(X = bulk, W = fra_melanoma$theta, C1 = NULL, C2 = NULL,
                parallel = TRUE,num_cores=3)
TCA_melanoma <- tensor(X = (as.matrix(bulk)), tca.mdl)
############################################
###We only look the gene level correlation to the samples which estimation sample fractions are not zero
res_alg_all_melanoma <- cell_deconvolve(X=bulk,
                    theta=fra_melanoma$theta,
					R=profileMat,
					epsilon=0.001,
					alpha=0.8,
					beta=0.1,tao_k=0.5,verbose=TRUE,Normalize=FALSE)
					
res_alg_all_melanoma2 <- cell_deconvolve_trace(O = as.matrix(bulk),
                                                  theta=fra_melanoma$theta,
                                                  R=profileMat,
                                                  epsilon=0.001,
                                                  alpha=0.8,beta=1,method="admm",gamma=1,
                                                  verbose=TRUE,max.iter = 100)
dimnames(deconv_melanoma$A)[[2]] <- colnames(profileMat)
HiDe <- NULL
HiDe2 <- NULL
bMIND <- NULL
TCA <- NULL
celltype <- NULL
bulk_exp <- NULL
reference <- NULL
for(ct in colnames(profileMat)){
CSE_pre2 <- deconv_melanoma$A[,ct,]
CSE <- P_melanoma[rownames(profileMat),,ct]
CSE_pre <- res_alg_all_melanoma$X_k[rownames(profileMat),,ct]
CSE_pre_trace <- res_alg_all_melanoma2[rownames(profileMat),,ct]
CSE_pre3 <- TCA_melanoma[[which(colnames(profileMat) == ct)]][rownames(profileMat),]
cor_pre <- NULL
cor_pre2 <- NULL
cor_pre3 <- NULL
cor_pre4 <- NULL
cor_bulk <- NULL
cor_ref <- NULL
for(i in which(fra_melanoma$theta[,ct]!=0)){
 cor_pre <- c(cor_pre, cor(log2(CSE[,i]+1),CSE_pre[,i],method="sp"))
 cor_pre2 <- c(cor_pre2, cor(log2(CSE[,i]+1),CSE_pre2[,i],method="sp"))
 cor_pre3 <- c(cor_pre3, cor(log2(CSE[,i]+1),CSE_pre3[,i],method="sp"))
 cor_pre4 <- c(cor_pre4, cor(log2(CSE[,i]+1),CSE_pre_trace[,i],method="sp"))
 cor_bulk <- c(cor_bulk, cor(log2(CSE[,i]+1),as.matrix(bulk)[rownames(CSE),i],method="sp"))
 cor_ref <- c(cor_ref, cor(log2(CSE[,i]+1),profileMat[rownames(CSE),ct],method="sp"))
}
HiDe <- c(HiDe,cor_pre)
bMIND <- c(bMIND,cor_pre2)
TCA <- c(TCA,cor_pre3)
HiDe2 <- c(HiDe2,cor_pre4)
bulk_exp <- c(bulk_exp, cor_bulk)
reference <- c(reference, cor_ref)
celltype <- c(celltype,rep(ct,length(cor_pre)))
}

dat <- data.frame(method=c(rep("ENIGMA",length(HiDe)),
                           rep("ENIGMA(trace)",length(HiDe2)),
                           rep("bMIND",length(bMIND)),
						   rep("TCA",length(TCA)),
						   rep("bulk",length(bulk_exp)),
						   rep("reference",length(reference))),
				  performance=c(HiDe,
				                HiDe2,
				                bMIND,
								TCA,
								bulk_exp,
								reference),
				  celltype=c(rep(celltype,6)))
dat$method <- factor(dat$method,levels=c("reference","bulk","TCA","bMIND","ENIGMA(trace)","ENIGMA"))
dat_sample <- dat
png("performance_compare_spearman_melanoma_sample.png",res=300,height=1000,width=2500)
p<-ggplot(dat, aes(x=celltype, y=performance, color=method)) +
    geom_boxplot()+theme_minimal()+labs(y="Correlation per sample")+facet_grid(~celltype, scales = "free_x") + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
    mytheme + theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))+ylim(c(0,1))
p
dev.off()


save(bulkSet_melanoma,profile_melanoma,P_melanoma,file="/Path/to/Data/melanoma_simulate_benchmark.Rdata")
###############################################################
###We only look the gene level correlation to the samples which estimation sample fractions are not zero
HiDe <- NULL
HiDe2 <- NULL
bMIND <- NULL
TCA <- NULL
celltype <- NULL
bulk_exp <- NULL
for(ct in colnames(profileMat)){
CSE_pre2 <- deconv_melanoma$A[,ct,]
CSE <- P_melanoma[rownames(profileMat),,ct]
CSE_pre <- res_alg_all_melanoma$X_k[rownames(profileMat),,ct]
CSE_pre_trace <- res_alg_all_melanoma2[rownames(profileMat),,ct]
CSE_pre3 <- TCA_melanoma[[which(colnames(profileMat) == ct)]][rownames(profileMat),]
cor_pre <- NULL
cor_pre2 <- NULL
cor_pre3 <- NULL
cor_pre4 <- NULL
cor_bulk <- NULL
for(i in rownames(profileMat)){
 cor_pre <- c(cor_pre, cor(log2(CSE[i,which(fra_melanoma$theta[,ct]!=0)]+1),CSE_pre[i,which(fra_melanoma$theta[,ct]!=0)],method="sp"))
 cor_pre2 <- c(cor_pre2, cor(log2(CSE[i,which(fra_melanoma$theta[,ct]!=0)]+1),CSE_pre2[i,which(fra_melanoma$theta[,ct]!=0)],method="sp"))
 cor_pre3 <- c(cor_pre3, cor(log2(CSE[i,which(fra_melanoma$theta[,ct]!=0)]+1),CSE_pre3[i,which(fra_melanoma$theta[,ct]!=0)],method="sp"))
 cor_pre4 <- c(cor_pre4, cor(log2(CSE[i,which(fra_melanoma$theta[,ct]!=0)]+1),CSE_pre_trace[i,which(fra_melanoma$theta[,ct]!=0)],method="sp"))
 cor_bulk <- c(cor_bulk, cor(log2(CSE[i,which(fra_melanoma$theta[,ct]!=0)]+1),as.matrix(bulk)[i,which(fra_melanoma$theta[,ct]!=0)],method="sp"))
}
HiDe <- c(HiDe,cor_pre)
bMIND <- c(bMIND,cor_pre2)
TCA <- c(TCA,cor_pre3)
HiDe2 <- c(HiDe2,cor_pre4)
bulk_exp <- c(bulk_exp, cor_bulk)
celltype <- c(celltype,rep(ct,length(cor_pre)))
}


dat <- data.frame(method=c(rep("ENIGMA",length(HiDe)),
                           rep("ENIGMA(trace)",length(HiDe2)),
                           rep("bMIND",length(bMIND)),
						   rep("TCA",length(TCA)),
						   rep("Bulk",length(bulk_exp))),
				  performance=c(HiDe,
				                HiDe2,
				                bMIND,
								TCA,
								bulk_exp),
				  celltype=c(rep(celltype,5)))
dat$method <- factor(dat$method,levels=c("Bulk","TCA","bMIND","ENIGMA(trace)","ENIGMA"))
dat_gene <- dat
png("performance_compare_spearman_melanoma(gene).png",res=300,height=1000,width=2000)
p<-ggplot(dat, aes(x=celltype, y=performance, color=method)) +
    geom_boxplot()+theme_minimal()+labs(y="Correlation per gene")
p
dev.off()








































