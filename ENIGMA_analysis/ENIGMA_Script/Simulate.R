GEP <- read.csv("/mnt/data1/weixu/HiDe/ST_data/GSE150580_Li_Brugge_10XscRNAseq_GeneCellMatrix_RNAcounts_GEO.csv",header=TRUE,row.names=1)
metadata <- read.csv("GSE150580_Li_Brugge_10XscRNAseq_Metadata_GEO.csv",header=TRUE,row.names=1)

metadata <- subset(metadata,m.CellTypes %in% c("Luminal-AV","Fibroblast","Myoepithelial"))
metadata$m.CellTypes <- as.character(as.matrix(metadata$m.CellTypes))
metadata$m.SampleName <- as.character(as.matrix(metadata$m.SampleName))

## Using Young Nulliparous B as the reference, A/C as the data for simulation
## Simulate ST datasets of 900 pixels, combined gene counts of a mixture of different cell-types up to 8 total individual cells.

sim_sc <- GEP[,metadata$m.SampleName %in% c("Young Nulliparous B","Young Nulliparous C")]
ref_sc <- GEP[,metadata$m.SampleName %in% c("Young Nulliparous A")]

## simulate cell type proportions according to coordinate

proLogistic <- function(coord){
   pro <- 1/(1+exp(-(coord-15)))
   pro
}

pro_fibro <- pro_lumi <- pro_myo <- matrix(0,nrow = 30, ncol = 30)
for(i in 1:30){
 for(j in 1:30){
   fibro_frac <- proLogistic(i)
   lumi_frac <- proLogistic(j)
   myo_frac <- proLogistic((i+j)/2)
   
   frac <- c(fibro_frac,lumi_frac,myo_frac)
   frac <- frac / sum(frac)
   fibro_frac <- frac[1]
   lumi_frac <- frac[2]
   myo_frac <- frac[3]
   
   pro_fibro[i,j] <- fibro_frac
   pro_lumi[i,j] <- lumi_frac
   pro_myo[i,j] <- myo_frac
   }
}

count_fibro <- floor(pro_fibro*8)
count_lumi <- floor(pro_lumi*8)
count_myo <- 8 - (count_fibro+count_lumi)



#########################################################
### random sampling and aggregate as the expression profile
metadata$m.CellTypes[metadata$m.CellTypes %in% c("Luminal-AV")] <- "Luminal cell"
sim_sc_metadata <- metadata[colnames(sim_sc),]
ref_sc_metadata <- metadata[colnames(ref_sc),]

fibro_gep <- NULL
lumi_gep <- NULL
myo_gep <- NULL
sample_id <- NULL
mixture <- NULL
for(i in 1:30){
  for(j in 1:30){
     num_fibro <- count_fibro[i,j]
	 num_lumi <- count_lumi[i,j]
	 num_myo <- count_myo[i,j]
	 
	 fibro_id <- rownames(sim_sc_metadata)[which(sim_sc_metadata$m.CellTypes == "Fibroblast")[sample(1:length(which(sim_sc_metadata$m.CellTypes == "Fibroblast")),num_fibro,replace=FALSE)]]
	 lumi_id <- rownames(sim_sc_metadata)[which(sim_sc_metadata$m.CellTypes == "Luminal cell")[sample(1:length(which(sim_sc_metadata$m.CellTypes == "Luminal cell")),num_lumi,replace=FALSE)]]
	 myo_id <- rownames(sim_sc_metadata)[which(sim_sc_metadata$m.CellTypes == "Myoepithelial")[sample(1:length(which(sim_sc_metadata$m.CellTypes == "Myoepithelial")),num_myo,replace=FALSE)]]
	 
	 fibro_gep <- cbind(fibro_gep, rowMeans(as.matrix(sim_sc[,fibro_id])))
	 lumi_gep <- cbind(lumi_gep, rowMeans(as.matrix(sim_sc[,lumi_id])))
	 myo_gep <- cbind(myo_gep, rowMeans(as.matrix(sim_sc[,myo_id])))
	 mixture <- cbind(mixture, rowSums(as.matrix(sim_sc[,fibro_id])) + rowSums(as.matrix(sim_sc[,lumi_id])) + rowSums(as.matrix(sim_sc[,myo_id])))
	 sample_id <- c(sample_id,paste("pixel_",i,"th_row_",j,"th_column",sep=""))
	}
}
fibro_gep[is.nan(fibro_gep)] <- 0
lumi_gep[is.nan(lumi_gep)] <- 0
myo_gep[is.nan(myo_gep)] <- 0

colnames(mixture) <- colnames(fibro_gep) <- colnames(lumi_gep) <- colnames(myo_gep) <- sample_id

## filtering genes
g <- which(rowSums(mixture) >= 10 & apply(mixture,1,function(x){sum(x!=0)}) >= 10)
mixture <- mixture[g,]
fibro_gep <- fibro_gep[g,]
lumi_gep <- lumi_gep[g,]
myo_gep <- myo_gep[g,]

#############
geneID <- rownames(mixture)
mixture[,colSums(mixture)!=0] <- mixture[,colSums(mixture)!=0] %*% diag(10^6/colSums(mixture[,colSums(mixture)!=0]))
fibro_gep[,colSums(fibro_gep)!=0] <- fibro_gep[,colSums(fibro_gep)!=0] %*% diag(10^6/colSums(fibro_gep[,colSums(fibro_gep)!=0]))
lumi_gep[,colSums(lumi_gep)!=0] <- lumi_gep[,colSums(lumi_gep)!=0] %*% diag(10^6/colSums(lumi_gep[,colSums(lumi_gep)!=0]))
myo_gep[,colSums(myo_gep)!=0] <- myo_gep[,colSums(myo_gep)!=0] %*% diag(10^6/colSums(myo_gep[,colSums(myo_gep)!=0]))

############
dataset <- list()
dataset$mixture <- mixture
dataset$fibro_gep <- fibro_gep
dataset$lumi_gep <- lumi_gep
dataset$myo_gep <- myo_gep
dataset$pro_fibro <- pro_fibro
dataset$pro_lumi <- pro_lumi
dataset$pro_myo <- pro_myo
dataset$ref_sc <- ref_sc
dataset$ref_sc_metadata <- ref_sc_metadata