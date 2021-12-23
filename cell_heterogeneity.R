library(splatter)
library(scatter)
library(cowplot)
library(TCA)
library(ENIGMA)
library(TED)

#Simulate profile with similar library size
params <- newSplatParams(seed=1002)
params <- setParams(params, update = list(nGenes = 8368, batchCells=200))
sim.groups.c2 <- splatSimulateGroups(params,
                                     group.prob = c(0.5,0.5),
                                     de.prob = c(0.01),
                                     verbose = FALSE)
sim.groups.c2 <- logNormCounts(sim.groups.c2)
sim.groups.c2 <- runPCA(sim.groups.c2)
sim.groups.c2 <- runTSNE(sim.groups.c2,dimred="PCA",n_dimred=5)
p_c2 <- plotTSNE(sim.groups.c2, colour_by = "Group",point_size=3)	
label.c2 <- sim.groups.c2$Group

params <- newSplatParams(seed=1004)
params <- setParams(params, update = list(nGenes = 8368, batchCells=200))
sim.groups.c3 <- splatSimulateGroups(params,
                                     group.prob = c(0.5,0.5),
                                     de.prob = c(0.1),
                                     verbose = FALSE)
sim.groups.c3 <- logNormCounts(sim.groups.c3)
sim.groups.c3 <- runPCA(sim.groups.c3)
sim.groups.c3 <- runTSNE(sim.groups.c3,dimred="PCA",n_dimred=5)
p_c3 <- plotTSNE(sim.groups.c3, colour_by = "Group",point_size=3)	
label.c3 <- sim.groups.c3$Group

params <- newSplatParams(seed=2022)
params <- setParams(params, update = list(nGenes = 8368, batchCells=200))
sim.groups.c4 <- splatSimulateGroups(params,
                                     group.prob = c(0.25,0.25,0.25,0.25),
                                     de.prob = c(0.3),
                                     verbose = FALSE)
sim.groups.c4 <- logNormCounts(sim.groups.c4)
sim.groups.c4 <- runPCA(sim.groups.c4)
sim.groups.c4 <- runTSNE(sim.groups.c4,dimred="PCA",n_dimred=5)
p_c4 <- plotTSNE(sim.groups.c4, colour_by = "Group",point_size=3)	
Senes_exp <- counts(sim.groups.c4)
label.c4 <- sim.groups.c4$Group

Senes_exp <- Senes_exp %*% diag(10000/colSums(Senes_exp))
sim.groups.c2 <- counts(sim.groups.c2) %*% diag(10000/colSums(counts(sim.groups.c2)))
sim.groups.c3 <- counts(sim.groups.c3) %*% diag(10000/colSums(counts(sim.groups.c3)))
	
rownames(sim.groups.c2) <- rownames(sim.groups.c3) <- rownames(Senes_exp)
colnames(sim.groups.c2) <- colnames(sim.groups.c3) <- colnames(Senes_exp) <- paste0("Sample-",1:ncol(Senes_exp))

mat <- list(Senes_exp,sim.groups.c2,sim.groups.c3)
names(mat) <- c("ct1","ct2","ct3")
H1_array <- array(0,
                  dim = c( 3,
                           8368,
                           200))
H1 <- array(0,
                  dim = c( 3,
                           8368))
for(i in 1:3){
  H1_array[i,,] <- mat[[i]]
  H1[i,] <- rowMeans(mat[[i]])
}

k <- 3 # number of cell types
ng <- 8368 # number of genes
p <- 200 # number of samples

cc <- matrix(runif(p*k), ncol=k)
cc <- t(scale(t(cc), center=FALSE, scale=rowSums(cc)))
colnames(cc) <- paste('cellType', 1:ncol(cc), sep="")


##evaluate differential expression genes
##calculate 
G <- NULL
for(i in 1:200){
    G <- cbind(G, t(as.matrix(t(as.matrix(cc[i,])) %*% H1_array[,,i])))
}
#noise <- t(matrix(rnorm(nrow(H1)*ng)*5, ncol=ng))
#noise[noise<0] <- 0
#H1 <- H1 + t(noise)
noise <- t(matrix(rnorm(p*ng)*15, ncol=ng))
noise[noise<0] <- 0
G <- G + noise
k <- apply(G,1,var) > 10^-8
G <- G[k,]
H1 <- H1[,k]
H1_array <- H1_array[,k,] 

rownames(G) <- colnames(H1) <- paste0("gene",c(1:nrow(G)))
colnames(G) <- paste0("Sample",1:ncol(G))
rownames(H1) <- colnames(cc)
##load ENIGMA object and perform deconvolution
library(scater)
egm = create_ENIGMA(bulk = G, ref = t(H1), ref_type = "bulk", meta_ref = as.matrix(colnames(t(H1))))
egm = batch_correct(egm)
egm@bulk <- G
egm = get_cell_proportion(egm)
plot_proportion(egm)
egm = ENIGMA_L2_max_norm(egm, epsilon=0.001, alpha=0.7,
                        beta=4500,tao_k=0.01,max.iter=1000,verbose=TRUE)

Fra_Simulate <- list(theta = egm@result_cell_proportion)

enigma <- egm@result_CSE_normalized[,egm@result_CSE_normalized$cell_type %in% "cellType1"]
enigma$Group <- label.c4				 
enigma <- enigma[Fra_Simulate$theta[,1]>0.05]
enigma <- runPCA(enigma)
p <- plotPCA(enigma, colour_by = "Group")
#enigma <- runUMAP(enigma,dimred="PCA",n_dimred=5)
#p_enigma_umap <- plotUMAP(enigma, colour_by = "Group",point_size=3)	
enigma <- runTSNE(enigma,dimred="PCA",n_dimred=5)
p_enigma_tsne1 <- plotTSNE(enigma, colour_by = "Group",point_size=3)	


enigma <- egm@result_CSE_normalized[,egm@result_CSE_normalized$cell_type %in% "cellType3"]
enigma$Group <- label.c3				 
enigma <- enigma[Fra_Simulate$theta[,3]>0.05]
enigma <- runPCA(enigma)
p <- plotPCA(enigma, colour_by = "Group")
#enigma <- runUMAP(enigma,dimred="PCA",n_dimred=5)
#p_enigma_umap <- plotUMAP(enigma, colour_by = "Group",point_size=3)	
enigma <- runTSNE(enigma,dimred="PCA",n_dimred=5)
p_enigma_tsne3 <- plotTSNE(enigma, colour_by = "Group",point_size=3)	


enigma <- egm@result_CSE_normalized[,egm@result_CSE_normalized$cell_type %in% "cellType2"]
enigma$Group <- label.c2				 
enigma <- enigma[Fra_Simulate$theta[,2]>0.05]
enigma <- runPCA(enigma)
p <- plotPCA(enigma, colour_by = "Group")
#enigma <- runUMAP(enigma,dimred="PCA",n_dimred=5)
#p_enigma_umap <- plotUMAP(enigma, colour_by = "Group",point_size=3)	
enigma <- runTSNE(enigma,dimred="PCA",n_dimred=5)
p_enigma_tsne2 <- plotTSNE(enigma, colour_by = "Group",point_size=3)	

########################################################
bulk <- G
time.tca <- system.time({tca.mdl <- tca(X = bulk, W = egm@result_cell_proportion, C1 = NULL, C2 = NULL,
                parallel = TRUE,num_cores=4,max_iters=5)
Z_hat_simulate <- tensor(X = (as.matrix(bulk)), tca.mdl)})

TCA <- SingleCellExperiment(assays=list(logcounts = Z_hat_simulate[[1]]))
TCA$Group <- label.c4				 
TCA <- TCA[Fra_Simulate$theta[,1]>0.05]
TCA <- runPCA(TCA)
p <- plotPCA(TCA, colour_by = "Group")
#TCA <- runUMAP(TCA,dimred="PCA",n_dimred=5)
#p_TCA_umap <- plotUMAP(TCA, colour_by = "Group",point_size=3)	
TCA <- runTSNE(TCA,dimred="PCA",n_dimred=5)
p_TCA_tsne1 <- plotTSNE(TCA, colour_by = "Group",point_size=3)	

TCA <- SingleCellExperiment(assays=list(logcounts = Z_hat_simulate[[2]]))
TCA$Group <- label.c2				 
TCA <- TCA[Fra_Simulate$theta[,2]>0.05]
TCA <- runPCA(TCA)
p <- plotPCA(TCA, colour_by = "Group")
#TCA <- runUMAP(TCA,dimred="PCA",n_dimred=5)
#p_TCA_umap <- plotUMAP(TCA, colour_by = "Group",point_size=3)	
TCA <- runTSNE(TCA,dimred="PCA",n_dimred=5)
p_TCA_tsne2 <- plotTSNE(TCA, colour_by = "Group",point_size=3)	


TCA <- SingleCellExperiment(assays=list(logcounts = Z_hat_simulate[[3]]))
TCA$Group <- label.c3				 
TCA <- TCA[Fra_Simulate$theta[,3]>0.05]
TCA <- runPCA(TCA)
p <- plotPCA(TCA, colour_by = "Group")
#TCA <- runUMAP(TCA,dimred="PCA",n_dimred=5)
#p_TCA_umap <- plotUMAP(TCA, colour_by = "Group",point_size=3)	
TCA <- runTSNE(TCA,dimred="PCA",n_dimred=5)
p_TCA_tsne3 <- plotTSNE(TCA, colour_by = "Group",point_size=3)	


#########################################################
rownames(H1) <- colnames(cc)
tcga.ted <- TED::run.Ted(ref.dat= H1,
X=t(G)*100,
cell.type.labels= rownames(H1),
tum.key="cellType1",
input.type="GEP",
n.cores=5,
pdf.name="Senescence")

exp <- t(norm.to.one(tcga.ted$res$first.gibbs.res$Znkg[,1,]))
ted <- SingleCellExperiment(assays=list(logcounts = exp))
ted$Group <- label.c4				 
ted <- ted
ted <- runPCA(ted)
p <- plotPCA(ted, colour_by = "Group")
#ted <- runUMAP(ted,dimred="PCA",n_dimred=5)
#p_ted_umap <- plotUMAP(ted, colour_by = "Group",point_size=3)	
ted <- runTSNE(ted,dimred="PCA",n_dimred=5)
p_ted_tsne1 <- plotTSNE(ted, colour_by = "Group",point_size=3)


exp <- t(norm.to.one(tcga.ted$res$first.gibbs.res$Znkg[,2,]))
ted <- SingleCellExperiment(assays=list(logcounts = exp))
ted$Group <- label.c2				 
ted <- ted
ted <- runPCA(ted)
p <- plotPCA(ted, colour_by = "Group")
#ted <- runUMAP(ted,dimred="PCA",n_dimred=5)
#p_ted_umap <- plotUMAP(ted, colour_by = "Group",point_size=3)	
ted <- runTSNE(ted,dimred="PCA",n_dimred=5)
p_ted_tsne2 <- plotTSNE(ted, colour_by = "Group",point_size=3)



exp <- t(norm.to.one(tcga.ted$res$first.gibbs.res$Znkg[,3,]))
ted <- SingleCellExperiment(assays=list(logcounts = exp))
ted$Group <- label.c3				 
ted <- ted
ted <- runPCA(ted)
p <- plotPCA(ted, colour_by = "Group")
#ted <- runUMAP(ted,dimred="PCA",n_dimred=5)
#p_ted_umap <- plotUMAP(ted, colour_by = "Group",point_size=3)	
ted <- runTSNE(ted,dimred="PCA",n_dimred=5)
p_ted_tsne3 <- plotTSNE(ted, colour_by = "Group",point_size=3)
###################################################
raw <- SingleCellExperiment(assays=list(logcounts = G))
raw$Group <- label.c4				 
raw <- raw
raw <- runPCA(raw)
p <- plotPCA(raw, colour_by = "Group")
#raw <- runUMAP(raw,dimred="PCA",n_dimred=5)
#p_raw_umap <- plotUMAP(raw, colour_by = "Group",point_size=3)	
raw <- runTSNE(raw,dimred="PCA",n_dimred=5)
p_raw1 <- plotTSNE(raw, colour_by = "Group",point_size=3)

raw$Group <- label.c2
p_raw2 <- plotTSNE(raw, colour_by = "Group",point_size=3)

raw$Group <- label.c3
p_raw3 <- plotTSNE(raw, colour_by = "Group",point_size=3)

p_true <- plot_grid(p_c4,p_c2,p_c3,nrow=3)
p_raw <- plot_grid(p_raw1,p_raw2,p_raw3,nrow=3)
p_enigma <- plot_grid(p_enigma_tsne1,p_enigma_tsne2,p_enigma_tsne3,nrow=3)
p_TCA <- plot_grid(p_TCA_tsne1,p_TCA_tsne2,p_TCA_tsne3,nrow=3)
p_ted <- plot_grid(p_ted_tsne1,p_ted_tsne2,p_ted_tsne3,nrow=3)

p_true <- p_true + ggtitle("Ground Truth") +
    theme(plot.title = element_text(hjust = 0.5))
p_raw <- p_raw + ggtitle("Pesudo-Bulk") +
    theme(plot.title = element_text(hjust = 0.5))	
p_TCA <- p_TCA + ggtitle("TCA") +
    theme(plot.title = element_text(hjust = 0.5))	
p_enigma <- p_enigma + ggtitle("ENIGMA") +
    theme(plot.title = element_text(hjust = 0.5))	
p_ted <- p_ted + ggtitle("BayesPrims") +
    theme(plot.title = element_text(hjust = 0.5))	

png("overall.png",res=300,height=3000,width=6000)
plot_grid(p_true,p_raw,p_TCA,p_enigma,p_ted,ncol=5)
dev.off()
