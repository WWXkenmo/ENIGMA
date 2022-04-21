###########################################################
#Simulate the cell type with the path
##############################################################################
###Using ESCO to simulate cell states
source("ENIGMA.R")
library(ESCO)
sim <- escoSimulateTraj(nGenes = 8368, nCells = 10000,paths.design = data.frame(
                          Path = c(1, 2, 3),
                          From = c(0, 1, 2),
                          Steps = c(1000, 1000, 1000)
                        ),
                        withcorr = TRUE,verbose=FALSE)

####Select one cell types
sim_traj <- sim[,sim$Path %in% 3]

###Class cell into 500 step bins
ind <- NULL
for(n in 1:200){
ind <- c(ind,(n-1)*5 + 1)
}
end <- NULL
for(n in 1:200){
end <- c(end,(n)*5 )
}
time <- (ind + end) / 2

step_aggre <- NULL
for(i in 1:length(sim_traj$Step)){
  id <- which(ind <= sim_traj$Step[i])
  id <- id[length(id)]
  if(sim_traj$Step[i] <= end[id]){
    step_aggre <- c(step_aggre,time[id])
  }
}

sim_traj <- assays(sim_traj)$TrueCounts


################
#Simulate other cells
params <- escoEstimate(sim_traj)

sim<-escoSimulateGroups(params=params)
sim <- assays(sim)$TrueCounts

###Simulate profile with two groups
sim2 <-escoSimulateGroups(params=params,
                        group.prob = c(0.5, 0.5), withcorr = TRUE, verbose =TRUE)
group <- sim2$Group
sim2 <- assays(sim2)$TrueCounts

rownames(sim_traj) <- rownames(sim) <-rownames(sim2)
#########################
##Randomly seperate the cells into two groups, once for generating CSE and Pseudo-bulk, the other generate the reference
sample.ref.ind = sample(1:ncol(sim2),0.1*ncol(sim2),replace=FALSE)
ref.sim1 = sim[,sample.ref.ind]
ref.sim2 = sim2[,sample.ref.ind];
ref.Traj = sim_traj[,sample.ref.ind];

##
bulk.sim1 = sim[,-sample.ref.ind]
bulk.sim2 = sim2[,-sample.ref.ind];bulk2.group = group[-sample.ref.ind]
bulk.Traj = sim_traj[,-sample.ref.ind];bulk.Traj.group = step_aggre[-sample.ref.ind]


##################
#Assign the CSE identity
idCell = matrix(NA,nrow=200,ncol=3)
for(i in 1:nrow(idCell)){
  idCell[i,] = c(NA,NA,sample(c("Group1","Group2"),1))
}
idCell[,1] = sample(names(table(bulk.Traj.group)),200,replace=FALSE)
idCell = as.data.frame(idCell)
colnames(idCell) = c("Cell.Trajectory","CellType2","CellType3")
rownames(idCell) = paste0("Sample-",1:200,sep="")

#################
####Generate CSE profile
H1_array <- array(0,
                  dim = c( 3,
                           8368,
                           200))
for(i in 1:dim(H1_array)[3]){
   ##For each CSE, we fixed sampled 50 cells to generate profile
   bulk.Traj.sub = bulk.Traj[,bulk.Traj.group %in% idCell$Cell.Trajectory[i]]
   bulk.Traj.sub = bulk.Traj.sub[,sample(1:ncol(bulk.Traj.sub),50,replace=TRUE)]
   ####Normalization and calculate means
   bulk.Traj.sub = bulk.Traj.sub %*% diag(10^4/colSums(bulk.Traj.sub))
   
   ##CellType2
   bulk.sim1.sub = bulk.sim1[,sample(1:ncol(bulk.sim1),50,replace=TRUE)]
   ####Normalization and calculate means
   bulk.sim1.sub = bulk.sim1.sub %*% diag(10^4/colSums(bulk.sim1.sub))
   
   ##CellType3
   bulk.sim2.sub = bulk.sim2[,bulk2.group %in% idCell$CellType3[i]]
   bulk.sim2.sub = bulk.sim2.sub[,sample(1:ncol(bulk.sim2.sub),50,replace=TRUE)]
   ####Normalization and calculate means
   bulk.sim2.sub = bulk.sim2.sub %*% diag(10^4/colSums(bulk.sim2.sub))
   
   profile = cbind(rowMeans(bulk.Traj.sub),rowMeans(bulk.sim1.sub),rowMeans(bulk.sim2.sub))
   H1_array[,,i] = t(profile)
}


###simulate the cell type fractions
k <- 3 # number of cell types
ng <- 8368 # number of genes
p <- 200 # number of samples

cc <- matrix(runif(p*k), ncol=k)
cc <- t(scale(t(cc), center=FALSE, scale=rowSums(cc)))
colnames(cc) <- colnames(idCell)
rownames(cc) <- rownames(idCell)

###Generate the pesudo-bulk
PseudoBulk <- NULL
for(i in 1:200){
    PseudoBulk <- cbind(PseudoBulk, t(as.matrix(t(as.matrix(cc[i,])) %*% H1_array[,,i])))
}
rownames(PseudoBulk) = rownames(bulk.Traj.sub)
colnames(PseudoBulk) = rownames(idCell)
###Attach Noise
noise <- t(matrix(rnorm(p*ng)*10, ncol=ng))
noise[noise<0] <- 0
PseudoBulk <- PseudoBulk + noise

#calculate reference
ref = cbind(rowMeans(ref.Traj %*% diag(10^4/colSums(ref.Traj))),rowMeans(ref.sim1 %*% diag(10^4/colSums(ref.sim1))),rowMeans(ref.sim2 %*% diag(10^4/colSums(ref.sim2))))
colnames(ref) = colnames(idCell)


########################################################################################
#Evaluation
##calculate cell fraction matrix
Frac <- get_proportion(PseudoBulk, ref)

##Running TCA
tca.mdl = tca(X = sqrt(PseudoBulk), W = Frac$theta, C1 = NULL, C2 = NULL,
                parallel = TRUE,num_cores=4,max_iters=20)
TCA_res2 <- tensor(X = (as.matrix(sqrt(PseudoBulk))), tca.mdl)


TCA <- SingleCellExperiment(assays=list(logcounts = (TCA_res2[[1]])))
TCA$Step <- as.numeric(as.matrix(idCell$Cell.Trajectory))				 
TCA <- TCA[,Frac$theta[,1]>0.05]
ENIGMA_diffusion <- DiffusionMap(TCA,n_local = 50, n_pcs = 5,verbose = TRUE)
DC_embedding <- cbind(ENIGMA_diffusion$DC1,ENIGMA_diffusion$DC2)
colnames(DC_embedding) <- c("DC1","DC2")
reducedDim(TCA, "DiffusionMap") <- DC_embedding

#TCA$Step <- as.numeric(as.matrix(idCell$Cell.Trajectory))				 
#TCA <- TCA[,Frac$theta[,1]>0.05]
#TCA <- runPCA(TCA)
#TCA <- runUMAP(TCA,dimred="PCA",n_dimred=5)

png("TCA_traj.png",res=300,height=1000,width=1500)
plotReducedDim(TCA,colour_by = "Step",dimred = "DiffusionMap",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))
dev.off()




#TCA$Step <- as.numeric(as.matrix(idCell$Cell.Trajectory))				 
#TCA <- TCA[,Frac$theta[,1]>0.05]
#TCA <- runPCA(TCA)
#TCA <- runUMAP(TCA,dimred="PCA",n_dimred=5)

exp <- TCA_res2[[1]]^2
exp.scale <- t(apply(exp,1,scale))
###chose the PC with the highest correlation with cell type fractions
d <- sqrt(svd(exp.scale)$d)
d <- d / sum(d)
prob_d <- NULL;for(i in 1:length(d)) prob_d <- c(prob_d, sum(d[1:i]))
PC <- svd(exp.scale)$v[,1:which(prob_d>0.8)[1]]
pc_cor <- apply(PC,2,function(x){cor(x,Frac$theta[,k],method="sp")})
PC <- PC[,which.max(abs(pc_cor))]
the <- (exp %*% as.matrix(PC) - length(PC) * mean(PC) * rowMeans(exp)) / (sum(PC^2) - length(PC)*mean(PC)^2)
exp.norm <- exp - as.matrix(the) %*% t(as.matrix(PC))

TCA2 <- SingleCellExperiment(assays=list(logcounts = exp.norm))
TCA2$Step <- as.numeric(as.matrix(idCell$Cell.Trajectory))				 
TCA2 <- TCA2[,Frac$theta[,1]>0.05]
ENIGMA_diffusion <- DiffusionMap(TCA2,n_local = 50, n_pcs = 5,verbose = TRUE)
DC_embedding <- cbind(ENIGMA_diffusion$DC1,ENIGMA_diffusion$DC2)
colnames(DC_embedding) <- c("DC1","DC2")
reducedDim(TCA2, "DiffusionMap") <- DC_embedding

png("TCA_traj2.png",res=300,height=1000,width=1500)
plotReducedDim(TCA2,colour_by = "Step",dimred = "DiffusionMap",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))
dev.off()
####Running ENIGMA
egm2 <- cell_deconvolve(X=as.matrix(sqrt(PseudoBulk)),
                                    theta=Frac$theta,
                                    R=sqrt(ref),
                                    epsilon=0.001,
                                    alpha=0.1,
                                    beta=0.5,tao_k=0.01,max.iter=1000,verbose=TRUE,pre.process = "sqrt",Norm.method = "PC")

###Running Trace Norm
egm_trace <- cell_deconvolve_trace(O = as.matrix(sqrt(PseudoBulk)),
                                                theta=Frac$theta,
                                                R=sqrt(ref),
                                                alpha=0.1,beta=1,solver="admm",
                                                verbose=TRUE,max.iter = 1000,pre.process = "sqrt",Norm.method = "PC")

#egm2$X_k_norm[egm2$X_k_norm<0] <- 0
enigma <- SingleCellExperiment(assays=list(logcounts = egm2$X_k_norm[,,1]))
enigma$Step <- as.numeric(as.matrix(idCell$Cell.Trajectory))				 
enigma <- enigma[,Frac$theta[,1]>0.05]
ENIGMA_diffusion <- DiffusionMap(enigma,n_local = 50, n_pcs = 5,verbose = TRUE)
DC_embedding <- cbind(ENIGMA_diffusion$DC1,ENIGMA_diffusion$DC2)
colnames(DC_embedding) <- c("DC1","DC2")
reducedDim(enigma, "DiffusionMap") <- DC_embedding

#enigma$Step <- as.numeric(as.matrix(idCell$Cell.Trajectory))				 
#enigma <- enigma[,Frac$theta[,1]>0.05]
#enigma <- runPCA(enigma)
#enigma <- runUMAP(enigma,dimred="PCA",n_dimred=5)

png("enigma_traj.png",res=300,height=1000,width=1500)
plotReducedDim(enigma,colour_by = "Step",dimred = "DiffusionMap",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))
dev.off()

enigma_trace <- SingleCellExperiment(assays=list(logcounts = egm_trace$X_k_norm[,,1]))
enigma_trace$Step <- as.numeric(as.matrix(idCell$Cell.Trajectory))				 
enigma_trace <- enigma_trace[,Frac$theta[,1]>0.05]
ENIGMA_diffusion <- DiffusionMap(enigma_trace,n_local = 50, n_pcs = 5,verbose = TRUE)
DC_embedding <- cbind(ENIGMA_diffusion$DC1,ENIGMA_diffusion$DC2)
colnames(DC_embedding) <- c("DC1","DC2")
reducedDim(enigma_trace, "DiffusionMap") <- DC_embedding

png("enigma_traj(trace).png",res=300,height=1000,width=1500)
plotReducedDim(enigma_trace,colour_by = "Step",dimred = "DiffusionMap",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))
dev.off()

#########################
##Using bMIND to evaluate
bmind_res = bMIND2(sqrt(PseudoBulk),ncore = 3, frac = Frac$theta, profile = sqrt(ref), noRE = FALSE)
####regress out the fraction effects
exp <- bmind_res$A[,1,]^2
exp.scale <- t(apply(exp,1,scale))
###chose the PC with the highest correlation with cell type fractions
d <- sqrt(svd(exp.scale)$d)
d <- d / sum(d)
prob_d <- NULL;for(i in 1:length(d)) prob_d <- c(prob_d, sum(d[1:i]))
PC <- svd(exp.scale)$v[,1:which(prob_d>0.8)[1]]
pc_cor <- apply(PC,2,function(x){cor(x,Frac$theta[,k],method="sp")})
PC <- PC[,which.max(abs(pc_cor))]
the <- (exp %*% as.matrix(PC) - length(PC) * mean(PC) * rowMeans(exp)) / (sum(PC^2) - length(PC)*mean(PC)^2)
exp.norm <- exp - as.matrix(the) %*% t(as.matrix(PC))

bmind2 <- SingleCellExperiment(assays=list(logcounts = exp.norm ))
bmind2$Step <- as.numeric(as.matrix(idCell$Cell.Trajectory))				 
bmind2 <- bmind2[,Frac$theta[,1]>0.05]
ENIGMA_diffusion <- DiffusionMap(bmind2,n_local = 50, n_pcs = 5,verbose = TRUE)
DC_embedding <- cbind(ENIGMA_diffusion$DC1,ENIGMA_diffusion$DC2)
colnames(DC_embedding) <- c("DC1","DC2")
reducedDim(bmind2, "DiffusionMap") <- DC_embedding


png("bmind_traj2.png",res=300,height=1000,width=1500)
plotReducedDim(bmind2,colour_by = "Step",dimred = "DiffusionMap",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))
dev.off()
########################

bmind <- SingleCellExperiment(assays=list(logcounts = bmind_res$A[,1,] ))
bmind$Step <- as.numeric(as.matrix(idCell$Cell.Trajectory))				 
bmind <- bmind[,Frac$theta[,1]>0.05]
ENIGMA_diffusion <- DiffusionMap(bmind,n_local = 50, n_pcs = 5,verbose = TRUE)
DC_embedding <- cbind(ENIGMA_diffusion$DC1,ENIGMA_diffusion$DC2)
colnames(DC_embedding) <- c("DC1","DC2")
reducedDim(bmind, "DiffusionMap") <- DC_embedding

png("bmind_traj.png",res=300,height=1000,width=1500)
plotReducedDim(bmind,colour_by = "Step",dimred = "DiffusionMap",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))
dev.off()


########################
##Using raw expression profile to compare
bulk_plot <- SingleCellExperiment(assays=list(logcounts = PseudoBulk))
bulk_plot$Step <- as.numeric(as.matrix(idCell$Cell.Trajectory))				 
bulk_plot <- bulk_plot[,Frac$theta[,1]>0.05]
ENIGMA_diffusion <- DiffusionMap(bulk_plot,n_local = 50, n_pcs = 5,verbose = TRUE)
DC_embedding <- cbind(ENIGMA_diffusion$DC1,ENIGMA_diffusion$DC2)
colnames(DC_embedding) <- c("DC1","DC2")
reducedDim(bulk_plot, "DiffusionMap") <- DC_embedding

png("bulk_plot_traj.png",res=300,height=1000,width=1500)
plotReducedDim(bulk_plot,colour_by = "Step",dimred = "DiffusionMap",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))
dev.off()

########################
##plot the ground truth
groundtruth <- SingleCellExperiment(assays=list(logcounts = H1_array[1,,]))
groundtruth$Step <- as.numeric(as.matrix(idCell$Cell.Trajectory))				 
groundtruth <- groundtruth[,Frac$theta[,1]>0.05]
ENIGMA_diffusion <- DiffusionMap(groundtruth,n_local = 50, n_pcs = 5,verbose = TRUE)
DC_embedding <- cbind(ENIGMA_diffusion$DC1,ENIGMA_diffusion$DC2)
colnames(DC_embedding) <- c("DC1","DC2")
reducedDim(groundtruth, "DiffusionMap") <- DC_embedding

png("groundtruth_traj.png",res=300,height=1000,width=1500)
plotReducedDim(groundtruth,colour_by = "Step",dimred = "DiffusionMap",point_size = 3,point_alpha=1)+theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),text = element_text(size=15))
dev.off()

#######################
#Evaluate through R2 with PseudoTime
dat = cbind(enigma$Step,reducedDim(enigma, "DiffusionMap"))
colnames(dat) = c("Step","DC1","DC2")
dat = as.data.frame(dat)
lm.model = lm(Step~., dat)
cor(lm.model$fitted.value, dat$Step,method="sp") # 0.4931314


dat = cbind(bulk_plot$Step,reducedDim(bulk_plot, "DiffusionMap"))
colnames(dat) = c("Step","DC1","DC2")
dat = as.data.frame(dat)
lm.model = lm(Step~., dat)
cor(lm.model$fitted.value, dat$Step) # 0.03373622


dat = cbind(TCA$Step,reducedDim(TCA, "DiffusionMap"))
colnames(dat) = c("Step","DC1","DC2")
dat = as.data.frame(dat)
lm.model = lm(Step~., dat)
cor(lm.model$fitted.value, dat$Step) # 0.06275233

dat = cbind(TCA$Step,reducedDim(bmind, "DiffusionMap"))
colnames(dat) = c("Step","DC1","DC2")
dat = as.data.frame(dat)
lm.model = lm(Step~., dat)
cor(lm.model$fitted.value, dat$Step) # 0.08520045

dat = cbind(enigma_trace$Step,reducedDim(enigma_trace, "DiffusionMap"))
colnames(dat) = c("Step","DC1","DC2")
dat = as.data.frame(dat)
lm.model = lm(Step~., dat)
cor(lm.model$fitted.value, dat$Step) # 0.02411792
















