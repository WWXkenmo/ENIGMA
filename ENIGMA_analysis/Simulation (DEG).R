############################################
##Simulation code for DEG detection
SNR <- 1.8 ## chose the signal to noise ratio level and generate expression profile for deconvolution
k <- 5 # number of cell types
ng <- 1000 # number of genes
p <- 100 # number of samples
ndiff <- 0.1*1000 # number of genes differentially expressed

H1 <- matrix(rnorm(5*ng), ncol=ng) ## generate the expression profile for control group
H2 <- H1 ## generate the expression profile for case group

# create differential expression for each cell type
DEG_list <- list()
seq <- 1:ncol(H2)
for(i in 1:nrow(H2)){
    DEG_id <- sample(1:ncol(H2),ndiff,replace=FALSE) ## randomly select the DEG
	add <- sample(c(SNR, -SNR),ndiff,replace=TRUE) ## randomly generate up-regulated or down-regulated genes
    H2[i,DEG_id] <- H2[i,DEG_id] + add ## add to the case group
    DEG_list[[i]] <- DEG_id * sign(add) ## generate the DEG list
    #seq <- seq[seq %in% DEG_id == FALSE]
}

# cell frequency matrix per sample
cc <- matrix(runif(p*k), ncol=k)
cc <- t(scale(t(cc), center=FALSE, scale=rowSums(cc)))
colnames(cc) <- paste('cellType', 1:ncol(cc), sep="")

###Add the noise to each cell type and generate a sample specific profile
H1_array <- array(0,
                  dim = c( nrow(H1),
                           ncol(H1),
                           (floor(100/2))))
for(i in 1:(floor(p/2))){
    H1_array[,,i] <- H1 + matrix(rnorm(5*ng), ncol=ng)
}
H2_array <- array(0,
                  dim = c( nrow(H2),
                           ncol(H2),
                           (floor(100/2))))
for(i in 1:(floor(p/2))){
    H2_array[,,i] <- H2 + matrix(rnorm(5*ng), ncol=ng)
}

##generate bulk GEP
G <- NULL
for(i in 1:50){
    G <- cbind(G, t(as.matrix(t(as.matrix(cc[i,])) %*% H1_array[,,i])))
}
for(i in 1:50){
    G <- cbind(G, t(as.matrix(t(as.matrix(cc[i+50,])) %*% H2_array[,,i])))
}
G <- G + t(matrix(rnorm(100*ng), ncol=ng))
# sample classes (2 groups)
y <- gl(2, p/2)

#################################
##return the bulk and reference profile

save(G,H1,DEG_list,file=paste("/Data/to/Path/Simulation/SNR_",SNR,"_Data.rdata",sep=""))
## our simulated datasets are provided in https://github.com/WWXkenmo/ENIGMA/tree/main/ENIGMA_analysis/Data/cts_DEG
