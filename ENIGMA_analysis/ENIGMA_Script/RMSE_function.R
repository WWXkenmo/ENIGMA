library(data.table)
## normalize result according to bulk read depth
calibrate1 <- function(X_k,Bulk){
  scale = 500*10^5
  index <- NULL
  for(j in 1:dim(X_k)[3]){#cell type
    for(i in 1:ncol(Bulk)){#sample
      if(sum(X_k[,i,j])!=0){
        X_k[,i,j] <- X_k[,i,j]*scale/ sum(X_k[,i,j])##norm
        index <- rbind(index, c(j,i))
      }
    }
  }
  return(
    X_k
  )
}
## normalize ref according to bulk read depth
calibrateRef1 <- function(ref,Bulk){
  scale = 500*10^5
  X_k = array(0,
              dim = c( nrow(Bulk),#gene
                       ncol(Bulk),#sample
                       ncol(ref)),#ct
              dimnames = list( rownames(Bulk),
                               colnames(Bulk),
                               colnames(ref))
  )
  for(j in 1:ncol(ref)){#ct
    for(i in 1:ncol(Bulk)){#sample
      X_k[,i,j] <- ref[,j]*scale/ sum(ref[,j])#norm
    }
  }
  X_k
}
###
RMSE_gene <- function(gt,X_k){#groundtrue
  rmse_m <- NULL
  for(i in 1:length(gt)){#ct
    rmse <- NULL
    for(j in intersect(rownames(gt[[i]]),rownames(X_k[,,i]))){#gene
      rmse <- c(rmse, norm(log2(gt[[i]][j,]*500+1) - log2(X_k[j,,i]+1),"2"))
    }
    rmse_m <- cbind(rmse_m,rmse)
  }
  colnames(rmse_m)=names(gt)
  rownames(rmse_m)=intersect(rownames(gt[[1]]),rownames(Bulk))
  rmse_m
}

RMSE_sample <- function(gt,X_k){
  rmse_m <- NULL
  for(i in 1:length(gt)){
    rmse <- NULL
    for(j in colnames(gt[[i]])){#sample
      rmse <- c(rmse, norm(log2(gt[[i]][,j]*500+1) - log2(X_k[,j,i]+1),"2"))
    }
    rmse_m <- cbind(rmse_m,rmse)
  }
  colnames(rmse_m)=names(gt)
  rownames(rmse_m)=colnames(gt[[1]])
  rmse_m
}

##################
ref_RMSE1<- function(ref,Bulk,gt,level){
  X_k <- calibrateRef1(ref,Bulk)
  if(level == "sample") res <- RMSE_sample(gt,X_k)
  if(level == "gene") res <- RMSE_gene(gt,X_k)
  res
}


bMIND_RMSE1 <- function(list,gt,Bulk,process,level){
  res_list <- list()
  
  for(model in 1:length(list)){
    nsample <- dim(list[[1]]$A)[3]
    ngene <- dim(list[[1]]$A)[1]
    ncelltype <- dim(list[[1]]$A)[2]
    name_s <- dimnames(list[[1]]$A)[[3]]
    name_g <- dimnames(list[[1]]$A)[[1]]
    name_c <- dimnames(list[[1]]$A)[[2]]
    deconv_X = array(0,
                     dim = c( ngene,
                              nsample,
                              ncelltype),
                     dimnames = list( name_g,
                                      name_s,
                                      name_c)
    )
    for(i in 1:ncelltype){
      deconv_X[,,i] = list[[model]]$A[,i,]
    }
    deconv_X[deconv_X<0] <- 0
    ##calibrate
    if(process[model]=="log") X_k <- calibrate1(2^deconv_X-1,Bulk)
    if(process[model]=="sqrt") X_k <- calibrate1(deconv_X^2,Bulk)
    if(process[model]=="none") X_k <- calibrate1(deconv_X,Bulk)
    
    if(level == "sample") res <- RMSE_sample(gt,X_k)
    if(level == "gene") res <- RMSE_gene(gt,X_k)
    res_list[[model]] <- res
  }
  
  res_list
}


TCA_RMSE1 <- function(list,gt,Bulk,process,level){
  res_list <- list()
  
  for(model in 1:length(list)){
    nsample <- dim(list[[1]][[1]])[2]
    ngene <- dim(list[[1]][[1]])[1]
    ncelltype <- length(list[[1]])
    name_s <- colnames(list[[1]][[1]])
    name_g <- rownames(list[[1]][[1]])
    name_c <- names(list[[1]])
    deconv_X = array(0,
                     dim = c( ngene,
                              nsample,
                              ncelltype),
                     dimnames = list( name_g,
                                      name_s,
                                      name_c)
    )
    for(i in 1:ncelltype){
      deconv_X[,,i] = list[[model]][[i]]
    }
    deconv_X[deconv_X<0] <- 0
    ##calibrate
    if(process[model]=="log") X_k <- calibrate1(2^deconv_X-1,Bulk)
    if(process[model]=="sqrt") X_k <- calibrate1(deconv_X^2,Bulk)
    
    if(level == "sample") res <- RMSE_sample(gt,X_k)
    if(level == "gene") res <- RMSE_gene(gt,X_k)
    res_list[[model]] <- res
  }
  names(res_list)=names(list)
  res_list
}

ENIGMA_RMSE1 <- function(list,gt,Bulk,process,level){
  res_list <- list()
  
  for(model in 1:length(list)){
    nsample <- dim(list[[1]]$X_k)[2]
    ngene <- dim(list[[1]]$X_k)[1]
    ncelltype <- dim(list[[1]]$X_k)[3]
    name_s <- dimnames(list[[1]]$X_k)[[2]]
    name_g <- dimnames(list[[1]]$X_k)[[1]]
    name_c <- dimnames(list[[1]]$X_k)[[3]]
    deconv_X = array(0,
                     dim = c( ngene,
                              nsample,
                              ncelltype),
                     dimnames = list( name_g,
                                      name_s,
                                      name_c)
    )
    for(i in 1:ncelltype){
      deconv_X[,,i] = list[[model]]$X_k[,,i]
    }
    deconv_X[deconv_X<0] <- 0
    ##calibrate
    if(process[model]=="log") X_k <- calibrate1(2^deconv_X-1,Bulk)
    if(process[model]=="sqrt") X_k <- calibrate1(deconv_X^2,Bulk)
    
    if(level == "sample") res <- RMSE_sample(gt,X_k)
    if(level == "gene") res <- RMSE_gene(gt,X_k)
    res_list[[model]] <- res
  }
  names(res_list)=names(list)
  res_list
}

