#L2-max Norm proximal point solve compare
cell_deconvolve <- function(X, theta, R, alpha=0.5, tao_k=0.005,beta=0.5,epsilon=0.001,max.iter=1000,verbose=FALSE,infer=FALSE,loss_his = TRUE,pos=TRUE,pre.process="log",Normalize=TRUE,Norm.method = "PC"){
    # unify geneid between X and R
    geneid = intersect( rownames(X), rownames(R) )
    X = X[geneid,]
    R = R[geneid,]
    
    # initialize
    P_old = array(0,
                  dim = c( nrow(X),
                           ncol(X),
                           ncol(theta)),
                  dimnames = list( rownames(X),
                                   colnames(X),
                                   colnames(theta))
    )
    for(i in 1:ncol(theta)){
        P_old[,,i] <- X
    }
    loss <- sub_loss(X, P_old, theta, alpha, beta, R)
    loss_new <- -1000
    delta <- abs(loss_new-loss$val)
    ###update iteractively
    P_old_new <- P_old
    P_old_pre <- P_old
    iter <- 1
    
    if(verbose) cat(date(), 'Optimizing cell type specific expression profile... \n')
    iter.exp <- 1
    loss <- NULL
	norm_all <- NULL
    repeat{
        ratio <- NULL
        dP <- derive_P2(X, theta,P_old,R,alpha)
		norm <- dP$norm
		dP <- dP$dP
        for(i in 1:ncol(theta)){
            P_hat <- proximalpoint(P_old[,,i], tao_k,dP[,,i],beta)
            P_old_new[,,i] <- P_hat
            
            ratio <- c(ratio, sum( (P_hat-P_old[,,i])^2 ))
        }
        if(verbose) writeLines( sprintf("   Ratio ranges from: %f - %f", min(ratio), max(ratio) ) )
        r <- sub_loss(X, P_old, theta, alpha,beta,R)
        if(loss_his) loss<- rbind(loss,c(r$part1,r$part2,r$part3))
        if(max(ratio) < epsilon||iter.exp >= max.iter){break}else{
            P_old <- P_old_new
            iter.exp <- iter.exp + 1
			norm_all <- rbind(norm_all,norm)
        }
    }
    
    
    ### optimized theta
    ### simply perform robust linear regression model
    if(verbose) cat('Optimizing cell type proportions... \n')
    if(infer){
        theta_new <- NULL
        for(j in 1:ncol(X)){
            Exp <- as.matrix(X[,j])
            rownames(Exp) <- rownames(P_old[,j,])
            colnames(Exp) <- colnames(X)[j]
            x <- P_old[,j,]
            x <- apply(x,2,scale)
            lm.o <- rlm(Exp ~ as.matrix(x),maxit=150)
            coef.v <- lm.o$coefficients[-1]
            coef.v[which(coef.v < 0)] <- 0
            total <- sum(coef.v)
            coef.v <- coef.v/total
            theta_new <- rbind(theta_new,coef.v)
        }
        colnames(theta_new) <- colnames(theta)
        rownames(theta_new) <- colnames(X);theta <- theta_new
    }
    ### optimize theta
    ### take the gradient of all theta and running gradient decent
    if(pos) P_old[P_old<0] <- 0
    loss_new.obj <- sub_loss(X, P_old, theta, alpha,beta,R)
    loss.obj <- sub_loss(X, P_old_pre, theta, alpha,beta,R)
    if(verbose) writeLines( sprintf("Total delta_loss: %f, %s", abs(loss_new.obj$val-loss.obj$val), date() ) )
    if(verbose) writeLines( paste("part1:",loss_new.obj$part1," part2:",loss_new.obj$part2," part3:",loss_new.obj$part3,sep="") )   
    
    if(Normalize){
        if(pre.process == "log") X_k_m <- 2^P_old - 1
        if(pre.process == "sqrt") X_k_m <- P_old^2
        if(pre.process == "none") X_k_m <- P_old
        pp_list = c("log","sqrt","none")
        msg = paste0("should be one of ", paste(pp_list, collapse = ", "), 
                     ".")
        if (!pre.process %in% pp_list) {
            stop("'preprocess method' ", msg)
        }
        if(verbose) cat("Perform Normalization...")
        X_k_norm <- X_k_m
        if(Norm.method == "PC"){
            for(k in 1:dim(X_k_m)[3]){
                exp <- X_k_m[,,k]
                exp.scale <- t(apply(exp,1,scale))
                ###chose the PC with the highest correlation with cell type fractions
                d <- sqrt(svd(exp.scale)$d)
                d <- d / sum(d)
                prob_d <- NULL;for(i in 1:length(d)) prob_d <- c(prob_d, sum(d[1:i]))
                PC <- svd(exp.scale)$v[,1:which(prob_d>0.8)[1]]
                pc_cor <- apply(PC,2,function(x){cor(x,theta[,k],method="sp")})
                PC <- PC[,which.max(abs(pc_cor))]
                the <- (exp %*% as.matrix(PC) - length(PC) * mean(PC) * rowMeans(exp)) / (sum(PC^2) - length(PC)*mean(PC)^2)
                exp.norm <- exp - as.matrix(the) %*% t(as.matrix(PC))
                X_k_norm[,,k] <- exp.norm
            }
        }
        
        if(Norm.method == "frac"){
            for(k in 1:dim(X_k_m)[3]){
                exp <- X_k_m[,,k]
                PC <- theta[,k]
                the <- (exp %*% as.matrix(PC) - length(PC) * mean(PC) * rowMeans(exp)) / (sum(PC^2) - length(PC)*mean(PC)^2)
                exp.norm <- exp - as.matrix(the) %*% t(as.matrix(PC))
                X_k_norm[,,k] <- exp.norm
            }
        }
    }
    writeLines( paste("Done! \nConverge in ",iter.exp," steps",sep="") )
    # return cell type specific gene expression matrix
    # return(P_old)
    res = list(
        X_k = P_old,
        theta = theta,
        loss_history = loss,
		gradient_norm = norm_all
    )
    if(Normalize) res$X_k_norm = X_k_norm
    return( res )
}
derive_P2 <- function(X, theta, P_old,R,alpha){
    ## P_old: a tensor variable with three dimensions
    ## theta: the cell type proportions variable
    ## cell_type_index: optimize which type of cells
    ## R: reference matrix
    dP1 <- dP2 <- array(0,
                        dim = c( nrow(X),
                                 ncol(X),
                                 ncol(theta)),
                        dimnames = list( rownames(X),
                                         colnames(X),
                                         colnames(theta))
    )
    for(cell_type_index in 1:ncol(theta)){
        R.m <- as.matrix(R[,cell_type_index])
        
        cell_type_seq <- c(1:ncol(theta))
        cell_type_seq <- cell_type_seq[cell_type_seq!=cell_type_index]
        
        X_summary = Reduce("+",
                           lapply(cell_type_seq, function(i) P_old[,,i]%*%diag(theta[,i]) )
        )
        X_summary <- X-X_summary
        
        dP1[,,cell_type_index] <- 2*(P_old[,,cell_type_index]%*%diag(theta[,cell_type_index]) - X_summary)%*%diag(theta[,cell_type_index])
        dP2[,,cell_type_index] <- 2*(as.matrix(rowMeans(P_old[,,cell_type_index]))-R.m)%*%t(as.matrix(rep((1/ncol(dP2[,,cell_type_index])),ncol(dP2[,,cell_type_index]))))
    }

    norm1 = sqrt( sum( dP1^2 ) )
	norm2 = sqrt( sum( dP2^2 ) )
    #calculate w1
    #if( crossprod(as.matrix(dP1), as.matrix(dP2)) >= crossprod(as.matrix(dP1)) ) {w1 = 1}
    #else if( crossprod(as.matrix(dP1), as.matrix(dP2)) >= crossprod(as.matrix(dP2)) ) {w1 = 0}
    #else {
    #    w1 = crossprod(as.matrix(dP2-dP1), as.matrix(dP2))/sum((dP1-dP2)^2)
    #}
    w1 <- alpha
    w2 <- 1-w1
    
    dP <- dP1*as.numeric(w1) + dP2*as.numeric(w2)
    res <- list()
	res$dP <- dP
	res$norm <- c(norm1,norm2)
	res
}