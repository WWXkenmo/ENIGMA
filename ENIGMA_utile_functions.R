##################utile functions#####################
squash <- function(V, beta){
    ## squash: calculate the optimal solution of the formula: X=argmin{ (||X-V||_F)^2 + beta*||X||_2_max }
    n <- NULL
    for(i in 1:nrow(V)){
        n <- c(n, sqrt( sum(V[i,]^2) ) )
    }
    pi <- order(n,decreasing=TRUE)
    s <- NULL
    for(i in 1:length(pi)){
        s <- c(s, sum(n[pi[1:i]]))
    }
    q <- max(which(n[pi]>=s/(c(1:length(s))+beta)))
    tao <- s[q]/(q+beta)
    
    for(i in 1:q){
        V[pi[i],] <- tao*V[pi[i],]/sqrt( sum(V[pi[i],]^2) )
    }
    
    V
}

proximalpoint <- function(P, tao_k,dP,miu){
    # X: Bulk gene expression dataset (g*n)
    # P_old: cell type specific gene expression profile (g*n*p)
    # theta: cell type ratio for each samples (n*p)
    # alpha: constraint parameters of the similarity between each estimated cell type specific expression and reference profile, constant
    # miu:  constraint parameters of the smoothness of gene expression, constant
    # R: reference profile (g*p)
    # P: the ith cell type specific gene expression profile needs to be undated
    # tao_k: gradient size
    # dP: gradient of matrix P
    # scale_alpha: the parameters for inequality decision
    # miu:  constraint parameters of the smoothness of gene expression, constant
    # cell_type_index: optimize which type of cells
    # gamma: the parameters for inequality decision
    
    P_hat <- t(squash(t(P-tao_k*dP),tao_k*miu))
    ##update P matrix
    return(P_hat)
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
    dP1 = dP1 / sqrt( sum( dP1^2 ) ) * 1e5
    dP2 = dP2 / sqrt( sum( dP2^2 ) ) * 1e5
    
    #calculate w1
    #if( crossprod(as.matrix(dP1), as.matrix(dP2)) >= crossprod(as.matrix(dP1)) ) {w1 = 1}
    #else if( crossprod(as.matrix(dP1), as.matrix(dP2)) >= crossprod(as.matrix(dP2)) ) {w1 = 0}
    #else {
    #    w1 = crossprod(as.matrix(dP2-dP1), as.matrix(dP2))/sum((dP1-dP2)^2)
    #}
    w1 <- alpha
    w2 <- 1-w1
    
    dP <- dP1*as.numeric(w1) + dP2*as.numeric(w2)
    return(dP)
}
