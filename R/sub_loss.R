loss <- function(X, P_old, theta, alpha,beta,R){
    # X: Bulk gene expression dataset (g*n)
    # P_old: cell type specific gene expression profile (g*n*p)
    # theta: cell type ratio for each samples (n*p)
    # alpha: constraint parameters of the similarity between each estimated cell type specific expression and reference profile, constant
    # miu:  constraint parameters of the smoothness of gene expression, constant
    # R: reference profile (g*p)

    part1 <- 0
    for(i in 1:ncol(theta)){
        part1 <- part1+P_old[,,i]%*%diag(theta[,i])
    }
    part1 <- part1
    # part1 <- norm((X-part1),"F")^2
    part1 <- sum( (X-part1)^2 )

    part2 <- 0
    for(i in 1:ncol(R)){
        ref <- matrix(rep(R[,i],ncol(X)),nrow=length(R[,i]))
        # part2 <- part2 + alpha*norm((P_old[,,i]-ref),"F")^2
        part2 <- part2 + sum( (rowMeans(P_old[,,i])-ref)^2 )
    }

    part3 <- 0
    for(i in 1:ncol(R)){
        part3 <- part3 + sum(svd(P_old[,,i])$d)
    }

    res <- NULL
    res$part1 <- part1*(alpha/2)
    res$part2 <- part2*(1-alpha)*(1/2)
    res$part3 <- beta*(part3)
    res
}
