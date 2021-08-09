#' @title Calculate the GEP of each cell types in bulk samples
#'
#' @param O
#' @param theta
#' @param R
#' @param alpha
#' @param beta
#' @param gamma
#' @param epsilon
#' @param max.iter
#' @param verbose
#'
#' @return
#' @export
#'
cell_deconvolve_trace <- function(O, theta, R, alpha=0.5,beta=5,gamma=1,epsilon=0.001,max.iter=100,verbose=FALSE){
    X = array(0,
              dim = c( nrow(O),
                       ncol(O),
                       ncol(theta)),
              dimnames = list( rownames(O),
                               colnames(O),
                               colnames(theta))
    )
    for(i in 1:ncol(theta)){
        X[,,i] <- O
    }

    ###
    A <- Y <- X
    A_k_1 <- A_k <- A
    Y_k_1 <- Y_k <- Y
    X_k_1 <- X_k <- X
    a <- as.matrix(rep(1/nrow(theta), nrow(theta)))

    #calculate F matrix
    F = array(0,
              dim = c( ncol(O),
                       ncol(O),
                       ncol(theta)),
              dimnames = list( colnames(O),
                               colnames(O),
                               colnames(theta))
    )
    for(i in 1:ncol(theta)){F[,,i] <- getF(theta[,i],alpha,gamma,a)}
    theta_hat <- colMeans(theta)
    k <- 0
    delta <- 10000
    repeat{
        if(abs(delta)<epsilon||k>max.iter){
            break;
        }else{
            # update
            X_k <- X_k_1
            Y_k <- Y_k_1
            A_k <- A_k_1

            ratio <- NULL
            for(j in 1:ncol(theta)){
                #a <- as.matrix(a.m[j,])
                T_k_j <- getT(j,X_k,theta,O,alpha)
                X_k_1[,,j] <- ((1-alpha)*as.matrix(R[,j])%*%t(a) - A_k[,,j] + gamma*Y_k[,,j] - T_k_j)%*%F[,,j]
                Y_k_1[,,j] <- SVT(((A_k[,,j]/gamma)+X_k_1[,,j]),(beta*theta_hat[j])/gamma)

                A_k_1[,,j] <- A_k[,,j] + gamma*(X_k_1[,,j]-Y_k_1[,,j])
                ratio <- c(ratio, sum( (X_k_1[,,j]-X_k[,,j])^2 )/sum(X_k[,,j]^2))
            }

            if(verbose) writeLines( sprintf("   Ratio ranges from: %f - %f", min(ratio), max(ratio) ) )
            r <- loss(O,X_k,theta,alpha,beta,R)
            if(verbose) writeLines( sprintf("   Loss: part1=%f , part2=%f , part3=%f", r$part1,r$part2,r$part3 ) )
            delta <- max(ratio)
            k <- k+1
        }
    }
    writeLines( paste("Converge in ",k," steps",sep="") )
    X_k
}

#Using nuclear norm to regularize the object function
getF <- function(theta,alpha,gamma,a){
    F_val <- alpha*diag(theta^2)+(1-alpha)*a%*%t(a)+gamma*diag(length(a))
    F_val <- solve(F_val)
    F_val
}

getT <- function(index,X,theta_m,O,alpha){
    X_summary <- 0;
    cell_type_seq <- c(1:ncol(theta_m))
    cell_type_seq <- cell_type_seq[cell_type_seq!=index]

    for(i in cell_type_seq){
        X_summary <- X_summary + X[,,i]%*%diag(theta_m[,i])
    }

    T_val <- alpha*(X_summary-O)%*%diag(theta_m[,index])
    T_val
}

SVT <- function(Mat,t){
    svd <- svd(Mat)
    d <- svd$d
    d <- d - t
    d[d<0] <- 0
    Mat_t <- svd$u %*% diag(d) %*% t(svd$v)
    Mat_t
}
