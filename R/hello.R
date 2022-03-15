X = t(matrix(c(1:12),nrow=3,ncol=4))
theta = matrix(c(0.5,0,1,0.5,1,0),nrow=3,ncol=2)
P = array(0,
              dim = c( nrow(X),
                       ncol(X),
                       ncol(theta)),
              dimnames = list( rownames(X),
                               colnames(X),
                               colnames(theta))
)
R <- matrix(c(2,4,4,6,3,2,5,4),nrow=4,ncol=2)
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
  print(dP1)
  print(dP2)
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
