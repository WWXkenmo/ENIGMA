cell_deconvolve <- function(X, theta, R, alpha=0.5, tao_k=0.005,miu=0.5,inner_epilson=0.001,outer_epilson=0.001,max.iter=1,max.iter.exp=1000,verbose=FALSE,infer=FALSE){
    # unify geneid between X and R
    geneid = intersect( rownames(X), rownames(R) )
    X = X[geneid,]
    R = R[geneid,]
    
    # 初始化各细胞类型表达谱矩阵
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
    loss <- sub_loss(X, P_old, theta, alpha, miu, R)
    loss_new <- -1000
    delta <- abs(loss_new-loss$val)
    ###update iteractively
    P_old_new <- P_old
    P_old_pre <- P_old
    iter <- 1
    
    repeat{
        if(delta<outer_epilson||iter>max.iter){
            break;
        }else{
            if(verbose) cat(date(), 'Optimizing cell type specific expression profile... \n')
            iter.exp <- 1
            repeat{
                ratio <- NULL
                ## 初始化更新细胞标签
                ## 更新第i个细胞类型表达谱梯度
                dP <- derive_P2(X, theta,P_old,R,alpha)
                ## 返回第i个细胞类型表达谱（P）
                ## 利用gradient decent对P进行迭代更新
                for(i in 1:ncol(theta)){
                    P_hat <- proximalpoint(P_old[,,i], tao_k,dP[,,i],miu)
                    P_old_new[,,i] <- P_hat
                    
                    ## check更新前后差异
                    ratio <- c(ratio, sum( (P_hat-P_old[,,i])^2 ))
                }
                if(verbose) writeLines( sprintf("   Ratio ranges from: %f - %f", min(ratio), max(ratio) ) )
                if(max(ratio) < inner_epilson||iter.exp >= max.iter.exp){break}else{
                    ## 更新P_old
                    P_old <- P_old_new
                    iter.exp <- iter.exp + 1
                }
            }
        }
        
        ### optimized theta
        ### simply perform robust linear regression model
        if(verbose) cat('Optimizing cell type proportions... \n')
        #theta_new <- NULL
        #C <- cbind(rep(1,ncol(theta)), diag(ncol(theta)))
        #b <- c(1,rep(0,ncol(theta)))
        
        #for(j in 1:ncol(X)){
        #   Exp = as.matrix(X[,j])
        #    x = P_old_new[,j,]
        
        # normalization
        #    Exp <- (Exp)/ sum(Exp)
        #    x = t( t(x)/colSums(x) )
        #    Rinv <- solve(chol( crossprod(x) ));
        #    d <- crossprod(Exp, x)
        #    coef.v <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)$solution
        #    theta_new <- rbind(theta_new,coef.v)
        #}
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
        rownames(theta_new) <- colnames(X)
        }
            ### optimize theta
            ### take the gradient of all theta and running gradient decent
            P_old[P_old<0] <- 0
            #loss_new.obj <- sub_loss(X, P_old, theta_new, alpha,miu,R)
            #loss.obj <- sub_loss(X, P_old_pre, theta, alpha,miu,R)
            #if(verbose) writeLines( sprintf("Total delta_loss: %f, %s", abs(loss_new.obj$val-loss.obj$val), date() ) )
            #if(verbose) writeLines( paste("part1:",loss_new.obj$part1," part2:",loss_new.obj$part2," part3:",loss_new.obj$part3,sep="") )
            #delta <- abs(loss_new.obj$val-loss.obj$val)
            P_old_pre <- P_old
            #theta <- theta_new
            iter <- iter + 1    
        
    }
    # 返回细胞特异表达谱矩阵
    # return(P_old)
    return( list(
        expr_array = P_old,
        theta_new = theta
    ) )
}


