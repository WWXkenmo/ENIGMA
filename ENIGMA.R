##ENIGMA
##################utile functions#####################
get_proportion <- function(X, ref) {
    cat( date(), "Calculating cell type proportion of bulk samples... \n" )
    gene_id = intersect( rownames(X), rownames(ref) )
    X_m = X[gene_id,]
    ref_m = ref[gene_id,]
    ref_m <- apply(ref_m,2,scale)
    theta <- NULL
    coefVec <- NULL
    for(i in 1:ncol(X_m)){
        Exp <- as.matrix(X_m[,i])
        rownames(Exp) <- rownames(ref_m)
        colnames(Exp) <- colnames(X_m)[i]
        Exp <- scale(Exp)
        
        rlm.o <- rlm(Exp ~ as.matrix(ref_m), maxit = 100)
        coef.v <- summary(rlm.o)$coef[2:(ncol(as.matrix(ref_m)) + 1), 1]
        coefVec <- rbind(coefVec,coef.v)
        coef.v[which(coef.v < 0)] <- 0
        total <- sum(coef.v)
        coef.v <- coef.v/total
        theta <- rbind(theta,coef.v)
    }
    colnames(theta) <- colnames(coefVec) <- colnames(ref_m)
    rownames(theta) <- rownames(coefVec) <- colnames(X_m)
    res <- list()
    res$theta <- theta
    res$coef <- coefVec
    return(res)
}

tsne_plot <- function(tsne,celltype){
    coord <- NULL
    for(i in names(table(celltype))){
        coord <- rbind(coord, colMeans(tsne[celltype == i,]))
    }
    plot(tsne[,1],tsne[,2],col=rainbow(length(table(celltype)))[as.integer(as.factor(celltype))],xlab="tSNE-1",ylab="tSNE-2",main="t-SNE plot")
         text(coord[,1], coord[,2],labels =names(table(celltype)))
}

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


###################################S_mode_correction##################################################
require("sva")
require("purrr")

remove_batch_effect <- function(bulk_eset, ref_eset, varname_main_ct, varnames_sub_ct=NULL, n_pseudo_bulk=1000) {
    cat(date(), "generating pseudo bulk... \n")
    pseudo_bulk_main_ct = generate_pseudo_bulk_from_scRNA( ref_eset, varname_main_ct, n = n_pseudo_bulk )
    colnames(pseudo_bulk_main_ct$mat) = paste(varname_main_ct, seq_len(ncol(pseudo_bulk_main_ct$mat)), sep = "_")
    
    if(is.null(varnames_sub_ct)) pseudo_bulk_subs = NULL
    else {
        pseudo_bulk_subs = lapply( varnames_sub_ct, function(x) {
            sample_ids = subset( pData(ref_eset), eval(parse(text = paste0("!is.na(", x, ")") )) ) %>% rownames
            tmp = generate_pseudo_bulk_from_scRNA( ref_eset[,sample_ids], x, n = n_pseudo_bulk )
            colnames(tmp$mat) = paste(x, seq_len(ncol(tmp$mat)), sep = "_")
            return(tmp)
        } )
    }
    
    cat(date(), "do ComBat... \n")
    M_list = c(
        list( exprs(bulk_eset), pseudo_bulk_main_ct$mat ),
        lapply( pseudo_bulk_subs, function(x) x$mat )
    )
    gene_id = M_list %>% lapply(rownames) %>% Reduce(intersect, .)
    M_mix = M_list %>% lapply(function(x) x[gene_id,] ) %>% do.call(cbind, .)
    M_mix_log2 = log2(M_mix+1)
    
    M_correct = ComBat(
        M_mix_log2[rowSums(M_mix_log2)>0,],
        rep( c("bulk", "pseudo_bulk"), c( ncol(bulk_eset), ncol(M_mix_log2)-ncol(bulk_eset) ) )
    )
    M_correct = 2^M_correct - 1
    
    cat(date(), "restore reference... \n")
    res = c(
        list(
            M_correct[,colnames(bulk_eset)],
            restore_ref( pseudo_bulk_main_ct$frac, M_correct[,colnames(pseudo_bulk_main_ct$mat)] )
        ),
        lapply( pseudo_bulk_subs, function(x) {
            restore_ref( x$frac, M_correct[,colnames(x$mat)] )
        } )
    )
    names(res) = c( "bulk", varname_main_ct, varnames_sub_ct )
    colnames(res[[2]]) <- colnames(pseudo_bulk_main_ct$frac)
    return(res)
    cat(date(), "done. \n")
}


generate_pseudo_bulk_from_scRNA <- function(ref_eset, ct_varname, n=1000) {
    frac_init = table( pData(ref_eset)[,ct_varname] )/nrow(pData(ref_eset))
    
    frac = map2_df(frac_init, frac_init*2, rnorm, n=n) %>% sapply(function(x) x)
    
    # set nagative equal to zero
    frac[frac<0] <- 0
    frac <- frac[rowSums(frac)>0,]
    # normalization
    frac <- t(t(frac) %*% diag(1/rowSums(frac)))
    
    M_star <- NULL
    for(i in 1:nrow(frac)){
        Ma = lapply( pData(ref_eset)[,ct_varname] %>% unique, function(ct) {
            sample_ids = subset( pData(ref_eset), eval(parse( text = paste0(ct_varname, "==\"", ct, "\"") )) ) %>%
                rownames %>%
                sample( 1000*frac[i,ct], replace = TRUE )
            exprs(ref_eset)[,colnames(ref_eset) %in% sample_ids]
        } ) %>%
            do.call(cbind, .)
        
        M_star <- cbind(M_star,rowSums(Ma))
    }
    
    return(list( frac=frac, mat=M_star ))
}


restore_ref <- function(frac, pseudo_bulk) {
    ref_exp <- NULL
    for(i in 1:nrow(pseudo_bulk)){
        coef <- nnls(frac, pseudo_bulk[i,])
        
        # pval = sapply( names(coef$coefficients), function(x) f.robftest(coef, x)$p.value )
        coef.v = coef$x
        # coef.v[which(pval>0.05)] = 0
        
        ref_exp <- rbind(ref_exp,coef.v)
    }
    rownames(ref_exp) <- rownames(pseudo_bulk)
    ref_exp[ref_exp<0] = 0
    
    return(ref_exp)
}


###########################B_mode_correction#######################################
B_mode_batch_effect_remove <- function(X,ref,fra){
    cat("Run B-mode to correct batch effect...")
    X_hat <- ref %*% t(fra)
    X_log <- log2(X+1)
    X_hat_log <- log2(X_hat+1)
    
    #
    cat("\n do Combat...")
    correct <- sva::ComBat(cbind(X_log,X_hat_log),batch=c(rep("1",ncol(X)),rep("2",nrow(fra))))
    X_correct <- 2^correct[,1:ncol(X)]-1
    
    correct <- 2^correct[,(ncol(X)+1):ncol(correct)]-1
    
    #
    #cat("/nRestore the gene specific expression...")
    ##Using nnls
    #gene_ref_exp <- NULL
    #for(i in 1:nrow(X_correct)){
    #    gene_exp <- rlm(fra, correct[i,])$coefficients
    #    gene_ref_exp <- rbind(gene_ref_exp,gene_exp)
    #}
    
    #rownames(gene_ref_exp) <- rownames(ref)
    #colnames(gene_ref_exp) <- colnames(ref)
    #gene_ref_exp[gene_ref_exp<0] <-0
    
    cat("\n Done")
    #res <- list()
    res <- X_correct
    res
}


##################Sub_loss for L2_max Norm###################
sub_loss <- function(X, P_old, theta, alpha,miu,R){
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
        part2 <- part2 + alpha*sum( (P_old[,,i]-ref)^2 )
    }
    
    part3 <- 0
    for(i in 1:ncol(R)){
        # norm <- apply(P_old[,,i],2,norm,"2")
        part3 <- part3 + max( colSums(P_old[,,i]^2) )
    }
    
    res <- list()
    val <- part1+part2+miu*part3
    res$val <- val
    res$part1 <- part1
    res$part2 <- part2/alpha
    res$part3 <- part3
    res
}


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


###############################################################
#Using nuclear norm to regularize the object function
getF <- function(theta,alpha,gamma,a){
  F <- alpha*diag(theta^2)+(1-alpha)*a%*%t(a)+gamma*diag(length(a))
  F <- solve(F)
  F
}

getT <- function(index,X,theta_m,O,alpha){
    X_summary <- 0;
    cell_type_seq <- c(1:ncol(theta_m))
    cell_type_seq <- cell_type_seq[cell_type_seq!=index]

    for(i in cell_type_seq){
        X_summary <- X_summary + X[,,i]%*%diag(theta_m[,i])
    }
	
	T <- alpha*(X_summary-O)%*%diag(theta_m[,index])
	T
}

SVT <- function(Mat,t){
    svd <- svd(Mat)
	d <- svd$d
	d <- d - t
	d[d<0] <- 0
	Mat_t <- svd$u %*% diag(d) %*% t(svd$v)
	Mat_t
}




###initialize the matrix X, Y, A
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

###calculate Loss function
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

