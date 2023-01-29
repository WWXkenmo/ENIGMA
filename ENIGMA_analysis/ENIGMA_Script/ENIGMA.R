##ENIGMA
##################utile functions#####################
####Loss curver plot
#Using nuclear norm to regularize the object function
plotMultiLossCurve <- function(loss_history_list,rlgType = "trace_norm",scale = "log",shape = FALSE,name_prefix,all=TRUE,plotTerm = 1){
    require("ggplot2")
    require("cowplot")
    
    if(rlgType == "trace_norm"){
        df_list <- list()
        for(l in 1:length(loss_history_list)){
            loss_history <- loss_history_list[[l]]
            loss_history_list[[l]] <- loss_history
            df <- data.frame(Iteration = rep(c(1:nrow(loss_history)),(ncol(loss_history)+1)),LossType = 
                                 c(rep("Distance to observed mixture",nrow(loss_history)),
                                   rep("Distance to reference",nrow(loss_history)),
                                   rep("Rank regularization",nrow(loss_history)),
                                   rep("Over All",nrow(loss_history))),Loss = c(as.numeric(as.matrix(loss_history)),rowSums(loss_history)))
            df_list[[l]] <- df
        }
        names(df_list) <- names(loss_history_list)
        name_vec <- NULL
        for(i in 1:length(loss_history_list)) name_vec <- c(name_vec,rep(names(loss_history_list)[i],nrow(loss_history_list[[i]])))
        
        df1 <- df2 <- df3 <- df4 <- NULL
        for(l in 1:length(loss_history_list)){
            df1 <- rbind(df1, subset(df_list[[l]],df_list[[l]]$LossType %in% "Distance to observed mixture"))
            df2 <- rbind(df2, subset(df_list[[l]],df_list[[l]]$LossType %in% "Distance to reference"))
            df3 <- rbind(df3, subset(df_list[[l]],df_list[[l]]$LossType %in% "Rank regularization"))
            df4 <- rbind(df4, subset(df_list[[l]],df_list[[l]]$LossType %in% "Over All"))
        }
        id <- c(colnames(df1),"Method")
        df1 <- as.data.frame(cbind(df1, name_vec))
        df2 <- as.data.frame(cbind(df2, name_vec))
        df3 <- as.data.frame(cbind(df3, name_vec))
        df4 <- as.data.frame(cbind(df4, name_vec))
        colnames(df1) <- colnames(df2) <- colnames(df3) <- colnames(df4) <- id
        
        
        if(scale == "log"){
            df1$Loss <- log(df1$Loss+1)
            df2$Loss <- log(df2$Loss+1)
            df3$Loss <- log(df3$Loss+1)
            df4$Loss <- log(df4$Loss+1)
        }
        
        if(!shape){
            g1 <- ggplot(data=df1, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Distance to observed mixture',colour = name_prefix)
            g2 <- ggplot(data=df2, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Distance to reference',colour = name_prefix)
            g3 <- ggplot(data=df3, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Rank regularization',colour = name_prefix)
            g4 <- ggplot(data=df4, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Overall',colour = name_prefix)
        }else{
            g1 <- ggplot(data=df1, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Distance to observed mixture',colour = name_prefix)
            g2 <- ggplot(data=df2, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Distance to reference',colour = name_prefix)
            g3 <- ggplot(data=df3, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Rank regularization',colour = name_prefix)
            g4 <- ggplot(data=df4, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Overall',colour = name_prefix)
        }
        if(all){
            if(scale != "log"){p1 <- plot_grid(g1+labs(y = 'Loss'),g2+labs(y = 'Loss'),g3+labs(y = 'Loss'),g4+labs(y = 'Loss'),nrow=2);print(p1)}
            if(scale == "log"){p1 <- plot_grid(g1,g2,g3,g4,nrow=2);print(p1)}
        }else{
            if(plotTerm == 1){
                if(scale != "log") print(g1+labs(y = 'Loss'))
                if(scale == "log") print(g1)}
            if(plotTerm == 2){
                if(scale != "log") print(g2+labs(y = 'Loss'))
                if(scale == "log") print(g2)}
            if(plotTerm == 3){
                if(scale != "log") print(g3+labs(y = 'Loss'))
                if(scale == "log") print(g3)}
            if(plotTerm == 4){
                if(scale != "log") print(g4+labs(y = 'Loss'))
                if(scale == "log") print(g4)}
        }
    }else if(rlgType == "L2_max_norm"){
        df_list <- list()
        for(l in 1:length(loss_history_list)){
            loss_history <- loss_history_list[[l]]
            loss_history_list[[l]] <- loss_history
            df <- data.frame(Iteration = rep(c(1:nrow(loss_history)),(ncol(loss_history)+1)),LossType = 
                                 c(rep("Distance to observed mixture",nrow(loss_history)),
                                   rep("Distance to reference",nrow(loss_history)),
                                   rep("Maximum L2 norm regularization",nrow(loss_history)),
                                   rep("Over All",nrow(loss_history))),Loss = c(as.numeric(as.matrix(loss_history)),rowSums(loss_history)))
            df_list[[l]] <- df
        }
        names(df_list) <- names(loss_history_list)
        name_vec <- NULL
        for(i in 1:length(loss_history_list)) name_vec <- c(name_vec,rep(names(loss_history_list)[i],nrow(loss_history_list[[i]])))
        
        df1 <- df2 <- df3 <- df4 <- NULL
        for(l in 1:length(loss_history_list)){
            df1 <- rbind(df1, subset(df_list[[l]],df_list[[l]]$LossType %in% "Distance to observed mixture"))
            df2 <- rbind(df2, subset(df_list[[l]],df_list[[l]]$LossType %in% "Distance to reference"))
            df3 <- rbind(df3, subset(df_list[[l]],df_list[[l]]$LossType %in% "Maximum L2 norm regularization"))
            df4 <- rbind(df4, subset(df_list[[l]],df_list[[l]]$LossType %in% "Over All"))
        }
        id <- c(colnames(df1),"Method")
        df1 <- as.data.frame(cbind(df1, name_vec))
        df2 <- as.data.frame(cbind(df2, name_vec))
        df3 <- as.data.frame(cbind(df3, name_vec))
        df4 <- as.data.frame(cbind(df4, name_vec))
        colnames(df1) <- colnames(df2) <- colnames(df3) <- colnames(df4) <- id
        
        
        if(scale == "log"){
            df1$Loss <- log(df1$Loss+1)
            df2$Loss <- log(df2$Loss+1)
            df3$Loss <- log(df3$Loss+1)
            df4$Loss <- log(df4$Loss+1)
        }
        
        if(!shape){
            g1 <- ggplot(data=df1, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Distance to observed mixture',colour = name_prefix)
            g2 <- ggplot(data=df2, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Distance to reference',colour = name_prefix)
            g3 <- ggplot(data=df3, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Maximum L2 norm regularization',colour = name_prefix)
            g4 <- ggplot(data=df4, aes(x=Iteration, y=Loss, group=Method, color=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                              y = 'log Loss',title = 'Overall',colour = name_prefix)
        }else{
            g1 <- ggplot(data=df1, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Distance to observed mixture',colour = name_prefix)
            g2 <- ggplot(data=df2, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Distance to reference',colour = name_prefix)
            g3 <- ggplot(data=df3, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Maximum L2 norm regularization',colour = name_prefix)
            g4 <- ggplot(data=df4, aes(x=Iteration, y=Loss, group=Method, color=Method,shape=Method)) + geom_line() + geom_point()+ theme_minimal() + labs(x = 'iteration',
                                                                                                                                                           y = 'log Loss',title = 'Overall',colour = name_prefix)
        }
        if(all){
            if(scale != "log"){p1 <- plot_grid(g1+labs(y = 'Loss'),g2+labs(y = 'Loss'),g3+labs(y = 'Loss'),g4+labs(y = 'Loss'),nrow=2);print(p1)}
            if(scale == "log"){p1 <- plot_grid(g1,g2,g3,g4,nrow=2);print(p1)}
        }else{
            if(plotTerm == 1){
                if(scale != "log") print(g1+labs(y = 'Loss'))
                if(scale == "log") print(g1)}
            if(plotTerm == 2){
                if(scale != "log") print(g2+labs(y = 'Loss'))
                if(scale == "log") print(g2)}
            if(plotTerm == 3){
                if(scale != "log") print(g3+labs(y = 'Loss'))
                if(scale == "log") print(g3)}
            if(plotTerm == 4){
                if(scale != "log") print(g4+labs(y = 'Loss'))
                if(scale == "log") print(g4)}
        }
    }
}


##########################


get_proportion <- function(X, ref) {
    require("MASS")
	
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

derive_P <- function(X, theta, P_old,R,alpha){
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
    
    w1 <- alpha
    w2 <- 1-w1
    
    dP <- dP1*as.numeric(w1) + dP2*as.numeric(w2)
    return(dP)
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
    
    w1 <- alpha
    w2 <- 1-w1
    
    dP <- dP1*as.numeric(w1) + dP2*as.numeric(w2)
    return(dP)
}

getX <- function(O,theta,R,A,Y,alpha,gamma){
    ##O: bulk data (m*n)
    ##theta: cell type fractions
    ##A: constraint variable (m*n*ct)
    ##Y: auxiliary variable (m*n*ct)
    ##R: reference (m*ct)
    ## reshape the input data
    theta_m <- A_m <- Y_m <- NULL
    for(i in 1:ncol(theta)){
        theta_m = cbind(theta_m, diag(theta[,i]))
        A_m = cbind(A_m, A[,,i])
        Y_m = cbind(Y_m, Y[,,i])
    }
    a <- matrix(0,nrow=(ncol(theta)*ncol(O)),ncol=ncol(theta))
    for(i in 1:ncol(theta)){
        a[((i-1)*ncol(O)+1):(i*ncol(O)),i] <- rep((1/ncol(O)),ncol(O))
    }
    ## update X
    F = solve(alpha*t(theta_m)%*%theta_m+(1-alpha)*a%*%t(a)+gamma*diag(nrow(a)))
    X = (alpha*O%*%theta_m+(1-alpha)*R%*%t(a)-A_m+gamma*Y_m)%*%F
    
    ## split X into CSE blocks
    X_k = array(0,
                dim = c( nrow(O),
                         ncol(O),
                         ncol(theta)),
                dimnames = list( rownames(O),
                                 colnames(O),
                                 colnames(theta))
    )
    
    for(i in 1:ncol(theta)){
        X_k[,,i] <- X[,((i-1)*ncol(O)+1):(i*ncol(O))]
    }
    
    ##return
    X_k
}

SVT <- function(Mat,t){
	require("corpcor")
	
    svd <- fast.svd(Mat)
	d <- svd$d
	d <- d - t
	d[d<0] <- 0
	Mat_t <- svd$u %*% diag(d) %*% t(svd$v)
	Mat_t
}

SVT_RM <- function(Mat){
    require("corpcor")
	#Estimate beta
	beta <- ncol(Mat)/nrow(Mat)

	#optimal threshold for hard thresholding
	lambda <- sqrt(2*(beta+1)+8*beta/((beta+1)+sqrt(beta^2+14*beta+1)))
	ratio <- (1+sqrt(beta))/lambda #(bulk edges)/(hard thresholding optimal)

	#practical value for soft thresholding
	omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43
	omega = omega * ratio / sqrt(nrow(Mat))

	#Calculate Singular Values
	svd <- fast.svd(Mat)
	d <- svd$d

	#Estimate t
	t <- median(d)*omega
	d <- d - t
	d[d<10^-6] <- 0
	Mat_t <- svd$u %*% diag(d) %*% t(svd$v)
	Mat_t
}

SVT_RM_value <- function(Mat){
    require("corpcor")
	#Estimate beta
	beta <- ncol(Mat)/nrow(Mat)

	#optimal threshold for hard thresholding
	lambda <- sqrt(2*(beta+1)+8*beta/((beta+1)+sqrt(beta^2+14*beta+1)))
	ratio <- (1+sqrt(beta))/lambda #(bulk edges)/(hard thresholding optimal)

	#practical value for soft thresholding
	omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43
	omega = omega * ratio / sqrt(nrow(Mat))

	#Calculate Singular Values
	svd <- fast.svd(Mat)
	d <- svd$d

	#Estimate t
	t <- median(d)*omega
	t
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
            repeat{
                ratio <- NULL
                dP <- derive_P2(X, theta,P_old,R,alpha)
                for(i in 1:ncol(theta)){
                    P_hat <- proximalpoint(P_old[,,i], tao_k,dP[,,i],beta*10^5)
					P_old_new[,,i] <- P_hat
                    
                    ratio <- c(ratio, sum( (P_hat-P_old[,,i])^2 ))
                }
                if(verbose) writeLines( sprintf("   Ratio ranges from: %f - %f", min(ratio), max(ratio) ) )
                r <- sub_loss(X, P_old, theta, alpha,beta,R)
                if(loss_his) loss<- rbind(loss,c(r$part1,r$part2,r$part3))
                if(max(ratio) < epsilon||iter.exp >= max.iter){break}else{
                    P_old <- P_old_new
                    iter.exp <- iter.exp + 1
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
				scale_x <- function(x){
				if(var(x)==0){x <- x - mean(x)}else{x <- scale(x)}
				x
				}
                exp.scale <- t(apply(exp,1,scale_x))
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
		
		if(Norm.method == "quantile"){
		     for(k in 1:dim(X_k_m)[3]){
		        X_k_norm[,,k] <- normalize.quantiles(X_k_norm[,,k])
		     }
		}
		
    }
	writeLines( paste("Done! \nConverge in ",iter.exp," steps",sep="") )
    # return cell type specific gene expression matrix
	# return(P_old)
	res = list(
	X_k = P_old,
	theta = theta,
	loss_history = loss
	)
	if(Normalize) res$X_k_norm = X_k_norm
	return( res )
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





###initialize the matrix X, Y, A
cell_deconvolve_trace <- function(O, theta, R, alpha=0.5,beta=5,tao_k=1,gamma=NULL,epsilon=NULL,max.iter=100,solver = "admm",verbose=FALSE,X_int=NULL,loss_his=TRUE,Normalize=TRUE,Norm.method = "PC",pre.process = "log",pos=TRUE,infer=FALSE){
    solver_list = c("admm","admm_fast","adaptive_admm","proximalpoint")
	msg = paste0("should be one of ", paste(solver_list, collapse = ", "), 
	".")
	if (!solver %in% solver_list) {
	stop("'solver' ", msg)
	}
	
    if(is.null(X_int)){
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
    }else{
        X <- X_int
    }
    
    ###
	A <- Y <- X
    A_k_1 <- A_k <- A
    Y_k_1 <- Y_k <- Y
    X_k_1 <- X_k <- X
    a <- as.matrix(rep(1/nrow(theta), nrow(theta)))
    
	if(is.null(epsilon)){
	   z1  <- format(median(abs(O)), scientific=TRUE, digit=7)
       z2  <- substring(z1, 1, 8)
       Power <- log( as.numeric( z1 ) / as.numeric( z2 ), 10 ) 
       epsilon <- 10^Power*(10^-4)
	}
	
	if(is.null(gamma)){
	   if(SVT_RM_value(O) < 1) gamma <- 2 
	   if(SVT_RM_value(O) > 1) gamma <- 0.5
	}
	
	#calculate F matrix 
	if(solver == "admm_fast"){
    F = array(0,
                  dim = c( ncol(O),
                           ncol(O),
                           ncol(theta)),
                  dimnames = list( colnames(O),
                                   colnames(O),
                                   colnames(theta))
    )
	for(i in 1:ncol(theta)){F[,,i] <- getF(theta[,i],alpha,gamma,a)}
	}
	
    theta_hat <- colMeans(theta)
    k <- 1
    delta <- 10000
    loss <- NULL
    if(verbose) cat(date(), 'Optimizing cell type specific expression profile... \n')
	if(verbose){if(solver == "adaptive_admm"){
	cat(paste('alpha: ',alpha, ' \n', 'gamma: ',gamma ,' \n','epsilon: ',epsilon,' \n',sep=""))}
	if(solver == "proximalpoint"){
	cat(paste('alpha: ',alpha, ' \n', 'step size: ',tao_k ,' \n','beta: ',beta ,' \n','gamma: ',gamma ,' \n','epsilon: ',epsilon,' \n',sep=""))}
	if(solver == "admm" || solver == "admm_fast"){
	cat(paste('alpha: ',alpha, ' \n', 'beta: ',beta ,' \n','gamma: ',gamma ,' \n','epsilon: ',epsilon,' \n',sep=""))
	}}
	
	writeLines(paste("Using ",solver," solver...",sep=""))
    repeat{
        if(abs(delta)<epsilon||k>max.iter){
            break;
        }else{
		
		
		###################################
		#using admm solver
		if(solver == "admm"){
		# update
		X_k <- X_k_1
		Y_k <- Y_k_1
		A_k <- A_k_1
		
		ratio <- NULL
		##update X
		updated_X <- getX(O,theta,R,A_k,Y_k,alpha,gamma)
		for(j in 1:ncol(theta)){
			#a <- as.matrix(a.m[j,])
			X_k_1[,,j] <- updated_X[,,j]
			Y_k_1[,,j] <- SVT(((A_k[,,j]/gamma)+X_k_1[,,j]),(beta*theta_hat[j])/gamma)
			
			A_k_1[,,j] <- A_k[,,j] + gamma*(X_k_1[,,j]-Y_k_1[,,j])
			ratio <- c(ratio, sum( (X_k_1[,,j]-X_k[,,j])^2 )/(nrow(X[,,j])*ncol(X[,,j])))
		}
		r1 <- loss(O,X_k,theta,alpha,beta,R) #calculating raw loss
		if(verbose){
		#print matrix ratio distance or absolute distance
		print <- paste("CSE inference step ",k," \n",sep="")
		for(c in 1:(length(ratio)-1)){
		print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
		}
		print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
		writeLines(print)
		}
		if(loss_his) loss<- rbind(loss,c(r1$part1,r1$part2,r1$part3))
		if(verbose) writeLines( sprintf("   Loss: part1=%f , part2=%f , part3=%f", r1$part1,r1$part2,r1$part3 ) )
		delta <- max(ratio)
		k <- k+1
		}
		
		###################################
		#using adaptive admm solver
		if(solver == "adaptive_admm"){
		# update
		X_k <- X_k_1
		Y_k <- Y_k_1
		A_k <- A_k_1

		ratio <- NULL
		##update X
		updated_X <- getX(O,theta,R,A_k,Y_k,alpha,gamma)
		for(j in 1:ncol(theta)){
			#a <- as.matrix(a.m[j,])
			X_k_1[,,j] <- updated_X[,,j]
			Y_k_1[,,j] <- SVT_RM(((A_k[,,j]/gamma)+X_k_1[,,j]))
			
			A_k_1[,,j] <- A_k[,,j] + gamma*(X_k_1[,,j]-Y_k_1[,,j])
			ratio <- c(ratio, sum( (X_k_1[,,j]-X_k[,,j])^2 )/(nrow(X[,,j])*ncol(X[,,j])))
		}
		r <- loss(O,X_k,theta,alpha,1,R)
		if(verbose){
		#print matrix ratio distance or absolute distance
		print <- paste("CSE inference step ",k," \n",sep="")
		for(c in 1:(length(ratio)-1)){
		  print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
		}
		print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
		writeLines(print)
		}
		if(loss_his) loss<- rbind(loss,c(r$part1,r$part2,r$part3))
		if(verbose) writeLines( sprintf("   Loss: part1=%f , part2=%f , part3=%f", r$part1,r$part2,r$part3 ) )
		delta <- max(ratio)
		k <- k+1
		}
		
		###################################
		#using admm_fast solver
		if(solver == "admm_fast"){
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
			ratio <- c(ratio, sum( (X_k_1[,,j]-X_k[,,j])^2 )/(nrow(X[,,j])*ncol(X[,,j])))
		}

		r <- loss(O,X_k,theta,alpha,beta,R)
		if(verbose){
		#print matrix ratio distance or absolute distance
		print <- paste("CSE inference step ",k," \n",sep="")
		for(c in 1:(length(ratio)-1)){
		  print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
		}
		print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
		writeLines(print)
		}
		if(loss_his) loss<- rbind(loss,c(r$part1,r$part2,r$part3))
		if(verbose) writeLines( sprintf("   Loss: part1=%f , part2=%f , part3=%f", r$part1,r$part2,r$part3 ) )
		delta <- max(ratio)
		k <- k+1
		}
		
		###################################
		#using admm_fast solver
		if(solver == "proximalpoint"){
		X_k <- X_k_1
		Y_k <- Y_k_1
		A_k <- A_k_1
		
		ratio <- NULL
        dP <- derive_P(O, theta,X_k,R,alpha)
        for(j in 1:ncol(theta)){
            X_i <- X_k[,,j]- tao_k*dP[,,j]
            X_i <- SVT(X_i,tao_k*theta_hat[j]*beta)
            X_k_1[,,j] <- X_i
            
            ratio <- c(ratio, sum( (X_k_1[,,j]-X_k[,,j])^2 )/(nrow(X[,,j])*ncol(X[,,j])))
        }
		r <- loss(O,X_k,theta,alpha,beta,R)
		if(verbose){
		#print matrix ratio distance or absolute distance
		print <- paste("CSE inference step ",k," \n",sep="")
		for(c in 1:(length(ratio)-1)){
		  print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
		}
		print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
		writeLines(print)
		}
		if(loss_his) loss<- rbind(loss,c(r$part1,r$part2,r$part3))
		if(verbose) writeLines( sprintf("   Loss: part1=%f , part2=%f , part3=%f", r$part1,r$part2,r$part3 ) )
		delta <- max(ratio)
		k <- k+1
		}
		
		}
		}
		steps <- k
		##########
		#Doing PC or Cell Fractions Normalization
		if(pos) X_k[X_k<0] <- 0
		if(Normalize){
		writeLines("Normalization...")
		if(pre.process == "log") X_k_m <- 2^X_k - 1
	    if(pre.process == "sqrt") X_k_m <- X_k^2
		if(pre.process == "none") X_k_m <- X_k
		pp_list = c("log","sqrt","none")
		msg = paste0("should be one of ", paste(pp_list, collapse = ", "), 
			".")
		if (!pre.process %in% pp_list) {
			stop("'preprocess method' ", msg)
		}
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
		
		if(Norm.method == "total"){
		     for(k in 1:dim(X_k_m)[3]){
		        X_k_norm[,,k] <- X_k_norm[,,k] %*% diag(10000/colSums(X_k_norm[,,k]))
			 }
		}
		
		if(Norm.method == "quantile"){
		     for(k in 1:dim(X_k_m)[3]){
		        X_k_norm[,,k] <- normalize.quantiles(X_k_norm[,,k])
			 }
		}
        }
    if(verbose) cat('Optimizing cell type proportions... \n')
    if(infer){
        if(Normalize) CSE = X_k_norm
        if(!Normalize) CSE = X_k_m
        theta_new <- NULL
        for(j in 1:ncol(O)){
            Exp <- as.matrix(O[,j])
            rownames(Exp) <- rownames(CSE[,j,])
            colnames(Exp) <- colnames(O)[j]
            x <- CSE[,j,]
            x <- apply(x,2,scale)
            lm.o <- rlm(Exp ~ as.matrix(x),maxit=150)
            coef.v <- lm.o$coefficients[-1]
            coef.v[which(coef.v < 0)] <- 0
            total <- sum(coef.v)
            coef.v <- coef.v/total
            theta_new <- rbind(theta_new,coef.v)
        }
        colnames(theta_new) <- colnames(theta)
        rownames(theta_new) <- colnames(O);theta <- theta_new
    }
	
    writeLines( paste("Done! \nConverge in ",steps," steps",sep="") )
    res <- list()
    res$X_k <- X_k
    res$loss_history <- loss
	res$theta <- theta
    if(Normalize) res$X_k_norm <- X_k_norm
    res
}

###calculate Loss function
loss <- function(X, P_old, theta, alpha,beta,R){
    require("corpcor")
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
        ref <- R[,i]
        # part2 <- part2 + alpha*norm((P_old[,,i]-ref),"F")^2
        part2 <- part2 + sum( (rowMeans(P_old[,,i])-ref)^2 )
    }
    
    part3 <- 0
    for(i in 1:ncol(R)){
        part3 <- part3 + colMeans(theta)[i]*sum(fast.svd(P_old[,,i])$d)
    }
    
    res <- NULL
    res$part1 <- part1*(alpha/2)
    res$part2 <- part2*(1-alpha)*(1/2)
    res$part3 <- beta*(part3)
    res
}
								  
############### ENIGMA code revision #############################

cell_deconvolve <- function(X, theta, R, alpha=0.5, tao_k=0.005,beta=0.5,epsilon=0.001,max.iter=1000,verbose=FALSE,infer=FALSE,loss_his = TRUE,pos=TRUE,pre.process="log",Normalize=TRUE,Norm.method = "PC"){
  # unify geneid between X and R
  geneid = intersect( rownames(X), rownames(R) )
  X = X[geneid,]
  R = R[geneid,]
  
  ## renormalization
	geneID <- rownames(X)
	sampleID <- colnames(X)
	ctID <- colnames(R)
	X <- X %*% diag(10^5/colSums(X))
	R <- R %*% diag(10^5/colSums(R))
	rownames(X) <- rownames(R) <- geneID
	colnames(X) <- sampleID
	colnames(R) <- ctID
	
	if(pre.process == "log"){
	 X <- log2(X+1)
	 R <- log2(R+1)
	}
	if(pre.process == "sqrt"){
	 X <- sqrt(X)
	 R <- sqrt(R)
	}
	
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
  mask_entry <- matrix(1,nrow = nrow(X), ncol = ncol(X)); mask_entry[X==0] <- 0
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
  repeat{
    ratio <- NULL
    dP <- derive_P2(X, theta,P_old,R,alpha,mask_entry)
    for(i in 1:ncol(theta)){
      P_hat <- proximalpoint(P_old[,,i], tao_k,dP[,,i],beta*10^5)
      P_old_new[,,i] <- P_hat

      ratio <- c(ratio, sum( (P_hat-P_old[,,i])^2 ))
    }
    if(verbose) writeLines( sprintf("   Ratio ranges from: %f - %f", min(ratio), max(ratio) ) )
    r <- sub_loss(X, P_old, theta, alpha,beta,R)
    if(loss_his) loss<- rbind(loss,c(r$part1,r$part2,r$part3))
    if(max(ratio) < epsilon||iter.exp >= max.iter){break}else{
      P_old <- P_old_new
      iter.exp <- iter.exp + 1
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
    loss_history = loss
  )
  if(Normalize) res$X_k_norm = X_k_norm
  return( res )
}


derive_P2 <- function(X, theta, P_old,R,alpha,mask_entry){
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
	dP1[,,cell_type_index] <- dP1[,,cell_type_index]*mask_entry
	dP2[,,cell_type_index] <- dP2[,,cell_type_index]*mask_entry
  }
  dP1 = dP1 / sqrt( sum( dP1^2 ) ) * 1e5
  dP2 = dP2 / sqrt( sum( dP2^2 ) ) * 1e5

  w1 <- alpha
  w2 <- 1-w1

  dP <- dP1*as.numeric(w1) + dP2*as.numeric(w2)
  return(dP)
}


cell_deconvolve_trace <- function(O, theta, R, alpha=0.5,beta=1,tao_k=2,gamma=1,epsilon=0.001,max.iter=100,solver = "proximalpoint",verbose=FALSE,X_int=NULL,loss_his=TRUE,Normalize=TRUE,Norm.method = "PC",pre.process = "log",pos=TRUE,infer=FALSE,epsilon_ks = 0.01,print_loss=FALSE,max_ks = 1){
    solver_list = c("admm","admm_fast","proximalpoint")
	msg = paste0("should be one of ", paste(solver_list, collapse = ", "), 
	".")
	if (!solver %in% solver_list) {
	stop("'solver' ", msg)
	}
	
	## renormalization
	geneID <- rownames(O)
	sampleID <- colnames(O)
	ctID <- colnames(R)
	O <- O %*% diag(10^5/colSums(O))
	R <- R %*% diag(10^5/colSums(R))
	rownames(O) <- rownames(R) <- geneID
	colnames(O) <- sampleID
	colnames(R) <- ctID
	
	if(pre.process == "log"){
	 O <- log2(O+1)
	 R <- log2(R+1)
	}
	if(pre.process == "sqrt"){
	 O <- sqrt(O)
	 R <- sqrt(R)
	}
	
	
	rm(geneID,sampleID,ctID);gc()
	
	## ref kappa score
	ref_kappa <- kappa(R)
	
    if(is.null(X_int)){
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
    }else{
        X <- X_int
    }
    
    ###
	mask_entry <- matrix(1,nrow = nrow(O), ncol = ncol(O)); mask_entry[O==0] <- 0
	A <- Y <- X
    A_k_1 <- A_k <- A
    Y_k_1 <- Y_k <- Y
    X_k_1 <- X_k <- X
    a <- as.matrix(rep(1/nrow(theta), nrow(theta)))
    
	
	#calculate F matrix 
	if(solver == "admm_fast"){
    F = array(0,
                  dim = c( ncol(O),
                           ncol(O),
                           ncol(theta)),
                  dimnames = list( colnames(O),
                                   colnames(O),
                                   colnames(theta))
    )
	for(i in 1:ncol(theta)){F[,,i] <- getF(theta[,i],alpha,gamma,a)}
	}
	
    theta_hat <- colMeans(theta)
    k <- 1
    delta <- delta_ks <- ks_new <- ks <- 10000
    loss <- NULL
    if(verbose) cat(date(), 'Optimizing cell type specific expression profile... \n')
	if(verbose){if(solver == "adaptive_admm"){
	cat(paste('alpha: ',alpha, ' \n', 'gamma: ',gamma ,' \n','epsilon: ',epsilon,' \n',sep=""))}
	if(solver == "proximalpoint"){
	cat(paste('alpha: ',alpha, ' \n', 'step size: ',tao_k ,' \n','beta: ',beta ,' \n','gamma: ',gamma ,' \n','epsilon: ',epsilon,' \n',sep=""))}
	if(solver == "admm" || solver == "admm_fast"){
	cat(paste('alpha: ',alpha, ' \n', 'beta: ',beta ,' \n','gamma: ',gamma ,' \n','epsilon: ',epsilon,' \n',sep=""))
	}}
	
	writeLines(paste("Using ",solver," solver...",sep=""))
    repeat{
	    cond1 <- abs(delta)<epsilon||k>max.iter
		cond2 <- abs(delta_ks) < epsilon_ks || ks_new < max_ks
        if(cond1&cond2){
            break;
        }else{
		
		
		###################################
		#using admm solver
		if(solver == "admm"){
		# update
		X_k <- X_k_1
		Y_k <- Y_k_1
		A_k <- A_k_1
		ks <- ks_new
		
		ratio <- NULL
		##update X
		updated_X <- getX(O,theta,R,A_k,Y_k,alpha,gamma)
		for(j in 1:ncol(theta)){
			#a <- as.matrix(a.m[j,])
			X_k_1[,,j] <- updated_X[,,j]*mask_entry
			Y_k_1[,,j] <- SVT(((A_k[,,j]/gamma)+X_k_1[,,j]),(beta*theta_hat[j])/gamma)*mask_entry
			
			A_k_1[,,j] <- A_k[,,j] + gamma*(X_k_1[,,j]-Y_k_1[,,j])
			ratio <- c(ratio, sum( (X_k_1[,,j]-X_k[,,j])^2 )/(nrow(X[,,j])*ncol(X[,,j])))
		}
		if(print_loss) r1 <- loss(O,X_k,theta,alpha,beta,R) #calculating raw loss
		if(verbose){
		#print matrix ratio distance or absolute distance
		print <- paste("CSE inference step ",k," \n",sep="")
		for(c in 1:(length(ratio)-1)){
		print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
		}
		print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
		writeLines(print)
		}
		if(loss_his&print_loss) loss<- rbind(loss,c(r1$part1,r1$part2,r1$part3))
		ks_new = abs(KappaScore(X_k) - ref_kappa)
		kss = format(ks_new, scientific=TRUE, digit=4)
		if(verbose) writeLines( paste(" Kappa Score: ",kss,sep="" ))
		delta <- max(ratio)
		delta_ks <- ks_new - ks
		k <- k+1
		}
		
		
		###################################
		#using admm_fast solver
		if(solver == "admm_fast"){
		X_k <- X_k_1
		Y_k <- Y_k_1
		A_k <- A_k_1
		ks <- ks_new
		
	    ratio <- NULL
		for(j in 1:ncol(theta)){
			#a <- as.matrix(a.m[j,])
			T_k_j <- getT(j,X_k,theta,O,alpha)
			X_k_1[,,j] <- (((1-alpha)*as.matrix(R[,j])%*%t(a) - A_k[,,j] + gamma*Y_k[,,j] - T_k_j)%*%F[,,j])*mask_entry
			Y_k_1[,,j] <- SVT(((A_k[,,j]/gamma)+X_k_1[,,j]),(beta*theta_hat[j])/gamma)*mask_entry

			A_k_1[,,j] <- A_k[,,j] + gamma*(X_k_1[,,j]-Y_k_1[,,j])
			ratio <- c(ratio, sum( (X_k_1[,,j]-X_k[,,j])^2 )/(nrow(X[,,j])*ncol(X[,,j])))
		}

		if(print_loss) r1 <- loss(O,X_k,theta,alpha,beta,R)
		if(verbose){
		#print matrix ratio distance or absolute distance
		print <- paste("CSE inference step ",k," \n",sep="")
		for(c in 1:(length(ratio)-1)){
		  print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
		}
		print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
		writeLines(print)
		}
		if(loss_his&print_loss) loss<- rbind(loss,c(r1$part1,r1$part2,r1$part3))
		ks_new = abs(KappaScore(X_k) - ref_kappa)
		kss = format(ks_new, scientific=TRUE, digit=4)
		if(verbose) writeLines( paste(" Kappa Score: ",kss,sep="" ))
		delta <- max(ratio)
		delta_ks <- ks_new - ks
		k <- k+1
		}
		
x = array(0,dim = c(10000,500,15))
for(i in 1:15){
  x[,,i] <- matrix(rnorm(10000*500),ncol=500)
}

reslist = foreach(
  i = 1:10
) %dopar% {
    svdd <- svd(x[,,i])
    res <- list(u = svdd$u)
  }

for(i in 1:15){
  svdd <- svd(x[,,i])
}



		###################################
		#using admm_fast solver
		if(solver == "proximalpoint"){
    ## setting multiple core
    if(ncore>1){
      doParallel::registerDoParallel(cl = ncore)
      foreach::getDoParRegistered()
      foreach::getDoParWorkers()
    }

		X_k <- X_k_1
		Y_k <- Y_k_1
		A_k <- A_k_1
		ks <- ks_new
		
		ratio <- NULL
        dP <- derive_P(O, theta,X_k,R,alpha)
    if(ncore==1){
        for(j in 1:ncol(theta)){
            X_i <- X_k[,,j]- tao_k*dP[,,j]
            X_i <- SVT(X_i,tao_k*theta_hat[j]*beta)
            X_k_1[,,j] <- X_i*mask_entry
            
            ratio <- c(ratio, sum( (X_k_1[,,j]-X_k[,,j])^2 )/(nrow(X[,,j])*ncol(X[,,j])))
        }
    }else{
        reslist <- foreach(
        j = 1:ncol(theta)
        ) %dopar% {
            X_i <- X_k[,,j]- tao_k*dP[,,j]
            X_i <- SVT(X_i,tao_k*theta_hat[j]*beta)
            X_k_1[,,j] <- X_i*mask_entry

            rr <- sum( (X_k_1[,,j]-X_k[,,j])^2 )/(nrow(X[,,j])*ncol(X[,,j]))
            list(paste(colnames(theta)[j],"_X_k_1",sep="") = X_k_1[,,j],paste(colnames(theta)[j],"_rr",sep="") = rr)
        }

    }
		if(print_loss) r1 <- loss(O,X_k,theta,alpha,beta,R)
		if(verbose){
		#print matrix ratio distance or absolute distance
		print <- paste("CSE inference step ",k," \n",sep="")
		for(c in 1:(length(ratio)-1)){
		  print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
		}
		print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
		writeLines(print)
		}
		if(loss_his&print_loss) loss<- rbind(loss,c(r1$part1,r1$part2,r1$part3))
		ks_new = abs(KappaScore(X_k) - ref_kappa)
		kss = format(ks_new, scientific=TRUE, digit=4)
		if(verbose) writeLines( paste(" Kappa Score: ",kss,sep="" ))
		delta <- max(ratio)
		delta_ks <- ks_new - ks
		k <- k+1
		}
		
		}
		}
		steps <- k
		##########
		#Doing PC or Cell Fractions Normalization
		if(pos) X_k[X_k<0] <- 0
		if(Normalize){
		writeLines("Normalization...")
		if(pre.process == "log") X_k_m <- 2^X_k - 1
	    if(pre.process == "sqrt") X_k_m <- X_k^2
		if(pre.process == "none") X_k_m <- X_k
		pp_list = c("log","sqrt","none")
		msg = paste0("should be one of ", paste(pp_list, collapse = ", "), 
			".")
		if (!pre.process %in% pp_list) {
			stop("'preprocess method' ", msg)
		}
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
		
		if(Norm.method == "total"){
		     for(k in 1:dim(X_k_m)[3]){
		        X_k_norm[,,k] <- X_k_norm[,,k] %*% diag(10000/colSums(X_k_norm[,,k]))
			 }
		}
		
		if(Norm.method == "quantile"){
		     for(k in 1:dim(X_k_m)[3]){
		        X_k_norm[,,k] <- normalize.quantiles(X_k_norm[,,k])
			 }
		}
        }
    if(verbose) cat('Optimizing cell type proportions... \n')
    if(infer){
        if(Normalize) CSE = X_k_norm
        if(!Normalize) CSE = X_k_m
        theta_new <- NULL
        for(j in 1:ncol(O)){
            Exp <- as.matrix(O[,j])
            rownames(Exp) <- rownames(CSE[,j,])
            colnames(Exp) <- colnames(O)[j]
            x <- CSE[,j,]
            x <- apply(x,2,scale)
            lm.o <- rlm(Exp ~ as.matrix(x),maxit=150)
            coef.v <- lm.o$coefficients[-1]
            coef.v[which(coef.v < 0)] <- 0
            total <- sum(coef.v)
            coef.v <- coef.v/total
            theta_new <- rbind(theta_new,coef.v)
        }
        colnames(theta_new) <- colnames(theta)
        rownames(theta_new) <- colnames(O);theta <- theta_new
    }
	
    writeLines( paste("Done! \nConverge in ",steps," steps",sep="") )
    res <- list()
    res$X_k <- X_k
    res$loss_history <- loss
	res$theta <- theta
    if(Normalize) res$X_k_norm <- X_k_norm
    res
}


KappaScore <- function(X_array){
   # first dim: gene
   # second dim: sample
   # third dim: cell type
   dim = dim(X_array)
   ks = 0
   for(i in 1:dim[2]){
   GEP <- X_array[,i,]
   ks <- ks + kappa(GEP)
   }
   ks = ks / dim[2]
   ks
}
KappaScore2 <- function(X_array){
   # first dim: gene
   # second dim: sample
   # third dim: cell type
   dim = dim(X_array)
   lapply(c(1:dim[2]),function(x){kappa(X_array[,x,])}) -> res
   ks = mean(unlist(res))
   ks
}

get_proportion <- function(X, ref) {
    require("MASS")
 
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
  if(total!=0){
        coef.v <- coef.v/total
  }
        theta <- rbind(theta,coef.v)
    }
    colnames(theta) <- colnames(coefVec) <- colnames(ref_m)
    rownames(theta) <- rownames(coefVec) <- colnames(X_m)
    res <- list()
    res$theta <- theta
    res$coef <- coefVec
    return(res)
}


GeneSigTest <- function(Bulk,frac,ref,nMC = 1000,p_threshold = 0.05,refine=FALSE,auto=TRUE,filtering=FALSE){
  if(auto){
    cor <- cor(ref)
	diag(cor) <- 0
	if(sum(cor>0.85) > 0){
	  refine = FALSE
	}else{
	  refine = TRUE
	}
  }
  
  writeLines("Using Linear Regression To Infer Probability of Expression...")
  p_tab <- NULL
  df <- (nrow(frac)-ncol(frac))
  for(i in 1:nrow(Bulk)){
   exp <- Bulk[i,] 
   rlm.o <- lm(exp ~ as.matrix(frac)-1)
   p <- pt(summary(rlm.o)$coef[,3], df,lower.tail=FALSE)
   p_tab <- rbind(p_tab,p)
  }
  p_tab[is.nan(p_tab)] <- 1
  pvalue.m <- p_tab
  for(i in 1:ncol(pvalue.m)) pvalue.m[,i] <- p.adjust(pvalue.m[,i],method="BH")
  colnames(pvalue.m) <- colnames(frac)
  rownames(pvalue.m) <- rownames(Bulk)
  ### soft max transformation
  if(refine){
  expp <- softmaxP(pvalue.m)
  rownames(expp) <- rownames(pvalue.m) <- rownames(Bulk)
  colnames(expp) <- colnames(pvalue.m) <- colnames(frac)
  
  score <- matrix(NA,nrow = nrow(expp),ncol = ncol(expp))
  rownames(score) <- rownames(expp)
  colnames(score) <- colnames(expp)
  
 
  writeLines("Bootstrapping...")
  for(i in 1:ncol(score)){
   ss <- expp[,i]/rowSums(expp[,-i])
   ss[is.nan(ss)] <- 0
   
   ### using bootstrapping to calculate pvalue
   tab <- list()
   nf <- nrow(p_tab)
   for(n in 1:nMC){
     ss_pe <- expp[,i]/(rowSums(expp[,-i])[sample(1:nf,nf,replace=FALSE)])
	 tab[[n]] <- ss_pe
   }
   tab <- t(matrix(unlist(tab), nrow=length(expp[,1])))
   
   vec <- NULL
   for(j in 1:nf){
      vec <- c(vec,sum(ss[j]<tab[,j]))
   }
   score[,i] <- vec/nMC
  }
  }else{
  score <- pvalue.m
  }
  ####### refine the results and only assign the gene to the most significant cell type
  if(refine) writeLines("Refining...")
  call <- matrix(0,nrow = nrow(score),ncol = ncol(score))
  rownames(call) <- rownames(score)
  colnames(call) <- colnames(score)
  
  for(ct in 1:ncol(score)){
   gene <- rownames(score)[score[,ct]<p_threshold]
   if(refine) order <- apply(pvalue.m[gene,],1,which.min)
   if(refine) gene <- gene[order == ct]
   call[gene,ct] <- 1
  }
  gene.call <- rowSums(call)
  
  ## return the results and filtered enigma object
  if(filtering){
  Bulk = Bulk[gene.call>0,]
  ref = ref[gene.call>0,]
  res <- list(
    call = call,
	pval = score,
	Bulk = Bulk,
	ref = ref
	)
  }else{
  res <- list(
    call = call,
	pval = score
	)
  }
  
  return(
   res
  )
}  


softmaxP <- function(p){
   expp <- exp(1-p)
   expp <- diag(1/rowSums(expp)) %*% expp
   expp
}



								  




