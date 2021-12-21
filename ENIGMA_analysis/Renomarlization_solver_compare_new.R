#####Compare ADMM, proximal point, and gradient renormalized proximal point method
##### implement ADMM based L2 max norm model (cell_deconvolve_admm_L2max). 
cell_deconvolve_admm_L2max <- function(O, theta, R, alpha=0.5,beta=5,gamma=1,epsilon=0.001,max.iter=100,verbose=FALSE,X_int=NULL,loss_his=TRUE){
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
    
    
    theta_hat <- colMeans(theta)
    k <- 0
    delta <- 10000
    loss <- NULL
    repeat{
        if(abs(delta)<epsilon||k>max.iter){
            break;
        }else{
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
                Y_k_1[,,j] <- t(squash(t(A_k[,,j]/gamma + X_k_1[,,j]),(beta)/gamma))
                
                A_k_1[,,j] <- A_k[,,j] + gamma*(X_k_1[,,j]-Y_k_1[,,j])
                ratio <- c(ratio, sum( (X_k_1[,,j]-X_k[,,j])^2 ))
            }
            
            if(verbose) writeLines( sprintf("   Ratio ranges from: %f - %f", min(ratio), max(ratio) ) )
            r <- sub_loss(O,X_k,theta,alpha,beta,R)
            if(loss_his) loss<- rbind(loss,c(r$part1,r$part2,r$part3))
            if(verbose) writeLines( sprintf("   Loss: part1=%f , part2=%f , part3=%f", r$part1,r$part2,r$part3 ) )
            delta <- max(ratio)
            k <- k+1
        }
    }
    writeLines( paste("Converge in ",k," steps",sep="") )
    res <- list()
    res$X_k <- X_k
    res$loss_history <- loss
    res
}
##################################################################################						   
#########################################################################################
#### To compare gradient renormalized model and non-renormalized model, we need to make sure the gradient size (dP1 and dP2) is comparable.
#### We chosed the gradient size for gradient renormalized model as 200, for more detail of parameter setting, please check gradient_size_setting.R

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
    dP1 = dP1 / sqrt( sum( dP1^2 ) ) * 2e2 ## renormalized gradient with new gradient size = 200
    dP2 = dP2 / sqrt( sum( dP2^2 ) ) * 2e2 ## renormalized gradient with new gradient size = 200 
    
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

ENIGMA_l2max_reG <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    epsilon=0.001,
                                    alpha=0.1,
                                    beta=50,tao_k=1,max.iter = 1000,verbose=TRUE,Normalize=FALSE,pos = FALSE)
									
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
    #dP1 = dP1 / sqrt( sum( dP1^2 ) ) * 1e2
    #dP2 = dP2 / sqrt( sum( dP2^2 ) ) * 1e2
    
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
							
ENIGMA_l2max <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    epsilon = 0.001,
                                    alpha=0.1,
                                    beta=50,tao_k=1,max.iter=1000,verbose=FALSE,pos=FALSE,Normalize=FALSE)

ENIGMA_l2max_miu5 <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    epsilon = 0.001,
                                    alpha=0.1,
                                    beta=5,tao_k=1,max.iter=1000,verbose=FALSE,pos=FALSE,Normalize=FALSE)

ENIGMA_l2max_miu1 <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    epsilon = 0.001,
                                    alpha=0.1,
                                    beta=1,tao_k=1,max.iter=1000,verbose=FALSE,pos=FALSE,Normalize=FALSE)	

ENIGMA_l2max_miu10 <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    epsilon = 0.001,
                                    alpha=0.1,
                                    beta=10,tao_k=1,max.iter=1000,verbose=FALSE,pos=FALSE,Normalize=FALSE)

ENIGMA_l2max_miu0 <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    epsilon = 0.001,
                                    alpha=0.1,
                                    beta=0,tao_k=1,max.iter=1000,verbose=FALSE,pos=FALSE,Normalize=FALSE)

ENIGMA_admm_l2max_miu1 <- cell_deconvolve_admm_L2max(O = as.matrix(Bulk),
                                                                   theta=Frac$theta,
                                                                   R=Reference,
                                                                   epsilon=0.001,
                                                                   alpha=0.1,beta=1,gamma=1,
                                                                   verbose=TRUE,max.iter = 1000)
ENIGMA_admm_l2max_miu5 <- cell_deconvolve_admm_L2max(O = as.matrix(Bulk),
                                                                   theta=Frac$theta,
                                                                   R=Reference,
                                                                   epsilon=0.001,
                                                                   alpha=0.1,beta=5,gamma=1,
                                                                   verbose=TRUE,max.iter = 1000)
ENIGMA_admm_l2max_miu10 <- cell_deconvolve_admm_L2max(O = as.matrix(Bulk),
                                                                   theta=Frac$theta,
                                                                   R=Reference,
                                                                   epsilon=0.001,
                                                                   alpha=0.1,beta=1,gamma=1,
                                                                   verbose=TRUE,max.iter = 1000)
ENIGMA_admm_l2max_miu50 <- cell_deconvolve_admm_L2max(O = as.matrix(Bulk),
                                                                   theta=Frac$theta,
                                                                   R=Reference,
                                                                   epsilon=0.001,
                                                                   alpha=0.1,beta=5,gamma=1,
                                                                   verbose=TRUE,max.iter = 1000)																   
LossList <- list(`Proximal Point(step size = 1, miu=50)` = ENIGMA_l2max$loss_history,
                 `Proximal Point(step size = 1, Gradient Renormalization,miu=50)` = ENIGMA_l2max_reG$loss_history,
				 `Proximal Point(step size = 1, miu=10)` = ENIGMA_l2max_miu10$loss_history,
				 `Proximal Point(step size = 1, miu=1)` = ENIGMA_l2max_miu1$loss_history,
				 `Proximal Point(step size = 1, miu=5)` = ENIGMA_l2max_miu5$loss_history)
plotMultiLossCurve(LossList,name_prefix = "Method",shape=TRUE,rlgType = "L2_max_norm")
png("/mnt/data1/weixu/HiDe/revised/Model_Compare/LossCurve/Renormalization_L2maxNorm.png",res=300,height=1500,width=3800)
plotMultiLossCurve(LossList,name_prefix = "Method",shape=TRUE,rlgType = "L2_max_norm")
dev.off()	



LossList <- list(`ADMM(step size = 1, miu=50)` = ENIGMA_admm_l2max_miu50$loss_history,
				 `ADMM(step size = 1, miu=10)` = ENIGMA_admm_l2max_miu10$loss_history,
				 `ADMM(step size = 1, miu=1)` = ENIGMA_admm_l2max_miu1$loss_history,
				 `ADMM(step size = 1, miu=5)` = ENIGMA_admm_l2max_miu5$loss_history)
plotMultiLossCurve(LossList,name_prefix = "Method",shape=TRUE,rlgType = "L2_max_norm")
png("ADMM_L2maxNorm.png",res=300,height=1500,width=3800)
plotMultiLossCurve(LossList,name_prefix = "Method",shape=TRUE,rlgType = "L2_max_norm")
dev.off()	
############################################
res0 <- DEG_test1(ENIGMA_l2max_miu0$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=TRUE)
res1 <- DEG_test1(ENIGMA_l2max_miu1$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=TRUE)
res2 <- DEG_test1(ENIGMA_l2max_miu5$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=TRUE)
res3 <- DEG_test1(ENIGMA_l2max_miu10$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=TRUE)
res4 <- DEG_test1(ENIGMA_l2max$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=TRUE)
res5 <- DEG_test1(ENIGMA_l2max_reG$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=TRUE)

#################################################
df <- data.frame(AUPRC = c(res0,res1,res2,res3,res4,res5),Method = c("Proximal Point(step size = 1, beta=0)","Proximal Point(step size = 1, beta=1)","Proximal Point(step size = 1, beta=5)","Proximal Point(step size = 1, beta=10)","Proximal Point(step size = 1, beta=50)","Proximal Point(step size = 1, Gradient Renormalization, beta=50)"))
df$Method <- factor(df$Method,level = c("Proximal Point(step size = 1, beta=0)","Proximal Point(step size = 1, beta=1)","Proximal Point(step size = 1, beta=5)","Proximal Point(step size = 1, beta=10)","Proximal Point(step size = 1, beta=50)","Proximal Point(step size = 1, Gradient Renormalization, beta=50)"))
p1<-ggplot(df, aes(x=Method, y=AUPRC, fill=Method)) +
    geom_bar(stat="identity")+theme_classic2()+theme(text = element_text(size=16),axis.text.x=element_blank ())+labs(x="")


###running steps
df <- data.frame(`Number of Iterations` = c(nrow(ENIGMA_l2max_miu0$loss_history),nrow(ENIGMA_l2max_miu1$loss_history),nrow(ENIGMA_l2max_miu5$loss_history),nrow(ENIGMA_l2max_miu10$loss_history),nrow(ENIGMA_l2max$loss_history),nrow(ENIGMA_l2max_reG$loss_history)),Method = c("Proximal Point(step size = 1, beta=0)","Proximal Point(step size = 1, beta=1)","Proximal Point(step size = 1, beta=5)","Proximal Point(step size = 1, beta=10)","Proximal Point(step size = 1, beta=50)","Proximal Point(step size = 1, Gradient Renormalization, beta=50)"))
df$Method <- factor(df$Method,level = c("Proximal Point(step size = 1, beta=0)","Proximal Point(step size = 1, beta=1)","Proximal Point(step size = 1, beta=5)","Proximal Point(step size = 1, beta=10)","Proximal Point(step size = 1, beta=50)","Proximal Point(step size = 1, Gradient Renormalization, beta=50)"))
p2<-ggplot(df, aes(x=Method, y=`Number.of.Iterations`, fill=Method)) +
    geom_bar(stat="identity")+theme_classic2()+theme(text = element_text(size=16),axis.text.x=element_blank ())+labs(x="",y="Number of Iterations")+NoLegend()
png("Barplot(L2maxNorm).png",res=300,height=2000,width=4000)
plot_grid(p2,p1+NoLegend(),nrow=1)
dev.off()
##############################################################
#Compare ADMM and proximal point learning performance
p <- 100
y <- gl(2, p/2)
res1_admm <- DEG_test1(ENIGMA_admm_l2max_miu1$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=TRUE)
res2_admm <- DEG_test1(ENIGMA_admm_l2max_miu5$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=TRUE)
res3_admm <- DEG_test1(ENIGMA_admm_l2max_miu10$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=TRUE)
res4_admm <- DEG_test1(ENIGMA_admm_l2max_miu50$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=TRUE)

df <- data.frame(AUPRC = c(res1,res2,res3,res4,res1_admm,res2_admm,res3_admm,res4_admm,res5),Method = c(rep(c("Proximal Point"),4),rep(c("ADMM"),4),"Proximal Point(Gradient Renormalization)"),beta = c(rep(c(1,5,10,50),2),50))
df$beta <- as.factor(as.character(df$beta))
df$beta <- factor(df$beta,levels=c("1","5","10","50"))
png("/mnt/data1/weixu/HiDe/revised/Model_Compare/ADMM_pps/Barplot(L2maxNorm).png",res=300,height=1500,width=3000)
pp1<-ggplot(df, aes(x=beta, y=AUPRC, fill=Method)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=16),axis.text.x=element_blank ())+labs(x="")
pp1
dev.off()

###running steps
df <- data.frame(`Number of Iterations` = c(nrow(ENIGMA_l2max_miu1$loss_history),nrow(ENIGMA_l2max_miu5$loss_history),nrow(ENIGMA_l2max_miu10$loss_history),nrow(ENIGMA_l2max$loss_history),nrow(ENIGMA_admm_l2max_miu1$loss_history),nrow(ENIGMA_admm_l2max_miu5$loss_history),nrow(ENIGMA_admm_l2max_miu10$loss_history),nrow(ENIGMA_admm_l2max_miu50$loss_history),nrow(ENIGMA_l2max_reG$loss_history)),Method = c(rep(c("Proximal Point"),4),rep(c("ADMM"),4),"Proximal Point(Gradient Renormalization)"),beta = c(rep(c(1,5,10,50),2),50))
df$beta <- as.factor(as.character(df$beta))
df$beta <- factor(df$beta,levels=c("1","5","10","50"))
pp2<-ggplot(df, aes(x=beta, y=`Number.of.Iterations`, fill=Method)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=16),axis.text.x=element_blank ())+labs(x="",y="Number of Iterations")+NoLegend()
png("Barplot(iteration,L2maxNorm).png",res=300,height=1500,width=3000)
pp2<-ggplot(df, aes(x=beta, y=`Number.of.Iterations`, fill=Method)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=16),axis.text.x=element_blank ())+labs(x="",y="Number of Iterations")
pp2
dev.off()
####################################################################################################
##prove L2-norm is working
model1 <- ENIGMA_l2max_miu0

model2 <- cell_deconvolve(X=as.matrix(Bulk),
                                theta=Frac$theta,
                                R=Reference,
                                inner_epilson=0.001,
                                outer_epilson=0.001,
                                alpha=0.3,
                                miu=0,tao_k=1,max.iter=1,max.iter.exp=1000,verbose=FALSE)

model3 <- cell_deconvolve(X=as.matrix(Bulk),
                                theta=Frac$theta,
                                R=Reference,
                                inner_epilson=0.001,
                                outer_epilson=0.001,
                                alpha=0.5,
                                miu=0,tao_k=1,max.iter=1,max.iter.exp=1000,verbose=FALSE)

model4 <- cell_deconvolve(X=as.matrix(Bulk),
                                theta=Frac$theta,
                                R=Reference,
                                inner_epilson=0.001,
                                outer_epilson=0.001,
                                alpha=0.7,
                                miu=0,tao_k=1,max.iter=1,max.iter.exp=1000,verbose=FALSE)

model5 <- cell_deconvolve(X=as.matrix(Bulk),
                                theta=Frac$theta,
                                R=Reference,
                                inner_epilson=0.001,
                                outer_epilson=0.001,
                                alpha=0.9,
                                miu=0,tao_k=1,max.iter=1,max.iter.exp=1000,verbose=FALSE)
								
res1 <- DEG_test1(model1$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE)
res2 <- DEG_test1(model2$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE)
res3 <- DEG_test1(model3$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE)
res4 <- DEG_test1(model4$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE)
res5 <- DEG_test1(model5$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE)
########################################################################
##
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
    dP1 = dP1 / sqrt( sum( dP1^2 ) ) * 1e2
    dP2 = dP2 / sqrt( sum( dP2^2 ) ) * 1e2
    
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

ENIGMA_l2max_reG1 <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    inner_epilson=0.001,
                                    outer_epilson=0.001,
                                    alpha=0.1,
                                    miu=50,tao_k=1,max.iter=1,max.iter.exp=150,verbose=TRUE)

ENIGMA_l2max_reG2 <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    inner_epilson=0.001,
                                    outer_epilson=0.001,
                                    alpha=0.3,
                                    miu=50,tao_k=1,max.iter=1,max.iter.exp=150,verbose=FALSE)

ENIGMA_l2max_reG3 <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    inner_epilson=0.001,
                                    outer_epilson=0.001,
                                    alpha=0.5,
                                    miu=50,tao_k=1,max.iter=1,max.iter.exp=150,verbose=FALSE)

ENIGMA_l2max_reG4 <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    inner_epilson=0.001,
                                    outer_epilson=0.001,
                                    alpha=0.7,
                                    miu=50,tao_k=1,max.iter=1,max.iter.exp=150,verbose=FALSE)

ENIGMA_l2max_reG6 <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    inner_epilson=0.001,
                                    outer_epilson=0.001,
                                    alpha=0.9,
                                    miu=50,tao_k=0.1,max.iter=1,max.iter.exp=1000,verbose=FALSE)

res1_reG <- DEG_test1(ENIGMA_l2max_reG1$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE)
res2_reG <- DEG_test1(ENIGMA_l2max_reG2$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE)
res3_reG <- DEG_test1(ENIGMA_l2max_reG3$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE)
res4_reG <- DEG_test1(ENIGMA_l2max_reG4$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE)
res5_reG <- DEG_test1(ENIGMA_l2max_reG5$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE)
########
perform_Rlg <- c(res1_reG,res2_reG,res3_reG,res4_reG,res5_reG)
perform_woRlg <- c(res1,res2,res3,res4,res5)
df2 <- data.frame(AUPRC = c(perform_Rlg,perform_woRlg),Method = c(rep("L2 MaxNorm",length(perform_Rlg)),rep("Without Regularization",length(perform_woRlg))),alpha = c(0.1,0.3,0.5,0.7,0.9))
p<-ggplot(df2, aes(x=alpha, y=AUPRC, group=Method)) +
  geom_line(aes(color=Method))+
  geom_point(aes(color=Method))+theme_classic()
png("DEG_Performance.png",res=300,width=1500,height=1000)
p
dev.off()
#######################################
##################
res <- NULL
ES <- c(1.8,2.4,3,3.6,4.2,4.8)
for(es in 1:length(ES)){
load(paste("/Path/to/save/Data/resCompare_bi_",ES[es],".Rdata",sep=""))
line <- NULL
for(k in 1:10){
Bulk <- testMatrixG[[k]]
Reference <- t(testMatrixH1[[k]])
DEG_list <- DEG_list_all[[k]]
Frac <- get_proportion(Bulk, Reference)

ENIGMA_l2max_reG <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    inner_epilson=0.001,
                                    outer_epilson=0.001,
                                    alpha=0.1,
                                    miu=50,tao_k=1,max.iter=1,max.iter.exp=150,verbose=FALSE)
									
reG <- DEG_test1(ENIGMA_l2max_reG$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE)
line <- c(line,reG)
}
res <- rbind(res, line)
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
    #dP1 = dP1 / sqrt( sum( dP1^2 ) ) * 1e2
    #dP2 = dP2 / sqrt( sum( dP2^2 ) ) * 1e2
    
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

res2 <- NULL
ES <- c(1.8,2.4,3,3.6,4.2,4.8)
for(es in 1:length(ES)){
load(paste("/Path/to/save/Data/resCompare_bi_",ES[es],".Rdata",sep=""))
line <- NULL
for(k in 1:10){
Bulk <- testMatrixG[[k]]
Reference <- t(testMatrixH1[[k]])
DEG_list <- DEG_list_all[[k]]
Frac <- get_proportion(Bulk, Reference)

ENIGMA_l2max_reG <- cell_deconvolve(X=as.matrix(Bulk),
                                    theta=Frac$theta,
                                    R=Reference,
                                    inner_epilson=0.001,
                                    outer_epilson=0.001,
                                    alpha=0.1,
                                    miu=0,tao_k=1,max.iter=1,max.iter.exp=1000,verbose=FALSE)
									
reG <- DEG_test1(ENIGMA_l2max_reG$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE)
line <- c(line,reG)
}
res2 <- rbind(res2, line)
}

##############################
#using boxplot to show the improvements

SNR <- rep(c(rep(ES[1],10),rep(ES[2],10),rep(ES[3],10),rep(ES[4],10),rep(ES[5],10),rep(ES[6],10)),2)
method <- c(rep("L2 MaxNorm",60),rep("Without Regularization",60))
performance <- c(as.numeric(t(res)),as.numeric(t(res2)))
sp.m <- data.frame(SNR = SNR, Method = method, AUPRC = performance)

mytheme <- readRDS("/Path/to/save/Data/mytheme.rds")
p_boxplot <- sp.m %>% 
    mutate(SNR=paste0("SNR=", SNR)) %>% 
    ggplot(aes(Method, AUPRC, color=Method)) + 
    geom_boxplot() + 
    facet_grid(~SNR) + 
    mytheme + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
    theme(legend.title = element_text(size = 12, color = "black", family = "Arial"), legend.text = element_text(size = 12, color = "black", family = "Arial")) + 
    theme(panel.border = element_rect(size = 0.3, linetype = "dashed", fill = NA))
png("/Path/to/save/Data/boxplot_DEG.png",res=300,height=1200,width=2400)
p_boxplot
dev.off()

###########################################################################
###compare the latent state identification
Bulk <- readRDS("/Path/to/save/Data/Bulk.rds")
Reference <- readRDS("/Path/to/save/Data/Reference.rds")
cellLabel <- readRDS("/Path/to/save/Data/CellLabel.rds")
Frac <- get_proportion(Bulk, Reference)

######
#raw model
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
    dP1 = dP1 / sqrt( sum( dP1^2 ) ) * 1e2
    dP2 = dP2 / sqrt( sum( dP2^2 ) ) * 1e2
    
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
ENIGMA_l2max <- cell_deconvolve(X=as.matrix(sqrt(Bulk)),
                                    theta=Frac$theta,
                                    R=sqrt(Reference),
                                    inner_epilson=0.001,
                                    outer_epilson=0.001,
                                    alpha=0.5,
                                    miu=5,tao_k=10,max.iter=1,max.iter.exp=1000,verbose=TRUE)

X_k <- X_k_norm <- ENIGMA_l2max$X_k
for(k in 1:dim(X_k)[3]){
	X_k[,,k][X_k[,,k]<0] <- 0
	exp <- X_k[,,k]^2
	
	exp.scale <- t(apply(exp,1,scale))
	PC <- svd(exp.scale)$v[,1]
	exp.norm <- NULL
	for(i in 1:nrow(exp)){
		lm.model <- lm(exp[i,]~PC)
		exp.norm <- rbind(exp.norm, (exp[i,] - lm.model$coefficients[2] * PC))
	}
	X_k_norm[,,k] <- exp.norm
}
        
sce <- SingleCellExperiment(assays = list(logcounts = X_k_norm[,,1]))
sce$cell_type <- cellLabel$c1
sce <- sce[,Frac$theta[,1]>0.05]
sce <- runPCA(sce)
sce <- runTSNE(sce,dimred="PCA",n_dimred=6)
label <- sce$cell_type
ARI_ct1 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct1 <- c(ARI_ct1,ARI)
}


sce2 <- SingleCellExperiment(assays = list(logcounts = X_k_norm[,,3]))
sce2$cell_type <- cellLabel$c3
sce2 <- sce2[,Frac$theta[,3]>0.05]
sce2 <- runPCA(sce2)
sce2 <- runTSNE(sce2,dimred="PCA",n_dimred=6)
label <- sce2$cell_type
ARI_ct2 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce2,"PCA")[,1:i],2)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct2 <- c(ARI_ct2,ARI)
}

ARI_Rlg <- list(ARI_ct1 = ARI_ct1,ARI_ct2 = ARI_ct2)

png("Rlg_tsne.png",res=300,height=1200,width=2600)
p1 <- plotTSNE(sce, colour_by="cell_type",point_size=3)
p2 <- plotTSNE(sce2, colour_by="cell_type",point_size=3)
cowplot::plot_grid(p1,p2,nrow=1)
dev.off()

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
    #dP1 = dP1 / sqrt( sum( dP1^2 ) ) * 1e2
    #dP2 = dP2 / sqrt( sum( dP2^2 ) ) * 1e2
    
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
ENIGMA_l2max_miu0 <- cell_deconvolve(X=as.matrix(sqrt(Bulk)),
                                    theta=Frac$theta,
                                    R=sqrt(Reference),
                                    inner_epilson=0.001,
                                    outer_epilson=0.001,
                                    alpha=0.5,
                                    miu=0,tao_k=0.01,max.iter=1,max.iter.exp=200,verbose=TRUE)
									
X_k <- X_k_norm <- ENIGMA_l2max_miu0$X_k
for(k in 1:dim(X_k)[3]){
	X_k[,,k][X_k[,,k]<0] <- 0
	exp <- X_k[,,k]^2
	
	exp.scale <- t(apply(exp,1,scale))
	PC <- svd(exp.scale)$v[,1]
	exp.norm <- NULL
	for(i in 1:nrow(exp)){
		lm.model <- lm(exp[i,]~PC)
		exp.norm <- rbind(exp.norm, (exp[i,] - lm.model$coefficients[2] * PC))
	}
	X_k_norm[,,k] <- exp.norm
}
        
sce <- SingleCellExperiment(assays = list(logcounts = X_k_norm[,,1]))
sce$cell_type <- cellLabel$c1
sce <- sce[,Frac$theta[,1]>0.05]
sce <- runPCA(sce)
sce <- runTSNE(sce,dimred="PCA",n_dimred=5)
label <- sce$cell_type
ARI_ct1 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce,"PCA")[,1:i],4)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct1 <- c(ARI_ct1,ARI)
}


sce2 <- SingleCellExperiment(assays = list(logcounts = X_k_norm[,,3]))
sce2$cell_type <- cellLabel$c3
sce2 <- sce2[,Frac$theta[,3]>0.05]
sce2 <- runPCA(sce2)
sce2 <- runTSNE(sce2,dimred="PCA",n_dimred=5)
label <- sce2$cell_type
ARI_ct2 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
clust <- cluster::pam(reducedDim(sce2,"PCA")[,1:i],2)
ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
print(ARI)
ARI_ct2 <- c(ARI_ct2,ARI)
}

ARI_WithOutRlg <- list(ARI_ct1 = ARI_ct1,ARI_ct2 = ARI_ct2)

png("WithOutRlg_tsne.png",res=300,height=1200,width=2600)
p1 <- plotTSNE(sce, colour_by="cell_type",point_size=3)
p2 <- plotTSNE(sce2, colour_by="cell_type",point_size=3)
cowplot::plot_grid(p1,p2,nrow=1)
dev.off()

#############
res1 <- max(ARI_Rlg$ARI_ct1)
res2 <- max(ARI_Rlg$ARI_ct2)
res3 <- max(ARI_WithOutRlg$ARI_ct1)
res4 <- max(ARI_WithOutRlg$ARI_ct2)
df <- data.frame(ARI = c(res1,res2,res3,res4),Method = c("L2 MaxNorm","L2 MaxNorm","Without Regularization","Without Regularization"),CellType = c("CellType1","CellType3","CellType1","CellType3"))
p<-ggplot(df, aes(x=CellType, y=ARI, fill=Method)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=20))+labs(x="")
png("BarplotARI.png",res=300,height=1800,width=2000)
p
dev.off()


df <- data.frame(`Number of Iterations` = c(nrow(ENIGMA_l2max$loss_history),nrow(ENIGMA_l2max_miu0$loss_history)),Method = c("L2 MaxNorm","Without Regularization"))
p2<-ggplot(df, aes(x=Method, y=`Number.of.Iterations`, fill=Method)) +
    geom_bar(stat="identity")+theme_classic2()+theme(text = element_text(size=16),axis.text.x=element_blank ())+labs(x="",y="Number of Iterations")
png("BarplotARI.png",res=300,height=1800,width=1700)
p2
dev.off()
