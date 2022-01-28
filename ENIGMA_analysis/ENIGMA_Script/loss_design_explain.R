########################################################
#explain why rank norm regularization should plugging in the average cell type proportion
Bulk <- readRDS("/mnt/data1/weixu/HiDe/revised/Model_Compare/CellStateIdentification3/Bulk.rds")
Reference <- readRDS("/mnt/data1/weixu/HiDe/revised/Model_Compare/CellStateIdentification3/Reference.rds")
cellLabel <- readRDS("/mnt/data1/weixu/HiDe/revised/Model_Compare/CellStateIdentification3/CellLabel.rds")
Frac <- get_proportion(Bulk, Reference)

############################################################
cell_deconvolve_trace_proximal_point_solver_final <- function(O, theta, R, alpha=0.5,beta=5,epsilon=0.001,tao_k=0.01,max.iter=100,X_int=NULL,verbose=FALSE,infer=FALSE,Normalize=TRUE,Norm.method = "PC",loss_his=TRUE){
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
    
    theta_hat <- colMeans(theta)
    k <- 0
    X_k_1 <- X_k <- X
    ##record loss
    if(verbose) cat(date(), 'Optimizing cell type specific expression profile... \n')
    loss <- NULL
    repeat{
        ratio <- NULL
        ## 更新第i个细胞类型表达谱梯度
        dP <- derive_P(O, theta,X_k,R,alpha)
        ## 返回第i个细胞类型表达谱（P）
        ## 利用gradient decent对P进行迭代更新
        for(i in 1:ncol(theta)){
            X_i <- X_k[,,i]- tao_k*dP[,,i]
            X_i <- SVT(X_i,tao_k*mean(theta_hat)*beta) ####without plugging cell type specific proportions
            X_k_1[,,i] <- X_i
            
            ## check更新前后差异
            ratio <- c(ratio, sum( (X_i-X_k[,,i])^2 )/(nrow(X_i)*ncol(X_i)))
        }
        r1 <- loss(O,X_k,theta,alpha,beta,R) #calculating raw loss
        #print updating
        if(verbose){
			   #print matrix ratio distance or absolute distance
               print <- paste("CSE inference step ",k," \n",sep="")
			   for(c in 1:(length(ratio)-1)){
			      print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
			   }
			   print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
			   writeLines(print)
		}
        if(verbose) writeLines( sprintf("   Loss: part1=%f , part2=%f , part3=%f", r1$part1,r1$part2,r1$part3 ) )
        if(max(ratio) < epsilon||k >= max.iter){break}else{
            ## 更新X
            X_k <- X_k_1
            k <- k + 1
            if(loss_his) loss<- rbind(loss,c(r1$part1,r1$part2,r1$part3))
        }
    }
    steps <- k
    ##########
    #Doing PC or Cell Fractions Normalization
    if(Normalize){
        if(verbose) cat("Perform Normalization...")
        X_k_norm <- X_k
        if(Norm.method == "PC"){
            for(k in 1:dim(X_k)[3]){
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
        }
        
        if(Norm.method == "frac"){
            for(k in 1:dim(X_k)[3]){
                exp <- X_k[,,k]^2
                exp.norm <- NULL
                for(i in 1:nrow(exp)){
                    lm.model <- lm(exp[i,]~theta[,k])
                    exp.norm <- rbind(exp.norm, (exp[i,] - lm.model$coefficients[2] * theta[,k]))
                }
                X_k_norm[,,k] <- exp.norm
            }
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
        rownames(theta_new) <- colnames(X)
    }
    
    writeLines( paste("Converge in ",steps," steps",sep="") )
    res <- list()
    res$X_k <- X_k
    res$loss_history <- loss
    if(Normalize) res$X_k_norm <- X_k_norm
    res
}

egm_without_theta <- cell_deconvolve_trace_proximal_point_solver_final(O = as.matrix(sqrt(Bulk)),
                                                            theta=Frac$theta,
                                                            R=sqrt(Reference),
                                                            epsilon=0.00001,
                                                            alpha=0.1,beta=1,tao_k = 0.5,
                                                            verbose=FALSE,max.iter = 1000)
															
														
#######################################
cell_deconvolve_trace_proximal_point_solver_final <- function(O, theta, R, alpha=0.5,beta=5,epsilon=0.001,tao_k=0.01,max.iter=100,X_int=NULL,verbose=FALSE,infer=FALSE,Normalize=TRUE,Norm.method = "PC",loss_his=TRUE){
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
    
    theta_hat <- colMeans(theta)
    k <- 0
    X_k_1 <- X_k <- X
    ##record loss
    if(verbose) cat(date(), 'Optimizing cell type specific expression profile... \n')
    loss <- NULL
    repeat{
        ratio <- NULL
        ## 更新第i个细胞类型表达谱梯度
        dP <- derive_P(O, theta,X_k,R,alpha)
        ## 返回第i个细胞类型表达谱（P）
        ## 利用gradient decent对P进行迭代更新
        for(i in 1:ncol(theta)){
            X_i <- X_k[,,i]- tao_k*dP[,,i]
            X_i <- SVT(X_i,tao_k*theta_hat[i]*beta) ####without plugging cell type specific proportions
            X_k_1[,,i] <- X_i
            
            ## check更新前后差异
            ratio <- c(ratio, sum( (X_i-X_k[,,i])^2 )/(nrow(X_i)*ncol(X_i)))
        }
        r1 <- loss(O,X_k,theta,alpha,beta,R) #calculating raw loss
        #print updating
        if(verbose){
			   #print matrix ratio distance or absolute distance
               print <- paste("CSE inference step ",k," \n",sep="")
			   for(c in 1:(length(ratio)-1)){
			      print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
			   }
			   print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
			   writeLines(print)
		}
        if(verbose) writeLines( sprintf("   Loss: part1=%f , part2=%f , part3=%f", r1$part1,r1$part2,r1$part3 ) )
        if(max(ratio) < epsilon||k >= max.iter){break}else{
            ## 更新X
            X_k <- X_k_1
            k <- k + 1
            if(loss_his) loss<- rbind(loss,c(r1$part1,r1$part2,r1$part3))
        }
    }
    steps <- k
    ##########
    #Doing PC or Cell Fractions Normalization
    if(Normalize){
        if(verbose) cat("Perform Normalization...")
        X_k_norm <- X_k
        if(Norm.method == "PC"){
            for(k in 1:dim(X_k)[3]){
                exp <- X_k[,,k]^2
                exp.scale <- t(apply(exp,1,scale))
                PC <- svd(exp.scale)$v[,1]
                the <- (exp %*% as.matrix(PC) - length(PC) * mean(PC) * rowMeans(exp)) / (sum(PC^2) - length(PC)*mean(PC)^2)
                exp.norm <- exp - as.matrix(the) %*% t(as.matrix(PC))
                X_k_norm[,,k] <- exp.norm
            }
        }
        
        if(Norm.method == "frac"){
            for(k in 1:dim(X_k)[3]){
                exp <- X_k[,,k]^2
                PC <- theta[,k]
                the <- (exp %*% as.matrix(PC) - length(PC) * mean(PC) * rowMeans(exp)) / (sum(PC^2) - length(PC)*mean(PC)^2)
                exp.norm <- exp - as.matrix(the) %*% t(as.matrix(PC))
                X_k_norm[,,k] <- exp.norm
            }
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
        rownames(theta_new) <- colnames(X)
    }
    
    writeLines( paste("Converge in ",steps," steps",sep="") )
    res <- list()
    res$X_k <- X_k
    res$loss_history <- loss
    if(Normalize) res$X_k_norm <- X_k_norm
    res
}

egm_with_theta <- cell_deconvolve_trace_proximal_point_solver_final(O = as.matrix(sqrt(Bulk)),
                                                            theta=Frac$theta,
                                                            R=sqrt(Reference),
                                                            epsilon=0.00001,
                                                            alpha=0.1,beta=1,tao_k = 0.5,
                                                            verbose=FALSE,max.iter = 1000)
##################################################################################

###############plot the ground truth eigenvalues and compare with the inferred
model1 <- egm_without_theta
model2 <- egm_with_theta

rm(egm_with_theta,egm_without_theta)
gc()

library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(ggpubr)

clustering <- function(array,ct,Fraction,num,clust_num){
sce <- SingleCellExperiment(assays = list(logcounts = array))
sce$cell_type <- ct
sce <- sce[,Fraction$theta[,num]>0.05]
sce <- runPCA(sce)
sce <- runTSNE(sce,dimred="PCA",n_dimred=5)
label <- sce$cell_type
ARI_ct1 <- NULL
for(i in c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)){
    clust <- cluster::pam(reducedDim(sce,"PCA")[,1:i],clust_num)
    ARI <- mclust::adjustedRandIndex(as.numeric(as.factor(label)),as.numeric(clust$cluster))
    ARI_ct1 <- c(ARI_ct1,ARI)
}
max(ARI_ct1)
}

ct1_ari <- clustering(model1$X_k_norm[,,1],cellLabel$c1,Frac,1,4)
ct1_ari2 <- clustering(model2$X_k_norm[,,1],cellLabel$c1,Frac,1,4)

ct3_ari <- clustering(model1$X_k_norm[,,3],cellLabel$c3,Frac,3,2)
ct3_ari2 <- clustering(model2$X_k_norm[,,3],cellLabel$c3,Frac,3,2)

df <- data.frame(ARI = c(ct1_ari,ct1_ari2,ct3_ari,ct3_ari2),Method = rep(c("model1","model2"),2),CellType = c(rep("cell type1",2),rep("cell type3",2)))
pg <- ggplot(df, aes(x=CellType, y=ARI, fill=Method)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=20))+labs(x="")
png("barplot_performance.png",res=300,height=1500,width=1700)
pg
dev.off()


##################
#DEG analysis
library(MIND)
library(TCA)
library(MASS)
library(ggplot2)
library(ggpubr)
###############################################
#####
source("DEG_analysis_uile_function.R")

load("DEG_example_data.Rdata")
Frac <- get_proportion(Bulk, Reference)


####################################
Tab <- NULL
alpha.set <- c(0.1,0.3,0.5,0.7,0.9)
for(as in alpha.set){

cell_deconvolve_trace_proximal_point_solver_final <- function(O, theta, R, alpha=0.5,beta=5,epsilon=0.001,tao_k=0.01,max.iter=100,X_int=NULL,verbose=FALSE,infer=FALSE,Normalize=TRUE,Norm.method = "PC",loss_his=TRUE){
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
    
    theta_hat <- colMeans(theta)
    k <- 0
    X_k_1 <- X_k <- X
    ##record loss
    if(verbose) cat(date(), 'Optimizing cell type specific expression profile... \n')
    loss <- NULL
    repeat{
        ratio <- NULL
        dP <- derive_P(O, theta,X_k,R,alpha)
      
        for(i in 1:ncol(theta)){
            X_i <- X_k[,,i]- tao_k*dP[,,i]
            X_i <- SVT(X_i,tao_k*mean(theta_hat)*beta) ####without plugging cell type specific proportions
            X_k_1[,,i] <- X_i
            
           
            ratio <- c(ratio, sum( (X_i-X_k[,,i])^2 )/(nrow(X_i)*ncol(X_i)))
        }
        r1 <- loss(O,X_k,theta,alpha,beta,R) #calculating raw loss
        #print updating
        if(verbose){
			   #print matrix ratio distance or absolute distance
               print <- paste("CSE inference step ",k," \n",sep="")
			   for(c in 1:(length(ratio)-1)){
			      print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
			   }
			   print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
			   writeLines(print)
		}
        if(verbose) writeLines( sprintf("   Loss: part1=%f , part2=%f , part3=%f", r1$part1,r1$part2,r1$part3 ) )
        if(max(ratio) < epsilon||k >= max.iter){break}else{
         
            X_k <- X_k_1
            k <- k + 1
            if(loss_his) loss<- rbind(loss,c(r1$part1,r1$part2,r1$part3))
        }
    }
    steps <- k
    ##########
    #Doing PC or Cell Fractions Normalization
    if(Normalize){
        if(verbose) cat("Perform Normalization...")
        X_k_norm <- X_k
        if(Norm.method == "PC"){
            for(k in 1:dim(X_k)[3]){
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
        }
        
        if(Norm.method == "frac"){
            for(k in 1:dim(X_k)[3]){
                exp <- X_k[,,k]^2
                exp.norm <- NULL
                for(i in 1:nrow(exp)){
                    lm.model <- lm(exp[i,]~theta[,k])
                    exp.norm <- rbind(exp.norm, (exp[i,] - lm.model$coefficients[2] * theta[,k]))
                }
                X_k_norm[,,k] <- exp.norm
            }
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
        rownames(theta_new) <- colnames(X)
    }
    
    writeLines( paste("Converge in ",steps," steps",sep="") )
    res <- list()
    res$X_k <- X_k
    res$loss_history <- loss
    if(Normalize) res$X_k_norm <- X_k_norm
    res
}

egm_without_theta <- cell_deconvolve_trace_proximal_point_solver_final(O = as.matrix(Bulk),
                                                            theta=Frac$theta,
                                                            R=Reference,
                                                            epsilon=0.00001,
                                                            alpha=as,beta=1,tao_k = 2,
                                                            verbose=FALSE,max.iter = 1000)
															
														

cell_deconvolve_trace_proximal_point_solver_final <- function(O, theta, R, alpha=0.5,beta=5,epsilon=0.001,tao_k=0.01,max.iter=100,X_int=NULL,verbose=FALSE,infer=FALSE,Normalize=TRUE,Norm.method = "PC",loss_his=TRUE){
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
    
    theta_hat <- colMeans(theta)
    k <- 0
    X_k_1 <- X_k <- X
    ##record loss
    if(verbose) cat(date(), 'Optimizing cell type specific expression profile... \n')
    loss <- NULL
    repeat{
        ratio <- NULL
        dP <- derive_P(O, theta,X_k,R,alpha)
        for(i in 1:ncol(theta)){
            X_i <- X_k[,,i]- tao_k*dP[,,i]
            X_i <- SVT(X_i,tao_k*theta_hat[i]*beta) ####without plugging cell type specific proportions
            X_k_1[,,i] <- X_i

            ratio <- c(ratio, sum( (X_i-X_k[,,i])^2 )/(nrow(X_i)*ncol(X_i)))
        }
        r1 <- loss(O,X_k,theta,alpha,beta,R) #calculating raw loss
        #print updating
        if(verbose){
			   #print matrix ratio distance or absolute distance
               print <- paste("CSE inference step ",k," \n",sep="")
			   for(c in 1:(length(ratio)-1)){
			      print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
			   }
			   print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
			   writeLines(print)
		}
        if(verbose) writeLines( sprintf("   Loss: part1=%f , part2=%f , part3=%f", r1$part1,r1$part2,r1$part3 ) )
        if(max(ratio) < epsilon||k >= max.iter){break}else{
            X_k <- X_k_1
            k <- k + 1
            if(loss_his) loss<- rbind(loss,c(r1$part1,r1$part2,r1$part3))
        }
    }
    steps <- k
    ##########
    #Doing PC or Cell Fractions Normalization
    if(Normalize){
        if(verbose) cat("Perform Normalization...")
        X_k_norm <- X_k
        if(Norm.method == "PC"){
            for(k in 1:dim(X_k)[3]){
                exp <- X_k[,,k]^2
                exp.scale <- t(apply(exp,1,scale))
                PC <- svd(exp.scale)$v[,1]
                the <- (exp %*% as.matrix(PC) - length(PC) * mean(PC) * rowMeans(exp)) / (sum(PC^2) - length(PC)*mean(PC)^2)
                exp.norm <- exp - as.matrix(the) %*% t(as.matrix(PC))
                X_k_norm[,,k] <- exp.norm
            }
        }
        
        if(Norm.method == "frac"){
            for(k in 1:dim(X_k)[3]){
                exp <- X_k[,,k]^2
                PC <- theta[,k]
                the <- (exp %*% as.matrix(PC) - length(PC) * mean(PC) * rowMeans(exp)) / (sum(PC^2) - length(PC)*mean(PC)^2)
                exp.norm <- exp - as.matrix(the) %*% t(as.matrix(PC))
                X_k_norm[,,k] <- exp.norm
            }
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
        rownames(theta_new) <- colnames(X)
    }
    
    writeLines( paste("Converge in ",steps," steps",sep="") )
    res <- list()
    res$X_k <- X_k
    res$loss_history <- loss
    if(Normalize) res$X_k_norm <- X_k_norm
    res
}

egm_with_theta <- cell_deconvolve_trace_proximal_point_solver_final(O = as.matrix(Bulk),
                                                            theta=Frac$theta,
                                                            R=Reference,
                                                            epsilon=0.00001,
                                                            alpha=as,beta=1,tao_k = 2,
                                                            verbose=FALSE,max.iter = 1000)

p <- 100
y <- gl(2, p/2)
res <- c(DEG_test1(egm_without_theta$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE),DEG_test1(egm_with_theta$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE))
Tab <- rbind(Tab,res)
}

####barplot perform
res <- cbind(as.numeric(Tab),rep(alpha.set,4),c(rep("model1",5),rep("model2",5)))
res <- as.data.frame(res)
colnames(res) <- c("AUPRC","alpha","Model")
res$AUPRC <- as.numeric(as.matrix(res$AUPRC))
res$alpha <- as.character(as.matrix(res$alpha))
barp <-ggplot(res, aes(x=alpha, y=AUPRC, fill=Model)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=20))+labs(x="")
png("barplot_DEG_performance(balance).png",res=300,height=1000,width=1800)
barp
dev.off()
Tab_balance <- Tab
#######################################################################################################
####generate the DEG with inbalanced cell type proportions
k <- 5 # number of cell types
ng <- 1000 # number of genes
p <- 100 # number of samples
ndiff <- 0.05*1000 # number of genes differentially expressed

H1 <- matrix(rnorm(5*ng), ncol=ng)
H2 <- H1
# create differential expression for 3rd cell type
DEG_list <- list()
seq <- 1:ncol(H2)
for(i in 1:nrow(H2)){
    DEG_id <- sample(1:ncol(H2),ndiff,replace=FALSE)
	add <- sample(c(4.8,-4.8),ndiff,replace=TRUE)
    H2[i,DEG_id] <- H2[i,DEG_id] + add
    DEG_list[[i]] <- DEG_id * sign(add)
    #seq <- seq[seq %in% DEG_id == FALSE]
}
# cell frequency matrix per sample
cc <- DirichletReg::rdirichlet(100, c(0.4,0.2,0.2,0.1,0.1))
colnames(cc) <- paste('cellType', 1:ncol(cc), sep="")


G <- rbind(cc[1:(floor(p/2)), ] %*% H1, cc[(floor(p/2)+1):p, ] %*%H2 ) + matrix(rnorm(p*ng), ncol=ng)
# sample classes (2 groups)
y <- gl(2, p/2)

###Add the noise to each cell type and generate a sample specific profile
H1_array <- array(0,
                  dim = c( nrow(H1),
                           ncol(H1),
                           (floor(100/2))))
for(i in 1:(floor(p/2))){
    H1_array[,,i] <- H1 + matrix(rnorm(5*ng), ncol=ng)
}
H2_array <- array(0,
                  dim = c( nrow(H2),
                           ncol(H2),
                           (floor(100/2))))
for(i in 1:(floor(p/2))){
    H2_array[,,i] <- H2 + matrix(rnorm(5*ng), ncol=ng)
}
##evaluate differential expression genes
##calculate 
G <- NULL
for(i in 1:50){
    G <- cbind(G, t(as.matrix(t(as.matrix(cc[i,])) %*% H1_array[,,i])))
}
for(i in 1:50){
    G <- cbind(G, t(as.matrix(t(as.matrix(cc[i+50,])) %*% H2_array[,,i])))
}
G <- G + t(matrix(rnorm(100*ng), ncol=ng))
# sample classes (2 groups)
y <- gl(2, p/2)
##########################################################
Bulk <- G
Reference <- t(H1)
rownames(Bulk) <- rownames(Reference) <- paste0("gene-",1:nrow(G),sep="")
colnames(Bulk) <- rownames(cc) <- paste("sample-",1:ncol(G),sep="")
colnames(Reference) <- colnames(cc)


Frac <- get_proportion(Bulk, Reference)
####################################
##Model W/o 1

Tab <- NULL
alpha.set <- c(0.1,0.3,0.5,0.7,0.9)
for(as in alpha.set){

cell_deconvolve_trace_proximal_point_solver_final <- function(O, theta, R, alpha=0.5,beta=5,epsilon=0.001,tao_k=0.01,max.iter=100,X_int=NULL,verbose=FALSE,infer=FALSE,Normalize=TRUE,Norm.method = "PC",loss_his=TRUE){
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
    
    theta_hat <- colMeans(theta)
    k <- 0
    X_k_1 <- X_k <- X
    ##record loss
    if(verbose) cat(date(), 'Optimizing cell type specific expression profile... \n')
    loss <- NULL
    repeat{
        ratio <- NULL
        dP <- derive_P(O, theta,X_k,R,alpha)
        
        for(i in 1:ncol(theta)){
            X_i <- X_k[,,i]- tao_k*dP[,,i]
            X_i <- SVT(X_i,tao_k*mean(theta_hat)*beta) ####without plugging cell type specific proportions
            X_k_1[,,i] <- X_i
            
            ratio <- c(ratio, sum( (X_i-X_k[,,i])^2 )/(nrow(X_i)*ncol(X_i)))
        }
        r1 <- loss(O,X_k,theta,alpha,beta,R) #calculating raw loss
        #print updating
        if(verbose){
			   #print matrix ratio distance or absolute distance
               print <- paste("CSE inference step ",k," \n",sep="")
			   for(c in 1:(length(ratio)-1)){
			      print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
			   }
			   print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
			   writeLines(print)
		}
        if(verbose) writeLines( sprintf("   Loss: part1=%f , part2=%f , part3=%f", r1$part1,r1$part2,r1$part3 ) )
        if(max(ratio) < epsilon||k >= max.iter){break}else{
            ## 更新X
            X_k <- X_k_1
            k <- k + 1
            if(loss_his) loss<- rbind(loss,c(r1$part1,r1$part2,r1$part3))
        }
    }
    steps <- k
    ##########
    #Doing PC or Cell Fractions Normalization
    if(Normalize){
        if(verbose) cat("Perform Normalization...")
        X_k_norm <- X_k
        if(Norm.method == "PC"){
            for(k in 1:dim(X_k)[3]){
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
        }
        
        if(Norm.method == "frac"){
            for(k in 1:dim(X_k)[3]){
                exp <- X_k[,,k]^2
                exp.norm <- NULL
                for(i in 1:nrow(exp)){
                    lm.model <- lm(exp[i,]~theta[,k])
                    exp.norm <- rbind(exp.norm, (exp[i,] - lm.model$coefficients[2] * theta[,k]))
                }
                X_k_norm[,,k] <- exp.norm
            }
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
        rownames(theta_new) <- colnames(X)
    }
    
    writeLines( paste("Converge in ",steps," steps",sep="") )
    res <- list()
    res$X_k <- X_k
    res$loss_history <- loss
    if(Normalize) res$X_k_norm <- X_k_norm
    res
}

egm_without_theta <- cell_deconvolve_trace_proximal_point_solver_final(O = as.matrix(Bulk),
                                                            theta=Frac$theta,
                                                            R=Reference,
                                                            epsilon=0.00001,
                                                            alpha=as,beta=3,tao_k = 1,
                                                            verbose=FALSE,max.iter = 1000,Normalize=FALSE)								
													
cell_deconvolve_trace_proximal_point_solver_final <- function(O, theta, R, alpha=0.5,beta=5,epsilon=0.001,tao_k=0.01,max.iter=100,X_int=NULL,verbose=FALSE,infer=FALSE,Normalize=TRUE,Norm.method = "PC",loss_his=TRUE){
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
    
    theta_hat <- colMeans(theta)
    k <- 0
    X_k_1 <- X_k <- X
    ##record loss
    if(verbose) cat(date(), 'Optimizing cell type specific expression profile... \n')
    loss <- NULL
    repeat{
        ratio <- NULL
        
        dP <- derive_P(O, theta,X_k,R,alpha)
        for(i in 1:ncol(theta)){
            X_i <- X_k[,,i]- tao_k*dP[,,i]
            X_i <- SVT(X_i,tao_k*theta_hat[i]*beta) ####without plugging cell type specific proportions
            X_k_1[,,i] <- X_i
            
            ratio <- c(ratio, sum( (X_i-X_k[,,i])^2 )/(nrow(X_i)*ncol(X_i)))
        }
        r1 <- loss(O,X_k,theta,alpha,beta,R) #calculating raw loss
        #print updating
        if(verbose){
			   #print matrix ratio distance or absolute distance
               print <- paste("CSE inference step ",k," \n",sep="")
			   for(c in 1:(length(ratio)-1)){
			      print <- paste(print, colnames(R)[c], ": ",ratio[c]," \n",sep="")
			   }
			   print <- paste(print, colnames(R)[length(ratio)], ": ",ratio[length(ratio)],sep="")
			   writeLines(print)
		}
        if(verbose) writeLines( sprintf("   Loss: part1=%f , part2=%f , part3=%f", r1$part1,r1$part2,r1$part3 ) )
        if(max(ratio) < epsilon||k >= max.iter){break}else{
            ## update X
            X_k <- X_k_1
            k <- k + 1
            if(loss_his) loss<- rbind(loss,c(r1$part1,r1$part2,r1$part3))
        }
    }
    steps <- k
    ##########
    #Doing PC or Cell Fractions Normalization
    if(Normalize){
        if(verbose) cat("Perform Normalization...")
        X_k_norm <- X_k
        if(Norm.method == "PC"){
            for(k in 1:dim(X_k)[3]){
                exp <- X_k[,,k]^2
                exp.scale <- t(apply(exp,1,scale))
                PC <- svd(exp.scale)$v[,1]
                the <- (exp %*% as.matrix(PC) - length(PC) * mean(PC) * rowMeans(exp)) / (sum(PC^2) - length(PC)*mean(PC)^2)
                exp.norm <- exp - as.matrix(the) %*% t(as.matrix(PC))
                X_k_norm[,,k] <- exp.norm
            }
        }
        
        if(Norm.method == "frac"){
            for(k in 1:dim(X_k)[3]){
                exp <- X_k[,,k]^2
                PC <- theta[,k]
                the <- (exp %*% as.matrix(PC) - length(PC) * mean(PC) * rowMeans(exp)) / (sum(PC^2) - length(PC)*mean(PC)^2)
                exp.norm <- exp - as.matrix(the) %*% t(as.matrix(PC))
                X_k_norm[,,k] <- exp.norm
            }
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
        rownames(theta_new) <- colnames(X)
    }
    
    writeLines( paste("Converge in ",steps," steps",sep="") )
    res <- list()
    res$X_k <- X_k
    res$loss_history <- loss
    if(Normalize) res$X_k_norm <- X_k_norm
    res
}

egm_with_theta <- cell_deconvolve_trace_proximal_point_solver_final(O = as.matrix(Bulk),
                                                            theta=Frac$theta,
                                                            R=Reference,
                                                            epsilon=0.00001,
                                                            alpha=as,beta=3,tao_k = 1,
                                                            verbose=FALSE,max.iter = 1000,Normalize=FALSE)
															
p <- 100
y <- gl(2, p/2)
res <- c(DEG_test1(egm_without_theta$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE),DEG_test1(egm_with_theta$X_k,y,Frac$theta,method = "enigma",10000,DEG_list,qval=FALSE))
print(res)
Tab <- rbind(Tab,res)
}

####barplot perform
res <- cbind(as.numeric(Tab),rep(alpha.set,4),c(rep("model1",5),rep("model2",5)))
res <- as.data.frame(res)
colnames(res) <- c("AUPRC","alpha","Model")
res$AUPRC <- as.numeric(as.matrix(res$AUPRC))
res$alpha <- as.character(as.matrix(res$alpha))
boxp <-ggplot(res, aes(x=alpha, y=AUPRC, fill=Model)) +
    geom_bar(stat="identity",position=position_dodge())+theme_classic2()+theme(text = element_text(size=20))+labs(x="")
png("barplot_DEGperformance(unbalance).png",res=300,height=1000,width=1800)
boxp
dev.off()
Tab_unbalance <- Tab




