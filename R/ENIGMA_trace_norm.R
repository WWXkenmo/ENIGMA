#' @title ENIGMA trace norm version
#'
#' @description trace norm version is the alternative version of ENIGMA, which is the regularized weighted matrix completion to constraint the trace norm of inferred cell type-specific gene expression matrix.
#' @param object ENIGMA object
#'
#' @param alpha
#' ENIGMA is a multi-objective optimization problem involve two object function,
#' the distance function between observed bulk RNA-seq and reconstitute RNA-seq
#' generated by weighted combination of CSE, and the distance function beween
#' average CSE expression and cell type reference matrix. The alpha is used to
#' determine weights of these two objects. If the alpha gets larger,
#' the optimization attach greater importance on the the first object.
#' Default: 0.5
#'
#' @param beta
#' This parameter is used to control the latent dimension of each CSE,
#' if this parameter gets larger, than the latent dimension of each CSE is smaller
#' (lower trace norm value), which means that each sample is more similar with
#' each others. The user need to tune this parameter based on the range of
#' the singular value of the bulk RNA-seq matrix.
#' Default: 1
#' 
#' @param tao_k
#' step size for proximal point method. Default: 1
#'
#' @param gamma
#' This parameter is used to control the distance between CSE (X) and auxiliary variable (Y). Default: 1
#'
#' @param epsilon
#' In trace norm based ENIGMA, the epsilon is not necessarily choose
#' a extremly small value, the number of iteration would influence
#' the latent dimensions of CSE, as each step is performing singular value thresholding.
#' Default: 0.0001
#'
#' @param max.iter
#' The maximum number of iterations. Default: 1000
#'
#' @param solver
#' The solver for solving trace norm model. method used: admm, admm_fast or proximal point method
#'
#'
#' @param verbose
#' Whether return the information after each step of processing. Default: TRUE
#'
#' @param Normalize
#' Whether perform normalization on resulted expression profile. Default: TRUE
#'
#' @param Norm.method
#' Method used to perform normalization. User could choose PC, frac or quantile. Default: frac
#' 
#' @param preprocess
#' Method used to perform variance stablization preprocessing. User could choose sqrt, log or none.sqrt: square root transformation; log: log2(*+1) transformation; none: no data transformation.
#'
#' @param loss_his
#' save the loss value of each round of iteration.
#'
#' @param model_tracker
#' save the model in returned object
#' 
#' @param model_name
#' name of the model
#'
#'
#' @param X_int
#' initialization for CSE profiles, an array object with three dimensions (the number of genes * the number of samples * the number of cell types), if user input a matrix (the number of genes * the number of samples), each cell type would be assigned the same start matrix.
#'
#' @return ENIGMA object where object@result_CSE contains the inferred CSE profile, object@result_CSE_normalized would contains normalized CSE profile if Normalize = TRUE, object@loss_his would contains the loss values of object functions during model training. If model_tracker = TRUE, then above results would be saved in the object@model.
#'
#'
#' @examples
#' \dontrun{
#' egm = ENIGMA_trace_norm(egm,model_tracker = TRUE, Norm.method="quantile")
#' egm@result_CSE
#' egm@result_CSE_normalized
#' }
#'
#' @export
#'
ENIGMA_trace_norm <- function(object, theta, R, alpha=0.5,beta=1,tao_k=1,gamma=NULL,epsilon=NULL,max.iter=1000,solver = "admm",verbose=FALSE,pos=TRUE,Normalize=TRUE,Norm.method = "frac",preprocess = "log",loss_his=TRUE,model_tracker=FALSE,model_name = NULL,X_int=FALSE){
    suppressPackageStartupMessages(require("scater"))
	suppressPackageStartupMessages(require("preprocessCore"))
	
	###Create a model assay
	if ( !(preprocess %in% c("none", "sqrt","log")) | (length(preprocess) != 1) ) {
        stop("Invalid data transformation method type. Please input 'none','sqrt' or 'log'. ")
    }
	if ( !(Norm.method %in% c("PC", "frac","quantile")) | (length(Norm.method) != 1) ) {
        stop("Invalid normalization method type. Please input 'PC','frac' or 'quantile'. ")
    }
	if ( !(solver %in% c("proximalpoint", "admm_fast","admm")) | (length(Norm.method) != 1) ) {
        stop("Invalid solver. Please input 'proximalpoint', 'admm_fast','admm'. ")
    }
	
    if(preprocess == "sqrt") O = sqrt(object@bulk)
	if(preprocess == "log") O = log2(object@bulk+1)
	
    theta = object@result_cell_proportion
	if(preprocess == "sqrt") R = sqrt(object@ref)
	if(preprocess == "log") R = log2(object@ref+1)
    # unify geneid between O and R
    geneid = intersect( rownames(O), rownames(R) )
    O = O[geneid,]
    R = R[geneid,]

    X = array(0,
              dim = c( nrow(O),
                       ncol(O),
                       ncol(theta)),
              dimnames = list( rownames(O),
                               colnames(O),
                               colnames(theta))
    )
    X_int_m = array(0,
              dim = c( nrow(O),
                       ncol(O),
                       ncol(theta)),
              dimnames = list( rownames(O),
                               colnames(O),
                               colnames(theta))
    )
    if(is.null(X_int) == FALSE){
       if(length(dim(X_int)) == 2){
           for(i in 1:ncol(theta)){
           X_int_m[,,i] = X_int
           }
        }
       if(length(dim(X_int)) == 3){
          for(i in 1:ncol(theta)){
           X_int_m[,,i] = X_int[,,i]
          }
        }
    }
    for(i in 1:ncol(theta)){
        if(is.null(X_int)){X[,,i] <- O}else{X[,,i] <- X_int_m[,,i]}
    }
    rm(X_int, X_int_m);gc()
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

    if(model_tracker){
	if(is.null(model_name)){
	  model_name = paste("trace_norm_model_",date(),"_trained",sep="")
	}
	  if(solver == "proximalpoint"){basic_infor = data.frame(alpha = alpha,beta = beta, step_size = tao_k, epsilon = epsilon,max_iter = max.iter,Normalize = Normalize,Normalize_method = Norm.method,preprocess = preprocess,solver = solver,pos=pos)}else{
	  basic_infor = data.frame(alpha = alpha,beta = beta, gamma=gamma, epsilon = epsilon,max_iter = max.iter,Normalize = Normalize,Normalize_method = Norm.method,preprocess = preprocess,solver = solver,pos=pos)}
	  object@model[[model_name]] <- list(basic_infor = basic_infor)
	}
	
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
    loss <- NULL
    if(verbose) cat(date(), 'Optimizing cell type specific expression profile... \n')
	if(verbose){
	if(solver == "proximalpoint"){
	cat(paste('alpha: ',alpha, ' \n', 'step size: ',tao_k ,' \n','beta: ',beta ,' \n','gamma: ',gamma ,' \n','epsilon: ',epsilon,' \n',sep=""))}
	if(solver == "admm" || solver == "admm_fast"){
	cat(paste('alpha: ',alpha, ' \n', 'beta: ',beta ,' \n','gamma: ',gamma ,' \n','epsilon: ',epsilon,' \n',sep=""))
	}}
	
	writeLines(paste("Using ",solver," solver...",sep=""))
    repeat{
        if(abs(delta)<epsilon||k>=max.iter){
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
		loss <- NULL
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
	if(preprocess == "log") X_k_m <- 2^X_k - 1
	if(preprocess == "sqrt") X_k_m <- X_k^2
	if(preprocess == "none") X_k_m <- X_k
	writeLines("Normalization...")
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
	
	###save
	object@result_CSE_normalized = res2sce(X_k_norm)
	if(model_tracker){
	  object@model[[model_name]]$result_CSE_normalized = res2sce(X_k_norm)
	}
	}
	writeLines( paste("Converge in ",steps," steps",sep="") )
	# return cell type specific gene expression profile
    if(preprocess == "sqrt") object@result_CSE = res2sce(X_k^2)
	if(preprocess == "log") object@result_CSE = res2sce(2^X_k - 1)
	if(preprocess == "none") object@result_CSE = res2sce(X_k)
	##loading loss history
	if(loss_his) object@loss_his = loss
	if(model_tracker){
	   if(loss_his){
	   object@model[[model_name]]$loss_his = loss
	   }else{
	   object@model[[model_name]]$loss_his = NULL
	   }
	   if(preprocess == "sqrt") object@model[[model_name]]$result_CSE = res2sce(X_k^2)
	   if(preprocess == "log") object@model[[model_name]]$result_CSE = res2sce(2^X_k - 1)
	   if(preprocess == "none") object@model[[model_name]]$result_CSE = res2sce(X_k)
	   
	   ### import the model name and model type
	   if(nrow(object@model_name)==0){
	   m = t(as.matrix(c(model_name, "trace norm model")))
	   m = as.data.frame(m)
	   colnames(m) = c("Model Name","Model Type")
	   rownames(m) = paste0("model",1:nrow(m))
	   object@model_name = m
	   }else{
	   m = rbind(as.matrix(object@model_name),t(as.matrix(c(model_name, "trace norm model"))))
	   m = as.data.frame(m)
	   colnames(m) = c("Model Name","Model Type")
	   rownames(m) = paste0("model",1:nrow(m))
	   object@model_name = m
	   }
	}
	if(verbose) cat(date(),'Done... \n')
    return(object)
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
	require("corpcor")
	
    svd <- fast.svd(Mat)
	d <- svd$d
	d <- d - t
	d[d<0] <- 0
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