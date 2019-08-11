########################################################################################################################
# Package: PQL Matrix Factorization
# Version: 1.0.1
# Date   : 2018-07-23
# Title  : Penalized Quasi-Likelihood Based Matrix Factorization on Single Cell RNAseq Data
# Authors: S.Q. Sun, Q.D. Feng, and X. Zhou
# Contact: shiquans@umich.edu and qidif@umich.edu 
#          University of Michigan, Department of Biostatistics
########################################################################################################################
#' @param RawCountDataSet count data.
#' @param Phenotypes the phenotype we are interested, such sex, age, disease status.
#' @param Covariates some fixed effect we are unwanted.
#' @param LibSize, total count for each sample
#' @param fit.model, two options, PMM.
#' @export
pqlmfmk <- function(RawCountDataSet, Covariates=NULL, RelatednessMatrix=NULL, LibSize=NULL, numPC=2, fit.model="PMM", fit.method="AI.REML", fit.maxiter=500, fit.tol=1e-5, numCore=1, filtering=TRUE, verbose=TRUE, ...) {
	if(numCore > 1){
		if(numCore>detectCores()){warning("PQLMFMK:: the number of cores you're setting is larger than detected cores!");numCore = detectCores()}
	}# end fi
	registerDoParallel(cores=numCore)
	
	# filtering counts
	if (filtering & fit.model == "PMM"){
	unfilterIdx <- apply( RawCountDataSet, 1, function(x) length(x[x>5])>=2 )
	CountData   <- RawCountDataSet[unfilterIdx,]
	}else{
		CountData <- RawCountDataSet
	}
	rm(RawCountDataSet)
    
	numVar <- dim(CountData)[1]
	numIDV <- dim(CountData)[2]
	
	# remove the intercept
	if(length(unique(Covariates[,1])) == 1){
		Covariates<-Covariates[,-1]
	}# end fi

	if(is.null(Covariates)){
		numCov <- 0
	}else{
		numCov     <- dim(Covariates)[2]
		Covariates <- as.matrix(Covariates)
	}# end fi
  
	cat(paste0("### the information of inputs:\n"))
	cat(paste0("## number of total individuals: ", numIDV,"\n"))
	cat(paste0("## number of total genes: ", numVar,"\n"))
	cat(paste0("## number of adjusted covariate",ifelse(numCov>0,"s",""),": ", numCov+1,"\n"))
	cat(paste0("## number of principal component",ifelse(numPC>0,"s",""),": ", numPC,"\n"))
  
	CountData  <- as.matrix(CountData)
  
  
	if(is.null(RelatednessMatrix)){
		stop("PQLMFMK::please input cell-relatedness matrix!")
	}else if(!class(RelatednessMatrix) %in% c("matrix", "list")){
		stop("PQLMFMK: \"RelatednessMatrix\" must be a matrix or a list.")
	}
	if(class(RelatednessMatrix) == "list" && length(table(sapply(RelatednessMatrix, dim))) != 1)
	stop("PQLMFMK: when \"RelatednessMatrix\" is a list, all its elements must be square matrices of the same size.")
	
	#####################
	# relatedness matrix should be positive definition
	if(is.matrix(RelatednessMatrix)){
		RelatednessMatrix <- as.matrix(RelatednessMatrix)
		eig               <- eigen(RelatednessMatrix)
		eigval            <- eig$value
		eigvector         <- eig$vectors
		if(any(eigval<1e-10)){ 
			warning("PQLMFMK::the relatedness matrix is singular, it has been modified!")
			# eigval[eigval<1e-10] <- runif(sum(eigval<1e-10), 0, 0.5)
			# RelatednessMatrix    <- eigvector%*%diag(eigval)%*%t(eigvector)
			#RelatednessMatrix <- as.matrix(nearPD(RelatednessMatrix, corr=T)$mat)	
			RelatednessMatrix <- as.matrix(nearPD(RelatednessMatrix)$mat)	
		}

		rm(eig)
		rm(eigval)
		rm(eigvector)
	}else if(is.list(RelatednessMatrix)){
		for(igrm in 1:length(RelatednessMatrix)){
			print(igrm)
			tmpGRM <- RelatednessMatrix[[igrm]]
			tmpGRM <- as.matrix(tmpGRM)
			eig               <- eigen(tmpGRM)
			eigval            <- eig$value
			eigvector         <- eig$vectors
			if(any(eigval<1e-10)){ 
				warning("PQLMFMK::the relatedness matrix is singular, it has been modified!")
				RelatednessMatrix[[igrm]] <- as.matrix(nearPD(tmpGRM)$mat)	
			}#end fi

			rm(eig)
			rm(eigval)
			rm(eigvector)
		}#end for
	}#end fi

	RelatednessMatrix <- c(RelatednessMatrix, list(diag(numIDV)))
	numKMat <- length(RelatednessMatrix)
	#***********************************#
	#       Poisson Mixed Model         #
	#***********************************#
	if(fit.model=="PMM"){
		
		if(is.null(LibSize)){
			LibSize <- apply(CountData, 2, sum)
			LibSize <- as.matrix(LibSize)
		}else{
			LibSize <- as.matrix(LibSize)
		}
	 	 
		#=============================================================
		# initialize W, alpha, Y, D using parallel
		y <- NULL  # dim = n x p
		Y <- NULL  # dim = n x p
		D <- NULL  # dim = n x p
		W <- NULL  # dim = n x c
		Alpha <- NULL # dim = p x c
		# do parallel using foreach function
		# resPMM <-foreach(iVar=1:numVar,.combine=cbind)%dopar%{
		for(iVar in 1:numVar){ 
			if(numCov==0){
				model0 <- glm(formula = CountData[iVar,]~1 + offset(log(LibSize)), family = poisson(link="log"))
				#idx   <- match(rownames(model.frame(formula = CountData[iVar,]~1 + offset(log(LibSize)), na.action = na.omit)), rownames(model.frame(formula = CountData[iVar,]~1 + offset(log(LibSize)), na.action = na.pass)))
			}else{
				model0 <- glm(formula = CountData[iVar,]~1 + Covariates + offset(log(LibSize)), family = poisson(link="log"))
				#idx   <- match(rownames(model.frame(formula = CountData[iVar,]~ 1 + Covariates + offset(log(LibSize)), na.action = na.omit)), rownames(model.frame(formula = CountData[iVar,]~ 1 + Covariates + offset(log(LibSize)), na.action = na.pass)))
			}# end fi
			
			res_init <- ComputeYAD(model0)
			
			# raw counts
			y <- cbind(y, res_init$y)
			# working vector
			Y <- cbind(Y, res_init$Y)
			# residual variance
			D <- cbind(D, res_init$D)	
			# covariate coefficient
			Alpha <- rbind(Alpha, res_init$alpha)
			
		}# end for iVar, parallel
		
		# covariates denoted by W in the model
		W <- model.matrix(model0)
		
		#=========================================================
		# main loop
		# relatedness matrix
		#cat(paste("length of relatedness matrix",length(RelatednessMatrix),"\n"))
		names(RelatednessMatrix) <- paste0("kins", 1:length(RelatednessMatrix))
		
		cat("# fitting Poisson mixed model ... \n")
		# iteration via AI to slove the variance components tau0, tau1, and tau2.
		if((fit.method=="AI.REML")){
			fixtau.old <- rep(0, length(RelatednessMatrix)+1)
			# to use average information method to fit alternative model
			model1 <- ComputeVC(y, Y, W, D, Alpha, RelatednessMatrix, model0, maxiter=fit.maxiter, tol=fit.tol)
			fixtau.new <- 1*(model1$tau < 1.01*fit.tol)

			while(any(fixtau.new != fixtau.old)){
				fixtau.old <- fixtau.new
				# to use average information method to fit alternative model
				model1 <- ComputeVC(y, Y, W, D, Alpha, RelatednessMatrix, model0, fixtau=fixtau.old, maxiter=fit.maxiter, tol=fit.tol)
				fixtau.new <- 1*(model1$tau < 1.01*fit.tol)
			}# end while
		}# end fi
		
		cat("# fitting generalized PCA model ... \n")
		X <- sweep(model1$Eta, 2, model1$Alpha, "-")	
		# generalized pca
		Sigma <- 0
		for(ik in 1:(numKMat-1) ){Sigma <- Sigma + model1$tau[ik]*RelatednessMatrix[[ik]]}
		#model2 <- AdaptiveGPCA(t(as.matrix(X)), as.matrix(Sigma), model1$tau, k=numPC)
		model2 <- AdaptiveGPCA(t(as.matrix(X)), as.matrix(Sigma), tau=c(1.0, model1$tau[numKMat]), k=numPC)
		
		#print(model1$tau)
		#print(numPC)
		# return the results
		#return(list(PQL=model1, PCA=model2, kmat=RelatednessMatrix$kins1))
		return(list(PQL=model1, PCA=model2))
	}# end PMM 
}# end function PQLMF

##########################
# initialization
ComputeYAD <- function(model0){
	if(model0$family$family %in% c("binomial")){
		y <- model0$numSucc
	}else{
		y <- model0$y
	}
	num_cell <- length(y)
	offset <- model0$offset
	if(is.null(offset)) {offset <- rep(0, num_cell)}
	
	family <- model0$family
	eta <- model0$linear.predictors
	mu <- model0$fitted.values
	mu.eta <- family$mu.eta(eta)
    D <- mu.eta/sqrt(model0$family$variance(mu))
  
	# if(family$family %in% c("binomial")){
	  # mu.eta <- model0$numTotal*mu.eta
	  # D <- mu.eta/sqrt(model0$numTotal*model0$family$variance(mu))
	  # mu <- model0$numTotal*mu
	# }

	# working vector for each gene
	Y <- eta - offset + (y - mu)/mu.eta	
	
	alpha <- model0$coef
	
	return(list(y=y, D=D, Y=Y, alpha=alpha))
}# end function ComputeYAD


##########
ComputeVC <- function(y, Y, W, D, Alpha, kmat, model0, tau=rep(0, length(kmat)+1), fixtau=rep(0, length(kmat)+1), maxiter=500, tol=1e-5){
	
	numIDV <- nrow(Y)
	offset <- model0$offset
	if(is.null(offset)) {offset <- rep(0, numIDV)}
	family <- model0$family
	# estimate the variance components
	if(family$family %in% c("poisson", "binomial")) {
		tau[1] <- 1
		fixtau[1] <- 1
	}# end fi
	
	#######################
	# intialize the tau
	idxtau <- which(fixtau == 0)
	numK2 <- sum(fixtau == 0)
	if(numK2 > 0){
	    tau[fixtau == 0] <- rep(min(0.9, median(apply(Y,2,var))/(length(kmat)+1)), numK2) 
	}
	# main loop function to estimate the variance components
	# the variable Y, D, and Alpha still need to be updated 
	for(iter in seq_len(maxiter)){
		# initialize the tau0
		tau0 <- tau
		# update alpha, tau, and working vector Y
		model1 <- AI_cpp(Y, W, length(kmat), kmat, D^2, tau, fixtau, tol)
		
		tau <- as.numeric(model1$tau)
		#cov <- as.matrix(model1$cov)
		Alpha <- as.numeric(model1$Alpha)
		#Eta <- as.numeric(model1$Eta) + offset
		Eta <- sweep(model1$Eta, 1, offset, "+")
		## working vector for all genes
		# for(i in 1:ncol(Eta)){
		# 	mu <- family$linkinv(Eta[,i])
		# 	mu.eta <- family$mu.eta(Eta[,i])
		#	D[,i] <- mu.eta/sqrt(family$variance(mu))
		## working vector
		# 	Y[,i] <- Eta[,i] - offset + (y[,i] - mu)/mu.eta
		# }# end for loop
		mu <- apply(Eta, 2, function(x){family$linkinv(x)})
		mu.eta <- apply(Eta, 2, function(x){family$mu.eta(x)})
		mu.var <- apply(mu, 2, function(x){family$variance(x)})
		D <- mu.eta/sqrt(mu.var)
		# working vectors
		Y <- Eta + (y - mu)/mu.eta
		Y <- sweep(Y, 1, offset,"-")
		
		# stop criterion
		if(abs(tau-tau0)/(abs(tau)+abs(tau0)+tol)<tol){break}
		
	}# end for loop
	converged <- ifelse(iter<maxiter, TRUE, FALSE)
	
	return(list(Y=Y, Alpha=Alpha, D=D, tau=tau, Eta=Eta, converged = converged))
}# end function ComputeVC


# InitialTau <- function(kmat, X, Y, D){
	# numIDV <- nrow(kmat)
	# numK <- length(kmat)
	# # the components need to be estimated
	# idxtau <- which(fixtau == 0)
	# numK2 <- sum(fixtau == 0)
	
	# # initialize the parameters tau1, and tau2
	# if(numK2 > 0){
	    # tau[fixtau == 0] <- rep(min(0.9, var(Y)/(numK+1)), numK2)  
		# # tau0 * V0
		# H <- tau[1]*diag(1/D^2)
		
		# # tau1*V1 + tau2*V2
		# for(ik in 1:numK) {H <- H + tau[ik+1]*kmat[[ik]]}
	
		# # compute the P matrix
		# Hinv <- chol2inv(chol(H))
		# HinvX <- crossprod(Hinv, X)
		# XHinvX <- crossprod(X, HinvX)
		# P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))
	
		# if(class(P) == "try-error"){
			# stop("Error in P matrix calculation!")
		# }

		# PY <- crossprod(P, Y)
		# tau0 <- tau
		# for(ik in 1:numK2){
		    # if(ik==1 && fixtau[1]==0){# initialize the tau0
				# tau[1] <- max(0, tau0[1] + tau0[1]^2 * (sum((PY/D)^2) - sum(diag(P)/D^2))/numIDV)
			# }else{# initialize tau1 and tau2
				# PAPY <- crossprod(P, crossprod(kmat[[idxtau[ik]-1]], PY))
				# tau[idxtau[ik]] <- max(0, tau0[idxtau[ik]] + tau0[idxtau[ik]]^2*(crossprod(Y, PAPY) - sum(P*kmat[[idxtau[ik]-1]]))/numIDV)
			# }# end fi
		# }# end for
	# }# end fi
# }# end function InitialTau
##########################################################################################
# test example 
# rm(list=ls())

# count_data <- read.csv("/net/mulan/shiquans/dataset/SingleCellRNAseq/Chuetal/GSE75748_sc_cell_type_ec.csv.gz",header=T,row.names=1)
# col_names <- colnames(count_data)
# split_col_names <- strsplit(col_names, "_")

# cell_type <- c()
# for(i in 1:length(col_names)){
        # cell_type <- c(cell_type, split_col_names[[i]][1])
# }
# cell_type <- as.numeric(as.factor(cell_type))

# test_idx1 <- which(cell_type==1)
# test_idx2 <- which(cell_type==2)
# sample_size <- length(test_idx1) + length(test_idx2)

# pheno <- rep(0,sample_size)
# pheno[1:length(test_idx1)] <- 1

# #########################
# count_data <- count_data[,c(test_idx1, test_idx2)]

# row_sum <- apply(count_data,1,sum)
# count_data <- count_data[which(row_sum>3000),]
# LibSize <- round(apply(count_data, 2, sum))
# normal_counts <- 1e+06*sweep(count_data, 2, LibSize, "/")
# log_normal_counts <- log(normal_counts + 0.1)

# RelatednessMatrix <- t(as.matrix(log_normal_counts))%*%as.matrix(log_normal_counts)/nrow(log_normal_counts)
# kmat <- RelatednessMatrix
# CountData <- round(as.matrix(count_data[1:10,]))
# numVar <- nrow(CountData)
# numIDV <- ncol(CountData)
# numCov <- 0

# RelatednessMatrix <- list(RelatednessMatrix, diag(numIDV))
# names(RelatednessMatrix) <- paste0("kins", 1:length(RelatednessMatrix))
# kmat <- RelatednessMatrix

# tau=rep(0, length(kmat)+1)
# tau[1] <-1
# fixtau=rep(0, length(kmat)+1)
# fixtau[1] <- 1
# tol=1e-5
# library(Rcpp)
# sourceCpp("/net/mulan/shiquans/Projects/SingleCellRNAseq/SingleCellRNAseq_PQL/scPQLcodeing/PQLMF.cpp")

#############################################
# original code is from "adaptiveGPCA" R package
#' variables. See \href{https://arxiv.org/abs/1702.00501}{Fukuyama,
#' J. (2017)} for more details.
#' @export
AdaptiveGPCA <- function(X, Q, tau, k=2, weights=rep(1, nrow(X))) {
    if(is.matrix(Q)){
        Qeig <- eigen(Q, symmetric = TRUE)
    }else if(is.list(Q) & !is.null(Q$vectors) & !is.null(Q$values)) {
        Qeig <- Q
    }else{
        stop("Q is not formatted correctly")
    }# end fi
    # normalize so that the trace of Q is the same as the trace of the identity matrix
    Qeig$values <- ncol(X) * Qeig$values / sum(Qeig$values)
    evecs <- Qeig$vectors
    #auto = estimateComponents(X, Q, Qeig = Qeig)
	#r = auto$r
    r <- tau[1]/(tau[1]+tau[2])
    #sigma2 = auto$sigma^2
    sigma2 <- tau[1]+tau[2]
    evals <- (rep(1/(sigma2 * (1-r)), ncol(X)) + (sigma2 * r)^(-1) * Qeig$values^(-1))^(-1)
    out.gpca <- GPCAEvecs(X, evecs, evals, weights, k)
    
	# return results
    return(list(V = out.gpca$V, U = out.gpca$U, QV = out.gpca$QV, lambda = out.gpca$lambda, vars = out.gpca$vars,
               r = r, evals = evals, sigma2 = sigma2))
}# end function AdaptiveGPCA

#' gPCA using pre-computed eidendecomposition
#' Performs gPCA with pre-computed eigenvectors and eigenvalues.
#' @param X Data matrix.
#' @param evecs Eigenvectors of \code{Q}, the inner product/similarity
#' matrix.
#' @param evals Eigenvalues of \code{Q}. 
#' @param D Sample weights
#' @param k The number of components to return.
#' @keywords internal
GPCAEvecs <- function(X, evecs, evals, D=rep(1, nrow(X)), k) {
    J <- X %*% evecs
    J <- sweep(J, 2, STATS = sqrt(evals), FUN = "*")
    J <- sweep(J, 1, STATS = D, FUN = "*")
    Jsvd <- svd(J, nv = k, nu = k)
    evalsginvsqrt <- evals^(-.5)
    evalsginvsqrt[evals == 0] <- 0
    V <- evecs %*% sweep(Jsvd$v, 1, STATS = evalsginvsqrt, FUN = "*")
    U <- sweep(Jsvd$u, 1, D^(-.5), FUN = "*")
    QV <- evecs %*% (sweep(t(evecs), 1, STATS = evals, FUN = "*") %*% as.matrix(V))
    colnames(U) = colnames(V) = colnames(QV) = paste("Axis", 1:ncol(U), sep = "")
	
	#return results
    return(list(V = V, U = U, QV = QV, lambda = Jsvd$d, vars = Jsvd$d^2 / sum(Jsvd$d^2)))
}# end function GPCAEvecs

####################################################################################################
# multiple kernels

############

#' @export
calculate_distance <- function(data, method){
    return(if (method == "spearman") {
        as.matrix(1 - cor(data, method = "spearman"))
    } else if (method == "pearson") {
        as.matrix(1 - cor(data, method = "pearson"))
    } else {
        ED2(data)
    })
}

#' Distance matrix transformation
#'
#' All distance matrices are transformed using either principal component 
#' analysis (PCA) or by calculating the 
#' eigenvectors of the graph Laplacian (Spectral). 
#' The columns of the resulting matrices are then sorted in 
#' descending order by their corresponding eigenvalues.
#'
#' @param dists distance matrix
#' @param method transformation method: either 'pca' or
#' 'laplacian'
#' @return transformed distance matrix
#'
#' @importFrom stats prcomp cmdscale
#' @export
transformation <- function(dists, method) {
    if (method == "pca") {
        t <- prcomp(dists, center = TRUE, scale. = TRUE)
        return(t$rotation)
    } else if (method == "laplacian") {
        L <- norm_laplacian(dists)
        l <- eigen(L)
        # sort eigenvectors by their eigenvalues
        return(l$vectors[, order(l$values)])
    }
}# end function


#' Calculate consensus matrix
#'
#' Consensus matrix is calculated using the Cluster-based Similarity 
#' Partitioning Algorithm (CSPA). For each clustering solution a binary 
#' similarity matrix is constructed from the corresponding cell labels: 
#' if two cells belong to the same cluster, their similarity is 1, otherwise 
#' the similarity is 0. A consensus matrix is calculated by averaging all 
#' similarity matrices.
#'
#' @param clusts a matrix containing clustering solutions in columns
#' @return consensus matrix
#' 
#' @useDynLib
#' @importFrom Rcpp sourceCpp
#' @export
ComputeConsensusMatrix <- function(clusts){
	# Rcpp function
    res <- consmx(clusts)
    colnames(res) <- as.character(c(1:nrow(clusts)))
    rownames(res) <- as.character(c(1:nrow(clusts)))
    return(res)
}# end function

#' Reindex cluster labels in ascending order
#' 
#' Given an \code{\link[stats]{hclust}} object and the number of clusters \code{k}
#' this function reindex the clusters inferred by \code{cutree(hc, k)[hc$order]}, so that
#' they appear in ascending order. This is particularly useful when plotting
#' heatmaps in which the clusters should be numbered from left to right.
#' 
#' @param hc an object of class hclust
#' @param k number of cluster to be inferred from hc
#'
#' @importFrom stats cutree
#' 
#' @examples
#' hc <- hclust(dist(USArrests), 'ave')
#' cutree(hc, 10)[hc$order]
#' reindex_clusters(hc, 10)[hc$order]
#' 
#' @export
reindex_clusters <- function(hc, k) {
    clusts <- stats::cutree(hc, k)
    labels <- names(clusts)
    names(clusts) <- 1:length(clusts)
    ordering <- clusts[hc$order]
    new.index <- NULL
    j <- 1
    for (i in unique(ordering)) {
        tmp <- rep(j, length(ordering[ordering == i]))
        names(tmp) <- names(ordering[ordering == i])
        new.index <- c(new.index, tmp)
        j <- j + 1
    }
    clusts <- new.index
    clusts <- clusts[order(as.numeric(names(clusts)))]
    names(clusts) <- labels
    return(clusts)
}


#' Calculate consensus matrix.
#' 
#' This function calculates consensus matrices based on the clustering solutions
#' @param object an object 
#' 
#' @return an list
#' 
ComputeConsens <- function(combined_idx, num_clusts=NULL){
    cluster_idx <- combined_idx
	object <- list()
    if(is.null(cluster_idx)) {
        stop(paste0("please run clustering method first!"))
        return(object)
    }else if(!is.matrix(cluster_idx)){
		cluster_idx <- as.matrix(cluster_idx)
	}else{# end fi
		# NULLing the variables to avoid notes in R CMD CHECK
		
        dat <- ComputeConsensusMatrix(cluster_idx)
        tmp <- ED2(dat)
        colnames(tmp) <- as.character(colnames(dat))
        rownames(tmp) <- as.character(colnames(dat))
        diss <- stats::as.dist(as.matrix(stats::as.dist(tmp)))
        hc <- stats::hclust(diss)
        if(!is.null(num_clusts)){
			clusts <- reindex_clusters(hc, num_clusts)
			silh <- cluster::silhouette(clusts, diss)
			return( list(dat=dat, hc=hc, cluster=clusts, silhouette=silh) )
		}else{
			return( list(dat=dat, hc=hc) )
		}# end fi
    }# end fi
}# end function


