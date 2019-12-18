
## extract the low-embedding
ExtractW <- function(res, method){
	if(method %in% c("ZIFA","scScope","DCA","umap","scvis")){
        W = res
    }else if (method == "FA"){
        W = res$res$loadings[, 1:ncol(res$res$loadings)]
    }else if(method %in% c("PCA", "MDS","DiffusionMap","tSNE")){
        W = res$res
    }else if(method == "pCMF"){
        W = res$res$factor$U
    }else if(method == "ZINBWaVE"){
        W = reducedDim(res$res)
    }else if (method == "ICA"){
        W = res$res$ica_comp
    }else if (method %in% c("PoissonNMF", "NMF")){
        W = res$res$W
    }else if (method %in% c("GLMPCA")){
        W = res$res$factors
    }else if (method %in% c("Isomap", "LLE", "LTSA")){
        W = res$res$Y
    }# end func
	return(W)
}# end func


##
ClusteringCells <- function(W, num_clust, method='all'){
    set.seed(0)
	res <- list()
	if(method=="kmeans"){
		res[[1]] <- kmeans(W, num_clust, nstart=1000, iter.max=100000)$cluster
		names(res) <- "kmeans"
		return(res)
	}else if(method=="gmm"){#Gaussian mixture model, GMM
		library(mclust)
		res[[1]] <- MUDAN::getComMembership(W, k=num_clust, method=igraph::cluster_infomap,verbose=FALSE)
		names(res) <- "commdetec"
		return(res)
	}else if(method=="hclust"){
		library(dynamicTreeCut)
		W <- scale(W)
		d <- as.dist(1-cor(t(W)))
		hc <- hclust(d, method='ward.D')
		res[[1]] <- cutree(hc, num_clust)
		names(res) <- "hclust"
		return(res)
	}else if(method=="snncliq"){
		## load the package and run main function snn_cliq
		##devtools::load_all("/net/mulan/disk2/shiquans/scRNAseqDRComparison/sctools/")
		dist_mat <- as.matrix(dist(W))
		res[[1]] <- sctools::snn_cliq(dist_mat, knn=30)$cell_labels
		names(res) <- "snncliq"
		return(res)
	}else if(method=="all"){
		res[[1]] = kmeans(W, num_clust, nstart=1000, iter.max=100000)$cluster
		library(dynamicTreeCut)
		d <- as.dist(1-cor(t(W)))
		hc <- hclust(d, method='ward.D')
		res[[2]] <- cutreeDynamic(hc, distM=as.matrix(d), deepSplit=num_clust, verbose=FALSE)
		names(res) <- c("kmeans","hclust")
		return(res)
	}# end fi
}# end function

EvaluationPerf <- function(pred_label, true_label, method="clustering"){
	clust_combined <- NULL
	for(iclust in 1:length(pred_label)){
		library(igraph)
		res <- rep(NA,2)
		if(method=="clustering"){
			# normalized mutual information
			nmi <- compare(true_label, pred_label[[iclust]], method = "nmi")
			# adjusted random index
			ari <- compare(true_label, pred_label[[iclust]], method = "adjusted.rand")
			res <- c(nmi, ari)
			names(res) <- c("nmi", "ari")
		}else if(method=="fmeasure"){
			## f-measure
			library("MLmetrics")
			f1 <- FMeasure(true_label, pred_label[[iclust]])
			res <- c(f1)
			names(res) <- c("fmeasure")
			return(res)
		}# end fi
		clust_combined <- cbind(clust_combined, res)
	}# end for
	colnames(clust_combined) <- names(pred_label)
	return(clust_combined)
}# end func




