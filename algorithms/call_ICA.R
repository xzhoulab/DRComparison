## Idenpendt Component Analysis,ICA
## 2019-4-2 21:08:49

## loading packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ica)
  #library(sctransform)
  #library(matrixStats)
  #library(magrittr)
  #library(ggplot2)
})

# main function
call_ICA <- function(sce, num_pc, params){
	# other parameter in FA method
	rotate <- params$rotate
	filtering_method <- params$filtering_method
	norm_method <- params$norm_method
	## normalized counts
	if(norm_method=="log"){
		norm_counts <- logcounts(sce)
	}else{
		norm_counts <- counts(sce)
	}#end func
	
	if(filtering_method=="nonzeros"){
		norm_counts <- norm_counts[which(rowSums(norm_counts>0)>5),]
		norm_counts <- norm_counts[, which(colSums(norm_counts>0)>10)]
	}
	if(norm_method=="vst"){
		norm_counts <- VST(norm_counts)
	}else if(norm_method=="lcpm"){
		norm_counts <- LCPM(norm_counts)
	}else if(norm_method=="lsd"){
		norm_counts <- LSD(norm_counts)
	}else if(norm_method=="seurat"){
		norm_counts <- SeuratN(norm_counts)
	}else if(norm_method=="fq"){
		norm_counts <- FQnorm(norm_counts)
	}
	#rm(sce)
	norm_counts <- t(norm_counts)
	tryCatch({
		ct <- system.time({
		#norm_counts <- selectFeatures_IQR(norm_counts, nb=1000)
		res_ica <- icafast(norm_counts, num_pc)
		#res_ica$pca_comp <- res_ica$X %*% res_ica$K
		res_ica$ica_comp <- res_ica$S
		})	
		# count time
		ct <- c(user.self = ct[["user.self"]], sys.self = ct[["sys.self"]], 
            user.child = ct[["user.child"]], sys.child = ct[["sys.child"]],
            elapsed = ct[["elapsed"]])
		return( list(res = res_ica, ctimes = ct) )
	},
	error = function(e) {
    list(res = structure(rep(NA, 1), ctimes = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,elapsed = NA), name_cell = colnames(sce)))
	})
}# end func

## example code
# S <- matrix(runif(10000), 500, 20)
# A <- matrix(c(1, 1, -1, 3), 20, 20, byrow = TRUE)
# X <- S %*% A
# a <- fastICA(X, 5, alg.typ = "parallel", fun = "logcosh", alpha = 1,
# method = "C", row.norm = FALSE, maxit = 200,
# tol = 0.0001, verbose = TRUE)
# par(mfrow = c(1, 3))
# plot(a$X, main = "Pre-processed data")
# plot(a$X %*% a$K, main = "PCA components")
# plot(a$S, main = "ICA components")
