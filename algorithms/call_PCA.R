## PCA
## 2019-4-2 21:08:49
## loading packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(BiocParallel)
})

# main function
call_PCA <- function(sce, num_pc, params){
	## other parameter in PCA method
	#num_clust <- params$num_clust
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
	}# end fi
	
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
	#
	#rm(sce)
	tryCatch({
		ct <- system.time({
		res_pca <- prcomp(norm_counts, center = TRUE, scale. = FALSE)
		res_pca <- res_pca$rotation[, seq_len(num_pc), drop = FALSE]
		})	
		# count time
		ct <- c(user.self = ct[["user.self"]], sys.self = ct[["sys.self"]], 
            user.child = ct[["user.child"]], sys.child = ct[["sys.child"]],
            elapsed = ct[["elapsed"]])
		list(res = res_pca, ctimes = ct)
	},
	error = function(e) {
    list(res = structure(rep(NA, 1), ctimes = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA), name_col = colnames(sce)))
	})
}# end func
	