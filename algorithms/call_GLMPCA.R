## GLM-PCA
## 2019-6-11 17:20:35

## loading packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(BiocParallel)
  library(matrixStats)
})

# ## sourcing all R files
# path.list = paste0(src.path, file.list)
# sapply(path.list, source)
## code from GitHub: https://github.com/willtownes/scrna2019
source("~/drmethods/scrna/algs/glmpca.R")
source("~/drmethods/scrna/algs/ortho.R")
source("~/drmethods/scrna/util/functions.R")
source("~/drmethods/scrna/util/functions_genefilter.R")
source("~/drmethods/scrna/util/txtparse.R")

call_GLMPCA <- function(sce, num_pc, params){
	# other parameter in pCMF method
	fam <- params$fam
	filtering_method <- params$filtering_method
	
	counts <- counts(sce)
	if(filtering_method=="nonzeros"){
		counts <- counts[which(rowSums(counts>0)>5),]
		counts <- counts[,which(colSums(counts>0)>10)]
	}
	#counts <- counts[rowSums(counts)>0,]

	tryCatch({
		ct <- system.time({ res_glmpca <- glmpca(counts, L=num_pc, fam=fam) })
		# extract the low dimension struct W
		ct <- c(user.self = ct[["user.self"]], sys.self = ct[["sys.self"]], 
            user.child = ct[["user.child"]], sys.child = ct[["sys.child"]],
            elapsed = ct[["elapsed"]])
		list(res = res_glmpca, ctimes = ct)
	},
	error = function(e) {
    list(res = structure(rep(NA, 1), ctimes = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA), name_col = colnames(sce)))
	})
}# end func
	