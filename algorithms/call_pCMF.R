## pCMF
## 2019-4-2 19:12:11

## loading packages
suppressPackageStartupMessages({
  library(pCMF)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(matrixStats)
})

# main function
call_pCMF <- function(sce, num_pc, params){
	# other parameter in pCMF method
	num_core <- params$num_core
	doParallel <- params$doParallel
	filtering_method <- params$filtering_method
	
	counts <- counts(sce)
	if(filtering_method=="nonzeros"){
		counts <- counts[which(rowSums(counts>0)>5),]
		counts <- counts[,which(colSums(counts>0)>10)]
	}
	#counts <- counts[rowSums(counts)>0,]
	#rm(sce)
	counts <- t(counts) ## for pCMF, the dimension of data should be n x p instead of p x n
	tryCatch({
		if(doParallel){
			# parallel to run
			ct1 <- system.time({ res_pcmf <- pCMF(counts, K=num_pc, verbose=FALSE, ncores=num_core) })
		}else{
			ct1 <- system.time({ res_pcmf <- pCMF(counts, K=num_pc, verbose=FALSE, ncores=1) })
		}# end fi
		
		# extract the low dimension struct W
		ct2 <- system.time({ 
			#W <- getU(res_pcmf) 
			W <- getV(res_pcmf) 
		})
		
		ct <- ct1 + ct2
		ct <- c(user.self = ct[["user.self"]], sys.self = ct[["sys.self"]], 
            user.child = ct[["user.child"]], sys.child = ct[["sys.child"]],
            elapsed = ct[["elapsed"]])
		list(res = res_pcmf, ctimes = ct)
	},
	error = function(e) {
    list(res = structure(rep(NA, 1), ctimes = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA), name_col = colnames(sce)))
	})
}# end func
	