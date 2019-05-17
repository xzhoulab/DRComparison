## ZINB-WaVE
## 2019-4-2 19:12:11

## loading packages
suppressPackageStartupMessages({
  library(zinbwave)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(scRNAseq)
  library(matrixStats)
})

# main function
call_ZINBWaVE <- function(sce, num_pc, params){
	
	# other parameter in zinbwave method
	num_core <- params$num_core
	doParallel <- params$doParallel
	filtering_method <- params$filtering_method
	X <- params$X
	V <- params$V
	
	counts <- counts(sce)
	counts <- round(counts)
	if(filtering_method=="nonzeros"){
		counts <- counts[which(rowSums(counts>0)>5),]
		counts <- counts[,which(colSums(counts>0)>0)]
	}

	rownames(counts) <- NULL
	colnames(counts) <- NULL
	colData <- DataFrame(Treatment=rep(c('scRNA'), ncol(counts)),
                     row.names=colnames(counts))
	colData <- colData(sce)
	zinbwave_data <- SummarizedExperiment(assays=list(counts=counts), colData=colData)
	rm(counts)
	#rm(sce)
	tryCatch({
		if(doParallel){
			# parallel to run
			if(is.null(X)&&is.null(V)){
				ct1 <- system.time({ res_zinb <- zinbwave(zinbwave_data, K=num_pc, BPPARAM=MulticoreParam(num_core)) })
			}else if(is.null(X) &!is.null(V)){
				ct1 <- system.time({ res_zinb <- zinbwave(zinbwave_data, K=num_pc, V=V, epsilon=1000, BPPARAM=MulticoreParam(num_core)) })
			}else if(is.null(V)&!is.null(X)){
				ct1 <- system.time({ res_zinb <- zinbwave(zinbwave_data, K=num_pc, X=X, epsilon=1000, BPPARAM=MulticoreParam(num_core)) })
			}# end fi	
		}else{
			if(is.null(X)&&is.null(V)){
				ct1 <- system.time({ res_zinb <- zinbwave(zinbwave_data, K=num_pc) })
			}else if(is.null(X)&!is.null(V)){
				ct1 <- system.time({ res_zinb <- zinbwave(zinbwave_data, K=num_pc, V=V) })
			}else if(is.null(V)&!is.null(X)){
				ct1 <- system.time({ res_zinb <- zinbwave(zinbwave_data, K=num_pc, X=X) })
			}# end fi
		}# end fi
		# extract the low dimension struct W
		ct2 <- system.time({ W <- reducedDim(res_zinb) })
		ct <- ct1 + ct2
		ct <- c(user.self = ct[["user.self"]], sys.self = ct[["sys.self"]], 
            user.child = ct[["user.child"]], sys.child = ct[["sys.child"]],
            elapsed = ct[["elapsed"]])
		list(res = res_zinb, ctimes = ct)
	},
	error = function(e) {
    list(res = structure(rep(NA, 1), ctimes = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,elapsed = NA), name_col = colnames(sce)))
	})
}# end func
	