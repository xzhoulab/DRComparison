## Poisson Nonegative Matrix Factorization, Poisson NMF
## 2019-4-2 21:08:49

## loading packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(psych)
  library(NNLM) #https://github.com/linxihui/NNLM/blob/master/R/nnmf.R
})

# main function
call_PoissonNMF <- function(sce, num_pc, params){
	## other parameter in NMF method
	filtering_method <- params$filtering_method
	num_core <- params$num_core
	
	## normalized counts
	counts <- counts(sce)
	
	if(filtering_method=="nonzeros"){
		counts <- counts[which(rowSums(counts>0)>5),]
		counts <- counts[, which(colSums(counts>0)>10)]
	}
	counts <- t(counts)
	tryCatch({
		ct <- system.time({
		res_nmf <- nnmf(counts, num_pc, method = "lee",loss = "mkl",n.threads = num_core)
		})	
		# count time
		ct <- c(user.self = ct[["user.self"]], sys.self = ct[["sys.self"]], 
            user.child = ct[["user.child"]], sys.child = ct[["sys.child"]],
            elapsed = ct[["elapsed"]])
		list(res = res_nmf, ctimes = ct)
	},
	error = function(e) {
    list(res = structure(rep(NA, 1), ctimes = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA, elapsed = NA), name_col = colnames(sce)))
	})
}# end func

