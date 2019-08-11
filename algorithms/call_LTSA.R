## LTSA: Local Tangent Space Alignment
## 2019-6-28 18:05:30

## loading packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(psych)
  library(Rdimtools)
})

FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
}


SeuratN <- function(counts){
library(Seurat)
library(cowplot)
# Set up control object
sce <- CreateSeuratObject(raw.data = counts, project = "IMMUNE_CTRL", min.cells = 5)
sce <- NormalizeData(sce)
sce <- ScaleData(sce)
return(sce@scale.data)
}
VST <- function(x,sv=1){
  varx  = apply(x,1,var)
  meanx   = apply(x,1,mean)
  phi   = coef(nls(varx~meanx+phi*meanx^2,start=list(phi=sv)))
  return(log(x+1/(2*phi)))
}

LCPM <- function(x){
	library(scater)
	return(log(calculateCPM(x) + 1))
}#end func

LSD <- function(x){
	lcounts <- log(x + 1)
	norm_counts <- t( scale(t(x)) ) # scale is a R generic func
	return(norm_counts)
}#end func

call_LTSA <- function(sce, num_pc, params){
	## other parameter in LLE method
	filtering_method <- params$filtering_method
	num_prop <- params$num_prop
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
	tryCatch({
		ct <- system.time({
		res_ltsa <- do.ltsa(t(norm_counts), ndim = num_pc, type = c("proportion", num_prop))
		})	
		# count time
		ct <- c(user.self = ct[["user.self"]], sys.self = ct[["sys.self"]], 
            user.child = ct[["user.child"]], sys.child = ct[["sys.child"]],
            elapsed = ct[["elapsed"]])
		list(res = res_ltsa, ctimes = ct)
	},
	error = function(e) {
    list(res = structure(rep(NA, 1), ctimes = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA), name_col = colnames(sce)))
	})
}# end func

