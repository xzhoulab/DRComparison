
# full quatile normalization
FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
}

# z-score normalization
SeuratN <- function(counts){
library(Seurat)
library(cowplot)
# Set up control object
sce <- CreateSeuratObject(raw.data = counts, project = "test", min.cells = 5)
sce <- NormalizeData(sce)
sce <- ScaleData(sce)
return(sce@scale.data)
}

# variance stabilizing transformation
VST <- function(x,sv=1){
  varx  = apply(x,1,var)
  meanx   = apply(x,1,mean)
  phi   = coef(nls(varx~meanx+phi*meanx^2,start=list(phi=sv)))
  return(log(x+1/(2*phi)))
}

# log-cpm
LCPM <- function(x){
	library(scater)
	return(log(calculateCPM(x) + 1))
}#end func


# log-standardized
LSD <- function(x){
	lcounts <- log(x + 1)
	norm_counts <- t( scale(t(x)) )
	return(norm_counts)
}#end func
