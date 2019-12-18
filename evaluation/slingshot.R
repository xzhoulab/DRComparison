rm(list=ls())
library(SingleCellExperiment)
library(slingshot)


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
}

datapath    = "~/datasets/gti_full/"
respath     = "~/scRNAseqDRComparison/results.TI/"
outpath     = "~/slingshot/"

imeth       <- 1
ipc         <- 2
irpt        <- 1



methods = c("PCA","DiffusionMap","ICA","FA","NMF","PoissonNMF","pCMF","ZINBWaVE","ZIFA","DCA","tSNE","umap","Isomap","MDS","LLE","LTSA","GLMPCA","scvis")


sce <- readRDS(paste0(datapath,"/gti_full_Hayashi.rds"))

countdata <- t(sce$counts)
countdata <- countdata[which(rowSums(countdata>0)>5),]
countdata <- countdata[, which(colSums(countdata>0)>10)]

sim <- SingleCellExperiment(assays = list(counts = countdata))

## obtain the dimension reduction object
if(methods[imeth]%in%c("DiffusionMap","FA","ICA","PoissonNMF","NMF","PCA","pCMF","ZINBWaVE","tSNE","Isomap","MDS","LLE","LTSA","GLMPCA")){
    est     <- readRDS(paste0(respath,"/Hayashi/res.Hayashi.nPC",ipc,".rpt",irpt,".",methods[imeth],".rds"))
    rd1     <- as.matrix(ExtractW(est,methods[imeth]))   

}else if(methods[imeth]%in%c("ZIFA","umap")){
    est     <- read.table(paste0(respath,"/Hayashi/res.Hayashi.nPC",ipc,".rpt",irpt,".",methods[imeth],".txt"))
    rownames(est) <- colnames(countdata)  
    rd1     <- as.matrix(ExtractW(est,methods[imeth]))  
}else if(methods[imeth]%in%c("DCA")){
    est     <- read.table(paste0(respath,"/Hayashi/res.nPC",ipc,".rpt",irpt,"/reduced.tsv"),row.names=1)
    rd1     <- as.matrix(ExtractW(est,methods[imeth]))  
}

reducedDims(sim) <- SimpleList(DRM = rd1)



## kmeans 
##-------------
set.seed(1)
colData(sim)$kmeans <- kmeans(rd1, centers = length(unique(sce$grouping)),iter.max=100000,nstart=5000)$cluster
slingout <- slingshot(sim, clusterLabels = 'kmeans', reducedDim = 'DRM')


## hclust
##-------------
rd_sc   <- scale(rd1)
d       <- as.dist(1-cor(t(rd_sc)))
set.seed(1)
hc      <- hclust(d, method='ward.D')
colData(sim)$hclust <- cutree(hc, length(unique(sce$grouping)))
slingout<- slingshot(sim, clusterLabels = 'hclust', reducedDim = 'DRM')

## louvain
##-------------
library("RANN")
library(igraph)
library(Rcpp)
sourceCpp("~/clustering.cpp")
set.seed(1)

num_neigh = 50
## Finding nearest neighbors...
neighborMatrix  <- try(nn2(rd1, rd1, num_neigh+1, searchtype = "standard")[[1]][,-1])

links           <- jaccard_coeff(neighborMatrix, FALSE)
links           <- links[links[,1]>0, ]
relations       <- as.data.frame(links)
colnames(relations) <- c("from","to","weight")
g               <- graph.data.frame(relations, directed=FALSE)

Qp              <- -1 # highest modularity score 
optim_res       <- NULL 
louvain_iter    <- 100
for(iter in 1:louvain_iter) {
    Q <- cluster_louvain(g)
    if(is.null(optim_res)) {
        Qp          <- max(Q$modularity)
        optim_res   <- Q
    }else{
        Qt          <- max(Q$modularity)
        if(Qt > Qp){ #use the highest values for clustering 
            optim_res   <- Q
            Qp          <- Qt
        }
    }# end fi
}# end for the cluster results
colData(sim)$louvain <- factor(membership(optim_res))    
slingout  <- try(slingshot(sim, clusterLabels = 'louvain', reducedDim = 'DRM'))

