rm(list=ls())
library(SingleCellExperiment)
library(monocle3)

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


## the one with the largest frequency in the start cell stage treated as the starter
get_earliest_principal_node <- function(input,startcell){
  cell_ids       <- which(colData(input)[, "milestone_id"] == startcell)
  closest_vertex <- input@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(input), ])
  root_pr_nodes  <- igraph::V(principal_graph(input)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}


methods = c("PCA","DiffusionMap","ICA","FA","NMF","PoissonNMF","pCMF","ZINBWaVE","ZIFA","DCA","tSNE","umap","Isomap","MDS","LLE","LTSA","GLMPCA","scvis")


datapath    = "~/scRNAseqDRComparison/datasets/gti_full/"
respath     = "~/scRNAseqDRComparison/results.TI/"

ipc = 2
irpt = 1
imeth = 1


for(ipc in c(2,6,14,20)){
    for(irpt in 1:5){
        for(imeth in 1:18){
            sce             <- readRDS(paste0(datapath,"/gti_full_Hayashi.rds"))
            countdata       <- t(sce$counts)
            countdata       <- countdata[which(rowSums(countdata>0)>5),]
            countdata       <- countdata[, which(colSums(countdata>0)>10)]

            cell_idx        <- which(colSums(countdata>0)>10)
            cell_metadata   <- as.data.frame(sce$cell_info[cell_idx,])
            rownames(cell_metadata) <- cell_metadata[,1]

            gene_annotation <- data.frame(gene_short_name=rownames(countdata))
            rownames(gene_annotation) <- rownames(countdata)


            cds <- new_cell_data_set(countdata,
                        cell_metadata = cell_metadata,
                        gene_metadata = gene_annotation)

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

            ## substitute the UMAP attributes with the DR object
            cds@reducedDims$UMAP <- rd1
            cds <- cluster_cells(cds,reduction_method="UMAP")
            cds <- learn_graph(cds)
            start_label <- as.data.frame(sce$milestone_network)[1,1]
            cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds,start_label))
        }
    }
}






