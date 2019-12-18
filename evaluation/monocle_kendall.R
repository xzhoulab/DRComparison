rm(list=ls())

datapath    = "~/scRNAseqDRComparison/datasets/gti_full/"
monocle_path = "~/output/monocle/"
outpath     = "~/kendall_monocle/"

kendall_func <- function(x,y,tot=NULL){
    if(length(x)>1&length(y)>1){
        P   <- 0
        Q   <- 0
        for(i in 1:(length(x)-1)){
            for(j in (i+1):length(x)){
                if(sign(x[i]-x[j])*sign(y[i]-y[j])<0){
                    Q = Q + 1
                }

                if(sign(x[i]-x[j])*sign(y[i]-y[j])>0){
                    P = P + 1
                }

            }
        }
        if(is.null(tot)){tot=length(x)}
        out <- (P-Q)/choose(tot,2) # option 1, slingshot, max
    }else{
        out <- 0
    }

    return(out)
}


methods = c("PCA","DiffusionMap","ICA","FA","NMF","PoissonNMF","pCMF","ZINBWaVE","ZIFA","DCA","tSNE","umap","Isomap","MDS","LLE","LTSA","GLMPCA","scvis")
imeth   = 1


sce         <- readRDS(paste0(datapath,"/gti_full_Hayashi.rds"))
celltype    <- c(sce$milestone_network$from,sce$milestone_network$to[nrow(sce$milestone_network)])
tti         <- 1:length(celltype)
names(tti)  <- celltype
tmp1        <- cbind.data.frame(cell=names(sce$grouping),tp=tti[sce$grouping])


combined_pc_kendall <- c()
icount <- 0
for(ipc in c(2,6,14,20)){
    single_pc   <- c()
    for(irpt in 1:2){
        load(paste0(monocle_path,"/res.Hayashi.nPC",ipc,".rpt",irpt,".",methods[imeth],".monocle.rds"))

        summary_monocle <- cbind.data.frame(closest_vertex=cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex,pseudotime=cds@principal_graph_aux$UMAP$pseudotime,clusters=cds@clusters$UMAP$clusters)
        cluster_summary <- aggregate(summary_monocle$pseudotime,by=list(summary_monocle$clusters),summary)[,2]
        rownames(cluster_summary)   <- paste0("c",1:nrow(cluster_summary))

        eti                         <- 1:length(unique(summary_monocle$clusters))
        cluster_order_by_median     <- rownames(cluster_summary[order(cluster_summary[,"Median"]),])

        names(eti)                  <- cluster_order_by_median
        summary_monocle$cell_order  <- eti[paste0("c",summary_monocle$clusters)]

        if(identical(as.character(tmp1$cell),rownames(summary_monocle))){
            single_rpt <- kendall_func(tmp1[,"tp"],summary_monocle$cell_order,tot=length(tmp1[,"tp"]))
        }

        rm(cds,summary_monocle)

        single_pc <- c(single_pc,single_rpt) 
        rm(single_rpt)
    }
        combined_pc_kendall <- rbind(combined_pc_kendall,single_pc) 
        rm(single_pc)
}

rownames(combined_pc_kendall) <- paste0("PC",c(2,6,14,20))
average_kendall <- apply(combined_pc_kendall,1,mean)
