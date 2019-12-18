##----------------------------------------- 
## assign name would help, numeric --> character
rm(list=ls())
library(slingshot)


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


imeth   <- 1


datapath    = "~/scRNAseqDRComparison/datasets/gti_full/"
outpath = "~/slingshot/"
respath     = "~/kendall_slingshot/"

fn <- list.files(datapath,"rds")


methods = c("PCA","DiffusionMap","ICA","FA","NMF","PoissonNMF","pCMF","ZINBWaVE","ZIFA","DCA","tSNE","umap","Isomap","MDS","LLE","LTSA","GLMPCA","scvis")

data_inventory <- c("Schlitzer","Petropoulos","Olsson","Hayashi","LiM","LiF",
                    "ZhangBeta" , "ZhangAlpha" ,"GuoF",  "GuoM","KowalczykYoung",
                    "KowalczykOld","ShalekLPS","Trapnell")



sce         <- readRDS(paste0(datapath,"/gti_full_Hayashi.rds"))
celltype    <- c(sce$milestone_network$from,sce$milestone_network$to[nrow(sce$milestone_network)])
tti         <- 1:length(celltype)
names(tti)  <- celltype
tmp1        <- cbind.data.frame(cell=names(sce$grouping),tp=tti[sce$grouping])


all_pc_rpt <- c()
for(irpt in 1:3){
combined_pc <- c()
for(ipc in c(2,6,14,20)){
    load(paste0(outpath,"/res.Hayashi.nPC",ipc,".rpt",irpt,".",methods[imeth],".SS.slingshot.rds"))
    pc_out <- c()
    for(ilin in 1:length(SlingshotDataSet(slingout)@lineages)){
        labels      <- SlingshotDataSet(slingout)@lineages[[ilin]]
        eti         <- 1:length(labels)
        names(eti)  <- labels
        if(is.null(names(slingout$kmeans))){
            names(slingout$kmeans) <- rownames(slingout@colData)
        }
        sub_kmeans  <- slingout$kmeans[which( slingout$kmeans %in% labels)]
        ep          <- eti[as.character(sub_kmeans)] # use charater to match the name of eti
        tmp2        <- cbind.data.frame(cell=names(sub_kmeans),ep=ep,cluster_id = sub_kmeans)
        res         <- merge(tmp1,tmp2,by="cell")
        union_cell  <- length(union(tmp1[,"cell"],tmp2[,"cell"]))
        cat("lineage",ilin,":\n")
        # print(kendall_func(res[,"tp"],res[,"ep"],tot=union_cell))
        single_lineage_out  <- kendall_func(res[,"tp"],res[,"ep"],tot=union_cell)
        pc_out              <- c(pc_out,single_lineage_out) 
    }
    combined_pc <- c(combined_pc, max(abs(pc_out)))
    rm(slingout)
}

all_pc_rpt <- cbind(all_pc_rpt,combined_pc)
}
average_kendall <- apply(all_pc_rpt,1,mean)

# save(combined_pc,file=paste0(respath,"/res.Hayashi.rpt",irpt,".",methods[imeth],".kmeans.SS.kendall.rds"))



