
## data format from transform


data_path <- "~/scRNAseqDRComparison/datasets/sce_full"

## load data set
ifiles <- "sce_Baron.rds"
sce <- readRDS(paste0(data_path,"/",ifile) )
sce <- scater::normalize(sce)

## ZIFA format
norm_counts <- logcounts(sce)
norm_counts <- norm_counts[which(rowSums(abs(norm_counts)>0)>5),]
write.table(norm_counts,file=paste0(data_path,"/sce_full_",idata,"_ZIFA.txt"),col.names=F,row.names=F,quote=F,sep="\t")

## DCA format
counts <- as.matrix(counts(sce))
counts <- counts[which(rowSums(abs(counts)>0)>5),]
write.table(round(counts),file=paste0(data_path,"/sce_full_",idata,"_DCA.tsv"),col.names=T,row.names=T,quote=F,sep="\t")

## scScope format
counts <- as.matrix(counts(sce))
counts <- counts[which(rowSums(abs(counts)>0)>5),]
counts <- t(counts)
write.csv(round(counts),file=paste0(data_path,"/sce_full_",idata,"_scScope.csv"),col.names=T,row.names=T)