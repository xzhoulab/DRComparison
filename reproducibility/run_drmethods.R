## run dimensionality reduction
## Modified by Shiquan Sun
## Date: 2019-4-3 07:56:04

args <- as.numeric(commandArgs(TRUE))
data_indx <- args[1]
method_indx <- args[2]
ipc <- args[3]
irpt <- args[4]

## library packages
library(pryr) ## to count memory usage

PCList <- list(
c(8, 20, 38, 58), #Baron
)
num_pc <- PCList[[data_indx]][ipc]
DATASET <- c("Baron")
METHODS <- c("FA", "PCA", "ZINBWaVE", "pCMF", "DiffusionMap", "ICA", "NMF", "PoissonNMF")

idata <- DATASET[data_indx]
imethod <- METHODS[method_indx]

data_path <- "./data"
method_path <- "./algorithms"
param_path <- "./parameter_settings"
res_path <- "./"
## load data set

## check folder exists, and create
if(!file.exists(paste0(res_path,"/",idata))){
    dir.create(file.path(res_path, idata))
}# end fi


if(!file.exists(paste0(res_path,"/",idata,"/res.",idata,".nPC",num_pc,".rpt",irpt,".",imethod,".rds")) ){
sce <- readRDS(paste0(data_path,"/sce_",idata,".rds") )
## source method
source(paste0(method_path,"/call_",imethod,".R") )
source(paste0(param_path,"/params_settings.R") )

params <- params_settings(imethod)

cat(paste0("run = ", irpt, "; data set: ",idata,"...\n"))
res_method <- get(paste0("call_",imethod))(sce=sce, num_pc=num_pc, params=params)
	
res <- list(dataset=idata, 
	method=imethod, 
	num_repeat=irpt, 
	res=res_method$res, 
	ctime=res_method$ctimes )

## save results as RDS file
saveRDS(res, file=paste0(res_path,"/",idata,"/res.",idata,".nPC",num_pc,".rpt",irpt,".",imethod,".rds") )
}# end fi
