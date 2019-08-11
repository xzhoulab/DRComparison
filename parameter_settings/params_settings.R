## Parameter Settings
## Modified by Shiquan Sun
## Date: 2019-6-11 17:28:17
## Normalization methods: log, lcpm, lsd, and vst
params_settings <- function(method, num_clust=NULL){
	
	if(method=="FA"){
		params <- list(rotate = "oblimin", filtering_method="nonzeros", norm_method="log")
	}else if(method=="PCA"){
		params <- list(filtering_method="nonzeros", norm_method="log")
	}else if(method=="ZINBWaVE"){
		params <- list(num_core=20, doParallel=TRUE, X=NULL, V=NULL, filtering_method="nonzeros")
	}else if(method=="PoissonNMF"){
		params <- list(num_core=5, doParallel=TRUE, algorithm="lee", loss="mkl", filtering_method="nonzeros")
	}else if(method=="NMF"){
		params <- list(num_core=5, doParallel=TRUE, algorithm="lee", loss="mse", filtering_method="nonzeros", norm_method="log")
	}else if(method=="pCMF"){
		params <- list(num_core=20, doParallel=TRUE, filtering_method="nonzeros")
	}else if(method=="DiffusionMap"){
		params <- list(filtering_method="nonzeros", norm_method="log")
	}else if(method=="ICA"){
		params <- list(filtering_method="nonzeros", norm_method="log")
	}else if(method=="GLMPCA"){
		params <- list(filtering_method="nonzeros", fam="nb")
	}else if(method=="tSNE"){
		params <- list(filtering_method="nonzeros", norm_method="log", perplexity = 30)
	}else if(method=="Isomap"){
		params <- list(filtering_method="nonzeros", norm_method="log", num_prop=0.1)
	}else if(method=="LLE"){
		params <- list(filtering_method="nonzeros", norm_method="log", num_prop=0.1)
	}else if(method=="MSD"){
		params <- list(filtering_method="nonzeros", norm_method="log", num_prop=0.1)
	}else if(method=="LLLE"){
		params <- list(filtering_method="nonzeros", norm_method="log", num_prop=0.1)
	}else if(method=="LTSA"){
		params <- list(filtering_method="nonzeros", norm_method="log", num_prop=0.3)
	}else{
		stop("no defined method!")
	}# end fi
	
	return(params)
}# end function parameter setting




