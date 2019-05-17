## Parameter Settings
## Modified by Shiquan Sun
## Date: 2019-4-3 08:26:59
## Normalization methods: log, lcpm, lsd, and vst
params_settings <- function(method, num_clust=NULL){
	
	if(method=="FA"){
		params <- list(rotate = "oblimin", filtering_method="nonzeros", norm_method="log")
	}else if(method=="PCA"){
		params <- list(filtering_method="nonzeros", norm_method="log")
	}else if(method=="ZINBWaVE"){
		params <- list(num_core=20, doParallel=TRUE, X=NULL, V=NULL, filtering_method="nonzeros")
	}else if(method=="PoissonNMF"){
		params <- list(num_core=5, doParallel=TRUE, filtering_method="nonzeros")
	}else if(method=="NMF"){
		params <- list(num_core=5, doParallel=TRUE, filtering_method="nonzeros", norm_method="log")
	}else if(method=="pCMF"){
		params <- list(num_core=20, doParallel=TRUE, filtering_method="nonzeros")
	}else if(method=="DiffusionMap"){
		params <- list(filtering_method="nonzeros", norm_method="log")
	}else if(method=="ICA"){
		params <- list(filtering_method="nonzeros", norm_method="log")
	}else{
		stop("no defined method!")
	}# end fi
	
	return(params)
}# end function parameter setting




