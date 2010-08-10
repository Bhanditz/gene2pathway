retrain = function(minnmap=30, level1Only="Metabolism", level2Only="Genetic Information Processing", organism="hsa", gene2Domains=NULL, remove.duplicates=FALSE, use.bagging=TRUE, nbag=11){		
	
	# 1. train codes
	# 2. used codes to train structure
	trainall = function(x, y){		
		ytmp = y
		ytmp[ytmp == 0] = -1	
		codes = gene2pathway:::code_train(x, ytmp, traindat$kegg_hierarchy$parentsLev1, traindat$kegg_hierarchy$parentsLev2, traindat$kegg_hierarchy$parentsLev12)	
		# fit model at all hierarchy levels
		fit = gene2pathway:::struct_train(codes$code, y, traindat$treesizes, traindat$kegg_hierarchy$parentPaths)
		model = list(W=fit$W, C=fit$C, detectors=codes$detectors, used_domains=codes$domains, alldomains=colnames(x), allpathways=colnames(fit$C), treesizes=traindat$treesizes, kegg_hierarchy=traindat$kegg_hierarchy)
		class(model) = "model"
		model
	}

	trainallBag = function(x, y){
		model = list()
		for(b in 1:nbag){
			boot = sample(1:nrow(y), nrow(y), replace=TRUE)
			model[[b]] = trainall(x[boot,], y[boot,])		
		}
		model
	}

	traindat = gene2pathway:::buildTrainingSet(minnmap=minnmap, level1Only=level1Only, level2Only=level2Only, organism=organism, remove.duplicates=remove.duplicates, gene2Domains=gene2Domains)	
	if(use.bagging)
		modelKEGG = trainallBag(traindat$features, traindat$labels)
	else
		modelKEGG = trainall(traindat$features, traindat$labels)	
	save(modelKEGG, file=paste("classificationModel_", organism,".rda",sep=""))
	modelKEGG
}

retrain.signaltrans = function(minnmap=10, organism="hsa", gene2Domains=NULL, remove.duplicates=FALSE, use.bagging=TRUE, nbag=11){
	# 1. train codes
	# 2. used codes to train structure
	trainall = function(x, y){
		ytmp = y
		ytmp[ytmp == 0] = -1	
		codes = gene2pathway:::code_train(x, ytmp, traindat$kegg_hierarchy$parentsLev1, traindat$kegg_hierarchy$parentsLev2, traindat$kegg_hierarchy$parentsLev12)	
		# fit model at all hierarchy levels
		fit = gene2pathway:::struct_train(codes$code, y, traindat$treesizes, traindat$kegg_hierarchy$parentPaths)
		model = list(W=fit$W, C=fit$C, detectors=codes$detectors, used_domains=codes$domains, alldomains=colnames(x), allpathways=colnames(fit$C), treesizes=traindat$treesizes, kegg_hierarchy=traindat$kegg_hierarchy, elemIDs=traindat$elemIDs)
		class(model) = "model"
		model
	}

	trainallBag = function(x, y){
		model = list()
		for(b in 1:nbag){
			boot = sample(1:nrow(y), nrow(y), replace=TRUE)
			model[[b]] = trainall(x[boot,], y[boot,])		
		}
		model
	}

	traindat = gene2pathway:::buildTrainingSet.signaltrans(minnmap=minnmap, remove.duplicates=remove.duplicates, organism=organism, gene2Domains=gene2Domains)
	if(use.bagging)
		modelSignalTrans = trainallBag(traindat$features, traindat$labels)
	else
		modelSignalTrans = trainall(traindat$features, traindat$labels)		
	save(modelSignalTrans, file=paste("classificationModelSignalTrans_",organism,".rda",sep=""))
	modelSignalTrans
}
