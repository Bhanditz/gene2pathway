predict.gene2pathway = function(object, newdata=NULL, diagnostics=FALSE, PDF=FALSE,...){		
	# 1. predict codes
	# 2. use codes to predict structures
	testall = function(fit, xtst){	
		cat(".")
		codes = gene2pathway:::code_test(fit$detectors, xtst)	
		# predict at all hierarchy levels
		gene2pathway:::struct_predict(fit, codes, diagnostics=diagnostics, PDF=PDF)		
	}
	
	testallBag = function(model, xtst){
		test = testall(model[[1]], xtst)$yhat[,,drop=FALSE]
		for(b in 2:length(model)){
			test = test + testall(model[[b]], xtst)$yhat[,,drop=FALSE]
		}		
		votes = test/length(model)
		test = (votes > 0.5)*1
		scores = apply(votes[,,drop=FALSE],1,function(v) 0.5*(mean(v[v > 0.5]) + 1-mean(v[v<=0.5])))
		list(yhat=test, scores=scores)	
	}

	if(missing(newdata))
		return(NULL)	
	if(class(object) == "model"){
		pred = testall(object, newdata)	
		scores = pred$diff	
		cnames = colnames(object$C)	
	}
	else{
		pred =  testallBag(object, newdata)	
		scores = pred$scores
		cnames = colnames(object[[1]]$C)
	}
	predcode = pred$yhat[,,drop=FALSE]
	dimnames(predcode) = list(rownames(newdata), cnames)	
	
# 	bereite Ergebnis in lesbarer Form auf
	gene2Path = lapply(1:nrow(predcode), function(p) colnames(predcode)[predcode[p,,drop=FALSE] == 1])	
	names(gene2Path) = rownames(predcode)	
	list(gene2Path=gene2Path, scores=scores)
}
