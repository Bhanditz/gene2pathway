# batch prediction
svmpredict = function(fit, xtst, type="response"){		
	ntst = nrow(xtst)
	if(is.null(ntst)){
		ntst = 1
		xtst = t(as.matrix(xtst))
	}
	yhat = double(ntst)
	if(ntst > 250){
		s = seq(1, ntst, by=250)
		for(ss in s){
			if(ss == s[length(s)])
				yhat[ss:ntst] = predict(fit, xtst[ss:ntst,], type=type)
			else
				yhat[ss:(ss+249)] = predict(fit, xtst[ss:(ss+249),], type=type)
		}
	}
	else
		yhat = predict(fit, xtst, type=type)
	yhat
}
