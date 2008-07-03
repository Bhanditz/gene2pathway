struct_predict = function(fit, ftst, diagnostics=FALSE, PDF=FALSE){			
	yhat = matrix(0, ncol=ncol(fit$C), nrow=nrow(ftst))
	diff = double(nrow(ftst))	
	for(i in 1:nrow(ftst)){	
		scores = (t(matrix(rep(ftst[i,],nrow(fit$C)),nrow=ncol(ftst)))*fit$C) %*% as.matrix(fit$W)
		maxscore = max(scores)
		mycodes = which(scores == maxscore)	
		sndbest = max(scores[scores < maxscore])		
		diff[i] = abs(maxscore - sndbest)		
		if(diagnostics){			
			idx = order(scores, decreasing=TRUE)	
			if(PDF)		
				pdf(file=paste(rownames(ftst)[i],".pdf",sep=""))
			plot(1:10, scores[idx[1:10]], xlab="solution", ylab="score", main=paste("score distribution (  1st-2nd = ", round(diff[i],2), ")"))
			points(1,maxscore,pch=21,cex=1.7,lwd=3)			
			if(PDF)
				dev.off()
			else
				readline("Press key to continue ...")
		}		
		if(length(mycodes) > 1)			
			yhat[i,] = pmin(colSums(fit$C[mycodes, ]),1)					
		else
			yhat[i,] = fit$C[mycodes, ]
	}
	list(yhat=yhat, diff=diff)
}
