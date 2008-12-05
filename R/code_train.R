# --> TODO:train output codes (ncol(labels) binary SVM, trained at each level of the hierarchy, exclude negative examples not belonging to the actual sub-hierarchy)
# y is a matrix of 1 and -1!
code_train = function(x, y, parentsLev1, parentsLev2, parentsLev12){
	code = matrix(0, ncol=ncol(y), nrow=nrow(y))
	detectors = list()
	domains = list()
	for(i in 1:ncol(y)){	
		cat("Level detector ",i, "\n\n")	
		ytmp = y[,i]
		pos = which(ytmp == 1)
		neg = which(ytmp == -1)	
# 		print(parentsLev1)
# 		print(parentsLev2)			
		# If possible, only take those negatives, which belong to the same super-level				
		if(i %in% parentsLev1)
			sel = 1:length(neg)
		else if(i %in% parentsLev2)
			sel = which(apply(y[neg,parentsLev1],1, function(xx) all(xx == y[pos[1],parentsLev1])))
		else # i is in level3
			sel = which(apply(y[neg,parentsLev12],1, function(xx) all(xx == y[pos[1],parentsLev12])))
		if(length(sel) < 20)
			sel = sample(1:length(neg), length(pos))		
		ytmp = ytmp[c(pos, neg[sel])]
		xtmp = x[c(pos, neg[sel]),]								
		cat("-->SVM training (#pos = ", length(pos), "#neg = ", length(neg[sel]), ", #features = ", ncol(xtmp), ")\n")	
# 		if(modsel){	
# 			C = -3:4
# 			bestC = 10^C[which.min(sapply(C, function(CC) errorbound(xtmp, ytmp, 10^CC)))]		
# 		}
# 		else
			bestC = 1	
		cat("-->SVM training with parameter C = ", bestC,"\n")		
		detectors[[i]] = gene2pathway:::svmlearn(xtmp, ytmp, bestC, prob.model=FALSE)		
		cat("-->Generating output codes\n")
		domains[[i]] = colnames(xtmp)				
		code[,i] = gene2pathway:::svmpredict(detectors[[i]], x[,domains[[i]]], type="decision")		
		wabs = abs(t(unlist(coef(detectors[[i]]))) %*% x[unlist(SVindex(detectors[[i]])),domains[[i]]])
		keep = (wabs > 0)	
		domains[[i]] = domains[[i]][keep]	
		cat(length(domains[[i]]), " domains important\n")		
	}
	list(code=code, detectors=detectors, domains=domains)
}
