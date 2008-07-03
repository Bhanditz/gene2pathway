# --> TODO Use output codes to train structure SVM
# C is a matrix of 1 and 0!
struct_train = function(f, C, treesizes, parentPaths, nu=0.1, m=2, T=1){		
	# ranking perceptron
	F = function(W, f, Cy){	
		crossprod(W,f*Cy)
	}


	# precalculate losses
	CC = unique(C)	
	W = double(ncol(f))
	codenrs = double(nrow(C))
	for(i in 1:nrow(CC)){		
		codenrs[which(apply(C, 1, function(x) all(x == CC[i,])))] = i
	}		
	codedists = matrix(0, nrow=nrow(CC), ncol=nrow(CC))
	dimnames(codedists) = list(1:nrow(CC), 1:nrow(CC))
	for(i in 1:nrow(CC)){
		if(i > 1){
			for(j in 1:(i-1)){
				codedists[i,j] = loss(CC[i,], CC[j,], treesizes, parentPaths)
				codedists[j,i] = codedists[i,j]
			}
		}
	}		
	maxFriendDist = 0
	for(t in 1:T){	
		rnd = sample(1:nrow(f), nrow(f))
		f = f[rnd,]
		C = C[rnd,]
		codenrs = codenrs[rnd]
		for(i in 1:nrow(f)){
			friends = which(codedists[codenrs[i],codenrs] <= maxFriendDist)			
			foes = setdiff(1:nrow(C), friends)			
			l = foes[which.max((t(matrix(rep(f[i,],length(foes)),nrow=ncol(f)))*C[foes,]) %*% as.matrix(W))]			
			if(F(W, f[i,], C[i,]) - m < F(W, f[i,], C[l,])){				
				lo = loss(C[i,],C[l,], treesizes, parentPaths)
				W = W + nu*lo  * (f[i,]*C[i,] - f[i,]*C[l,])
			}
			cat(".")
		}		
		codepred = gene2pathway:::struct_predict(list(W=W, C=CC), f)$yhat
		err = sapply(1:nrow(f), function(j) loss(C[j,], codepred[j,], treesizes, parentPaths))
		cat("\ntraining error = ", median(err)," (min = ", min(err), "max = ", max(err), "avg = ", mean(err), ")\n")
	}
	fit = list(W = W, C = CC)
}
