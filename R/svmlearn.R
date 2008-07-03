svmlearn = function(x, y, C){
		cat("Training SVM ( C = ", C, ") ... ")
		weights = table(y)	
		weights = c(1,1)
		names(weights) = c("1","-1") 			
		y = as.factor(y)
		fit = ksvm(x, y, type="C-svc", kernel="vanilladot", kpar=list(), scaled=FALSE, class.weights=weights, C=C)
		cat("done.\n")
		fit
}
