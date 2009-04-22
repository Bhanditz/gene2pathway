svmlearn = function(x, y, C, prob.model=FALSE, cross=0){
		cat("Training SVM ( C = ", C, ") ... ")
# 		weights = table(y)	
# 		weights = c(1,1)
# 		names(weights) = c("1","-1") 			
		y = as.factor(y)
		if(class(x) == "kernelMatrix")
			fit = ksvm(x, y, type="C-svc", C=C, prob.model=prob.model, cross=cross)
		else
			fit = ksvm(x, y, type="C-svc", kernel="vanilladot", kpar=list(), scaled=FALSE, C=C, prob.model=prob.model, cross=cross)
		cat("done.\n")		
		fit
}
