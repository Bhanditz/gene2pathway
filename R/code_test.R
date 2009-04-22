code_test = function(detectors, xtst){
	code = matrix(0, ncol=length(detectors), nrow=nrow(xtst))	
	for(i in 1:length(detectors)){
		code[,i] = gene2pathway:::svmpredict(detectors[[i]], xtst, type="decision")
	}
	rownames(code) = rownames(xtst)
	code
}
