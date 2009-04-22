getKEGGHierarchy = function(level1Only="Metabolism", level2Only="Genetic Information Processing"){
	kegg = RCurl:::getURL("ftp://ftp.genome.jp/pub/kegg/brite/ko/ko00001.keg")
	kegg = unlist(strsplit(kegg,"\n"))
	level1 = grep("A<B>", kegg)
	level2 = grep("B  <B>", kegg)
	level3 = grep("PATH", kegg)
	pathIDsLev1 = sapply(kegg[level1], substr, 5, 9)
	pathIDsLev2 = sapply(kegg[level2], substr, 7, 11)
	pathIDsLev3 = sapply(kegg[level3], substr, 6, 10)
	pathNamesLev1 = sapply(kegg[level1], function(k) substr(k, 11, gregexpr("<",k)[[1]][2]-1))
	pathNamesLev2 = sapply(kegg[level2], function(k) substr(k, 13, gregexpr("<",k)[[1]][2]-1))
	pathNamesLev3 = sapply(kegg[level3], function(k) substr(k, 12, gregexpr("\\[",k)[[1]][1]-2))
	names(pathNamesLev1) = pathIDsLev1
	names(pathNamesLev2) = pathIDsLev2
	names(pathNamesLev3) = pathIDsLev3	

	parentPaths = lapply(level2, function(p) pathIDsLev1[which.max(level1[level1 < p])])	
	parentPaths = c(parentPaths, lapply(level3, function(p) c(pathIDsLev1[which.max(level1[level1 < p])], pathIDsLev2[which.max(level2[level2 < p])])))
	names(parentPaths) = c(pathIDsLev2, pathIDsLev3)
	
	code_vector = matrix(0, ncol=length(level1)+length(level2)+length(level3), nrow=length(level3))
	dimnames(code_vector) = list(pathIDsLev3, c(pathIDsLev1,pathIDsLev2,pathIDsLev3))
	if(length(level1Only) > 0) level1Only = pathIDsLev1[sapply(level1Only, function(l) grep(l, names(pathIDsLev1)))] 
	if(length(level2Only) > 0) level2Only = pathIDsLev1[sapply(level2Only, function(l) grep(l, names(pathIDsLev1)))]
	for(p in pathIDsLev3){
		pa = unlist(parentPaths[p])		
		code_vector[p, pa[1]] = 1	
		if(!(pa[1] %in% level1Only))
			code_vector[p, pa[2]] = 1
		if(!(pa[1] %in% c(level1Only, level2Only)))
			code_vector[p, p] = 1
	}
	code_vector = code_vector[(code_vector[,grep("Human Diseases", names(pathIDsLev1))] != 1),] 	
	cs = colSums(code_vector)
	code_vector = code_vector[,(cs > 0) & (cs < nrow(code_vector))]
	list(code_vector=code_vector, parentPaths=parentPaths, pathIDsLev1=pathIDsLev1, pathIDsLev2=pathIDsLev2, pathIDsLev3=pathIDsLev3, pathNamesLev1=pathNamesLev1, pathNamesLev2=pathNamesLev2, pathNamesLev3=pathNamesLev3)
}

