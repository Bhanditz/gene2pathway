test.overrepresentation = function(genesOfInterest, predpath, KEGGonly=FALSE, cutoff=0.1){
	others = setdiff(names(predpath$gene2Path), genesOfInterest)
	if(KEGGonly){
		KEGGgenes = names(predpath$byKEGG[predpath$byKEGG])
		gene2Path = predpath$gene2Path[predpath$byKEGG]
		genesOfInterest = intersect(genesOfInterest, KEGGgenes)
		others = intersect(others, KEGGgenes)
	}
	else{
		gene2Path = predpath$gene2Path
	}
	pathways = setdiff(unlist(gene2Path), NA)	
	freqsig = table(unlist(gene2Path[genesOfInterest]))
	freqothers = table(unlist(gene2Path[others]))
	p.values = sapply(pathways, function(p){
		conf.tab = matrix(c(freqsig[p], freqothers[p], length(genesOfInterest) - freqsig[p], length(others) - freqothers[p]),nrow=2, dimnames=list(c("differential", "not differential"),c(p,"others")))
		conf.tab[is.na(conf.tab)] = 0
		fisher.test(conf.tab, alternative="g")$p.value
	})
	p.values = p.adjust(p.values, method="BY")
	p.values = p.values[p.values < cutoff]
	if(exists("kegg_hierarchy", envir=gene2pathwayEnv))
		kegg_hierarchy = get("kegg_hierarchy", envir=gene2pathwayEnv)
	else{
		kegg_hierarchy = gene2pathway:::getKEGGHierarchy(level1Only=c(), level2Only=c())
		assign("kegg_hierarchy", kegg_hierarchy, envir=gene2pathwayEnv)		
	}
	pathnames = c(kegg_hierarchy$pathNamesLev1, kegg_hierarchy$pathNamesLev2, kegg_hierarchy$pathNamesLev3)	
	p.values = as.data.frame(cbind(pathnames[names(p.values)],p.values))
	colnames(p.values) = c("pathname", "p.value")
	p.values
}
