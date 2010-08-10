buildTrainingSet = function(minnmap=30, level1Only="Metabolism", level2Only="Genetic Information Processing", organism="hsa", remove.duplicates=FALSE, gene2Domains=NULL){		
#	if(KEGG.package){ # this is fast			
#		cat("Retrieving KEGG information via KEGG.db package ...\n")
#		organisms = unique(gsub("[0-9]*","",AnnotationDbi::ls(KEGGPATHID2EXTID)))		
#		if(!(organism %in% organisms))
#			stop(paste("Organism '", organism, "' is unknown in KEGG package! Please retry with KEGG.package=FALSE (slow). Please refer also to <URL:http://www.genome.jp/kegg-bin/create_kegg_menu> for a complete list of organisms supported by KEGG.", sep=""))
#		hKEGGids <- grep(paste("^",organism,sep=""), AnnotationDbi::ls(KEGGPATHID2EXTID), value=TRUE)
#		path2Genes <- AnnotationDbi::mget(hKEGGids, KEGGPATHID2EXTID)
#		hKEGGgenes <- unique(unlist(path2Genes, use.names=FALSE))
#		hKEGGgenes <-  hKEGGgenes[!is.na(hKEGGgenes)]
#		genes2Path = AnnotationDbi::mget(hKEGGgenes,KEGGEXTID2PATHID)
#		genes2Path = sapply(genes2Path, function(g) sapply(g, function(gg) sub(organism,"", gg)))
#		if(organism == "dme")
#			flyBase = TRUE		
#		else
#			flyBase = FALSE
#	}
#	else{ # this is slow
		cat("Retrieving KEGG information via SOAP ...\n")
		organisms=list.organisms()
		if(!(organism %in% names(organisms)))
			stop(paste("Organism '", organism, "' is unknown in KEGG! Please refer to <URL:http://www.genome.jp/kegg-bin/create_kegg_menu> for a complete list of supported organisms.", sep=""))
		pathways = names(list.pathways(organism))
		path2Genes = try(lapply(pathways, get.genes.by.pathway))
		names(path2Genes) = pathways
		hKEGGgenes <- unique(unlist(path2Genes, use.names=FALSE))		
		hKEGGgenes <-  hKEGGgenes[!is.na(hKEGGgenes)]				
		genes2Path = sapply(hKEGGgenes, function(g) names(path2Genes)[sapply(path2Genes, function(p) any(g == p))]) # this is probably faster than calling get.pathways.by.genes		
		genes2Path = sapply(genes2Path, function(g) sapply(g, function(gg) sub(paste("path:",organism,sep=""),"", gg)))
		flyBase = FALSE
#	}
	kegg_hierarchy = gene2pathway:::getKEGGHierarchy(level1Only, level2Only)
	parentPaths = kegg_hierarchy$parentPaths
	code_vector = kegg_hierarchy$code_vector	
	genes2Path = sapply(genes2Path, function(g) intersect(g, rownames(code_vector)))
	genes2Path = genes2Path[sapply(genes2Path, length) > 0]
	hKEGGgenes = names(genes2Path)		
	if(length(grep(":", hKEGGgenes)) > 0){
		if(organism == "hsa")
			hKEGGgenes = sub(paste(organism,":",sep=""), "", hKEGGgenes)
		else
			hKEGGgenes = gene2pathway:::KEGG2Entrez(hKEGGgenes, organism=organism)	
		names(genes2Path) = hKEGGgenes
	}
	cat("done \n")
	
	if(!is.null(gene2Domains) && length(intersect(hKEGGgenes, names(gene2Domains))) == 0){
		stop("There is a conflict between KEGG Entrez gene / open reading frame IDs and the gene identifiers in list 'gene2Domains'.\n 'gene2Domains' should be a list of Entrez gene IDs and their correponding InterPro domains.")		
	}
	features = gene2pathway:::getInterProDomains(hKEGGgenes, gene2Domains=gene2Domains, organism=organism, flyBase=flyBase)	
	if(remove.duplicates){
		cat("Removing genes with duplicate feature vectors \n")		
		features = unique(features)
	}		
	common = intersect(hKEGGgenes, rownames(features))
	features = features[common,]
	features = features[,colSums(features) > 0]
	
	genes2Path = genes2Path[common]
	labels = matrix(0, ncol=ncol(code_vector), nrow=length(genes2Path))
	colnames(labels) = colnames(code_vector)
	rownames(labels) = names(genes2Path)
	for(g in 1:length(genes2Path)){
		pathways = intersect(genes2Path[[g]], rownames(code_vector)) # we only want pathways at level 3
		if(length(pathways) > 1){
			labels[g, ] = pmin(colSums(code_vector[pathways,]),1) 
		}
		else
			labels[g, ] = code_vector[pathways,]	
	}	
					
	# require a minimum of mapping genes per branch	
	nmaps = colSums(labels)	
	labels = labels[,nmaps > minnmap]
	freq = table(unlist(kegg_hierarchy$parentPaths[colnames(labels)]))
	redchild = intersect(names(kegg_hierarchy$parentPaths)[sapply(kegg_hierarchy$parentPaths, function(pa) any(pa %in% names(freq[freq == 1])))], colnames(labels))
	labels = labels[,!(colnames(labels) %in% redchild)]
	cat("Removing hierarchy branches with < ", minnmap, " mapping genes\n--->", ncol(labels), "hierarchy branches left\n")
	if(ncol(labels) == 0)
		stop("No sufficient number of mapping genes! Impossible to train model.")

	code_vector = code_vector[,colnames(labels)]
	code_vector = code_vector[rowSums(code_vector) > 0,]
	
	kegg_hierarchy$parentsLev1 = which(colnames(labels) %in% kegg_hierarchy$pathIDsLev1)
	kegg_hierarchy$parentsLev2 = which(colnames(labels) %in% kegg_hierarchy$pathIDsLev2)
	kegg_hierarchy$parentsLev12 = which(colnames(labels) %in% c(kegg_hierarchy$pathIDsLev1, kegg_hierarchy$pathIDsLev2))
	
	treesizes = colMeans(code_vector)
	
	list(features=features, labels=labels, treesizes=treesizes, kegg_hierarchy=kegg_hierarchy)
}