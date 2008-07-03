buildTrainingSet.signaltrans = function(minnmap=10, organism="hsa", remove.duplicates=FALSE, gene2Domains=NULL){
	organisms=list.organisms()
	if(!(organism %in% names(organisms)))
		stop(paste("Organism '", organism, "' unknown in KEGG! Please refer to <URL:http://www.genome.jp/kegg-bin/create_kegg_menu> for a complete list of supported organisms."))
	cat("Retrieving KEGG information via SOAP ...\n", sep="")
	kegg_hierarchy = gene2pathway:::getKEGGHierarchy()
	pathways = names(which(sapply(kegg_hierarchy$parentPaths, function(p) "01320" %in% unlist(p))))
	pathways = paste("path:",organism, pathways,sep="")
	genes = try(lapply(pathways, get.genes.by.pathway))
	names(genes) = pathways
	genes = lapply(genes, function(g) sub(paste(organism,"\\:",sep=""),"", g))
	pathways = pathways[sapply(genes, length) > 0]

	KEGG2Entrez.table = gene2pathway:::KEGG2Entrez(organism=organism)
	clustall = try(sapply(pathways, getComponents.internal, KEGG2Entrez.table))
	clust = clustall["geneIDs",]
	elems = clustall["elemIDs",]
	names(clust) = sub(paste("path:",organism,sep=""),"",names(clust))
	names(elems) = sub(paste("path:",organism,sep=""),"",names(elems))
	clust = clust[!sapply(clust, function(c) is.null(c))]		
	allen = sapply(sapply(clust, function(c) if(class(c) == "list") sapply(c,length) else length(c)), length)	
	roots = names(clust)
	childs = roots[allen > 1]
	childs = unlist(sapply(childs, function(c) paste(c, seq(1:allen[c]),sep=".")))
	levels = c(roots, childs)	
	hKEGGgenes = unique(unlist(clust))		
	cat("done\n")
	
	if(!is.null(gene2Domains) && length(intersect(hKEGGgenes, names(gene2Domains))) == 0)
		stop("There is a conflict between KEGG Entrez gene / open reading frame IDs and the gene identifiers in list 'gene2Domains'.\n 'gene2Domains' should be a list of Entrez gene IDs and their correponding InterPro domains.")
	features = gene2pathway:::getInterProDomains(hKEGGgenes, gene2Domains=gene2Domains, organism=organism)
	if(remove.duplicates)
		features = unique(features)
	common = intersect(hKEGGgenes, rownames(features))
	features = features[common,]
	features = features[,colSums(features) > 0]
	
	parentPaths = vector("list", length=length(childs))
	names(parentPaths) = childs
	labels = matrix(0, ncol=length(levels), nrow=nrow(features))
	colnames(labels) = levels 
	rownames(labels) = rownames(features)
	for(l in 1:length(clust)){		
		myparent = names(clust)[l]
		if(class(clust[[l]]) == "list"){
			for(ll in 1:length(clust[[l]])){ # put labels for each sub-class
				genesInClust = intersect(clust[[l]][[ll]], rownames(labels))		
				mypath = paste(myparent,ll,sep=".")
				labels[genesInClust, mypath] = 1		
				parentPaths[mypath] = myparent
			}
		}	
		genesInClust = intersect(unlist(clust[[l]]), rownames(labels))
		labels[genesInClust, myparent] = 1			
	}		
	labels = labels[,colSums(labels) > minnmap]	
	freq = table(unlist(parentPaths[colnames(labels)]))
	redchild = paste(names(freq[freq == 1]),"[1-9]",sep=".")
	labels = labels[,setdiff(1:ncol(labels), unlist(sapply(redchild, grep, colnames(labels))))]
	cat("Removing hierarchy branches with < ", minnmap, " mapping genes\n--->", ncol(labels), "hierarchy branches left\n")
	if(ncol(labels) == 0)
		stop("No sufficient number of mapping genes! Impossible to train model.")
		
	treesizes = double(ncol(labels))+1
	names(treesizes) = colnames(labels)
	freq = table(unlist(parentPaths[colnames(labels)]))
	treesizes[names(freq)] = freq
	treesizes = treesizes / sum(treesizes)
	
	parentsLev1 = which(colnames(labels) %in% names(clust))
	parentsLev2 = setdiff(1:ncol(labels), parentsLev1)
	parentsLev12 = c()
	elems = elems[colnames(labels)[parentsLev1]]
	
	list(features=features, labels=labels, treesizes=treesizes, parentPaths=parentPaths, parentsLev1=parentsLev1, parentsLev2=parentsLev2, parentsLev12=parentsLev12, elemIDs=elems)
}

