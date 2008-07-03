KEGG2Entrez = function(KEGGIDs=NULL, geneID.list=NULL, organism="hsa"){
	if(is.null(geneID.list))
		geneID.list = getURL(paste("ftp://ftp.genome.jp/pub/kegg/genes/organisms/", organism, "/", organism, "_ncbi-geneid.list",sep=""))
	geneID.list = unlist(strsplit(geneID.list,"\n"))
	geneID.list = unlist(strsplit(geneID.list,"\t"))
	geneID.list = matrix(geneID.list, ncol=2, byrow=TRUE)
	if(!is.null(KEGGIDs)){
		entrez.ids = geneID.list[geneID.list[,1] %in% KEGGIDs,2]	
		entrez.ids = sapply(entrez.ids, function(x) strsplit(x, ":")[[1]][2])
		return(entrez.ids)
	}
	else{
		entrez.ids = geneID.list[,2]	
		entrez.ids = sapply(entrez.ids, function(x) strsplit(x, ":")[[1]][2])
		return(cbind(geneID.list[,1], entrez.ids))
	}
	
}

Entrez2ORF.internal = function(entrezIDs, KEGG2Entrez.tab, organism="dme"){	
	ORFIDs = KEGG2Entrez.tab[KEGG2Entrez.tab[,2] %in% entrezIDs,1]
	ORFIDs = sapply(ORFIDs, function(x) strsplit(x, ":")[[1]][2])
	return(ORFIDs)
}


ORF2Entrez = function(ORFIDs, organism="dme"){
	gene2pathway:::KEGG2Entrez(paste(organism,":",ORFIDs,sep=""), organism=organism)
}

Entrez2ORF = function(entrezIDs, organism="dme"){
	KEGG2Entrez.tab = gene2pathway:::KEGG2Entrez(organism=organism)
	ORFIDs = KEGG2Entrez.tab[KEGG2Entrez.tab[,2] %in% entrezIDs,1]
	ORFIDs = sapply(ORFIDs, function(x) strsplit(x, ":")[[1]][2])
	return(ORFIDs)
}
