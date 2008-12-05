getComponents = function(pathway.id){	
	cat("--> connected components for", pathway.id,"\n")
	rels = gene2pathway:::get.element.relations.by.pathway(pathway.id)		
	geneIDs=NULL
	comp=NULL
	if(nrow(rels) > 1){
		elems = gene2pathway:::get.elements.by.pathway(pathway.id)
		myelems = unique(c(rels[,1],rels[,2]))
		which.gene = sapply(myelems, function(e) elems[[e]]$type == "gene")
		myelems = myelems[which.gene]
		r = (rels[,1] %in% myelems) & (rels[,2] %in% myelems)
		rels = rels[r,,drop=FALSE]		
		if(nrow(rels) > 0){
			G = matrix(0, ncol=length(myelems), nrow=length(myelems))	
			dimnames(G) = list(myelems, myelems)		
			for(r in 1:nrow(rels)){
				G[as.character(rels[r,1]),as.character(rels[r,2])] = 1
				G[as.character(rels[r,2]),as.character(rels[r,1])] = 1
			}
			diag(G) = 1		
			gr = as(G, "graphNEL")
			comp = RBGL::connectedComp(gr)
			compGenes = sapply(comp, function(c) unlist(sapply(elems[c], function(e) e$name)))
# 			geneIDs = sapply(compGenes, function(g) sub(paste(organism,":",sep=""),"", g))
			organism = gsub("[0-9]","",strsplit(pathway.id, ":")[[1]][2])
			geneID.list = RCurl:::getURL(paste("ftp://ftp.genome.jp/pub/kegg/genes/organisms/", organism, "/", organism, "_ncbi-geneid.list",sep=""))
			geneIDs = sapply(compGenes, function(g) gene2pathway:::KEGG2Entrez(g, geneID.list=geneID.list, organism=organism))
		}	
	}			
	list(geneIDs=geneIDs, elemIDs=comp)
}

getComponents.internal = function(pathway.id, KEGG2Entrez.table){	
	cat("--> connected components for", pathway.id,"\n")
	rels = gene2pathway:::get.element.relations.by.pathway(pathway.id)		
	geneIDs=NULL
	comp=NULL
	if(nrow(rels) > 1){
		elems = gene2pathway:::get.elements.by.pathway(pathway.id)
		myelems = unique(c(rels[,1],rels[,2]))
		which.gene = sapply(myelems, function(e) elems[[e]]$type == "gene")
		myelems = myelems[which.gene]
		r = (rels[,1] %in% myelems) & (rels[,2] %in% myelems)
		rels = rels[r,,drop=FALSE]		
		if(nrow(rels) > 0){
			G = matrix(0, ncol=length(myelems), nrow=length(myelems))	
			dimnames(G) = list(myelems, myelems)		
			for(r in 1:nrow(rels)){
				G[as.character(rels[r,1]),as.character(rels[r,2])] = 1
				G[as.character(rels[r,2]),as.character(rels[r,1])] = 1
			}
			diag(G) = 1		
			gr = as(G, "graphNEL")
			comp = RBGL::connectedComp(gr)
			compGenes = sapply(comp, function(c) unlist(sapply(elems[c], function(e) e$name)))
			geneIDs = sapply(compGenes, function(g) KEGG2Entrez.table[KEGG2Entrez.table[,1] %in% g, 2])
		}	
	}			
	list(geneIDs=geneIDs, elemIDs=comp)
}
