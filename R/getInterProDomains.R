getInterProDomains = function(geneIDs, gene2Domains=NULL, organism="hsa", alldoms=NULL, flyBase=FALSE){	
	geneIDs = as.character(geneIDs)
	if(!is.null(gene2Domains)){
		gene2Domains = gene2Domains[geneIDs]
	}
	else{		
		cat("Retrieving information from InterPro database for organism '", organism, "' via Ensembl ...\n")
		if(organism == "hsa")
			ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
		else if(organism == "dme")
			ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
		else if(organism == "mmu")
			ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
		else if(organism == "rno")	
			ensembl <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
		else if(organism == "sce")	
			ensembl <- useMart("ensembl", dataset = "scerevisiae_gene_ensembl")
		else if(organism == "oaa")	
			ensembl <- useMart("ensembl", dataset = "oanatinus_gene_ensembl")
		else if(organism == "ptr")	
			ensembl <- useMart("ensembl", dataset = "ptroglodytes_gene_ensembl")
		else if(organism == "gga")	
			ensembl <- useMart("ensembl", dataset = "ggallus_gene_ensembl")
		else if(organism == "mcc")	
			ensembl <- useMart("ensembl", dataset = "mmulatta_gene_ensembl")
		else if(organism == "bta")	
			ensembl <- useMart("ensembl", dataset = "btaurus_gene_ensembl")
		else if(organism == "cel")	
			ensembl <- useMart("ensembl", dataset = "celegans_gene_ensembl")
		else if(organism == "aga")	
			ensembl <- useMart("ensembl", dataset = "agambiae_gene_ensembl")
		else if(organism == "xtr")	
			ensembl <- useMart("ensembl", dataset = "xtropicalis_gene_ensembl")
		else if(organism == "dre")	
			ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")
		else if(organism == "mdo")	
			ensembl <- useMart("ensembl", dataset = "mdomestica_gene_ensembl")
		else if(organism == "cfa")	
			ensembl <- useMart("ensembl", dataset = "cfamiliaris_gene_ensembl")
		else{
			ensembl = useMart("ensembl")
			organisms = listDatasets(ensembl)
			orgs = sub("_gene_ensembl","",organisms[,1])
			if(organism %in% orgs){
				o = organisms[orgs == organism,1]
				class(o) = "character"
				ensembl = useDataset(o, ensembl)
			}
			else
				stop(paste("Genes for organism '", organism, "' cannot be mapped to InterPro via Ensembl!\nPlease provide your own mapping from genes to domains.\nPlease have also a look at <URL:http://www.ebi.ac.uk/ensembl/> for a list of organisms supported by Ensembl.",sep=""))
		}					
		if(!flyBase)
			tmp <- getBM(attributes=c("entrezgene", "interpro"), filters="entrezgene", values=geneIDs, mart = ensembl)
		else
			tmp <- getBM(attributes=c("flybase_gene_id", "interpro"), filters="flybase_gene_id", values=geneIDs, mart = ensembl)		
		if(is.null(tmp) | all(is.na(tmp[,"interpro"])))
			stop("No mapping Entrez gene ID -> InterPro found!")
		if(!flyBase)
			gene2Domains <- split(tmp$interpro, tmp$entrezgene, drop=FALSE)		
		else
			gene2Domains <- split(tmp$interpro, tmp$flybase_gene_id, drop=FALSE)		
		gene2Domains = gene2Domains[geneIDs]	
		empty = sapply(gene2Domains, function(x) all(x == ""))	
		gene2Domains = gene2Domains[!empty]
	}		
	
	interpro_hierarchy = getURL("ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt")
	interpro_hierarchy = unlist(strsplit(interpro_hierarchy,"\n"))
	doms = c()
	oldlevel = -2
	code = c()
	parentlist = list()
	for(l in 1:length(interpro_hierarchy)){
		s = unlist(strsplit(as.character(interpro_hierarchy[l]),"IPR",fixed=TRUE))	
		level = nchar(s[1])/2	
		dom = paste("IPR", unlist(strsplit(s[2], "::", fixed=TRUE))[1], sep="")
		doms = c(doms, dom)
		if(level > oldlevel)
			code = c(code, dom)
		else if(level < oldlevel){
			if(level > 0)
				code = c(code[1:level], dom)
			else
				code = dom
		}
		else
			code[length(code)] = dom
		parentlist[[l]] = code
		oldlevel = level
	}
	names(parentlist) = doms

	if(is.null(alldoms))
		alldoms = unique(c(unlist(gene2Domains),doms))
	features = matrix(0, ncol=length(alldoms), nrow=length(gene2Domains))
	dimnames(features) = list(names(gene2Domains), alldoms)
	for(g in 1:length(gene2Domains)){
		do = setdiff(gene2Domains[[g]],"")
		features[g,intersect(unique(c(do, unlist(parentlist[do]))), alldoms)] = 1		
	}	
	cat("done: Information found for ", length(gene2Domains), "out of ", length(geneIDs), " genes\n")	
	features
}
