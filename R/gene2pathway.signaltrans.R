gene2pathway.signaltrans = function(geneIDs=NULL, flyBase=FALSE, gene2Domains=NULL, organism="hsa", useKEGG=TRUE){	
	if(is.null(geneIDs) & is.null(gene2Domains))
		stop("You have to provide either a list of Entrez gene IDs or a mapping of genes to InterPro domains")	
	if(!exists("gene2pathwayEnv"))
		assign("gene2pathwayEnv",new.env(parent=globalenv()),envir=.GlobalEnv)
	if(exists("organismSignalTrans", envir=gene2pathwayEnv))
		old.organism = get("organismSignalTrans", envir=gene2pathwayEnv)
	else
		old.organism = organism
	assign("organismSignalTrans", organism, envir=gene2pathwayEnv)
	if(!is.null(gene2Domains))	
		geneIDs = names(gene2Domains)
	else
		geneIDs = as.character(geneIDs) 
	if(!exists("modelSignalTrans", envir=gene2pathwayEnv) | (organism != old.organism)){
		cat("Loading classification model ...\n")
		myfile = paste("classificationModelSignalTrans_",organism,sep="")
		tryCatch(data(list=myfile,package="gene2pathway", envir=gene2pathwayEnv), warning=function(w) stop("No Model for organism '", organism, "' available.\nPlease invoke 'retrain.signaltrans' to generate one."))
# 		data(list=myfile,package="gene2pathway", envir=gene2pathwayEnv)
# 		load(paste("classificationModelSignalTrans_",organism,".rda",sep=""), envir=gene2pathwayEnv)
	}
	model = get("modelSignalTrans", envir=gene2pathwayEnv)		
	if(class(model) == "model"){
		alldomains = model$alldomains
		elemIDs = model$elemIDs
		pathways = model$allpathways		
		parentPaths = model$parentPaths
	}
	else{
		alldomains = model[[1]]$alldomains		
		elemIDs = model[[1]]$elemIDs
		pathways = model[[1]]$allpathways
		parentPaths = model[[1]]$parentPaths
	}
	roots = unique(unlist(parentPaths))
	if(useKEGG){			
		cat("Using KEGG information from SOAP service ...\n")			
		organisms=list.organisms()
		if(!(organism %in% names(organisms)))
			stop(paste("Organism '", organism, "' unknown in KEGG! Please refer to <URL:http://www.genome.jp/kegg-bin/create_kegg_menu> for a complete list of supported organisms."))
		cat("Mapping to signal transduction pathway components via KEGG database ...\n",sep="")
# 		if(!(organism %in% c("hsa", "mmu", "rno"))){
			if(flyBase){
				geneIDs = unlist(AnnotationDbi::mget(geneIDs, org.Dm.egFLYBASE2EG, ifnotfound=NA))
				geneIDs = geneIDs[!is.na(geneIDs)]
				if(length(geneIDs) == 0)
					stop("No mapping FlyBase -> Entrez gene ID found!")
				flyBase = FALSE			
			}
			if(organism != "hsa"){
				KEGG2Entrez.tab = gene2pathway:::KEGG2Entrez(organism=organism)
				geneIDs.conv = gene2pathway:::Entrez2ORF.internal(geneIDs, KEGG2Entrez.tab, organism=organism)	
			}
			else{
				geneIDs.conv = geneIDs
				KEGG2Entrez.tab = NULL
			}
		comp2gene = sapply(roots, function(p){
			elems = gene2pathway:::get.elements.by.pathway(paste("path:",organism,p,sep=""))
			compGenes = sapply(elemIDs[[p]], function(c) unlist(sapply(elems[c], function(e) e$names)))
			sapply(compGenes, function(g) sub(paste(organism,":",sep=""),"", g))
		})		
		comp2gene = comp2gene[roots]		
		anno.genes = geneIDs.conv[geneIDs.conv %in% unique(unlist(comp2gene))]
		cat("---> Information found for ", length(anno.genes), "genes\n")
		KEGGgenes = lapply(anno.genes, function(kg){
			unlist(sapply(1:length(comp2gene), function(i){
				comp = sapply(comp2gene[[i]], function(cc) any(cc %in% kg))
				if(any(comp)){
					compname = intersect(paste(roots[i],which(comp),sep="."), pathways)
					if(length(compname) > 0)
						return(c(roots[i], compname))
					else
						return(roots[i])
				}				
			}))				
		})		
		if(length(grep(":", names(KEGGgenes))) > 0){
			if(is.null(KEGG2Entrez.tab))
				anno.genes = sub(paste(organism,":",sep=""), "", names(KEGGgenes))
			else
				anno.genes = gene2pathway:::KEGG2Entrez(names(KEGGgenes), geneID.list=KEGG2Entrez.tab, organism=organism)
			names(KEGGgenes) = anno.genes			
		}
		cat("done.\n")		
		if(!is.null(gene2Domains)  && length(intersect(names(KEGGgenes), geneIDs)) == 0)
			warning("There may be a conflict between KEGG Entrez gene / open reading frame IDs and the gene identifiers in list 'gene2Domains'.\n Gene identifiers in 'gene2Domains' should be Entrez gene IDs, if useKEGG=TRUE.")
	}
	else
		KEGGgenes = list()		
	topredict = setdiff(geneIDs, names(KEGGgenes))
	if(length(topredict) > 0){		
		cat(length(topredict), " genes to predict\n")		
		features = gene2pathway:::getInterProDomains(topredict, gene2Domains=gene2Domains, alldoms=alldomains, organism=organism, flyBase=flyBase)	
		
		cat("Model prediction possible for ", nrow(features), " genes ")
		gene2Pathall = gene2pathway:::predict.gene2pathway(model, features)
		cat("done\nPreparing output\n")
		gene2Path = gene2Pathall$gene2Path	
		scorestmp = gene2Pathall$scores
		names(scorestmp) = rownames(features)
	}
	else{
		scorestmp = NA
		gene2Path = list()
	}
	
	totallist = vector("list", length=length(geneIDs))	
	scores = vector("list", length=length(geneIDs))	
	components = vector("list", length=length(geneIDs))
	names(totallist) = geneIDs
	names(scores) = geneIDs
	names(components) = geneIDs	
	totallist[names(gene2Path)] = gene2Path	
	scores[names(gene2Path)] = scorestmp	
	byKEGG = rep(FALSE, length(totallist))
	names(byKEGG) = names(totallist)
	if(useKEGG){
		totallist[names(KEGGgenes)] = KEGGgenes
		byKEGG[names(KEGGgenes)] = TRUE
	}	
	components[names(gene2Path)] = sapply(gene2Path, function(gp){			
		 lapply(intersect(gp, names(parentPaths)), function(gpp){			
			p = unlist(strsplit(gpp, "\\."))
			unlist(elemIDs[[p[1]]][p[2]])			
		})				
	})

	if(exists("kegg_hierarchy", envir=gene2pathwayEnv))
		kegg_hierarchy = get("kegg_hierarchy", envir=gene2pathwayEnv)
	else{
		kegg_hierarchy = gene2pathway:::getKEGGHierarchy(level1Only=c(), level2Only=c())
		assign("kegg_hierarchy", kegg_hierarchy, envir=gene2pathwayEnv)		
	}
	pathnames = c(kegg_hierarchy$pathNamesLev1, kegg_hierarchy$pathNamesLev2, kegg_hierarchy$pathNamesLev3)		
	totallist[totallist == ""] = NA
	totallist[is.na(names(totallist))] = NA
	gene2Pathname = sapply(totallist, function(gp) setdiff(pathnames[as.character(gp)], NA))	
	names(gene2Pathname) = names(totallist)
	
	cat("finished\n")
	list(gene2Path=totallist, gene2Pathname=gene2Pathname, scores=scores, byKEGG=byKEGG, elemIDs=components)
}