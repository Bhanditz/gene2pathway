gene2pathway = function(geneIDs=NULL, flyBase=FALSE, gene2Domains=NULL, organism="hsa", useKEGG=TRUE, KEGG.package=TRUE){
	if(is.null(geneIDs) & is.null(gene2Domains))
		stop("You have to provide either a list of Entrez gene IDs or a mapping of genes to InterPro domains")
	if(!exists("gene2pathwayEnv"))
		assign("gene2pathwayEnv",new.env(parent=globalenv()),envir=.GlobalEnv) 
	if(exists("organismKEGG", envir=gene2pathwayEnv))
		old.organism = get("organismKEGG", envir=gene2pathwayEnv)
	else
		old.organism = organism	
	assign("organismKEGG", organism, envir=gene2pathwayEnv)
	if(!is.null(gene2Domains))	
		geneIDs = names(gene2Domains)
	else
		geneIDs = as.character(geneIDs)		
	if(!exists("modelKEGG", envir=gene2pathwayEnv) | (organism != old.organism)){
		cat("Loading classification model ...\n")
		myfile = paste("classificationModel_",organism,sep="")		
		tryCatch(data(list=myfile,package="gene2pathway", envir=gene2pathwayEnv), warning=function(w) stop("No Model for organism '", organism, "' available.\nPlease invoke 'retrain' to generate one."))
# 		data(list=myfile,package="gene2pathway", envir=gene2pathwayEnv)
# 		load(paste("classificationModel_",organism,".rda",sep=""), envir=gene2pathwayEnv)		
	}
	model = get("modelKEGG", envir=gene2pathwayEnv)	
	if(class(model) == "model"){
		alldomains = model$alldomains		
		pathways = model$allpathways		
		parentPaths = model$parentPaths		
	}
	else{
		alldomains = model[[1]]$alldomains				
		pathways = model[[1]]$allpathways
		parentPaths = model[[1]]$parentPaths
	}
	if(useKEGG){		
		if(KEGG.package){
			cat("Using KEGG information from KEGG.db package ...\n")			
			organisms = unique(gsub("[0-9]*","",AnnotationDbi::ls(KEGGPATHID2EXTID)))
			if(!(organism %in% organisms))
				stop(paste("Organism '", organism, "' is unknown in KEGG package! Please retry with KEGG.package=FALSE (slow). Please refer also to <URL:http://www.genome.jp/kegg-bin/create_kegg_menu> for a complete list of organisms supported by KEGG.", sep=""))	
			if((organism == "dme") & !flyBase){
				geneIDs = unlist(AnnotationDbi::mget(geneIDs, org.Dm.egFLYBASE, ifnotfound=NA))
				geneIDs = geneIDs[!is.na(geneIDs)]
				if(length(geneIDs) == 0)
					stop("No mapping Entrez gene ID -> FlyBase found!")
				flyBase = TRUE
			}
		}		
		else{
			cat("Using KEGG information from SOAP service ...\n")
			organisms = list.organisms()
			if(!(organism %in% names(organisms)))
				stop(paste("Organism '", organism, "' is unknown in KEGG! Please refer to <URL:http://www.genome.jp/kegg-bin/create_kegg_menu> for a complete list of supported organisms.",sep=""))	
			if(flyBase){
				geneIDs = unlist(AnnotationDbi::mget(geneIDs, org.Dm.egFLYBASE2EG, ifnotfound=NA))
				geneIDs = geneIDs[!is.na(geneIDs)]
				if(length(geneIDs) == 0)
					stop("No mapping FlyBase -> Entrez gene ID found!")				
				flyBase = FALSE
			}
		}
		if(!(organism %in% c("hsa", "mmu", "rno")) & !flyBase){			
			KEGG2Entrez.tab = gene2pathway:::KEGG2Entrez(organism=organism)
			geneIDs.conv = gene2pathway:::Entrez2ORF.internal(geneIDs, KEGG2Entrez.tab, organism=organism)	
		}
		else
			geneIDs.conv = geneIDs
		if(KEGG.package){ # fast						
			KEGGgenes = AnnotationDbi::mget(geneIDs.conv, KEGGEXTID2PATHID, ifnotfound=NA)		
			KEGGgenes = lapply(KEGGgenes, function(kg) sub(organism,"",kg))
			KEGGgenes = KEGGgenes[!is.na(KEGGgenes)]
			cat("Information from KEGG package available for ", length(KEGGgenes), " genes ...\n")
			KEGGgenes = lapply(KEGGgenes, function(kg) unique(c(kg,  unlist(parentPaths[kg]))))
		}
		else{ # slow			
			pathways = names(list.pathways(organism))
			path2Genes = try(lapply(pathways, get.genes.by.pathway))
			names(path2Genes) = pathways
			hKEGGgenes <- unique(unlist(path2Genes, use.names=FALSE))		
			hKEGGgenes <-  hKEGGgenes[!is.na(hKEGGgenes)]					
			genes2Path = sapply(hKEGGgenes, function(g) names(path2Genes)[sapply(path2Genes, function(p) any(g == p))]) # this is probably faster than calling get.pathways.by.genes
			genes2Path = sapply(genes2Path, function(g) sapply(g, function(gg) sub(paste("path:",organism,sep=""),"", gg)))
					
			KEGGgenes = genes2Path[paste(organism,":",geneIDs.conv,sep="")]
			KEGGgenes = KEGGgenes[!is.na(KEGGgenes)]
			cat("Information from KEGG package available for ", length(KEGGgenes), " genes ...\n")
			KEGGgenes = lapply(KEGGgenes, function(kg) unique(c(kg,  unlist(parentPaths[kg]))))		
		}
		if(!(organism %in% c("hsa", "mmu", "rno")) & !flyBase){
			anno.genes = gene2pathway:::KEGG2Entrez(names(KEGGgenes), geneID.list=KEGG2Entrez.tab, organism=organism)
			names(KEGGgenes) = anno.genes			
		}
		if(!is.null(gene2Domains)  && length(intersect(names(KEGGgenes), geneIDs)) == 0)
			warning("There may be a conflict between KEGG Entrez gene IDs and the gene identifiers in list 'gene2Domains'.\n Gene identifiers in 'gene2Domains' should be Entrez gene IDs, if useKEGG=TRUE.")
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
	names(totallist) = geneIDs
	names(scores) = geneIDs
	totallist[names(gene2Path)] = gene2Path	
	scores[names(gene2Path)] = scorestmp
	byKEGG = rep(FALSE, length(totallist))
	names(byKEGG) = names(totallist)
	if(useKEGG){
		totallist[names(KEGGgenes)] = KEGGgenes
		byKEGG[names(KEGGgenes)] = TRUE
	}
	if(exists("kegg_hierarchy", envir=gene2pathwayEnv))
		kegg_hierarchy = get("kegg_hierarchy", envir=gene2pathwayEnv)
	else{
		kegg_hierarchy = gene2pathway:::getKEGGHierarchy(level1Only=c(), level2Only=c())
		assign("kegg_hierarchy", kegg_hierarchy, envir=gene2pathwayEnv)		
	}
	pathnames = c(kegg_hierarchy$pathNamesLev1, kegg_hierarchy$pathNamesLev2, kegg_hierarchy$pathNamesLev3)		
	totallist[totallist == ""] = NA
	totallist[is.na(names(totallist))] = NA
	gene2Pathname = sapply(totallist, function(gp) pathnames[as.character(gp)])	
	names(gene2Pathname) = names(totallist)
	cat("finished\n")
	list(gene2Path=totallist, gene2Pathname=gene2Pathname, byKEGG=byKEGG, scores=scores)
}
