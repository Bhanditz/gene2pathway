loss = function(y, u, treesizes, parentPaths){
	err = which(y != u)
	if(length(err) == 0)
		return(0)			
	roots = setdiff(names(err), names(parentPaths))
	ls = 0
	for(e in names(err)){		
		if(e %in% roots)
			ls = ls + treesizes[e]
		else{
			par = intersect(unlist(parentPaths[e]), names(u))
			if(all(y[par] == u[par]))
				ls = ls + treesizes[e]
		}
	}	
	ls
}
