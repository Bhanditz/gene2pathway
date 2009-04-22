.onLoad <- function(lib, pkgname, where){
#   library.dynam(pkgname, pkgname, lib)	
	assign("gene2pathwayEnv",new.env(parent=globalenv()),envir=.GlobalEnv) 	
}
