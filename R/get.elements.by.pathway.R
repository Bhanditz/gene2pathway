get.elements.by.pathway = function(pathway.id){    		
	KEGGserver <-  SOAPServer("http://soap.genome.jp/keggapi/request_v6.0.cgi")
	KEGGxmlns = "SOAP/KEGG"
	KEGGaction = "SOAP/KEGG"
	KEGGsoapns = "1.1"
	e = .SOAP(KEGGserver, "get_elements_by_pathway", .soapArgs = list(pathway_id = pathway.id),
	action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version = KEGGsoapns))
	names(e) = sapply(e, function(ee) ee$element_id)
	e
}
	