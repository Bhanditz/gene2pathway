get.elements.by.pathway = function(pathway.id){    		
	KEGGserver <-  SSOAP:::SOAPServer("http://soap.genome.jp/keggapi/request_v6.2.cgi")
	KEGGxmlns = "SOAP/KEGG"
	KEGGaction = "SOAP/KEGG"
	KEGGsoapns = "1.1"
	e = SSOAP:::.SOAP(KEGGserver, "get_elements_by_pathway", .soapArgs = list(pathway_id = pathway.id),
	action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SSOAP:::SOAPNameSpaces(version = KEGGsoapns))
	names(e) = sapply(e, function(ee) ee$element_id)
	e
}
	