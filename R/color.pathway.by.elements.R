color.pathway.by.elements = function(pathway.id, elements){    		
	KEGGserver <-  SSOAP:::SOAPServer("http://soap.genome.jp/keggapi/request_v6.2.cgi")
	KEGGxmlns = "SOAP/KEGG"
	KEGGaction = "SOAP/KEGG"
	KEGGsoapns = "1.1"
	names(elements) = NULL
	SSOAP:::.SOAP(KEGGserver, "get_html_of_colored_pathway_by_elements", .soapArgs = list(pathway_id = pathway.id, element_id_list=elements, fg_list=rep("red", length(elements)), bg_list=rep("yellow", length(elements))),
	action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SSOAP:::SOAPNameSpaces(version = KEGGsoapns))	
}
