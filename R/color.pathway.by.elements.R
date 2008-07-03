color.pathway.by.elements = function(pathway.id, elements){    		
	KEGGserver <-  SOAPServer("http://soap.genome.jp/keggapi/request_v6.0.cgi")
	KEGGxmlns = "SOAP/KEGG"
	KEGGaction = "SOAP/KEGG"
	KEGGsoapns = "1.1"
	.SOAP(KEGGserver, "color_pathway_by_elements", .soapArgs = list(pathway_id = pathway.id, element_id_list=elements, rep("red", length(elements)), rep("yellow", length(elements))),
	action = KEGGaction, xmlns = KEGGxmlns, nameSpaces = SOAPNameSpaces(version = KEGGsoapns))	
}
