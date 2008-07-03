\name{getComponents}
\alias{getComponents}
\alias{get.elements.by.pathway}
\alias{color.pathway.by.elements}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{KEGG pathway information}
\description{
1. get connected pathway components; 2. get all elements of a given pathway; 3. color certain elements in a pathway.
}
\usage{
getComponents(pathway.id)

get.elements.by.pathway(pathway.id)

color.pathway.by.elements(pathway.id, elements)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{pathway.id}{KEGG pathway ID, e.g. "path:hsa04012"}
\item{elements}{KEGG element IDs: character vector of numbers}
}
\details{
All functions use the KEGG SOAP service.
}
\value{
getComponents: a list with the entries
\item{geneIDs}{Entrez gene IDs mapping to each pathway component}
\item{elemIDs}{KEGG element IDs mapping to each pathway component}

get.elements.by.pathway: list, see <URL http://www.genome.jp/kegg/soap/doc/keggapi\_manual.html> for details

color.pathway.by.elements: an URL of a colored gif file, see <URL http://www.genome.jp/kegg/soap/doc/keggapi\_manual.html> for details
}
\author{ Holger Froehlich }
\examples{
\dontrun{
  comp = getComponents("path:hsa04020") # get all connected components
  color.pathway.by.elements("path:hsa04020", comp$elemIDs[[1]]) # mark first component
}
}
\keyword{ file }% at least one, from doc/KEYWORDS
