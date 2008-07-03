\name{ORF2Entrez}
\alias{ORF2Entrez}
\alias{Entrez2ORF}

\title{ Conversion locus tag -> Entrez ID and Entrez ID -> locus tag according to KEGG }
\description{
  Conversion locus tag -> Entrez ID and Entrez ID -> locus tag according to KEGG (see <URL:ftp://ftp.genome.jp/pub/kegg/genes/organisms/>)
}
\usage{
	ORF2Entrez(ORFIDs, organism="dme")
	Entrez2ORF(entrezIDs, organism="dme")
}
\arguments{
	\item{ORFIDs}{character vector of locus tags / ORF identifiers}
	\item{organism}{KEGG letter code describing an organism.  Please refer to <URL:http://www.genome.jp/kegg-bin/create\_kegg\_menu> for a complete list of organisms (and their letter codes) supported by KEGG.}
	\item{entrezIDs}{character vector of Entrez gene IDs}
}
\value{
	character vector
}
\examples{
	entrez=ORF2Entrez("Dmel_CG4942", organism="dme")
	Entrez2ORF(entrez, organism="dme") # should be "Dmel_CG4942"
}
