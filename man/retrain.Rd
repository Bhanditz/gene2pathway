\name{retrain}
\alias{retrain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Retrain classification model}
\description{
Retrains the hierarchical classification model. This way new information from InterPro and KEGG databases can be incorporated to give better predictions. Retraining should be done on a regular basis from time to time.
}
\usage{
retrain(minnmap=30, level1Only="Metabolism", level2Only="Genetic Information Processing", organism="hsa", gene2Domains=NULL, KEGG.package=FALSE, remove.duplicates=FALSE, use.bagging=TRUE, nbag=11)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{minnmap}{ prune hierarchy branches with < minnmap mapping genes}
  \item{level1Only}{ for these hierarchy branches only the first level is used }
  \item{level2Only}{ for these hierarchy branches only the first and the second levels are used }
  \item{organism}{KEGG letter code describing an organism.  Please refer to <URL:http://www.genome.jp/kegg-bin/create\_kegg\_menu> for a complete list of organisms (and their letter codes) supported by KEGG.}
  \item{gene2Domains}{By default associations between genes and InterPro domains are retrieved via biomaRt from Ensembl. Alternatively, the user can provide its own mapping of genes to InterPro domains in form of a list here (see details).}
  \item{KEGG.package}{Instead of retrieving information directly from KEGG, one can use the KEGG.db package instead, which is significantly faster. However, the KEGG.db package only supports a fraction of organisms so far. Please refer to the manual pages of the KEGG.db package for further information. Default: Don't use KEGG.db package}
  \item{remove.duplicates}{ remove genes having the same InterPro domains prior training. Default: Don't do this }
  \item{use.bagging}{ use bagging }
  \item{nbag}{number of models to average over}
}
\details{
A hierarchical classification model based on SVMs and a ranking perceptron algorithm is trained. This model is usually additionally bagged to improve prediction qualitiy. The method produces a "classificationModel\_[organism].rda" (e.g. "classificationModel\_hsa.rda") file, which should be stored in the package data directory. Once a new model has been trained, the complete package should be reloaded.

The current version of the KEGG hierarchy is always retrieved directly from KEGG via FTP. By default associations between genes and InterPro domains are retrieved automatically via biomaRt from Ensembl. Please refer to <URL:http://www.ebi.ac.uk/ensembl/> for a list of organisms supported by Ensembl. Alternatively to using Ensembl and biomaRt, the user can provide its own mapping of genes to InterPro domains in form of a list. This especially allows for using organisms, which are supported by KEGG, but not by Ensembl so far. The list has the form genes -> InterPro domains, and each list entry is named by the Entrez gene ID of the corresponding gene. This is, because KEGG uses Entrez gene IDs for the mapping genes -> KEGG pathways.
}
\value{
The model structure. See \code{\link{classificationModel}} for details.
}
\author{ Holger Froehlich }
\seealso{ \code{\link{gene2pathway}}, \code{\link{classificationModel}}}
\examples{
\dontrun{
	retrain(KEGG.package=TRUE, organism="dme") # retrain classification model for drosophila using information from the KEGG.db package
}
}
\keyword{ file }% at least one, from doc/KEYWORDS
