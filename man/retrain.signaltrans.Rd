\name{retrain.signaltrans}
\alias{retrain.signaltrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Retrain classification model for signaling pathways}
\description{
Retrains the hierarchical classification model for signaling pathway components. This way new information from InterPro and KEGG databases can be incorporated to give better predictions. Retraining should be done on a regular basis from time to time.
}
\usage{
retrain.signaltrans(minnmap=10, organism="hsa", gene2Domains=NULL, remove.duplicates=FALSE, use.bagging=TRUE, nbag=11)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{minnmap}{ prune hierarchy branches with < minnmap mapping genes}
  \item{organism}{KEGG letter code describing an organism.  Please refer to <URL:http://www.genome.jp/kegg-bin/create\_kegg\_menu> for a complete list of organisms (and their letter codes) supported by KEGG.}
  \item{gene2Domains}{By default associations between genes and InterPro domains are retrieved via biomaRt from Ensembl. Alternatively, the user can provide its own mapping of genes to InterPro domains in form of a list here (see details).}
  \item{remove.duplicates}{ remove genes having the same InterPro domains prior training }
  \item{use.bagging}{ use bagging }
  \item{nbag}{number of models to average over}
}
\details{
A hierarchical classification model based on SVMs and a ranking perceptron algorithm is trained. This model is usually additionally bagged to improve prediction qualitiy. The method produces a "classificationModelSignalTrans\_[organism].rda" (e.g. "classificationModelSignalTrans\_hsa.rda") file, which should be stored in the package data directory. Once a new model has been trained, the complete package should be reloaded.

The current version of the KEGG hierarchy is always retrieved directly from KEGG via FTP. Labels for the training set are obtained via the function \code{\link{getComponents}}, which uses the KEGG SOAP service. By default associations between genes and InterPro domains are retrieved automatically via biomaRt from Ensembl. Please refer to <URL:http://www.ebi.ac.uk/ensembl/> for a list of organisms supported by Ensembl. Alternatively to using Ensembl and biomaRt, the user can provide its own mapping of genes to InterPro domains in form of a list. This especially allows for using organisms, which are supported by KEGG, but not by Ensembl so far. The list has the form genes -> InterPro domains, and each list entry is named by the Entrez gene ID of the corresponding gene. This is, because KEGG uses Entrez gene IDs for the mapping genes -> KEGG pathways.
}
\value{
The model structure. See \code{\link{classificationModelSignalTrans}} for details.
}
\author{ Holger Froehlich }
\seealso{ \code{\link{gene2pathway.signaltrans}}, \code{\link{classificationModelSignalTrans}}}
\keyword{ file }% at least one, from doc/KEYWORDS