\name{gene2pathway.signaltrans}
\alias{gene2pathway.signaltrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pathway membership prediction}
\description{
Predicts a gene's membership to a KEGG signaling pathway and/or pathway component via the contained InterPro domains.
}
\usage{
gene2pathway.signaltrans(geneIDs=NULL, flyBase=FALSE, gene2Domains=NULL, organism="hsa", useKEGG=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{geneIDs}{ a character vector of Entrez gene IDs or FlyBase identifiers (not necessary, if the argument gene2Domains is provided)}
  \item{flyBase}{ Are FlyBase identifiers provided? Default: No}
  \item{gene2Domains}{By default associations between genes and InterPro domains are retrieved via biomaRt from Ensembl. Alternatively, the user can provide its own mapping of genes to InterPro domains in form of a list here (see details).}
  \item{organism}{KEGG letter code describing an organism.  Please refer to <URL:http://www.genome.jp/kegg-bin/create\_kegg\_menu> for a complete list of organisms (and their letter codes) supported by KEGG.}
  \item{useKEGG}{ Should KEGG information instead of a prediction be used when possible?}  
}
\details{
A hierarchical classification model based on SVMs and a ranking perceptron is used. This model is usually additionally bagged to improve prediction qualitiy. The model is stored in the package data directory and is recommended to be retrained from time to time.

The current version of the KEGG hierarchy is always retrieved directly from KEGG via FTP. By default associations between genes and InterPro domains are retrieved automatically via biomaRt from Ensembl. Please refer to <URL:http://www.ebi.ac.uk/ensembl/> for a list of organisms supported by Ensembl. Alternatively to using Ensembl and biomaRt, the user can provide its own mapping of genes to InterPro domains in form of a list. This especially allows for using organisms, which are supported by KEGG, but not by Ensembl so far. The list has the form genes -> InterPro domains, and each list entry is named by a gene identifier of the corresponding gene. If useKEGG=TRUE, Entrez gene IDs or FlyBase identifiers have to be used. Otherwise, arbitrary identifiers can are allowed.
}
\value{
\item{gene2Path}{mapping of gene IDs to corresponding KEGG pathway IDs}
\item{gene2Pathname}{mapping of gene IDs to corresponding KEGG pathway names}
\item{byKEGG}{inticates by TRUE/FALSE for each gene whether the mapping information was obtained directly from KEGG or whether it was predicted}
\item{scores}{confindence scores for the prediction (0, if no prediction was performed): see notes for details}
\item{elemIDs}{KEGG elements mapping to the corresponding predicted pathway components, if there are any, otherwise NULL. May be used to highlight pathway components with \code{\link{color.pathway.by.elements}}.}
}
\author{ Holger Froehlich }
\note{ By default a bagged model prediction is used, i.e. each of the individual sub-models is giving a vote for a specific output. The final output is determined by the majority of the votes for each hierarchy branch separately. The corresponding fraction voting for a specific branch may be interpreted as its probability. In the ideal case all individual branch probabilites should always be close to 1, if the gene maps to that part of the KEGG hierarchy, and close to 0 otherwise. A cumulative measure of confindence is thus the average over all probabilities > 0.5 and one minus the average over all probabilites < 0.5. We combine both measure by taking the average of both and report it as a reliability score.

If the user decides to retrain a model WITHOUT using bagging, then the reliability score is simply the margin between the highest and the second highest ranked solution. This margin should be larger 2 for good confindence.
}
\seealso{ \code{\link{retrain.signaltrans}}, \code{\link{classificationModelSignalTrans}}}
\examples{
\dontrun{
 gene2pathway.signaltrans("1443") 
}
}
\keyword{ file }% at least one, from doc/KEYWORDS
