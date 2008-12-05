\name{classificationModelSignalTrans}
\alias{classificationModelSignalTrans}
\alias{modelSignalTrans}
\docType{data}	
\title{Hierarchical Classification Model for Signaling Transduction Pathways and Pathway Components}
\description{
  This file contains the hierarchical classification model to predict KEGG signaling pathways and pathway components for genes. The model contains only pathway components, to which a specified minimum number of genes could be mapped in the training phase (see \code{\link{retrain.signaltrans}}). Important: There exists one separate model file for each organism.
}	
\format{List of class "model", where each model has the following entries:
 \describe{
  \item{W}{learned decision hyperplane normal vector}
  \item{C}{dictionary of label vectors, which can be predicted individually or which can be used to predict combinations of them}
  \item{detectors}{SVM models trained to separate one specific pathway branch from the rest of the hierarchy}
  \item{used\_domains}{InterPro domains used by the classifier to separate the specific branch from the rest of the hierarchy}
  \item{alldomains}{all InterPro domains used to build feature vectors}
  \item{allpathways}{hierarchy branches, which can be predicted}
  \item{treesizes"}{relative size of hierarchy below the corresponding branch}
  \item{kegg\_hierarchy}{a nested list with information (parent branches, pathway names, pathway IDs, hierarchy level) on all higher hierarchy branches for each pathway}
  \item{elemIDs}{a list of KEGG element IDs mapping to each pathway component -
may be used to highlight pathway components with
\code{\link{color.pathway.by.elements}}.}
 }
}		
\author{Holger Froehlich}
\seealso{\code{\link{classificationModel}}}
\keyword{datasets}
