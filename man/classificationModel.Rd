\name{classificationModel}
\alias{classificationModel}
\alias{classificationModel_hsa}
\alias{classificationModel_dme}
\alias{modelKEGG}
\docType{data}	
\title{Hierarchical Classification Model}
\description{
  This file contains the hierarchical classification model to predict KEGG pathway branches for genes. The model uses a pruned KEGG hierarchy, where metabolic pathways are not distinguished further, and the KEGG hierarchy for "cellular processes" and "genetic information processing" is pruned at the 2nd level. By default the model uses bagging to improve prediction accuracy. Important: There exists one separate model file for each organism.
}	
\format{List of class "model", where each model has the following entries:
 \describe{
  \item{W}{learned decision hyperplane normal vector}
  \item{C}{dictionary of label vectors, which can be predicted individually or which can be used to predict combinations of them}
  \item{detectors}{SVM models trained to separate one specific pathway branch from the rest of the hierarchy}
  \item{used\_domains}{InterPro domains used by the classifier to separate the specific branch from the rest of the hierarchy}
  \item{alldomains}{all InterPro domains used to build feature vectors}
  \item{allpathways}{hierarchy branches, which can be predicted}
  \item{treesizes}{relative size of hierarchy below the corresponding branch}
  \item{kegg\_hierarchy}{a nested list with information (parent branches, pathway names, pathway IDs, hierarchy level) on all higher hierarchy branches for each pathway}
 }
}		
\author{Holger Froehlich}
\seealso{\code{\link{classificationModelSignalTrans}}}
\keyword{datasets}
