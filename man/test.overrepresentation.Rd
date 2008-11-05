\name{test.overrepresentation}
\alias{test.overrepresentation}

\title{ Test statistical overrepresentation of KEGG pathways in a list of genes }
\description{
Test the statistical overrepresentation of KEGG pathways in a group of genes using Fisher's exact test. The analysis can either be based on all KEGG pathways predicted by \code{\link{gene2pathway}}/\code{\link{gene2pathway.signaltrans}} or on original KEGG annotation only.
}
\usage{
test.overrepresentation(genesOfInterest, predpath, KEGGonly=FALSE, cutoff=0.1)
}
\arguments{
\item{genesOfInterest}{a character vector of gene identifiers (see \code{\link{gene2pathway}}, \code{\link{gene2pathway.signaltrans}}) for a gene list of interest}
\item{predpath}{predictions of \code{\link{gene2pathway}} or \code{\link{gene2pathway.signaltrans}}}
\item{KEGGonly}{use KEGG annotation only}
\item{cutoff}{p-value significance cutoff}
}
\value{
	Table with two columns: KEGG pathway and adjusted p-value (adjustment according to Benjamini-Yekutieli)
}
