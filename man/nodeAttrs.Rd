\name{nodeAttrs}
\alias{nodeAttrs}
\title{Generate node attributes}
\usage{
  nodeAttrs(nbAn, nbComp, labAn, labComp, file)
}
\arguments{
  \item{nbAn}{Number of analyses being considered, i.e
  number of IcaSet objects}

  \item{nbComp}{Number of components by analysis, if of
  length 1 then it is assumed that each analysis has the
  same number of components.}

  \item{labAn}{Labels of the analysis, if missing it will
  be generated as an1, an2, ...}

  \item{labComp}{List containing the component labels
  indexed by analysis, if missing will be generated as
  comp1, comp2, ...}

  \item{file}{File where the description of the node
  attributes will be written}
}
\value{
  A data.frame describing each node/component
}
\description{
  This function builds a data.frame describing for each
  node of the graph its ID and which analysis/data it comes
  from.
}
\details{
  The created file is used in Cytoscape.
}
\examples{
## 4 datasets, 20 components calculated in each dataset, labAn
nodeAttrs(nbAn=4, nbComp=20, labAn=c("tutu","titi","toto","tata"))
}
\author{
  Anne Biton
}

