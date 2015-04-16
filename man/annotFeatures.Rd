\name{annotFeatures}
\alias{annotFeatures}
\title{Annotation of features using an annotation package}
\usage{
  annotFeatures(features, type, annotation)
}
\arguments{
  \item{features}{Feature IDs to be annotated}

  \item{type}{The object from the package used to annotate
  the features, must be available in
  \code{ls("package:package_name")}}

  \item{annotation}{An annotation package}
}
\value{
  A vector of gene/object IDs indexed by the feature IDs.
}
\description{
  This function annotates a set of features
}
\examples{
library(hgu133a.db)
annotFeatures(features = c("1007_s_at", "1053_at", "117_at", "121_at", "1255_g_at"),
              type="SYMBOL", annotation="hgu133a.db")
}
\author{
  Anne Biton
}

