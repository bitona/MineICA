\name{relativePath}
\alias{relativePath}
\title{Relative path}
\usage{
  relativePath(path1, path2)
}
\arguments{
  \item{path1}{The first path}

  \item{path2}{The second path}
}
\value{
  The relative path between path1 and path2
}
\description{
  Computes the relative path between two imbricated paths
}
\details{
  \code{path1} and \code{path2} must be imbricated.
}
\examples{
path1 <- "home/lulu/res/gene2comp/"
path2 <- "home/lulu/res/comp2annot/invasive/"
relativePath(path1,path2)
}
\author{
  Anne Biton
}

