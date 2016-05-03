\name{annotReciprocal}
\alias{annotReciprocal}
\title{annotReciprocal}
\usage{
  annotReciprocal(dataGraph, file,
    keepOnlyReciprocal = FALSE)
}
\arguments{
  \item{dataGraph}{data.frame which contains the graph
  description, must have two columns \code{n1} and
  \code{n2} filled with node IDs, each row denoting there
  is an edge from \code{n1} to \code{n2}.}

  \item{file}{file where the graph description is written}

  \item{keepOnlyReciprocal}{if TRUE \code{dataGraph} is
  restricted to reciprocal edges, else all edges are kept
  (default).}
}
\value{
  This function returns the argument \code{dataGraph} with
  an additional column named 'reciprocal' which contains
  TRUE if the edge described by the row is reciprocal, and
  FALSE if it is not reciprocal.
}
\description{
  This function notes edges of a graph as reciprocal or
  not.
}
\examples{
dg <- data.frame(n1=c("A","B","B","C","C","D","E","F"),n2=c("B","G","A","B","D","C","F","E"))
annotReciprocal(dataGraph=dg)
}
\author{
  Anne Biton
}

