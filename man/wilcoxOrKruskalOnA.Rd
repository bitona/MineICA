\name{wilcoxOrKruskalOnA}
\alias{wilcoxOrKruskalOnA}
\title{Comparison of distributions of sample groups}
\usage{
  wilcoxOrKruskalOnA(A, colAnnot, annot)
}
\arguments{
  \item{A}{A matrix of dimensions 'samples x components'
  containing the sample contributions}

  \item{annot}{A matrix of dimensions 'samples x variables'
  containing the sample annotations}

  \item{colAnnot}{The name of the column of \code{annot} to
  be considered}
}
\value{
  A vector of p-values
}
\description{
  Compare the sample contributions according to their
  annotation level across the components.
}
\details{
  Wilcoxon or Kruskal-Wallis tests are performed depending
  on the number of levels in the considered annotation.
}
\author{
  Anne Biton
}
\seealso{
  \code{wilcox.test}, \code{kruskal.test}
}
\keyword{internal}

