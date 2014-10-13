\name{readA}
\alias{readA}
\title{read A}
\usage{
  readA(Afile, datfile, dat, annot = TRUE)
}
\arguments{
  \item{Afile}{The file which contains the matrix of sample
  contributions. It must be a txt file where the separator
  is \code{white space}, that is one or more spaces, tabs,
  newlines or carriage returns}

  \item{datfile}{The file which contains the matrix (of
  dimension features x samples) based on which the matrix A
  was calculated}

  \item{dat}{The data based on which the matrix A was
  calculated (features x samples)}

  \item{annot}{TRUE (default) if the Afile contains
  rownames of matrix A, FALSE if the rownames has to be
  extracted from dat}
}
\value{
  This function returns a matrix of dimension samples x
  components with rownames filled with sample IDs.
}
\description{
  readA
}
\details{
  This function reads and annotates matrix A.

  The matrix dat must be the one on which the matrix A was
  calculated. It is assumed that the number of components
  is lower than the number of samples, the matrix will be
  transposed to have dimension 'samples x components'
  according to this assumption. If \code{annot} is FALSE,
  colnames of dat are used to annotate rownames of A.
}
\author{
  Anne Biton
}
\keyword{internal}

