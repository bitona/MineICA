\name{getSdExpr}
\alias{getSdExpr}
\title{getSdExpr}
\usage{
  getSdExpr(features, dat)
}
\arguments{
  \item{features}{IDs}

  \item{dat}{Expression data indexed by IDs}
}
\value{
  Returns a vector
}
\description{
  Compute standard deviation of the gene expression
}
\examples{
dat <- matrix(rnorm(1000),ncol=10,nrow=100)
rownames(dat) <- 1:100
MineICA:::getSdExpr(features = 2:20, dat = dat)
}
\author{
  Anne Biton
}
\keyword{internal}

