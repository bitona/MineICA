\name{plotMclust}
\alias{plotMclust}
\title{Plots the Gaussian fitted  by \code{\link[mclust]{Mclust}}}
\usage{
  plotMclust(mc, data)
}
\arguments{
  \item{mc}{The result of Mclust function applied to
  argument \code{data}}

  \item{data}{The vector of numeric values on which was
  applied Mclust}
}
\value{
  NULL
}
\description{
  Given a result of function Mclust applied on a numeric
  vector, this function add the fitted Gaussian to a
  previous plot. This is an internal function called by
  \code{plotPosSamplesInComp}.
}
\details{
  This function can only deal with at the most three
  Gaussian.
}
\examples{
## create a mix of two Gaussian
v <-c(rnorm(80,mean=-0.5,sd=1),rnorm(80,mean=1,sd=0.2))
## apply Mclust
mc <- Mclust(v)
## plot fitted Gaussian on histogram of v
hist(v, freq=FALSE)
MineICA:::plotMclust(mc=mc,data=v)
}
\author{
  Anne Biton
}
\keyword{internal}

