\name{selectFeatures_IQR}
\alias{selectFeatures_IQR}
\title{Selection of features based on their IQR}
\usage{
  selectFeatures_IQR(data, nb)
}
\arguments{
  \item{data}{Measured data of dimension features x samples
  (e.g, gene expression data)}

  \item{nb}{The number of features to be selected}
}
\value{
  A subset of \code{data} restricted to the features having
  the \code{nb} highest IQR value
}
\description{
  This function selects the features having the largest
  Inter Quartile Range (IQR).
}
\examples{
dat <- matrix(rnorm(10000),ncol=10,nrow=1000)
rownames(dat) <- 1:1000
selectFeatures_IQR(data=dat, nb=500)
}
\author{
  Pierre Gestraud
}

