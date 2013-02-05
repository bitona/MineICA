\name{plotAllMix}
\alias{plotAllMix}
\title{Plots the Gaussian fitted  by \code{Mclust} on several numeric vectors}
\usage{
  plotAllMix(mc, A, nbMix = NULL, pdf, nbBreaks = 20,
    xlim = NULL)
}
\arguments{
  \item{mc}{A list consisting of outputs of function
  \code{Mclust} applied to each column of \code{A}, if this
  argument is missing \code{Mclust} is applied by the
  function.}

  \item{A}{A data.frame of dimensions 'samples x
  components'.}

  \item{nbMix}{The number of Gaussian to be fitted.}

  \item{nbBreaks}{The number of breaks for the histogram.}

  \item{xlim}{x-axis limits to be used in the plot.}

  \item{pdf}{A pdf file.}
}
\value{
  A list of \code{Mclust} results.
}
\description{
  Given a result of function \code{Mclust} applied on
  several numeric vectors, this function plots the fitted
  Gaussian on their histograms.
}
\details{
  This function can only deal with at the most three
  Gaussian
}
\examples{
A <-matrix(c(c(rnorm(80,mean=-0.5,sd=1),rnorm(80,mean=1,sd=0.2)),rnorm(160,mean=0.5,sd=1),
           c(rnorm(80,mean=-1,sd=0.3),rnorm(80,mean=0,sd=0.2))),ncol=3)
## apply function Mclust to each column of A
mc <- apply(A,2,Mclust)
## plot the corresponding Gaussians on the histogram of each column
plotAllMix(mc=mc,A=A)
## apply function Mclust to each column of A, and impose the fit of two Gaussian (G=2)
mc <- apply(A,2,Mclust,G=2)
## plot the corresponding Gaussians on the histogram of each column
plotAllMix(mc=mc,A=A)
## When arg 'mc' is missing, Mclust is applied by the function
plotAllMix(A=A)
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{plotMix}}, \code{\link{hist}},
  \code{\link[mclust]{Mclust}}
}

