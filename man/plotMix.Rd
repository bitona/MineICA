\name{plotMix}
\alias{plotMix}
\title{Plots an histogram and Gaussian fitted  by \code{\link[mclust]{Mclust}}}
\usage{
  plotMix(mc, data, nbBreaks, traceDensity = TRUE,
    title = "", xlim, ylim, ...)
}
\arguments{
  \item{mc}{The result of Mclust function applied to
  argument \code{data}}

  \item{data}{A vector of numeric values}

  \item{nbBreaks}{The number of breaks for the histogram}

  \item{traceDensity}{If TRUE (default) density are
  displayed on the y-axis, else if FALSE counts are
  displayed on the y-acis}

  \item{title}{A title for the plot}

  \item{xlim}{x-axis limits to be used in the plot}

  \item{ylim}{y-axis limits to be used in the plot}

  \item{...}{additional arguments for hist}
}
\value{
  NULL
}
\description{
  Given a result of function \code{Mclust} applied to a
  numeric vector, this function draws the fitted Gaussian
  on the histogram of the data values.
}
\details{
  A shapiro test p-value is added to the plot title. This
  function can only deal with at the most three Gaussian.
}
\examples{
## create a mix of two Gaussian
v <-c(rnorm(80,mean=-0.5,sd=1),rnorm(80,mean=1,sd=0.2))
## apply Mclust
mc <- Mclust(v)
## plot fitted Gaussian on histogram of v
plotMix(mc=mc,data=v,nbBreaks=30)
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{hist}}, \code{\link[mclust]{Mclust}}
}

