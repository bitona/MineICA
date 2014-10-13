\name{runICA}
\alias{runICA}
\title{Run of fastICA and JADE algorithms}
\usage{
  runICA(method = c("fastICA", "JADE"), X, nbComp,
    alg.type = c("deflation", "parallel"),
    fun = c("logcosh", "exp"), maxit = 500, tol = 10^-6,
    ...)
}
\arguments{
  \item{method}{The ICA method to use, either "JADE" (the
  default) or "fastICA".}

  \item{X}{A data matrix with n rows representing
  observations (e.g genes) and p columns representing
  variables (e.g samples).}

  \item{nbComp}{The number of components to be extracted.}

  \item{alg.type}{If \code{alg.type="parallel"} the
  components are extracted simultaneously (the default), if
  \code{alg.type="deflation"} the components are extracted
  one at a time, see \code{\link[fastICA]{fastICA}}.}

  \item{fun}{The functional form of the G function used in
  the approximation to neg-entropy (see 'details' of the
  help of function \code{\link[fastICA]{fastICA}}).}

  \item{maxit}{The maximum number of iterations to
  perform.}

  \item{tol}{A positive scalar giving the tolerance at
  which the un-mixing matrix is considered to have
  converged.}

  \item{...}{Additional parameters for \code{fastICA} and
  \code{JADE}}
}
\value{
  A list, see outputs of \code{\link[fastICA]{fastICA}} and
  \code{\link[JADE]{JADE}}. This list includes at least
  three elements: \describe{\item{A}{the estimated mixing
  matrix} \item{S}{the estimated source matrix},
  item{W}{the estimated unmixing matrix}}
}
\description{
  This function performs ICA decomposition of a matrix
  using functions \code{\link[fastICA]{fastICA}} and
  \code{\link[JADE]{JADE}}.
}
\details{
  See details of the functions
  \code{\link[fastICA]{fastICA}} and
  \code{\link[JADE]{JADE}}.
}
\examples{
set.seed(2004);
M <- matrix(rnorm(5000*6,sd=0.3),ncol=10)
M[1:10,1:3] <- M[1:10,1:3] + 2
M[1:100,1:3] <- M[1:100,1:3] +1
resJade <- runICA(X=M, nbComp=2, method = "JADE", maxit=10000)
}
\author{
  Anne Biton
}

