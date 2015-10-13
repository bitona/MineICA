\name{clusterFastICARuns}
\alias{clusterFastICARuns}
\title{Run of fastICA and JADE algorithms}
\usage{
  clusterFastICARuns(X, nbComp, nbIt = 100,
    alg.type = c("deflation", "parallel"),
    fun = c("logcosh", "exp"), maxit = 500, tol = 10^-6,
    funClus = c("hclust", "agnes", "pam", "kmeans"),
    row.norm = FALSE, bootstrap = FALSE, ...)
}
\arguments{
  \item{X}{A data matrix with n rows representing
  observations (e.g genes) and p columns representing
  variables (e.g samples).}

  \item{nbComp}{The number of components to be extracted.}

  \item{nbIt}{The number of iterations of FastICA}

  \item{alg.type}{If \code{alg.type="parallel"} the
  components are extracted simultaneously (the default), if
  \code{alg.type="deflation"} the components are extracted
  one at a time, see \code{\link[fastICA]{fastICA}}.}

  \item{fun}{The functional form of the G function used in
  the approximation to neg-entropy (see 'details' of the
  help of function \code{\link[fastICA]{fastICA}}).}

  \item{row.norm}{a logical value indicating whether rows
  of the data matrix \code{X} should be standardized
  beforehand (see help of function \code{fastICA})}

  \item{maxit}{The maximum number of iterations to
  perform.}

  \item{tol}{A positive scalar giving the tolerance at
  which the un-mixing matrix is considered to have
  converged.}

  \item{funClus}{The clustering function to be used to
  cluster the estimates}

  \item{bootstrap}{if TRUE the data is bootstraped before
  each fastICA iteration, else (default) only random
  initializations are done}

  \item{...}{Additional parameters for code{funClus}}
}
\value{
  A list consisting of: \describe{\item{A}{the estimated
  mixing matrix} \item{S}{the estimated source matrix},
  item{W}{the estimated unmixing matrix}, \item{Iq}{Iq
  indices.}}
}
\description{
  This function runs the fastICA algorithm several times
  with random initializations. The obtained components are
  clustered and the medoids of these clusters are used as
  the final estimates. The returned estimates are ordered
  by decreasing Iq values which measure the compactness of
  the clusters (see details).
}
\details{
  This function implements in R fastICA iterations followed
  by a clustering step, as defined in the matlab package
  'icasso'. Among the indices computed by icasso, only the
  Iq index is currently computed. As defined in 'icasso',
  the Iq index measures the difference between the
  intra-cluster similarity and the extra-cluster
  similiarity. No visualization of the clusters is yet
  available.

  If \code{bootstrap=TRUE} a bootstrap (applied to the
  observations) is used to perturb the data before each
  iteration, then function \code{fastICA} is applied with
  random initializations.

  By default, in 'icasso', agglomerative hierarchical
  clustering with average linkage is performed. To use the
  same clustering, please use \code{funClus="hclust"} and
  \code{method="average"}. But this function also allows
  you to apply the clustering of your choice among
  \code{kmeans, pam, hclust, agnes} by specifying
  \code{funClus} and adding the adequat additional
  parameters.

  See details of the functions
  \code{\link[fastICA]{fastICA}}.
}
\examples{
## generate a data
set.seed(2004);
M <- matrix(rnorm(5000*6,sd=0.3),ncol=10)
M[1:100,1:3] <- M[1:100,1:3] + 2
M[1:200,1:3] <- M[1:200,4:6] +1

## Random initializations are used for each iteration of FastICA
## Estimates are clustered using hierarchical clustering with average linkage
res <- clusterFastICARuns(X=M, nbComp=2, alg.type="deflation",
                          nbIt=3, funClus="hclust", method="average")

## Data are boostraped before each iteration and random initializations
## are used for each iteration of FastICA
## Estimates are clustered using hierarchical clustering with ward
res <- clusterFastICARuns(X=M, nbComp=2, alg.type="deflation",
                          nbIt=3, funClus="hclust", method="ward")
}
\author{
  Anne Biton
}

