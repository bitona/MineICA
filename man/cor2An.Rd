\name{cor2An}
\alias{cor2An}
\title{Correlation between two matrices}
\usage{
  cor2An(mat1, mat2, lab,
    type.corr = c("pearson", "spearman"), cutoff_zval = 0)
}
\arguments{
  \item{mat1}{matrix of dimension features/genes x number
  of components, e.g the results of an ICA decomposition}

  \item{mat2}{matrix of dimension features/genes x number
  of components, e.g the results of an ICA decomposition}

  \item{lab}{The vector of labels for mat1 and mat2, e.g
  the the names of the two datasets on which were
  calculated the two decompositions}

  \item{type.corr}{Type of correlation, either
  \code{'pearson'} or \code{'spearman'}}

  \item{cutoff_zval}{cutoff_zval: 0 (default) if all genes
  are used to compute the correlation between the
  components, or a threshold to compute the correlation on
  the genes that have at least a scaled projection higher
  than cutoff_zval.}
}
\value{
  This function returns a list consisting of: \item{cor}{a
  matrix of dimensions '(nbcomp1+nbcomp2) x
  (nbcomp1*nbcomp2)', containing the correlation values
  between each pair of components,} \item{pval}{ a matrix
  of dimension '(nbcomp1+nbcomp2) x (nbcomp1*nbcomp2)',
  containing the p-value of the correlation tests for each
  pair of components,} \item{inter}{ the intersection
  between the features/genes of \code{mat1} and
  \code{mat2},} \item{labAn}{ the labels of the compared
  matrices.}
}
\description{
  This function measures the correlation between two
  matrices containing the results of two decompositions.
}
\details{
  Before computing the correlations, the components are
  scaled and restricted to common row names.

  It must be taken into account by the user that if
  \code{cutoff_zval} is different from NULL or zero, the
  computation will be slowler since each pair of component
  is treated individually.

  When \code{cutoff_zval} is specified, for each pair of
  components, genes that are included in the circle of
  center 0 and radius \code{cutoff_zval} are excluded from
  the computation of the correlation between the gene
  projection of the two components.
}
\examples{
cor2An(mat1=matrix(rnorm(10000),nrow=1000,ncol=10), mat2=matrix(rnorm(10000),nrow=1000,ncol=10),
       lab=c("An1","An2"), type.corr="pearson")
}
\author{
  Anne Biton
}
\seealso{
  \code{rcorr}, \code{cor.test}, \code{\link{compareAn}}
}

