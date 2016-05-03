\name{correl2Comp}
\alias{correl2Comp}
\title{correl2Comp}
\description{This function computes the correlation between two components.}
\usage{correl2Comp(comp1, comp2, type.corr = "pearson", plot = FALSE, 
    cutoff_zval = 0, test = FALSE, alreadyTreat = FALSE)}
\arguments{
  \item{comp1}{The first component, a vector of projections or contributions indexed by labels}
  \item{comp2}{The second component, a vector of projections or contributions indexed by labels}
  \item{type.corr}{Type of correlation to be computed, either \code{'pearson'} or \code{'spearman'}}
  \item{plot}{if \code{TRUE}, plot \code{comp1} vs \code{comp2}}
  \item{cutoff_zval}{either NULL or 0 (default) if all genes are used to compute the correlation between the components, or a threshold to compute the correlation on the genes that have at least a scaled projection higher than cutoff_zval. }
  \item{test}{if TRUE the correlation test p-value is returned instead of the correlation value}
  \item{alreadyTreat}{if TRUE comp1 and comp2 are considered as being already treated (i.e scaled and restricted to common elements) }
}
\details{Before computing the correlation, the components are scaled and restricted to common labels.  
When \code{cutoff_zval} is different from \code{0}, the elements that are included in the circle of center 0 and radius \code{cutoff_zval} are not taken into account during the computation of the correlation.}
\value{This function returns either the correlation value or the p-value of the correlation test.}

\author{Anne Biton}





