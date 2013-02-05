\name{compareAn2graphfile}
\alias{compareAn2graphfile}
\title{compareAn2graphfile}
\usage{
  compareAn2graphfile(listPairCor, useMax = TRUE,
    cutoff = NULL, useVal = c("cor", "pval"), file = NULL)
}
\arguments{
  \item{listPairCor}{The output of the function
  \code{\link{compareAn}}, containing the correlation
  between several pairs of objects of class
  \code{\link{IcaSet}}.}

  \item{useMax}{If TRUE, the graph is restricted to edges
  that correspond to maximum score, see details}

  \item{cutoff}{Cutoff used to select pairs that will be
  included in the graph.}

  \item{useVal}{The value on which is based the graph,
  either \code{"cor"} for correlation or \code{"pval"} for
  p-values of correlation tests.}

  \item{file}{File name.}
}
\value{
  A data.frame with the graph description, has two columns
  \code{n1} and \code{n2} filled with node IDs, each row
  denotes that there is an edge from \code{n1} to
  \code{n2}. Additional columns quantify the strength of
  association: correlation (\code{cor}), p-value
  (\code{pval}), (\code{1-abs(cor)}) (\code{distcor}),
  log10-pvalue (\code{logpval}).
}
\description{
  This function builds a correlation graph from the outputs
  of function \code{\link{compareAn}}.
}
\details{
  When correlations are considered (\code{useVal}="cor"),
  absolute values are used since the components have no
  direction.

  If \code{useMax} is \code{TRUE} each component is linked
  to the most correlated component of each different
  \code{IcaSet}.

  If \code{cutoff} is specified, only correlations
  exceeding this value are taken into account during the
  graph construction. For example, if \code{cutoff} is 1,
  only relationships between components that correspond to
  a correlation value larger than 1 will be included.

  When \code{useVal="pval"} and \code{useMax=TRUE}, the
  minimum value is taken instead of the maximum.
}
\examples{
dat1 <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
rownames(dat1) <- paste("g", 1:1000, sep="")
colnames(dat1) <- paste("s", 1:10, sep="")
dat2 <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
rownames(dat2) <- paste("g", 1:1000, sep="")
colnames(dat2) <- paste("s", 1:10, sep="")

## run ICA
resJade1 <- runICA(X=dat1, nbComp=3, method = "JADE")
resJade2 <- runICA(X=dat2, nbComp=3, method = "JADE")

## build params
params <- buildMineICAParams(resPath="toy/")

## build IcaSet object
icaSettoy1 <- buildIcaSet(params=params, A=data.frame(resJade1$A), S=data.frame(resJade1$S),
                          dat=dat1, alreadyAnnot=TRUE)$icaSet
icaSettoy2 <- buildIcaSet(params=params, A=data.frame(resJade2$A), S=data.frame(resJade2$S),
                          dat=dat2, alreadyAnnot=TRUE)$icaSet

resCompareAn <- compareAn(icaSets=list(icaSettoy1,icaSettoy2), labAn=c("toy1","toy2"),
                         type.corr="pearson", level="genes", cutoff_zval=0)

## Build a graph where edges correspond to maximal correlation value (useVal="cor"),
compareAn2graphfile(listPairCor=resCompareAn, useMax=TRUE, useVal="cor", file="myGraph.txt")


\dontrun{
#### Comparison of 2 ICA decompositions obtained on 2 different gene expression datasets.
## load the two datasets
library(breastCancerMAINZ)
library(breastCancerVDX)
data(mainz)
data(vdx)

## Define a function used to build two examples of IcaSet objects
treat <- function(es, annot="hgu133a.db") {
   es <- selectFeatures_IQR(es,10000)
   exprs(es) <- t(apply(exprs(es),1,scale,scale=FALSE))
   colnames(exprs(es)) <- sampleNames(es)
   resJade <- runICA(X=exprs(es), nbComp=10, method = "JADE", maxit=10000)
   resBuild <- buildIcaSet(params=buildMineICAParams(), A=data.frame(resJade$A), S=data.frame(resJade$S),
                        dat=exprs(es), pData=pData(es), refSamples=character(0),
                        annotation=annot, typeID= typeIDmainz,
                        chipManu = "affymetrix", mart=mart)
   icaSet <- resBuild$icaSet
}
## Build the two IcaSet objects
icaSetMainz <- treat(mainz)
icaSetVdx <- treat(vdx)

## Compute correlation between every pair of IcaSet objects.
resCompareAn <- compareAn(icaSets=list(icaSetMainz,icaSetVdx),
labAn=c("Mainz","Vdx"), type.corr="pearson", level="genes", cutoff_zval=0)

## Same thing but adding a selection of genes on which the correlation between two components is computed:
# when considering pairs of components, only projections whose scaled values are not located within
# the circle of radius 1 are used to compute the correlation (cutoff_zval=1).
resCompareAn <-  compareAn(icaSets=list(icaSetMainz,icaSetVdx),
labAn=c("Mainz","Vdx"), type.corr="pearson", cutoff_zval=1, level="genes")

## Build a graph where edges correspond to maximal correlation value (useVal="cor"),
## i.e, component A of analysis i is linked to component B of analysis j,
## only if component B is the most correlated component to A amongst all component of analysis j.
compareAn2graphfile(listPairCor=resCompareAn, useMax=TRUE, useVal="cor", file="myGraph.txt")

## Restrict the graph to correlation values exceeding 0.4
compareAn2graphfile(listPairCor=resCompareAn, useMax=FALSE, cutoff=0.4,  useVal="cor", file="myGraph.txt")

}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{compareAn}}, \code{\link{cor2An}}
}

