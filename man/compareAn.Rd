\name{compareAn}
\alias{compareAn}
\title{Comparison of IcaSet objects using correlation}
\usage{
  compareAn(icaSets, labAn,
    type.corr = c("pearson", "spearman"), cutoff_zval = 0,
    level = c("samples", "features", "genes"))
}
\arguments{
  \item{icaSets}{list of IcaSet objects, e.g results of ICA
  decompositions obtained on several datasets.}

  \item{labAn}{vector of names for each icaSet, e.g the the
  names of the datasets on which were calculated the
  decompositions.}

  \item{type.corr}{Type of correlation to compute, either
  \code{'pearson'} or \code{'spearman'}.}

  \item{cutoff_zval}{either NULL or 0 (default) if all
  genes are used to compute the correlation between the
  components, or a threshold to compute the correlation on
  the genes that have at least a scaled projection higher
  than cutoff_zval. Will be used only when correlations are
  calculated on S or SByGene.}

  \item{level}{Data level of the \code{IcaSet} objects on
  which is applied the correlation. It must correspond to a
  feature shared by the IcaSet objects: \code{'samples'} if
  they were applied to common samples (correlations are
  computed between matrix \code{A}), \code{'features'} if
  they were applied to common features (correlations are
  computed between matrix \code{S}), \code{'genes'} if they
  share gene IDs after annotation into genes (correlations
  are computed between matrix \code{SByGene}).}
}
\value{
  A list whose length equals the number of pairs of
  \code{IcaSet} and whose elements are outputs of function
  \code{\link{cor2An}}.
}
\description{
  Compare \code{\link{IcaSet}} objects by computing the
  correlation between either projection values of common
  features or genes, or contributions of common samples.
}
\details{
  The user must carefully choose the object on which the
  correlation will be computed. If \code{level='samples'},
  the correlations are based on the mixing matrices of the
  ICA decompositions (of dimension samples x components).
  \code{'A'} will be typically chosen when the ICA
  decompositions were computed on the same dataset, or on
  datasets that include the same samples. If
  \code{level='features'} is chosen, the correlation is
  calculated between the source matrices (of dimension
  features x components) of the ICA decompositions.
  \code{'S'} will be typically used when the ICA
  decompositions share common features (e.g same
  microarrays). If \code{level='genes'}, the correlations
  are calculated on the attributes \code{'SByGene'} which
  store the projections of the annotated features.
  \code{'SByGene'} will be typically chosen when ICA were
  computed on datasets from different technologies, for
  which comparison is possible only after annotation into a
  common ID, like genes.

  \code{cutoff_zval} is only used when \code{level} is one
  of \code{c('genes','features')}, in order to restrict the
  correlation to the contributing features or genes.

  When \code{cutoff_zval} is specified, for each pair of
  components, genes or features that are included in the
  circle of center 0 and radius \code{cutoff_zval} are
  excluded from the computation of the correlation.

  It must be taken into account by the user that if
  \code{cutoff_zval} is different from \code{NULL} or
  \code{0}, the computation will be much slowler since each
  pair of component is treated individually.
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

listPairCor <- compareAn(icaSets=list(icaSettoy1,icaSettoy2), labAn=c("toy1","toy2"),
                         type.corr="pearson", level="genes", cutoff_zval=0)


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

## The pearson correlation is used as a measure of association between the gene projections
# on the different components (type.corr="pearson").
listPairCor <- compareAn(icaSets=list(icaSetMainz,icaSetVdx),
labAn=c("Mainz","Vdx"), type.corr="pearson", level="genes", cutoff_zval=0)

## Same thing but adding a selection of genes on which the correlation between two components is computed:
# when considering pairs of components, only projections whose scaled values are not located within
# the circle of radius 1 are used to compute the correlation (cutoff_zval=1).
listPairCor <-  compareAn(icaSets=list(icaSetMainz,icaSetVdx),
labAn=c("Mainz","Vdx"), type.corr="pearson", cutoff_zval=1, level="genes")
}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{cor2An}}
}

