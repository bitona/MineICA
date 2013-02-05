\name{runCompareIcaSets}
\alias{runCompareIcaSets}
\title{runCompareIcaSets}
\usage{
  runCompareIcaSets(icaSets, labAn,
    type.corr = c("pearson", "spearman"), cutoff_zval = 0,
    level = c("genes", "features", "samples"),
    fileNodeDescr = NULL, fileDataGraph = NULL,
    plot = TRUE, title = "", col, cutoff_graph = NULL,
    useMax = TRUE, tkplot = FALSE)
}
\arguments{
  \item{icaSets}{List of \code{\link{IcaSet}} objects, e.g
  results of ICA decompositions obtained on several
  datasets.}

  \item{labAn}{Vector of names for each icaSet, e.g the the
  names of the datasets on which were calculated the
  decompositions.}

  \item{type.corr}{Type of correlation to compute, either
  \code{'pearson'} or \code{'spearman'}.}

  \item{cutoff_zval}{Either NULL or 0 (default) if all
  genes are used to compute the correlation between the
  components, or a threshold to compute the correlation
  using the genes that have at least a scaled projection
  higher than cutoff_zval. Will be used only when
  \code{level} is one of \code{c("features","genes")}.}

  \item{level}{Data level of the \code{IcaSet} objects on
  which is applied the correlation. It must correspond to a
  data level shared by the IcaSet objects: \code{'samples'}
  if they were applied to common samples (correlations are
  computed between matrix \code{A}), \code{'features'} if
  they were applied to common features (correlations are
  computed between matrix \code{S}), \code{'genes'} if they
  share gene IDs after annotation into genes (correlations
  are computed between matrix \code{SByGene}).}

  \item{fileNodeDescr}{File where node descriptions are
  saved (useful when the user wants to visualize the graph
  using Cytoscape).}

  \item{fileDataGraph}{File where graph description is
  saved (useful when the user wants to visualize the graph
  using Cytoscape).}

  \item{plot}{if \code{TRUE} (default) plot the correlation
  graph}

  \item{title}{title of the graph}

  \item{col}{vector of colors indexed by elements of labAn;
  if missing, colors will be automatically attributed}

  \item{cutoff_graph}{the cutoff used to select pairs that
  will be included in the graph}

  \item{useMax}{if \code{TRUE}, the graph is restricted to
  edges that correspond to maximum correlation between
  components, see details}

  \item{tkplot}{If TRUE, performs interactive plot with
  function \code{tkplot}, else uses \code{plot.igraph}}
}
\value{
  A list consisting of \describe{ \item{dataGraph:}{a
  data.frame defining the correlation graph}
  \item{nodeAttrs:}{a data.frame describing the node of the
  graph,} \item{graph:}{the graph as an object of class
  \code{igraph},} \item{graphid}{the id of the graph
  plotted with \code{tkplot}}. }
}
\description{
  This function encompasses the comparison of several
  IcaSet objects using correlations and the plot of the
  corresponding correlation graph. The IcaSet objects are
  compared by calculating the correlation between either
  projection values of common features or genes, or
  contributions of common samples.
}
\details{
  This function calls four functions:
  \code{\link{compareAn}} which computes the correlations,
  \code{\link{compareAn2graphfile}} which builds the graph,
  \code{\link{nodeAttrs}} which builds the node description
  data, and \code{\link{plotCorGraph}} which uses tkplot to
  plot the graph in an interactive device.

  If the user wants to see the correlation graph in
  Cytoscape, he must fill the arguments
  \code{fileDataGraph} and \code{fileNodeDescr}, in order
  to import the graph and its node descriptions as a .txt
  file in Cytoscape.

  When \code{labAn} is missing, each element i of
  \code{icaSets} is labeled as 'Ani'.

  The user must carefully choose the data level used in the
  comparison: If \code{level='samples'}, the correlations
  are based on the mixing matrices of the ICA
  decompositions (of dimension samples x components).
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
  of \code{c('features','genes')}, in order to restrict the
  correlation to the contributing features or genes.

  When \code{cutoff_zval} is specified, for each pair of
  components, genes or features that are included in the
  circle of center 0 and radius \code{cutoff_zval} are
  excluded from the computation of the correlation.

  It must be taken into account by the user that if
  cutoff_zval is different from NULL or zero, the
  computation will be much slowler since each pair of
  component is treated individually.

  Edges of the graph are built based on the correlation
  values between the components. Absolute values of
  correlations are used since components have no direction.

  If \code{useMax} is \code{TRUE} each component will be
  linked to only one component of each other IcaSet that
  corresponds to the most correlated component among all
  components of the same IcaSet. If \code{cutoff_graph} is
  specified, only correlations exceeding this value are
  taken into account to build the graph. For example, if
  \code{cutoff} is 1, only relationships between components
  that correspond to a correlation value higher than 1 will
  be included. Absolute correlation values are used since
  the components have no direction.

  The contents of the returned list are \describe{
  \item{dataGraph:}{\code{dataGraph} data.frame that
  describes the correlation graph,}
  \item{nodeAttrs:}{\code{nodeAttrs} data.frame that
  describes the node of the graph}
  \item{graph}{\code{graph} the graph as an igraph-object,}
  \item{graphid:}{\code{graphid} the id of the graph
  plotted using tkplot.} }
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

## build IcaSet objects
icaSettoy1 <- buildIcaSet(params=params, A=data.frame(resJade1$A), S=data.frame(resJade1$S),
                          dat=dat1, alreadyAnnot=TRUE)$icaSet
icaSettoy2 <- buildIcaSet(params=params, A=data.frame(resJade2$A), S=data.frame(resJade2$S),
                          dat=dat2, alreadyAnnot=TRUE)$icaSet

## compare IcaSet objects
rescomp <- runCompareIcaSets(icaSets=list(icaSettoy1, icaSettoy2), labAn=c("toy1","toy2"),
                             type.corr="pearson", level="genes")


\dontrun{
## load the microarray-based gene expression datasets
## of breast tumors
library(breastCancerMAINZ)
library(breastCancerVDX)
data(mainz)
data(vdx)

## Define a function used to build two examples of IcaSet objects
## and annotate the probe sets into gene Symbols
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

## compare the IcaSets
runCompareIcaSets(icaSets=list(icaSetMainz, icaSetVdx), labAn=c("Mainz","Vdx"), type.corr="pearson", level="genes")
}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{compareAn2graphfile}},
  \code{\link{compareAn}}, \code{\link{cor2An}},
  \code{\link{plotCorGraph}}
}

