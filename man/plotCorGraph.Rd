\name{plotCorGraph}
\alias{plotCorGraph}
\title{Plots graph using}
\usage{
  plotCorGraph(dataGraph, edgeWeight = "cor", nodeAttrs,
    nodeShape, nodeCol = "labAn", nodeName = "indComp",
    col, shape, title = "", reciproCol = "reciprocal",
    tkplot = FALSE, ...)
}
\arguments{
  \item{dataGraph}{A data.frame containing the graph
  description. It must have two columns \code{n1} and
  \code{n2}, each row denoting that there is an edge from
  n1 to n2.  Node labels in columns \code{n1} and \code{n2}
  of \code{dataGraph} must correspond to node IDs in column
  \code{id} of \code{nodeAttrs}.}

  \item{edgeWeight}{The column of dataGraph used to weight
  edges.}

  \item{nodeAttrs}{A data.frame with node description, see
  function \code{nodeAttrs}.}

  \item{nodeShape}{Denotes the column of \code{nodeAttrs}
  used to attribute the node shapes.}

  \item{nodeCol}{Denotes the column of \code{nodeAttrs}
  used to color the nodes in the graph.}

  \item{nodeName}{Denotes the column of \code{nodeAttrs}
  used as labels for the nodes in the graph.}

  \item{col}{A vector of colors, for the nodes, indexed by
  the unique elements of \code{nodeCol} column from
  \code{nodeAttrs}. If missing, colors will be
  automatically attributed.}

  \item{shape}{A vector of shapes indexed by the unique
  elements of column \code{nodeShape} from
  \code{nodeAttrs}. If missing, shapes will be
  automatically attributed.}

  \item{title}{Title for the plot}

  \item{reciproCol}{Denotes the column of \code{dataGraph}
  containing \code{TRUE} if the row defines a reciprocal
  node, else \code{FALSE}. See
  \code{\link{annotReciprocal}}.}

  \item{tkplot}{If TRUE, performs interactive plot with
  function \code{tkplot}, else uses \code{plot.igraph}.}

  \item{\dots}{Additional parameters as required by
  \code{tkplot}.}
}
\value{
  A list consisting of \describe{ \item{dataGraph}{a
  data.frame defining the correlation graph}
  \item{nodeAttrs}{a data.frame describing the node of the
  graph} \item{graph}{the graph as an object of class
  \code{igraph}} \item{graphid}{the id of the graph plotted
  using \code{tkplot}} }
}
\description{
  This function plots the correlation graph in an
  interactive device using function \code{tkplot}.
}
\details{
  You have to slighly move the nodes to see cliques because
  strongly related nodes are often superimposed. The
  \code{edgeWeight} column is used to weight the edges
  within the fruchterman.reingold layout available in the
  package \code{igraph}.

  The argument \code{nodeCol} typically denotes the column
  containing the names of the datasets. Colors are
  automatically attributed to the nodes using palette Set3
  of package \code{RColorBrewer}. The corresponding colors
  can be directly specified in the 'col' argument. In that
  case, 'col' must be a vector of colors indexed by the
  unique elements contained in \code{nodeCol} column (e.g
  dataset ids).

  As for colors, one can define the column of
  \code{nodeAttrs} that is used to define the node shapes.
  The corresponding shapes can be directly specified in the
  \code{shape} argument. In that case, \code{shape} must be
  one of \code{c("circle","square", " vcsquare",
  "rectangle", "crectangle", "vrectangle")} and must be
  indexed by the unique elements of \code{nodeShape}
  column.

  Unfortunately, shapes can't be taken into account when
  tkplot is TRUE (interactive plot).

  If \code{reciproCol} is not missing, it is used to color
  the edges, either in grey if the edge is not reciprocal
  or in black if the edge is reciprocal.
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
icaSets <- list(icaSettoy1, icaSettoy2)

resCompareAn <- compareAn(icaSets=list(icaSettoy1,icaSettoy2), labAn=c("toy1","toy2"),
                         type.corr="pearson", level="genes", cutoff_zval=0)

## Build a graph where edges correspond to maximal correlation value (useVal="cor"),
dataGraph <- compareAn2graphfile(listPairCor=resCompareAn, useMax=TRUE, useVal="cor", file="myGraph.txt")

## construction of the data.frame with the node description
nbComp <- rep(3,2) #each IcaSet contains 3 components
nbAn <- 2 # we are comparing 2 IcaSets
# labels of components created as comp*i*
labComp <- foreach(icaSet=icaSets, nb=nbComp, an=1:nbAn) \%do\% {
                  paste(rep("comp",sum(nb)),1:nbComp(icaSet),sep = "")}

# creation of the data.frame with the node description
nodeDescr <- nodeAttrs(nbAn = nbAn, nbComp = nbComp, labComp = labComp,
                       labAn = c("toy1","toy2"), file = "nodeInfo.txt")

## Plot correlation graph, slightly move the attached nodes to make the cliques visible
## use tkplot=TRUE to have an interactive graph
res <- plotCorGraph(title = "Compare toy 1 and 2", dataGraph = dataGraph, nodeName = "indComp", tkplot = FALSE,
                 nodeAttrs = nodeDescr, edgeWeight = "cor", nodeShape = "labAn", reciproCol = "reciprocal")


\dontrun{
## load two microarray datasets
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

icaSets <- list(icaSetMainz, icaSetVdx)
labAn <- c("Mainz", "Vdx")

## correlations between gene projections of each pair of IcaSet
resCompareAn <- compareAn(icaSets = icaSets, level = "genes", type.corr= "pearson",
                          labAn = labAn, cutoff_zval=0)

## construction of the correlation graph using previous output
dataGraph <- compareAn2graphfile(listPairCor=resCompareAn, useMax=TRUE, file="corGraph.txt")

## construction of the data.frame with the node description
nbComp <- rep(10,2) #each IcaSet contains 10 components
nbAn <- 2 # we are comparing 2 IcaSets
# labels of components created as comp*i*
labComp <- foreach(icaSet=icaSets, nb=nbComp, an=1:nbAn) \%do\% {
                  paste(rep("comp",sum(nb)),1:nbComp(icaSet),sep = "")}

# creation of the data.frame with the node description
nodeDescr <- nodeAttrs(nbAn = nbAn, nbComp = nbComp, labComp = labComp,
    labAn = labAn, file = "nodeInfo.txt")

## Plot correlation graph, slightly move the attached nodes to make the cliques visible
res <- plotCorGraph(title = "Compare two ICA decomsitions obtained on \\n two
                 microarray-based data of breast tumors", dataGraph = dataGraph, nodeName = "indComp",
                 nodeAttrs = nodeDescr, edgeWeight = "cor", nodeShape = "labAn", reciproCol = "reciprocal")

}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{compareAn}}, \code{\link{nodeAttrs}},
  \code{\link{compareAn2graphfile}},
  \code{\link{runCompareIcaSets}}
}

