\name{buildIcaSet}
\alias{buildIcaSet}
\title{This function builds an object of class \code{\link{IcaSet}}.}
\usage{
  buildIcaSet(params, A, S, dat, pData = new("data.frame"),
    fData = new("data.frame"), witGenes = new("character"),
    compNames = new("character"),
    refSamples = new("character"),
    annotation = new("character"),
    chipManu = new("character"),
    chipVersion = new("character"), alreadyAnnot = FALSE,
    typeID = c(geneID_annotation = "SYMBOL", geneID_biomart = "hgnc_symbol", featureID_biomart = ""),
    runAnnot = TRUE, organism = "Human",
    mart = new("Mart"))
}
\arguments{
  \item{params}{An object of class
  \code{\link{MineICAParams}} containing the parameters of
  the analysis}

  \item{A}{The mixing matrix of the ICA decomposition (of
  dimension samples x components).}

  \item{S}{The source matrix of the ICA decomposition (of
  dimension features x components).}

  \item{dat}{The data matrix the ICA was applied to (of
  dimension features x samples).}

  \item{pData}{Phenotype data, a data.frame which contains
  the sample informations of dimension samples x
  annotations.}

  \item{fData}{Feature data, a data.frame which contrains
  the feature descriptions of dimensions features x
  annotations.}

  \item{witGenes}{A vector of witness genes. They are
  representative of the expression behavior of the
  contributing genes of each component. If missing or NULL,
  they will be automatically attributed using function
  \code{\link{selectWitnessGenes}}.}

  \item{compNames}{A vector of component labels.}

  \item{refSamples}{A vector of reference sample IDs (e.g
  the "normal" samples).}

  \item{annotation}{An annotation package (e.g a ".db"
  package specific to the microarray used to generate
  \code{dat})}

  \item{chipManu}{If microarray data, the manufacturer:
  either 'affymetrix' or 'illumina'.}

  \item{chipVersion}{For illumina microarrays: the version
  of the microarray.}

  \item{alreadyAnnot}{TRUE if the feature IDs contained in
  the row names of \code{dat} and \code{S} already
  correspond to the final level of annotation (e.g if they
  are already gene IDs). In that case, no annotation is
  performed.}

  \item{typeID}{A character vector specifying the
  annotation IDs, it includes three elements : \describe{
  \item{geneID_annotation}{the IDs from the package to be
  used to annotate the features into genes. It will be used
  to fill the attributes \code{datByGene} and
  \code{SByGene} of the \code{icaSet}. It must match one of
  the objects the corresponding package supports (you can
  access the list of objects by typing
  ls("package:packagename")). If no annotation package is
  provided, this element is not useful.}
  \item{geneID_biomart}{the type of gene IDs, as available
  in \code{listFilters(mart)}; where mart is specified as
  described in \code{\link[biomaRt]{useMart}}. If you have
  directly built the IcaSet at the gene level (i.e if no
  annotation package is used), \code{featureID_biomart} and
  \code{geneID_biomart} will be identical.}
  \item{featureID_biomart}{the type of feature IDs, as
  available in \code{listFilters(mart)}; where \code{mart}
  is specified as described in function
  \code{\link[biomaRt]{useMart}}. Not useful if you work at
  the gene level.} }}

  \item{runAnnot}{If TRUE, \code{icaSet} is annotated with
  function \code{annotInGene}.}

  \item{organism}{The organism the data correspond to.}

  \item{mart}{The mart object (database and dataset) used
  for annotation, see function \code{useMart} of package
  \code{biomaRt}}
}
\value{
  An object of class IcaSet
}
\description{
  This function builds an object of class
  \code{\link{IcaSet}}.
}
\examples{
dat <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
rownames(dat) <- paste("g", 1:1000, sep="")
colnames(dat) <- paste("s", 1:10, sep="")

## build a data.frame containing sample annotations
annot <- data.frame(type=c(rep("a",5),rep("b",5)))
rownames(annot) <- colnames(dat)

## run ICA
resJade <- runICA(X=dat, nbComp=3, method = "JADE")

## build params
params <- buildMineICAParams(resPath="toy/")

## build IcaSet object
icaSettoy <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S),
                         dat=dat, pData=annot, alreadyAnnot=TRUE)
params <- icaSettoy$params
icaSettoy <- icaSettoy$icaSet

\dontrun{
## load data
library(breastCancerMAINZ)
data(mainz)

## run ICA
resJade <- runICA(X=dataMainz, nbComp=10, method = "JADE", maxit=10000)

## build params
params <- buildMineICAParams(resPath="mainz/")

## build IcaSet object

# fill typeID, Mainz data originate from affymetrix HG-U133a  microarray and are indexed by probe sets
# we want to annotate the probe sets into Gene Symbols
typeIDmainz <-  c(geneID_annotation="SYMBOL", geneID_biomart="hgnc_symbol", featureID_biomart="affy_hg_u133a")

icaSetMainz <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S),
                             dat=exprs(mainz), pData=pData(mainz),
                             annotation="hgu133a.db", typeID= c(geneID_annotation = "SYMBOL",
                             geneID_biomart = "hgnc_symbol", featureID_biomart = "affy_hg_u133a"),
                             chipManu = "affymetrix", runAnnot=TRUE,
                             mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl"))
}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{selectWitnessGenes}},
  \code{\link{annotInGene}}
}

