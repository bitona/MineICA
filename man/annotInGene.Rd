\name{annotInGene}
\alias{annotInGene}
\title{Features annotation of an object of class IcaSet.}
\usage{
  annotInGene(icaSet, params, annot = TRUE)
}
\arguments{
  \item{icaSet}{An object of class \code{\link{IcaSet}} to
  be annotated, must contain a valid \code{annotation}
  attribute.}

  \item{params}{An object of class
  \code{\link{MineICAParams}} containing the parameters of
  the analysis.}

  \item{annot}{TRUE (default) if the IcaSet object must
  indeed be annotated}
}
\value{
  The modified argument \code{icaSet}, with filled
  attributes \code{SByGene} and \code{datByGene}.
}
\description{
  This function annotates the features of an
  \code{\link{IcaSet}} object and fills its attributes
  \code{SByGene} and \code{datByGene}.
}
\details{
  When attribute \code{annotation} of \code{icaSet} is not
  specified (of length \code{0}), \code{biomaRt} is used to
  annotate the features through function
  \code{\link{annotFeaturesWithBiomaRt}}.

  When specified, attribute \code{annotation} of argument
  \code{icaSet} must be an annotation package and will be
  used to annotate the \code{featureNames} of
  \code{icaSet}. In addition, the attribute \code{typeID}
  (a vector) of argument \code{icaSet} must contain a valid
  element \code{geneID_annotation} that determines the
  object of the package to be used for the annotation, see
  \code{\link{IcaSet}}.

  When argument \code{annot} is TRUE, this function fills
  the attributes \code{SByGene} and \code{datByGene} of
  \code{icaSet}. When several feature IDs are available for
  a same gene ID, the median value of the corresponding
  features IDs is attributed to the gene (the median of the
  projection values is used for attribute \code{SByGene},
  and the median of the expression values is used for
  attribute \code{datByGene}).

  When attribute \code{chipManu} of the argument
  \code{icaSet} is "illumina", the features are first
  converted into nuID using the package 'lumi*Mapping' and
  then annotated into genes. In that case, features can
  only be annotated in ENTREZID or SYMBOL. It means that
  \code{typeID(icaSet)['geneID_annotation']} must be either
  'ENTREZID' or 'SYMBOL'. You will need to annotate
  yourself the IcaSet object if you want to use different
  IDs.
}
\examples{
#load data
data(icaSetCarbayo)
require(hgu133a.db)

# run annotation of the features into gene Symbols as specified in 'typeID(icaSetCarbayo)["geneID_annotation"]',
# using package hgu133a.db as defined in 'annotation(icaSetMainz)'
icaSetCarbayo <- annotInGene(icaSet=icaSetCarbayo, params=buildMineICAParams())

\dontrun{
#load data
library(breastCancerMAINZ)
data(mainz)
#run ICA
resJade <- runICA(X=exprs(mainz), nbComp=5, method = "JADE", maxit=10000)

#build params
params <- buildMineICAParams(resPath="mainz/")

#build a new IcaSet object, omitting annotation of the features (runAnnot=FALSE)
#but specifying the element "geneID_annotation" of argument 'typeID'
icaSetMainz <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S),
                             dat=exprs(mainz), pData=pData(mainz),
                             annotation="hgu133a.db", typeID= c(geneID_annotation = "SYMBOL",
                             geneID_biomart = "hgnc_symbol", featureID_biomart = "affy_hg_u133a"),
                             chipManu = "affymetrix", runAnnot=FALSE,
                             mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl"))

#Attributes SByGene is empty and attribute datByGene refers to assayData
SByGene(icaSetMainz)
head(datByGene(icaSetMainz))

# run annotation of the features into gene Symbols as specified in 'typeID(icaSetMainz)["geneID_annotation"]',
# using package hgu133a.db as defined in 'annotation(icaSetMainz)'
icaSetMainz <- annotInGene(icaSet=icaSetMainz, params=params)
}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{annotFeaturesComp}}
}

