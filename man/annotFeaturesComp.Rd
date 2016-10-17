\name{annotFeaturesComp}
\alias{annotFeaturesComp}
\title{Features annotation}
\usage{
  annotFeaturesComp(icaSet, params,
    type = toupper(typeID(icaSet)["geneID_annotation"]),
    featureId = typeID(icaSet)["featureID_biomart"],
    geneId = typeID(icaSet)["geneID_biomart"])
}
\arguments{
  \item{icaSet}{An object of class \code{\link{IcaSet}}
  whose features have to be annotated. The attribute
  \code{annotation} of this object contains the annotation
  package to be used.}

  \item{params}{An object of class
  \code{\link{MineICAParams}} containing the parameters of
  the analysis.}

  \item{type}{The ID of the object of the annotation
  package to be used for the annotation, must be available
  in \code{ls("package:package_name")}}

  \item{featureId}{The type of the feature IDs, in the
  \code{biomaRt} way (type \code{listFilters(mart)} to
  choose one). Used when \code{annotation(icaSet)} is of
  length 0.}

  \item{geneId}{The type of the gene IDs, in the
  \code{biomaRt} way (type \code{listAttributes(mart)} to
  choose one). Used when \code{annotation(icaSet)} is of
  length 0.}
}
\value{
  This function returns the argument \code{icaSet} with
  attributes \code{SByGene} and \code{datByGene} filled.
}
\description{
  ##' This function annotates the features of an object of
  class \code{\link{IcaSet}}, and fills its attributes
  \code{SByGene} and \code{datByGene}.
}
\details{
  This function is called by function
  \code{\link{annotInGene}} which will check the validity
  of the attributes \code{annotation, typeID, chipManu} and
  eventually \code{chipVersion} of \code{icaSet}. If
  available, the attribute \code{annotation} of argument
  \code{icaSet} must be an annotation package and will be
  used to annotate the \code{featureNames} of
  \code{icaSet}. If attribute \code{annotation} of argument
  \code{icaSet} is not available (of length 0),
  \code{biomaRt} is used to annotate the features.

  This function fills the attributes \code{SByGene} and
  \code{datByGene} of the argument \code{icaSet}. When
  several feature IDs are available for a same gene ID, the
  median value of the corresponding features IDs is
  attributed to the gene (the median of projection values
  is used for attribute \code{SByGene}, and the median of
  expression values is used for attribute
  \code{datByGene}).

  When attribute \code{chipManu} of the argument
  \code{icaSet} is "illumina", the features are first
  converted into nuID using the package 'lumi*Mapping' and
  then annotated into genes. In that case, features can
  only be annotated in ENTREZID or SYMBOL. It means that
  \code{typeID(icaSet)['geneID_annotation']} must be either
  'ENTREZID' or 'SYMBOL'. You will need to annotate
  yourself the \code{\link{IcaSet}} object if you want to
  use different IDs.
}
\examples{
## load an example of IcaSet
data(icaSetCarbayo)
params <- buildMineICAParams()
require(hgu133a.db)
####===================================================
## Use of annotation package contained in annotation(icaSet)
####====================================================
## annotation in SYMBOL
icaSetCarbayo_annot <- annotFeaturesComp(icaSet=icaSetCarbayo, params=params, type="SYMBOL")
# arg 'type' is optional since the function uses contents of typeID(icaSet) as the defaults,
# it is specified in these examples for pedagogy views

## annotation in Entrez Gene
icaSetCarbayo_annot <- annotFeaturesComp(icaSet=icaSetCarbayo, params=params, type="ENTREZID")

\dontrun{
####===================================================
## Use of biomaRt, when annotation(icaSet) is of length 0
####====================================================
## empty attribute 'annotation' of the IcaSet object
# when this attribute is not specified, biomaRt is used for annotation
annotation(icaSetCarbayo) <- character()

# make sure the mart attribute is correctly defined
mart(icaSetCarbayo) <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

## make sure elements "featureID_biomaRt" and "geneID_biomaRt" of typeID(icaSet) are correctly filled
# they will be used by function 'annotFeaturesComp' through biomaRt to query the database
typeID(icaSetCarbayo)

## run annotation of HG-U133A probe set IDs into Gene Symbols using biomaRt
icaSetCarbayo_annot <- annotFeaturesComp(icaSet=icaSetCarbayo, params=params)

}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{annotFeatures}},
  \code{\link{annotFeaturesWithBiomaRt}},
  \code{\link{annotInGene}}
}

