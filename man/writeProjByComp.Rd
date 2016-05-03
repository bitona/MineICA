\name{writeProjByComp}
\alias{writeProjByComp}
\title{writeProjByComp}
\description{This function writes in an html file the description of the
  features, or genes, that contribute to each component. It also writes
  an html file containing, for each feature or gene, its projection value on every component.}
\usage{writeProjByComp(icaSet, params, mart = useMart(biomart = "ensembl", 
    dataset = "hsapiens_gene_ensembl"), typeRetrieved = NULL, addNbOcc =
    TRUE, selectionByComp = NULL, level = c("features", "genes"), typeId, selCutoffWrite=2.5)}
\arguments{
  \item{icaSet}{An object of class \code{\link{IcaSet}}}
  \item{params}{An object of class \code{\link{MineICAParams}} containing the
    parameters of the analysis. The files are written in the path
    \code{genesPath(params)}. \code{selCutoff(params)} is used to select the features or genes
    by component.}
  \item{mart}{An output of function \code{\link[biomaRt]{useMart}}
    containing the database used for annotation.}
  \item{typeRetrieved}{The annotations biomaRt is queried about. They
    describe the feature or gene IDs of the argument \code{icaSet}, see \code{\link[biomaRt]{listFilters}}.}
  \item{addNbOcc}{If TRUE, the number of components the features/genes
    contribute to is added to the output. A gene/feature is considered
    as a contributor of a component if its absolute scaled projection
    value is higher than \code{selCutoff(icaSet)}.}
  \item{selectionByComp}{A list containing the feature/gene projections
    on each component, already restricted to the ones considered as contributors.}
  \item{level}{The data level of \code{icaSet} that will be annotated:
    either the feature projections ("features"), or the gene projections ("genes").}
  \item{typeId}{The type of ID the features or the genes of
    \code{icaSet} correspond to. By default \code{typeID(icaSet)} is used.  It must be provided in the biomaRt way
    (type \code{listFilters(mart)} to choose the appropriate value).}
  \item{selCutoffWrite}{The cutoff applied to the absolute projection values to select the features/genes that will be annotated using package \code{biomaRt}, default is 2.5.}
}
\details{One file is created by component, each file is named by the
  index of the components (\code{indComp(icaSet)}) and located in the
  path \code{genePath(params)}.

In case you are interested in writing the description of features and
their annotations, please remember to modify code{genesPath(params)}, or
the previous files will be overwritten.

The genes are ranked according to their absolute projection values.

This function also writes an html file named "genes2comp"
providing, for each feature or gene, the number of components it
contributes to (according to the threshold \code{cutoffSel(params)}),
and its projection value on all the components.  
The projection values are scaled.

See function \code{\link{writeGenes}} for details.
}
\value{This function returns a list of two elements:\describe{
    \item{listAnnotComp:}{a list with the output of
      \code{\link{writeGenes}} for each component} 
    \item{nbOccInComp:}{a data.frame storing the projection
  values of each feature/gene (row) across all the components
  (columns).}}
}
\examples{
\dontrun{
## load IcaSet object
## We will use 'icaSetCarbayo', whose features are hgu133a probe sets
## and feature annotations are Gene Symbols. 
data(icaSetCarbayo)

## define database to be used by biomaRt
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

## define the parameters of the analysis
params <- buildMineICAParams(resPath="~/resMineICACarbayo/", selCutoff=0)

## Make sure the elements "_biomaRt" of attribute 'typeID' are defined
typeID(icaSetCarbayo) 

### Query biomaRt and write gene descriptions in HTML files
### The files will be located in the directory 'genesPath(params)'

## 1. Write description of genes 
res <- writeProjByComp(icaSet=icaSetCarbayo, params=params, mart=mart,
           level="genes") #, typeId="hgnc_symbol")

## 2. Write description of features 
# change attribute 'genesPath' of params to preserve the gene descriptions
genesPath(params) <- paste(resPath(params),"comp2features/",sep="")
res <- writeProjByComp(icaSet=icaSetCarbayo, params=params, mart=mart,
           level="features") #, typeId="affy_hg_u133a")

}
}

\author{Anne Biton}

\seealso{\code{\link{writeGenes}}, \code{getBM}, \code{listFilters}, \code{listAttributes}, \code{useMart}, \code{\link{selectContrib}}, \code{\link{nbOccInComp}}}

