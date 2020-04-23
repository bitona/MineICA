\name{runAn}
\alias{runAn}
\title{Run analysis of an IcaSet object}
\usage{
  runAn(params, icaSet, keepVar,
    heatmapCutoff = params["selCutoff"],
    funClus = c("Mclust", "kmeans"), nbClus,
    clusterOn = "A", keepComp, keepSamples,
    adjustBy = c("none", "component", "variable"),
    typePlot = c("boxplot", "density"),
    mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"),
    dbGOstats = c("KEGG", "GO"), ontoGOstats = "BP",
    condGOstats = TRUE,
    cutoffGOstats = params["pvalCutoff"],
    writeGenesByComp = TRUE, writeFeaturesByComp = FALSE,
    selCutoffWrite = 2.5, runVarAnalysis = TRUE,
    onlySign = T, runClustering = FALSE, runGOstats = TRUE,
    plotHist = TRUE, plotHeatmap = TRUE)
}
\arguments{
  \item{params}{An object of class
  \code{\link[MineICA:class-MineICAParams]{MineICAParams}}
  containing the parameters of the analysis.}

  \item{icaSet}{An object of class
  \code{\link[MineICA:class-IcaSet]{IcaSet}}.}

  \item{keepVar}{The variable labels to be considered, i.e
  a subset of the annotation variables available in
  (\code{varLabels(icaSet)}).}

  \item{keepSamples}{The samples to be considered, i.e a
  subset of (\code{sampleNames(icaSet)}).}

  \item{heatmapCutoff}{The cutoff (applied to the scaled
  feature/gene projections contained in S/SByGene) used to
  select the contributing features/genes.}

  \item{funClus}{The function to be used to cluster the
  samples, must be one of
  \code{c("Mclust","kmeans","pam","pamk","hclust","agnes")}.
  Default is \code{"Mclust"}.}

  \item{nbClus}{The number of clusters to be computed when
  applying \code{funClus}. Can be missing (default) if
  \code{funClus="Mclust"} or \code{funClus="pamk"}.}

  \item{keepComp}{The indices of the components to be
  analyzed, must be included in \code{indComp(icaSet)}. If
  missing, all components are treated.}

  \item{adjustBy}{The way the p-values of the Wilcoxon and
  Kruskal-Wallis tests should be corrected for multiple
  testing: \code{"none"} if no p-value correction has to be
  done, \code{"component"} if the p-values have to be
  corrected by component, \code{"annotation"} if the
  p-values have to be corrected by variable}

  \item{typePlot}{The type of plot used to show
  distribution of sample-groups contributions, either
  "density" or "boxplot"}

  \item{mart}{A mart object used for annotation, see
  function \code{\link[biomaRt]{useMart}}}

  \item{dbGOstats}{The used database to use ('GO' and/or
  'KEGG'), default is both.}

  \item{ontoGOstats}{A string specifying the GO ontology to
  use. Must be one of 'BP', 'CC', or 'MF', see
  \code{\link[Category:class-GOHyperGParams]{GOHyperGParams}}.
  Only used when argument \code{dbGOstats} is 'GO'.}

  \item{condGOstats}{A logical indicating whether the
  calculation should conditioned on the GO structure, see
  \code{\link[Category:class-GOHyperGParams]{GOHyperGParams}}.}

  \item{cutoffGOstats}{The p-value threshold used for
  selecting enriched gene sets, default is
  params["pvalCutoff"]}

  \item{writeGenesByComp}{If TRUE (default) the gene
  projections (\code{SByGene(icaSet)}) are written in an
  html file and annotated using \code{biomaRt} for each
  component.}

  \item{writeFeaturesByComp}{If TRUE (default) the feature
  projections (\code{S(icaSet)}) are written in an html
  file and annotated using \code{biomaRt} for each
  component.}

  \item{runGOstats}{If TRUE the enrichment analysis of the
  contributing genes is run for each component using
  package \code{GOstats} (default is TRUE).}

  \item{plotHist}{If TRUE the position of the sample
  annotations within the histograms of the sample
  contributions are plotted.}

  \item{plotHeatmap}{If TRUE the heatmap of the
  contributing features/genes are plotted for each
  component.}

  \item{runClustering}{If TRUE the potential associations
  between a clustering of the samples (performed according
  to the components), and the sample annotations, are
  tested using chi-squared tests.}

  \item{runVarAnalysis}{If TRUE the potential associations
  between sample contributions (contained in
  \code{A(icaSet)}) are tested using Wilcoxon or
  Kruskal-Wallis tests.}

  \item{onlySign}{If TRUE (default), only the significant
  results are plotted in functions \code{qualVarAnalysis,
  quantVarAnalysis, clusVarAnalysis}, else all plots are
  done.}

  \item{selCutoffWrite}{The cutoff applied to the absolute
  feature/gene projection values to select the
  features/genes that will be annotated using package
  \code{biomaRt}, default is 2.5.}

  \item{clusterOn}{Specifies the matrix used to apply
  clustering if \code{runClustering=TRUE}: \describe{
  \item{\code{"A"}:}{the clustering is performed in one
  dimension, on the vector of sample contributions,}
  \item{"S":}{the clustering is performed on the original
  data restricted to the contributing individuals,}
  \item{"AS":}{the clustering is performed on the matrix
  formed by the product of the column of A and the row of
  S.}}}
}
\value{
  NULL
}
\description{
  This function runs the analysis of an ICA decomposition
  contained in an IcaSet object, according to the
  parameters entered by the user and contained in a
  MineICAParams.
}
\details{
  This function calls functions of the MineICA package
  depending on the arguments: \describe{
  \item{\code{\link{writeProjByComp}} (if
  \code{writeGenesByComp=TRUE} or
  \code{writeFeaturesByComp})}{which writes in html files
  the description of the features/genes contributing to
  each component, and their projection values on all the
  components.} \item{\code{\link{plot_heatmapsOnSel}} (if
  \code{plotHeatmap=TRUE})}{which plots heatmaps of the
  data restricted to the contributing features/genes of
  each component.} \item{\code{\link{plotPosAnnotInComp}}
  (if \code{plotHist=TRUE})}{which plots, within the
  histogram of the sample contribution values of every
  component, the position of groups of samples formed
  according to the sample annotations contained in
  \code{pData(icaSet)}.}
  \item{\code{\link{clusterSamplesByComp}} (if
  \code{runClustering=TRUE})}{which clusters the samples
  according to each component.}
  \item{\code{\link{clusVarAnalysis}} (if
  \code{runClustering=TRUE})}{which computes the
  chi-squared test of association between a given
  clustering of the samples and each annotation level
  contained in \code{pData(icaSet)}, and summarizes the
  results in an HTML file. } \item{\code{\link{runEnrich}}
  (if \code{runGOstats=TRUE})}{which perforns enrichment
  analysis of the contributing genes of the components
  using package \link{GOstats}.}
  \item{\code{\link{qualVarAnalysis}} and
  \code{\link{quantVarAnalysis}} (if
  \code{varAnalysis=TRUE})}{which tests if the groups of
  samples formed according to sample annotations contained
  in \code{pData(icaSet)} are differently distributed on
  the components, in terms of contribution value. } }

  Several directories containing the results of each
  analysis are created by the function: \describe{
  \item{ProjByComp:}{contains the annotations of the
  features or genes, one file per component;}
  \item{varAnalysisOnA:}{contains two directories: 'qual/'
  and 'quant/' which respectively contain the results of
  the association between components qualitative and
  quantitative variables;} \item{Heatmaps:}{contains the
  heatmaps (one pdf file per component) of contributing
  genes by component;} \item{varOnSampleHist:}{contains
  athe histograms of sample contributions superimposed with
  the histograms of the samples grouped by variable;}
  \item{cluster2var:}{contains the association between a
  clustering of the samples performed on the mixing matrix
  \code{A} and the variables.} }
}
\examples{
\dontrun{

## load an example of IcaSet
data(icaSetCarbayo)
## make sure the 'mart' attribute is correctly defined
mart(icaSetCarbayo) <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

## creation of an object of class MineICAParams
## here we use a low threshold because 'icaSetCarbayo' is already
# restricted to the contributing features/genes
params <- buildMineICAParams(resPath="~/resMineICACarbayotestRunAn/", selCutoff=2, pvalCutoff=0.05)
require(hgu133a.db)

runAn(params=params, icaSet=icaSetCarbayo)
}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{writeProjByComp}},
}

