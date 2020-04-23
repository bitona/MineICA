\name{buildMineICAParams}
\alias{buildMineICAParams}
\title{Creates an object of class MineICAParams}
\usage{
  buildMineICAParams(Sfile = new("character"),
    Afile = new("character"), datfile = new("character"),
    annotfile = new("character"), resPath = "", genesPath,
    annot2col = new("character"), pvalCutoff = 0.05,
    selCutoff = 3)
}
\arguments{
  \item{Sfile}{A txt file containing the Source matrix S.}

  \item{Afile}{A txt file containing the Mixing matrix A.}

  \item{datfile}{A txt file containing the data (e.g
  expression data) on which the decomposition was
  calculated.}

  \item{annotfile}{Either a "rda" or "txt" file containing
  the annotation data for the samples (must be of
  dimensions samples x annotations).}

  \item{resPath}{The path where the outputs of the analysis
  will be written, default is the current directory.}

  \item{genesPath}{The path _within_ the resPath where the
  gene projections will be written. If missing, will be
  automatically attributed as \code{resPath}/ProjByComp/.}

  \item{annot2col}{A vector of colors indexed by annotation
  levels. If missing, will be automatically attributed
  using function \code{annot2Color}.}

  \item{pvalCutoff}{The cutoff used to consider a p-value
  significant, default is 0.05.}

  \item{selCutoff}{The cutoff applied to the absolute
  feature/gene projection values to consider them as
  contributors, default is 3. Must be either of length 1
  and the same treshold is applied to all components, or of
  length equal to the number of components in order to a
  specific threshold is for each component.}
}
\value{
  An object of class \code{\link{MineICAParams}}
}
\description{
  This function builds an object of class
  \code{\link[MineICA:class-MineICAParams]{MineICAParams}}.
  It contains the parameters that will be used by function
  \code{\link{runAn}} to analyze the ICA decomposition
  contained in an object of class
  \code{\link[MineICA:class-IcaSet]{IcaSet}}.
}
\examples{
## define default parameters and fill resPath
params <- buildMineICAParams(resPath="resMineICACarbayo/")

## change the default cutoff for selection of contribugint genes/features
params <- buildMineICAParams(resPath="resMineICACarbayo/", selCutoff=4)
}
\author{
  Anne Biton
}
\seealso{
  \code{\link[MineICA:class-MineICAParams]{MineICAParams}},
  \code{\link[MineICA]{runAn}}
}

