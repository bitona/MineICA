\name{MineICAParams}
\docType{class}
\alias{class:MineICAParams}
\alias{MineICAParams}
\alias{MineICAParams-class}
\alias{Sfile}
\alias{getSfile}
\alias{setSfile}
\alias{Sfile<-}
\alias{Sfile,MineICAParams-method}
\alias{Sfile<-,MineICAParams,character-method}
\alias{datfile}
\alias{getdatfile}
\alias{setdatfile}
\alias{datfile<-}
\alias{annotfile<-,MineICAParams-method}
\alias{genesPath<-,ANY-method}
\alias{resPath<-,ANY-method}
\alias{selCutoff<-,MineICAParams-method}
\alias{pvalCutoff<-,MineICAParams-method}
\alias{Afile<-,MineICAParams-method}
\alias{Sfile<-,MineICAParams-method}
\alias{annot2col<-,MineICAParams-method}
\alias{datfile<-,MineICAParams-method}
\alias{datfile,MineICAParams-method}
\alias{datfile<-,MineICAParams,character-method}
\alias{Afile}
\alias{getAfile}
\alias{setAfile}
\alias{Afile<-}
\alias{Afile,MineICAParams-method}
\alias{Afile<-,MineICAParams,character-method}
\alias{annotfile}
\alias{getAnnotfile}
\alias{setAnnotfile}
\alias{annotfile<-}
\alias{annotfile,MineICAParams-method}
\alias{annotfile<-,MineICAParams,character-method}
\alias{resPath}
\alias{getResPath}
\alias{setResPath}
\alias{resPath<-}
\alias{resPath,MineICAParams-method}
\alias{resPath<-,MineICAParams,character-method}
\alias{genesPath}
\alias{getGenesPath}
\alias{setGenesPath}
\alias{genesPath<-}
\alias{genesPath,MineICAParams-method}
\alias{genesPath<-,MineICAParams,character-method}
\alias{annot2col}
\alias{getAnnot2col}
\alias{setAnnot2col}
\alias{annot2col<-}
\alias{annot2col,MineICAParams-method}
\alias{annot2col<-,MineICAParams,character-method}
\alias{pvalCutoff}
\alias{getPvalCutoff}
\alias{setPvalCutoff}
\alias{pvalCutoff<-}
\alias{pvalCutoff,MineICAParams-method}
\alias{pvalCutoff<-,MineICAParams,numeric-method}
\alias{selCutoff}
\alias{getSelCutoff}
\alias{setSelCutoff}
\alias{selCutoff<-}
\alias{selCutoff,MineICAParams-method}
\alias{selCutoff<-,MineICAParams,numeric-method}

%\alias{[}
\alias{[,ANY,MineICAParams-method} 
\alias{[,ANY,ANY,MineICAParams-method}
\alias{[,ANY,ANY,ANY,MineICAParams-method}
\alias{[,MineICAParams,ANY-method} 
\alias{[,MineICAParams,ANY,ANY-method} 
\alias{[,MineICAParams,ANY,ANY,ANY-method} 
%\alias{[<-}
\alias{[<-,MineICAParams,ANY,ANY,ANY,ANY-method} 
\alias{[<-,MineICAParams,ANY,ANY,ANY-method} 
\alias{[<-,MineICAParams,ANY,ANY-method} 
%\alias{[<-,MineICAparams,ANY-method} 


\title{
  Class to contain parameters for the analysis of an ICA decomposition.
}

\description{
  Container for parameters used during the analysis of an ICA
  decomposition obtained on genomics data.
}
\section{Creating Objects}{

  \code{new("MineICAParams")}

  \code{new("MineICAParams",
    resPath="",
    genesPath="ProjByComp",
    pvalCutoff=0.05,
    selCutoff=3)
  }
}

\section{Slots}{
  \describe{
\item{\code{Sfile}}{ A txt file containing the Source matrix S.} 
\item{\code{Afile}}{ A txt file containing the Mixing matrix A.} 
\item{\code{datfile}}{ A txt file containing the data (typically expression data) on which the decomposition was calculated.}
\item{\code{annotfile}}{ Either a RData or txt file containing the annotation data for the samples (must be of dimensions samples*annotations).} 
\item{\code{resPath}}{ The path where the outputs of the analysis will be written.} 
\item{\code{genesPath}}{ The path _within_ the resPath where the gene projections will be written. If missing, will be automatically attributed as resPath/gene2components/.}
\item{\code{annot2col}}{ A vector of colors indexed by annotation levels. If missing, will be automatically attributed using function \code{annot2Color}.}
\item{\code{pvalCutoff}}{ The cutoff used to consider a p-value significant, default is 0.05.}
\item{\code{selCutoff}}{ The cutoff applied on the absolute feature/gene
  projection values to consider gene as contributing to a component,
  default is 3. Must be either of length 1 and the same treshold is
  applied to all components, or of length equal to the number of
  components in order to use a specific threshold for each component.} 
  }
  
}
\section{Methods}{
  \describe{
     For any slot:
     \item{\code{slotName(MineICAParams)} and
       \code{slotName(MineICAParams)<-}:}{Accessing and setting any slot
       of name \code{slotName} contained in an MineICAParams object.}   
      \item{\code{MineICAParams["slotName"]} and
	\code{MineICAParams["slotName"]<-}:}{Accessing and setting any
	slot of name \code{slotName} contained in an MineICAParams object.}   
   }
}
   

  



\author{Anne Biton}

\seealso{
  \code{\link{MineICAParams-class}}, \code{\link{runAn}}.
}

\examples{
# create an instance of LocSet
new("MineICAParams")


}

\keyword{classes}
