\name{IcaSet}
\docType{class}
\alias{class:IcaSet}
\alias{IcaSet}
\alias{IcaSet-class}

%\alias{dat}  % pour la generique
%\alias{dat,IcaSet-method}% # avec la signature de ma classe
%\alias{dat<-} %# la generique de remplacement
%\alias{dat<-,IcaSet,matrix-method} %# avec la signature de la classe 

\alias{[}
\alias{[,ANY,ANY,IcaSet-method}
\alias{[,IcaSet,ANY-method} 
\alias{[,IcaSet,ANY,ANY-method} 
\alias{[,IcaSet,ANY,ANY,ANY-method} 
\alias{[<-}
\alias{[<-,IcaSet,ANY,ANY,ANY,ANY-method} 
\alias{[<-,IcaSet,ANY,ANY,ANY-method} 
\alias{[<-,IcaSet,ANY,ANY-method} 

\alias{organism}
\alias{organism<-}
\alias{organism,IcaSet-method}
\alias{organism<-,IcaSet-method}

\alias{mart}
\alias{mart<-}
\alias{mart,IcaSet-method}
\alias{mart<-,IcaSet-method}
\alias{mart<-,IcaSet,character-method}
\alias{setMart,IcaSet-method}
\alias{getMart,IcaSet-method}

\alias{chipVersion}
\alias{chipVersion,IcaSet-method}
\alias{chipVersion<-}
\alias{chipVersion<-,IcaSet-method}
\alias{chipVersion<-,IcaSet,character-method}

%\alias{indComp}
%\alias{indComp<-}
\alias{indComp,IcaSet-method}
\alias{setIndComp,IcaSet-method}
\alias{getIndComp,IcaSet-method}
\alias{indComp<-,IcaSet,character-method}
%\alias{compNames}
%\alias{compNames<-}
\alias{compNames,IcaSet-method}
\alias{setLabelsComp,IcaSet-method}
\alias{getLabelsComp,IcaSet-method}
\alias{compNames<-,IcaSet,character-method}
%\alias{witGenes}
%\alias{witGenes<-}
\alias{witGenes,IcaSet-method}
%\alias{setWitGenes,IcaSet-method}
%\alias{getWitGenes,IcaSet-method}
\alias{witGenes<-,IcaSet,character-method}
\alias{refSamples}
\alias{refSamples<-}
\alias{refSamples,IcaSet-method}
\alias{setRefSamples,IcaSet-method}
\alias{getRefSamples,IcaSet-method}
\alias{refSamples<-,IcaSet-method}
\alias{refSamples<-,IcaSet,character-method}
\alias{typeID}
\alias{typeID<-}
\alias{typeID,IcaSet-method}
\alias{typeID<-,IcaSet-method}
\alias{setTypeID,IcaSet-method}
\alias{getTypeID,IcaSet-method}
\alias{typeID<-,IcaSet,list-method}
\alias{chipManu}
\alias{chipManu<-}
\alias{chipManu<-,IcaSet-method}
\alias{chipManu<-,IcaSet,character-method}
\alias{chipManu,IcaSet-method}
\alias{setChipManu,IcaSet-method}
\alias{getChipManu,IcaSet-method}


\alias{Slist,IcaSet-method}
\alias{SlistByGene,IcaSet-method}
\alias{Alist,IcaSet-method}



\title{
  Class to Contain and Describe an ICA decomposition of High-Throughput Data.
}

\description{
  Container for high-throughput data and results of ICA decomposition
  obtained on these data. \code{IcaSet} class is derived from
  \code{\link{eSet}}, and requires a matrix named \code{dat} as
  \code{assayData} member.
}
\section{Extends}{
   Directly extends class \code{\link{eSet}}.
}
\section{Creating Objects}{

  \code{new("IcaSet")}

  \code{new("IcaSet",
    annotation = character(0),
    experimentData = new("MIAME"),
    featureData = new("AnnotatedDataFrame"),
    phenoData = new("AnnotatedDataFrame"),
    protocolData = phenoData[,integer(0)],
    dat = new("matrix"),
    A=new("data.frame"),
    S=new("data.frame"), ...)
    }

  This creates an \code{IcaSet} with \code{assayData} implicitly
  created to contain \code{dat}. %Additional named matrix arguments
  %with the same dimensions as \code{dat} are added to
  %\code{assayData}; the row and column names of these additional
  %matrices should match those of \code{dat}. 

  \code{new("IcaSet",
    annotation = character(0),
    assayData = assayDataNew(dat=new("matrix")),
    experimentData = new("MIAME"),
    featureData = new("AnnotatedDataFrame"),
    phenoData = new("AnnotatedDataFrame"),
    protocolData = phenoData[,integer(0)],
    A=new("data.frame"),
    S=new("data.frame"), ...)
  }

  This creates an \code{IcaSet} with \code{assayData} provided
  explicitly. %In this form, the only required named arguments are \code{assayData}.


  \code{IcaSet} instances are usually created through
  \code{new("IcaSet", ...)}. Usually the arguments to \code{new}
  include \code{dat} ('features x samples', e.g a matrix of expression
  data), \code{phenoData} ('samples x annotations', a
  matrix of sample annotations), \code{S} the Source
  matrix of the ICA decomposition ('features x comp'), \code{A} the Mixing matrix of the ICA
  decomposition ('samples x comp'), \code{annotation} the annotation
  package, \code{typeID} the description of the feature and gene IDs.  

  The other attributes can be missing, in which case
  they are assigned default values.

  The function \code{\link{buildIcaSet}} is a more convenient way to
  create \code{IcaSet} instances, and allows to automatically annotate
  the features.
}
\section{Slots}{
  Inherited from \code{eSet}:
  \describe{
    \item{\code{annotation}:}{See \code{\link{eSet}}}
    \item{\code{assayData}:}{Contains matrices with equal
      dimensions, and with column number equal to
      \code{nrow(phenoData)}. \code{assayData} must contain a matrix
      \code{dat} with rows representing features (e.g., reporters)
      and columns representing samples. Class:\code{\link{AssayData-class}}}
    \item{\code{experimentData}:}{See \code{\link{eSet}}}
    \item{\code{featureData}:}{See \code{\link{eSet}}}
    \item{\code{phenoData}:}{See \code{\link{eSet}}}
    \item{\code{protocolData}:}{See \code{\link{eSet}}}
  
  Specific slot:

    \item{\code{organism}:}{Contains the name of the species. Currently
  only Human ("Human" or "Homo sapiens") and Mouse ("Mouse" or "Mus
  musculus") are supported. Only used when \code{chipManu}="illumina"
  }
    \item{\code{mart}:}{An output of \code{\link[biomaRt]{useMart}} of package \code{biomaRt}. Only useful if no annotation package is available for argument \code{icaSet}.
  }
  
  \item{\code{datByGene}:}{Data.frame containing the data \code{dat} where
  features have been replaced by their annotations (e.g, gene IDs). Rows 
  represent annotations of the features (e.g., gene IDs) and
  columns represent samples.}
  \item{\code{A}:}{The mixing matrix of the ICA decomposition, contained
  in a data.frame whose
  column number equals the number of components and row number equals
  \code{nrow(phenoData)} (dimension: 'samples x comp').}
  \item{\code{S}:}{The source matrix of the ICA decomposition, contained
  in a data.frame whose
  column number equals the number of components and row number equals
  \code{nrow(assayData)} (dimension: 'features x comp').}
  \item{\code{SByGene}:}{The matrix Source of the ICA decomposition, contained
  in a data.frame whose
  column number equals the number of components and row number equals
  \code{nrow(datByGene)} (dimension: 'annotatedFeatures x comp').}
  \item{\code{compNames}:}{A vector of component labels with length
  equal to the number of component.}
  \item{\code{indComp}:}{A vector of component indices with length
  equal to the number of component.}
  \item{\code{witGenes}:}{A vector of gene IDs with length
  equal to the number of component.}
  \item{\code{chipManu}:}{The manufacturer of the technology the
  data originates from. Useful for the annotation of the features when
  data originates from an _illumina_ microarray.}
  \item{\code{chipVersion}:}{The version of the chip, only useful for
  when \code{chipManu}="illumina"}
   \item{\code{refSamples}:}{A vector of sample IDs including the reference samples, e.g the "normal" samples.
  Must be included in \code{sampleNames(object)}, i.e in \code{colnames(dat)}.}
   \item{\code{typeID}:}{A vector of characters providing the annotation IDs. It includes
     three elements:
 \describe{ 
 \item{geneID_annotation}{the IDs from the
 package to be used to annotate the features into genes. It will be used to
 fill the attributes \code{datByGene} and \code{SByGene} of the \code{icaSet}.
 It must match one of the objects the corresponding package supports
 (you can access the list of objects by typing ls("package:packagename")). If
 no annotation package is provided, this element is not useful.}
 \item{geneID_biomart}{the type of gene IDs, as available in
 \code{listFilters(mart)}; where mart is specified as described in \code{\link[biomaRt]{useMart}}.
 If you have directly built the IcaSet at the
 gene level (i.e if no annotation package is used), \code{featureID_biomart} and
 \code{geneID_biomart} will be identical.} 
 \item{featureID_biomart}{the
 type of feature IDs, as available in \code{listFilters(mart)}; where
 \code{mart} is specified as described in function \code{\link[biomaRt]{useMart}}.
 Not useful if you work at the gene level.} }}

}
}
\section{Methods}{
  
  Class-specific methods.
  \describe{
   
      \item{\code{getComp(IcaSet, ind,
      level=c("features","genes"))}}{Given a component index, extract
      the corresponding sample contribution values from A, and the
      feature (\code{level}="features") or gene (\code{level}="genes")
      projections from S. Returns a list with two elements:
      \code{contrib} the sample contributions and \code{proj} the
      feature or gene projections.}
     
     Access and set any slot specific to IcaSet:
     \item{\code{slotName(IcaSet)}, and
       \code{slotName(IcaSet)<-}:}{Accessing and setting any slot
       of name \code{slotName} contained in an IcaSet object.}
     
    \item{\code{IcaSet["slotName"]}, and
       \code{IcaSet["slotName"]<-}:}{Accessing and setting any slot
       of name \code{slotName} contained in an IcaSet object.}
     
     Most used accessors and settors:
     \item{\code{A(IcaSet)}, and
       \code{A(IcaSet)<-}:}{Accessing and setting Mixing matrix \code{A}.}
     \item{\code{S(IcaSet)}, and
       \code{S(IcaSet)<-}:}{Accessing and setting
       the data.frame Source \code{S}.}
     \item{\code{Slist(IcaSet)}:}{Accessing
       the data.frame Source as a list where names are preserved.}
     \item{\code{SByGene(IcaSet)}, and
       \code{SByGene(IcaSet)<-}:}{Accessing
       and setting the _annotated_ data.frame Source \code{SByGene}.}
     \item{\code{SlistByGene(IcaSet)}:}{Accessing
       the _annotated_ Source matrix as a list where names are preserved.}
      \item{\code{organism(IcaSet)}, \code{organism(IcaSet,characte)<-}}{Access and
	set value in the \code{organism} slot.}
     \item{\code{dat(IcaSet)}, \code{dat(IcaSet,matrix)<-}}{Access and
       set elements named \code{dat} in the \code{AssayData-class}
       slot.}

    }

   Derived from \code{\link{eSet}}:
   \describe{
     \item{\code{pData(IcaSet)}, \code{pData(IcaSet,value)<-}:}{See \code{\link{eSet}}}
     \item{\code{assayData(IcaSet)}:}{See \code{\link{eSet}}}
     \item{\code{sampleNames(IcaSet)} and \code{sampleNames(IcaSet)<-}:}{See \code{\link{eSet}}}
     \item{\code{featureNames(IcaSet)}, \code{featureNames(IcaSet, value)<-}:}{See \code{\link{eSet}}}
     \item{\code{dims(IcaSet)}:}{See \code{\link{eSet}}}
     \item{\code{phenoData(IcaSet)}, \code{phenoData(IcaSet,value)<-}:}{See \code{\link{eSet}}}
     \item{\code{varLabels(IcaSet)}, \code{varLabels(IcaSet, value)<-}:}{See \code{\link{eSet}}}
     \item{\code{varMetadata(IcaSet)}, \code{varMetadata(IcaSet,value)<-}:}{See \code{\link{eSet}}}
     \item{\code{varMetadata(IcaSet)}, \code{varMetadata(IcaSet,value)}}{See \code{\link{eSet}}}
     \item{\code{experimentData(IcaSet)},\code{experimentData(IcaSet,value)<-}:}{See \code{\link{eSet}}}
     \item{\code{pubMedIds(IcaSet)}, \code{pubMedIds(IcaSet,value)}}{See \code{\link{eSet}}}
     \item{\code{abstract(IcaSet)}:}{See \code{\link{eSet}}}
     \item{\code{annotation(IcaSet)}, \code{annotation(IcaSet,value)<-}}{See \code{\link{eSet}}}
     \item{\code{protocolData(IcaSet)}, \code{protocolData(IcaSet,value)<-}}{See \code{\link{eSet}}}
     \item{\code{combine(IcaSet,IcaSet)}:}{See \code{\link{eSet}}}
     \item{\code{storageMode(IcaSet)}, \code{storageMode(IcaSet,character)<-}:}{See \code{\link{eSet}}}
   }
  
  Standard generic methods:
  \describe{
    \item{\code{initialize(IcaSet)}:}{Object instantiation, used
      by \code{new}; not to be called directly by the user.}
    
    \item{\code{validObject(IcaSet)}:}{Validity-checking method, ensuring
      that \code{dat} is a member of
      \code{assayData}, and that the number of features, genes, samples,
      and components are consistent across all the attributes of the
      IcaSet object. \code{checkValidity(IcaSet)} imposes this
      validity check, and the validity checks of \code{eSet}.}
    \item{\code{IcaSet[slotName]}, \code{IcaSet[slotName]<-}:}{Accessing
      and setting any slot of name \code{slotName} contained in an
      IcaSet object.}
    \item{\code{IcaSet[i, j, k]}:}{Extract object of class "IcaSet"
          for features or genes with names i, samples with names or
      indices j, and components with names or indices k.}
       \item{\code{makeDataPackage(object, author, email, packageName, packageVersion, license, biocViews, filePath, description=paste(abstract(object), collapse="\n\n"), ...)}}{
      Create a data package based on an IcaSet object. See
      \code{\link{makeDataPackage}}.}
     \item{\code{show(IcaSet)}:}{See \code{\link{eSet}}}
     \item{\code{dim(IcaSet)}, \code{ncol}:}{See \code{\link{eSet}}}
     \item{\code{IcaSet[(index)]}:}{See \code{\link{eSet}}}
     \item{\code{IcaSet$}, \code{IcaSet$<-}:}{See \code{\link{eSet}}}
     \item{\code{IcaSet[[i]]}, \code{IcaSet[[i]]<-}:}{See \code{\link{eSet}}}
  }
}

\author{Anne Biton}

\seealso{
  \code{\link{eSet-class}}, \code{\link{buildIcaSet}},
  \code{\link{class-IcaSet}}, \code{\link{class-MineICAParams}}.
}

\examples{
# create an instance of IcaSet
new("IcaSet")
dat <- matrix(runif(100000), nrow=1000, ncol=100)
rownames(dat) <- 1:nrow(dat)
new("IcaSet",
    dat=dat, 
    A=as.data.frame(matrix(runif(1000), nrow=100, ncol=10)),
    S=as.data.frame(matrix(runif(10000), nrow=1000, ncol=10), row.names = 1:nrow(dat)))


}

\keyword{classes}
