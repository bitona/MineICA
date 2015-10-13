setMethod("initialize",
          signature="IcaSet",
          definition = function(.Object,
                   annotation      = character(),
                   assayData       = assayDataNew(dat=dat, ...),
                   experimentData  = new("MIAME"),
                   featureData     = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   phenoData       = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   A=new("data.frame"),
                   S=new("data.frame"),
                   dat=new("matrix"),   
                   SByGene=new("data.frame"),
                   datByGene=new("data.frame"),
                   compNames=character(), 
                   indComp=numeric(), 
                   witGenes=character(),
                   refSamples=character(),
                   chipManu=character(),
                   chipVersion=character(),
                   typeID=character(),
                   organism=character(),
                   mart=useMart("ensembl"),...){

              .Object@A <- A
              .Object@S <- S

              .Object@SByGene <- SByGene
              if (length(compNames)==0 & ncol(A)>0) 
                  compNames <- paste("Component",as.character(c(1:ncol(A))))
              if (length(indComp)==0 & ncol(A)>0)
                  indComp <- 1:ncol(A)

              colnames(.Object@A) <- colnames(.Object@S) <-indComp

              .Object@compNames <- compNames
              .Object@indComp <- indComp
              .Object@refSamples <- refSamples
              .Object@witGenes <- witGenes
              

              .Object@datByGene <- datByGene
              .Object@chipManu <- chipManu
              .Object@chipVersion <- chipVersion

              .Object@typeID <- typeID

              .Object@organism <- organism


              if (length(annotation)>0)
                  library(annotation,character.only=TRUE)
              
              .Object@mart <- mart
              
              callNextMethod(.Object,
                             assayData       = assayData,
                             phenoData       = phenoData,
                             featureData     = featureData,
                             experimentData  = experimentData,
                             annotation      = annotation,
                             organism = organism)

              
          })


validIcaSet <- function(object) {
    msg <- NULL
    msg <- Biobase:::validMsg(msg, Biobase:::isValidVersion(object, "IcaSet"))

    if (dim(assayDataElement(object,"dat"))[[1]] != dim(S(object))[[1]])
        msg <- paste(msg, "Assay data and matrix S must have the same number of rows (i.e the same number of genes or probe sets)", sep="\n")
    if (dim(datByGene(object))[[2]] != dim(assayDataElement(object,"dat"))[[2]])
        msg <- paste(msg, "expression data by gene (datByGene) and assayData must have the same number of columns", sep="\n")
    if (dim(assayDataElement(object,"dat"))[[2]] != dim(A(object))[[1]])
        msg <- paste(msg, "Assay data and matrix A must include the same number of samples",sep = "\n")
    if (dim(A(object))[[2]] != dim(S(object))[[2]])
        msg <- paste(msg, "Matrix A and S must include the same number of components", sep = "\n")
    if (!all(rownames(A(object)) == sampleNames(assayData(object))))
        msg <- paste(msg, "Row names of A must match the sampleNames of assayData", sep = "\n")
    if (!all(rownames(S(object)) == featureNames(object)))
        msg <- paste(msg, "Rownames of S must match the rownames of assayData.", sep = "\n")
    if (nrow(SByGene(object))>0 & nrow(datByGene(object))>0) {
        if (!all(rownames(SByGene(object)) == rownames(datByGene(object))))
            msg <- paste(msg, "Names of each SByGene element must match the rownames of datByGene.", sep = "\n")
        if (ncol(SByGene(object)) != nbComp(object))
            msg <- paste(msg, "The number of elements in SByGene must match the number of components", sep = "\n")
 
    }
    if (length(compNames(object)) != nbComp(object))
        msg <- paste(msg, "The number of labels must match the number of components", sep = "\n")
    if (length(indComp(object)) != nbComp(object))
        msg <- paste(msg, "The number of indices must match the number of components", sep = "\n")
    if (length(witGenes(object))>0) {
        if (length(witGenes(object)) != nbComp(object))
            msg <- paste(msg, "The number of 'witness' genes must match the number of components", sep = "\n")
        if (length(setdiff(witGenes(object),rownames(assayDataElement(object,"dat"))))>0) {
            if (length(setdiff(witGenes(object),geneNames(object)))>0)
                msg <- paste(msg, "The 'witness' genes must be included in the rownames of the assay data.", sep = "\n")
        }
    }
   if (length(refSamples(object))>0)
        if (length(setdiff(refSamples(object),sampleNames(object)))>1)
            msg <- paste(msg, "All the reference samples are not included within the sample names.")       


    if (is.null(msg))
        TRUE
    else 
        msg
    
}

setValidity("IcaSet", function(object) {
     validIcaSet(object)
 })

     setMethod(
      f ="[",
      signature(x = "IcaSet", "ANY","ANY","ANY"),
      definition = function (x, i, j, k, ..., drop=FALSE){

      if (!missing(i)) {    
          if (is.character(i)) {
              ifeatures <- intersect(i,featureNames(x))
              if (length(intersect(i,featureNames(x)))>0) {
                  x@S <- x@S[ifeatures,,drop=FALSE]
                  ff <- matrix(nrow=length(ifeatures),ncol=0)
                  rownames(ff) <- ifeatures
                  x@featureData <- annotatedDataFrameFrom(ff, byrow=TRUE)
                  dat(x) <- dat(x)[ifeatures,,drop=FALSE]
              }
              else {
                  if (length(intersect(i,geneNames(x)))>0) {
                      igenes <- intersect(i,geneNames(x))
                      x@SByGene <- x@S[igenes,,drop=FALSE]
                      datByGene(x) <- datByGene(x)[igenes,,drop=FALSE]
                  }
                  else {
                      switch (EXPR =i ,
                              "dat" ={ return ( assayDataElement(x,"dat") )} ,
                              "datByGene" ={ if (nrow(x@datByGene)!=0) return ( x@datByGene ) else return(assayDataElement(x,"dat"))} ,
                              "pData" ={ return ( pData(x) )} ,
                              "featureData" ={ return ( featureData(x) )} ,
                              "annotation" ={ return ( annotation(x) )} ,
                              "package" ={ return ( annotation(x) )} ,
                              "A" ={ return ( x@A )} ,
                              "S" ={ return ( x@S )} ,                   
                              "SByGene" ={ return ( x@SByGene )} ,
                              "compNames" ={ return ( x@compNames )} ,
                              "indComp" ={ return ( x@indComp )} ,
                              "witGenes" ={ return ( x@witGenes )} ,
                              "refSamples" ={ return ( x@refSamples)} ,
                              "typeID" ={ return ( x@typeID)} ,
                              "organism" ={ return ( x@organism)} ,
                              "mart" ={ return ( x@mart)} ,
                              stop ( "This attribute is not valid! " )
                              )
                  }
              }
          }
      }

          if (!missing(j)) {

              if (is.numeric(j))
                  j <- sampleNames(x)[j]
          
                  keepSamples <- j
                  diffSamples <- setdiff(keepSamples,sampleNames(x))
                  if (length(diffSamples) == length(keepSamples))
                      stop("The sample ids are not available in the object.")
                  else if (length(diffSamples) > 0)
                      warning(paste("The samples:",paste(diffSamples,collapse=", "),"are not available in the object."))
                  keepSamples <- intersect(keepSamples,sampleNames(x)) 
                  
                  x@A <- A(x)[keepSamples,,drop=FALSE]
                  varLabelsor <- varLabels(x)
                  pData(x) <- pData(x)[keepSamples,,drop=FALSE]
                  varLabels(x) <-  varLabelsor
                  pData(protocolData(x))<- pData(protocolData(x))[keepSamples,,drop=FALSE]
                  assayDataElement(x,"dat") <- assayDataElement(x,"dat")[,keepSamples,drop=FALSE]
                  x@datByGene <- datByGene(x)[,keepSamples, drop=FALSE]
                  refSamples(x) <- intersect(refSamples(x),keepSamples)
                  

          }

          if (!missing(k)) {
              keepComp <- k
              
              if (is.character(keepComp)) {
                  diffComp <- setdiff(keepComp,compNames(x))
                  if (length(diffComp) == length(keepComp))
                      stop("The component ids are not available in the object.")
                  else if (length(diffComp) > 0)
                      warning(paste("The components:",paste(diffComp,collapse=", "),"are not available in the object."))
                  keepComp <- match(keepComp,compNames(x))
              }

              keepComp <- intersect(keepComp,1:nbComp(x))
            
              diffcomp <- setdiff(keepComp,indComp(x))
              if (length(diffcomp) == length(keepComp))
                  stop("The component numbers provided in 'keepComp' are not included in 'indComp(x)'.")
              else if (length(diffcomp) > 0)
                  warning(paste("The component(s):",paste(diffcomp,collapse=", ",sep=""),"are not available in 'x'."))
              
              keepComp <- intersect(keepComp,indComp(x))
              keepCompor <- keepComp
              keepComp <- match(keepComp, indComp(x))
              
              x@compNames <- compNames(x)[keepComp]
              
              if (length(witGenes(x))>0)
                  x@witGenes <- witGenes(x)[keepComp]
              
              x@A <- A(x)[,keepComp, drop=FALSE]
              x@S <- S(x)[,keepComp, drop=FALSE]
              x@SByGene <- SByGene(x)[,keepComp, drop=FALSE]
                
              indComp(x) <- keepCompor
          }
          
          return(x)
              
     }
)



setMethod( "show" ,"IcaSet" ,
       function (object){
           cat("Number of components:",nbComp(object),"\n")
           cat("Component labels:",compNames(object),"\n")
           callNextMethod()
           
      }
)




setMethod( "selectContrib" ,signature("IcaSet", "numeric", "character") ,
       function (object, cutoff=3, level=c("features","genes")){

           
           level <- level[1]
           level <- match.arg(level)
           switch(level,
                  features={Sl=Slist(object)},
                  genes={Sl=SlistByGene(object)}
                  )

           if (length(cutoff)==1)
               cutoff <- rep(cutoff[1],nbComp(object))
           else
               if ((length(cutoff)>1) & (length(cutoff)<nbComp(object))) {
                   cutoff <- rep(cutoff[1],nbComp(object))
                   warning("The length of arg 'cutoff' or 'selCutoff' must be either 1 or equal to the number of components")    
               }

           
           sel <- foreach (comp.proj=Sl, cutt=cutoff) %do% {
               sources.zval <- (comp.proj-mean(comp.proj))/sd(comp.proj)
               sources.keep <- sources.zval[abs(sources.zval) >= cutt]
               return(sources.keep)
           }

           names(sel) <- compNames(object)

           return(sel)
       }
          )


setMethod( "selectContrib" ,signature("list", "numeric") ,
       function (object, cutoff=3, ...){

           if (is.numeric(object))
           {
               comp <- object
               sel <- (comp-mean(comp))/sd(comp)
               sel <- sel[abs(sel) >= cutoff[1]]
               
           }
           else {

           if (length(cutoff)==1)
               cutoff <- rep(cutoff[1],length(object))
           else
               if ((length(cutoff)>1) & (length(cutoff)<length(object))) {
                   cutoff <- rep(cutoff[1],length(object))
                   warning("The length of arg 'cutoff' or 'selCutoff' must be either 1 or equal to the number of components")    
               }
                   

               comp.proj <- cutt <- NULL
               sel <- foreach (comp.proj=object, cutt=cutoff) %do% {
                   sources.zval <- (comp.proj-mean(comp.proj))/sd(comp.proj)
                   sources.keep <- sources.zval[abs(sources.zval) >= cutt]
                   return(sources.keep)
               }
               names(sel) <- names(object)
           }
           

           return(sel)
       }
          )


setMethod( "getComp" ,signature("IcaSet","character","numeric") ,
       function (object, level=c("features","genes"), ind){
           ac <- Alist(object)
           level <- level[1]
           callf <- c("features"="Slist","genes"="SlistByGene")[level]
           sc <- do.call(callf,list(object))
           
           if (length(ind)==1) {
               sc <- sc[[ind]]
               ac <- ac[[ind]]
           }
           else {
               sc <- sc[ind]
               ac <- ac[ind]
           }
           
           return (list(contrib=ac,proj=sc))
      }
)




setMethod( "geneNames" ,"IcaSet" ,
       function (object){
               return(rownames(datByGene(object)))
      }
)

setMethod( "nbComp" ,"IcaSet" ,
       function (object){
               return(ncol(A(object)))
      }
)


setMethod( "datByGene" ,"IcaSet" ,
          function (object){
           if (nrow(object@datByGene)>0)
               return ( object@datByGene)
           else
               return(object@assayData[["dat"]])

          }

)


setMethod(f = "A",
          signature = "IcaSet" ,
          definition = function (object){
              return (object@A)
          }
          )

setMethod(f = "dat",
          signature = "IcaSet" ,
          definition = function (object){
              return (assayDataElement(object,"dat") )
          }
          )


setMethod( "getA" ,"IcaSet" ,
       function (object){
           return (object@A)
      }
)


setMethod( "getS" ,"IcaSet" ,
       function (object){
           return (object@S)
      }
)

setMethod( "S" ,"IcaSet" ,
       function (object){
           return (object@S)
      }
)


setMethod( "Slist" ,"IcaSet" ,
       function (object){
           S <- S(object)
           Slist <- lapply(as.list(S),
                           function(x,n) {
                               names(x) <- n
                               return(x)
                           },
                           n = rownames(S)
                           )
           return(Slist)

      }
)

setMethod( "Alist" ,"IcaSet" ,
       function (object){
           A <- A(object)
           Alist <- lapply(as.list(A),
                           function(x,n) {
                               names(x) <- n
                               return(x)
                           },
                           n = rownames(A)
                           )
           return(Alist)

      }
)



setMethod( "getSByGene" ,"IcaSet" ,
       function (object){
           return(object@SByGene)
      }
)
setMethod( "SByGene" ,"IcaSet" ,
       function (object){
           return(object@SByGene)
      }
)

setMethod( "SlistByGene" ,"IcaSet" ,
       function (object){
           S <- SByGene(object)
           Slist <- lapply(as.list(S),
                           function(x,n) {
                               names(x) <- n
                               return(x)
                           },
                           n = rownames(S)
                           )
           return(Slist)
      }
)

setMethod( "compNames" ,"IcaSet" ,
       function (object){
           return (object@compNames)
      }
)

setMethod( "getLabelsComp" ,"IcaSet" ,
       function (object){
           return (object@compNames)
      }
)

setMethod( "indComp" ,"IcaSet" ,
       function (object){
           return (object@indComp)
      }
)

setMethod( "getIndComp" ,"IcaSet" ,
       function (object){
           return (object@indComp)
      }
)
        
 setMethod( "witGenes" ,"IcaSet" ,
       function (object){
           return (object@witGenes)
      }
)

 setMethod( "getWitGenes" ,"IcaSet" ,
       function (object){
           return (object@witGenes)
      }
)
        

setMethod( "getPackage" ,"IcaSet" ,
       function (object){
           return (object@annotation)
      }
)

setMethod( "package" ,"IcaSet" ,
       function (object){
           return (object@annotation)
      }
)


setMethod( "getChipVersion" ,"IcaSet" ,
       function (object){
           return (object@chipVersion)
      }
)

setMethod( "chipVersion" ,"IcaSet" ,
       function (object){
           return (object@chipVersion)
      }
)

setMethod( "getChipManu" ,"IcaSet" ,
       function (object){
           return (object@chipManu)
      }
)

setMethod( "chipManu" ,"IcaSet" ,
       function (object){
           return (object@chipManu)
      }
)

setMethod( "refSamples" ,"IcaSet" ,
       function (object){
           return (object@refSamples)
      }
)

setMethod( "getRefSamples" ,"IcaSet" ,
       function (object){
           return (object@refSamples)
      }
)

setMethod( "typeID" ,"IcaSet" ,
       function (object){
           return (object@typeID)
      }
)

setMethod( "getTypeID" ,"IcaSet" ,
       function (object){
           return (object@typeID)
      }
)
setMethod( "organism" ,"IcaSet" ,
       function (object){
           return (object@organism)
      }
)

setMethod( "mart" ,"IcaSet" ,
       function (object){
           return (object@mart)
      }
)

setMethod( "getMart" ,"IcaSet" ,
       function (object){
           return (object@mart)
      }
)


setReplaceMethod(f ="[",
                 signature = signature(x="IcaSet", i="ANY", j="ANY", value="ANY"),
                 definition = function (x, i, j, value ){
                     switch ( EXPR =i ,
                             "dat" ={ assayDataElement(x,"dat") <- value} ,
                             "datByGene" ={ x@datByGene <- value} ,
                             "pData" ={pData(x) <- value} ,
                             "featureData" ={ fData(x) <- value} ,
                             "fData" ={ fData(x) <- value} ,
                             "annotation" ={ x@annotation <- value} ,
                             "package" ={ x@annotation <- value} ,
                             "A" ={ x@A <- value} ,
                             "S" ={ x@S <- value} ,                             
                             "SByGene" = {x@SByGene <- value},
                             "compNames" = {x@compNames <- value},
                             "indComp" = {x@indComp <- value},
                             "witGenes" ={ x@witGenes <- value} ,
                             "chipManu" ={ x@chipManu <- value} ,
                             "chipVersion" ={ x@chipVersion <- value} ,
                             "refSamples" ={ x@refSamples <- value} ,
                             "typeID" ={ x@typeID <- value} ,
                             "organism" ={ x@organism <- value} ,
                             "mart" ={ x@mart <- value} ,
                             stop ( " This attribute doesn't exist " )
                             )
                     validObject ( x )
                     return ( x )
                 }
                 )



setReplaceMethod(
      f = "datByGene" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@datByGene <- value
        validObject (object)
        return (object)
     }
)

setReplaceMethod(f = "sampleNames" ,
                 signature = "IcaSet" ,
                 definition = function (object, value){
                     if (nrow(A(object))>0 & ncol(A(object))>0) {
                         newA <- A(object)
                         rownames(newA) <- value
                         object@A <- newA
                     }
                     if (nrow(datByGene(object))>0 & ncol(datByGene(object))>0){
                         newd <- datByGene(object)
                         colnames(newd) <- value
                         object@datByGene <- newd
                     }
                    callNextMethod()

                      
                 }
)


setReplaceMethod(
      f = "setA" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@A <- value
        validObject(object)
        return (object)
     }
)

setReplaceMethod(
      f = "A" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@A <- value
        validObject(object)
        return (object)
     }
)

setReplaceMethod(
      f = "dat" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        assayDataElement(object,"dat") <- value
        validObject(object)
        return (object)
     }
)


setReplaceMethod(
      f = "setS" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@S <- value
        return (object)
     }
)

setReplaceMethod(
      f = "S" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@S <- value
        validObject(object)
        return (object)
     }
)


setReplaceMethod(
      f = "setSByGene" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@SByGene <- value
        validObject(object)
        return (object)
     }
)



setReplaceMethod(
      f = "SByGene" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@SByGene <- value
        validObject(object)
        return (object)
     }
)

setReplaceMethod(
      f = "setLabelsComp" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@compNames <- value
        validObject(object)
        return (object)
     }
)

setReplaceMethod(
      f = "compNames" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@compNames <- value
        validObject(object)
        return (object)
     }
)

setReplaceMethod(
      f = "setIndComp" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@indComp <- value
        validObject(object)
        return (object)
     }
)

setReplaceMethod(
      f = "indComp" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@indComp <- value
        validObject(object)
        return (object)
     }
)


setReplaceMethod(
      f = "setWitGenes" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@witGenes <- value
        validObject (object)
        return (object)
     }
)

setReplaceMethod(
      f = "witGenes" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@witGenes <- value
        validObject (object)
        return (object)
     }
)


setReplaceMethod(
      f = "chipVersion" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@chipVersion <- value
        validObject (object)
        return (object)
     }
)

setReplaceMethod(
      f = "setChipVersion" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@chipVersion <- value
        validObject (object)
        return (object)
     }
)

setReplaceMethod(
      f = "chipManu" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@chipManu <- value
        validObject (object)
        return (object)
     }
)

setReplaceMethod(
      f = "setChipManu" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@chipManu <- value
        validObject (object)
        return (object)
     }
)


setReplaceMethod( f="refSamples",
      signature = "IcaSet" ,
       definition = function (object, value){
           object@refSamples <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="setRefSamples",
      signature = "IcaSet" ,
       definition = function (object, value){
           object@refSamples <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="typeID",
      signature = "IcaSet" ,
       definition = function (object, value){
           object@typeID <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="setTypeID",
      signature = "IcaSet" ,
       definition = function (object, value){
           object@typeID <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="organism",
      signature = "IcaSet" ,
       definition = function (object, value){
           object@organism <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="setMart",
      signature = "IcaSet" ,
       definition = function (object, value){
           object@mart <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="mart",
      signature = "IcaSet" ,
       definition = function (object, value){
           object@mart <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod(
      f = "package" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@annotation <- value
        validObject (object)
        return (object)
     }
)

setReplaceMethod(
      f = "setPackage" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@annotation <- value
        validObject (object)
        return (object)
     }
)

setReplaceMethod(
      f = "setAnnotation" ,
      signature = "IcaSet" ,
      definition = function (object, value){
        object@annotation <- value
        validObject (object)
        return (object)
     }
)
