
# Extract attributes of an object of the class MineICAParams
#
# @name [
# @aliases [,MineICAParams-method
# @docType methods
# @rdname extract-methods
#
setMethod(
      f ="[",
      signature(x = "MineICAParams", i = "ANY", j="ANY"),
      definition = function (x, i, j, ..., drop){
              switch ( EXPR =i ,
                      "Afile" ={ return ( x@Afile )} ,
                      "Sfile" ={ return ( x@Sfile )} ,
                      "datfile" ={ return ( x@datfile )} ,
                      "annotfile" ={ return ( x@annotfile )} ,
                      "resPath" ={ return ( x@resPath )} ,
                      "genesPath" ={ return ( x@genesPath )} ,
                      "annot2col" ={ return ( x@annot2col )} ,
                      "selCutoff" ={ return ( x@selCutoff )} ,
                      "pvalCutoff" ={ return ( x@pvalCutoff )} ,
                      stop ( "This attribute is not valid! " )
             )
     }
)





setMethod( "Afile" ,"MineICAParams" ,
       function (object){
           return (object@Afile)
      }
)

setMethod( "getAfile" ,"MineICAParams" ,
       function (object){
           return (object@Afile)
      }
)


setMethod( "Sfile" ,"MineICAParams" ,
       function (object){
           return (object@Sfile)
      }
)

setMethod( "getSfile" ,"MineICAParams" ,
       function (object){
           return (object@Sfile)
      }
)

setMethod( "datfile" ,"MineICAParams" ,
       function (object){
           return (object@datfile)
      }
)

setMethod( "getDatfile" ,"MineICAParams" ,
       function (object){
           return (object@datfile)
      }
)

setMethod( "annotfile" ,"MineICAParams" ,
       function (object){
           return (object@annotfile)
      }
)

setMethod( "getAnnotfile" ,"MineICAParams" ,
       function (object){
           return (object@annotfile)
      }
)

setMethod( "genesPath" ,"MineICAParams" ,
       function (object){
           return (object@genesPath)
      }
)

setMethod( "getGenesPath" ,"MineICAParams" ,
       function (object){
           return (object@genesPath)
      }
)

setMethod( "resPath" ,"MineICAParams" ,
       function (object){
           return (object@resPath)
      }
)

setMethod( "getResPath" ,"MineICAParams" ,
       function (object){
           return (object@resPath)
      }
)



setMethod( "getAnnot2col" ,"MineICAParams" ,
       function (object){
           return (object@annot2col)
      }
)

setMethod( "annot2col" ,"MineICAParams" ,
       function (object){
           return (object@annot2col)
      }
)

setMethod( "pvalCutoff" ,"MineICAParams" ,
       function (object){
           return (object@pvalCutoff)
      }
)

setMethod( "getPvalCutoff" ,"MineICAParams" ,
       function (object){
           return (object@pvalCutoff)
      }
)

setMethod( "selCutoff" ,"MineICAParams" ,
       function (object){
           return (object@selCutoff)
      }
)

setMethod( "getSelCutoff" ,"MineICAParams" ,
       function (object){
           return (object@selCutoff)
      }
)


## Replace attributes of an object of the class MineICAParams
#
## @name [
## @aliases [<-,MineICAParams-method
## @docType methods
## @rdname extract-methods
setReplaceMethod(
      f ="[",
      signature = "MineICAParams" ,
      definition = function (x, i, j, value){
              switch ( EXPR =i ,
                      "Afile" ={ x@Afile <- value} ,
                      "Sfile" ={ x@Sfile <- value} ,
                      "datfile" ={ x@datfile <- value} ,
                      "annotfile" ={ x@annotfile <- value} ,
                      "resPath" ={ x@resPath <- value;
                                   home <- system("echo $HOME",intern=TRUE);
                                   value <- gsub(value,pattern="~",replacement=home);
                                   system(paste("mkdir",value)) } ,
                      "genesPath" ={ x@genesPath <- value},#; system(paste("mkdir",value)) } ,
                      "annot2col" ={ x@annot2col <- value} ,
                      "selCutoff" ={ x@selCutoff <- value} ,
                      "pvalCutoff" ={ x@pvalCutoff <- value} ,
                      stop ( "This attribute desn't exist" )
             )
             validObject ( x )
             return ( x )
              
     }
)


setReplaceMethod(
      f = "Afile" ,
      signature = "MineICAParams" ,
      definition = function (object, value){
           object@Afile <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( 
      f = "setAfile" ,
      signature = "MineICAParams" ,
       definition = function (object, value){
           object@Afile <- value           
           validObject (object)
           return (object)
      }
)


setReplaceMethod( f = "Sfile",
      signature = "MineICAParams" ,
       definition = function (object, value){
           object@Sfile <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f = "setSfile" ,
      signature = "MineICAParams" ,
       definition = function (object, value){
           object@Sfile <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="datfile" ,
      signature = "MineICAParams" ,
       definition = function (object, value){
           object@datfile <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="setDatfile" ,
      signature = "MineICAParams" ,
       definition = function (object, value){
           object@datfile <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="annotfile" ,
      signature = "MineICAParams" ,
       definition = function (object, value){
           object@annotfile <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="setAnnotfile" ,
      signature = "MineICAParams" ,
       definition = function (object, value){
           object@annotfile <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="resPath" ,
       definition = function (object, value){
           home <- system("echo $HOME",intern=TRUE)
           value <- gsub(value,pattern="~",replacement=home)
            object@resPath <- value
            validObject (object)
            system(paste("mkdir",value))
           return (object)
      }
)

setReplaceMethod( f="setGenesPath" ,
       definition = function (object, value){
            object@genesPath <- value
            validObject (object)
            #system(paste("mkdir",value)) 
           return (object)
      }
)

setReplaceMethod( f="genesPath" ,
       definition = function (object, value){
            object@genesPath <- value
            validObject (object)
           return (object)
      }
)

setReplaceMethod( f="setResPath" ,
       definition = function (object, value){
           home <- system("echo $HOME",intern=TRUE)
           value <- gsub(value,pattern="~",replacement=home)
            object@resPath <- value
            validObject (object)
           return (object)
      }
)


setReplaceMethod( f="setAnnot2col" ,
      signature = "MineICAParams" ,
       definition = function (object, value){
           object@annot2col <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="annot2col" ,
       signature = "MineICAParams" ,
       definition =
          function (object, value){
           object@annot2col <- value
           validObject (object)
           return (object)
          }
)

setReplaceMethod( f="pvalCutoff" ,
      signature = "MineICAParams" ,
       definition = function (object, value){
           object@pvalCutoff <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="setPvalCutoff" ,
      signature = "MineICAParams" ,
       definition = function (object, value){
           object@pvalCutoff <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="selCutoff" ,
      signature = "MineICAParams" ,
       definition = function (object, value){
           object@selCutoff <- value
           validObject (object)
           return (object)
      }
)

setReplaceMethod( f="setSelCutoff" ,
      signature = "MineICAParams" ,
       definition = function (object, value){
           object@selCutoff <- value
           validObject (object)
           return (object)
      }
)

