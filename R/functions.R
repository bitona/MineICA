##' Object of class \code{\link{IcaSet}} containing an ICA decomposition calculated by the FastICA algorithm
##' (through matlab function "icasso") on
##' bladder cancer expression data measured on HG-U133A Affymetrix microarrays.
##' The original expression data  were normalized with MAS5 by the authors of the paper
##' followed by log2-transformation.
##' ICA was run on the dataset restricted to the 10000 most variable probe sets (based on IQR values).
##' 10 components were computed.
##' Only probe sets/genes having an absolute projection higher than 3 are stored in this object.
##' @title IcaSet-object containing a FastICA decomposition of gene expression microarrray-based data of bladder cancer samples.
##' @name icaSetCarbayo
##' @docType data
##' @author Anne Biton
##' @references \url{http://jco.ascopubs.org/content/24/5/778/suppl/DC1}
##' @keywords datasets
NULL

##' readA
##' 
##' This function reads and annotates matrix A.
##' 
##' The matrix dat must be the one on which the matrix A was calculated.
##' It is assumed that the number of components is lower than the number of
##' samples, the matrix will be transposed to have dimension 'samples x components' according to this assumption. If \code{annot} is FALSE, colnames of dat are used to annotate
##' rownames of A.
##' @title read A
##' @param Afile The file which contains the matrix of sample contributions. It
##' must be a txt file where the separator is \code{white space}, that is one or
##' more spaces, tabs, newlines or carriage returns
##' @param datfile The file which contains the matrix (of dimension features x
##' samples) based on which the matrix A was calculated
##' @param dat The data based on which the matrix A was calculated (features x
##' samples)
##' @param annot TRUE (default) if the Afile contains rownames of matrix A,
##' FALSE if the rownames has to be extracted from dat
##' @return This function returns a matrix of dimension samples x components
##' with rownames filled with sample IDs.
##' @keywords internal
##' @author Anne Biton
readA <- function (
                   Afile
                   ,datfile
                   ,dat
                   ,annot = TRUE
                   )
{

	message("..Reading matrix A.. ")
	A <- read.table(Afile, header = if (!annot) TRUE else FALSE) #sep = sep, 

	nr <- nrow(A)
	nc <- ncol(A)
	if (nr < nc)
            A <- t(A) 
	A <- as.data.frame(A, check.names = FALSE, stringsAsFactors = FALSE)

        if (annot) {
            
            if (missing(dat)) {
                if (missing(datfile))
                    stop("Either the expression matrix or the corresponding file must be provided.")
                dat <- read.table(datfile, header = TRUE,  nrows = 1, stringsAsFactors= FALSE, check.names = FALSE)
                if (ncol(dat) != nrow(A))
                    stop("The  number of rows of A must equal the number of columns of dat.")
            }
        
            rownames(A) = colnames(dat)
        }
        
	gc()
	return(A)
}

##' This function reads and annotates matrix S.
##' 
##' The matrix dat must be the one on which the matrix S was calculated.
##' It is assumed that the number of components is lower than the number of features, the matrix will be transposed to have dimension 'features x components' according to this assumption. 
##' If \code{annot} is FALSE, rownames of dat are used to annotate rownames of S.
##' @title read S
##' @param Sfile The file which contains the matrix of feature projections. It
##' must be a txt file where the separator is \code{white space}, that is one or
##' more spaces, tabs, newlines or carriage returns.
##' @param datfile The file which contains the matrix (of dimension features x
##' samples) based on which the matrix S was calculated. It must be a txt file
##' where the separator is \code{white space}, that is one or more spaces, tabs,
##' newlines or carriage returns.
##' @param dat The data based on which the matrix A was calculated (features x
##' samples)
##' @param annot TRUE (default) if the Afile contains rownames of matrix A,
##' FALSE if the rownames has to be extracted from dat
##' @return This function returns a matrix of dimension features x components
##' with rownames filled with feature IDs.
##' @keywords internal
##' @author Anne Biton
readS <- function (Sfile,
                   datfile,
                   dat,
                   annot = TRUE)
{

   message("..Reading matrix S..")
	S <- read.table(Sfile, header = if (!annot) TRUE else FALSE)# sep = sep, 
	nr <- nrow(S)
	nc <- ncol(S)
        if (nr < nc)
            S <- t(S) 
        S <- as.data.frame(S, check.names = FALSE, stringsAsFactors = FALSE)


    if (annot) {
    
        message("...Put feature labels on S matrix...")
	if (missing(dat)) {
            if (missing(datfile))
                stop("Either the data matrix or the corresponding file must be provided.")
            dat <- read.table(datfile, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
        }
	rownames(S) <- rownames(dat)
     }

	gc()
	return(S)
}



## 
## Builds a subset of an IcaSet object
##
## \code{keepSamples} must be available in \code{sampleNames(icaSet)}. 
## \code{keepComp} must be available in indComp(icaSet), i.e if \code{icaSet} has been previously restricted to a few components by subIcaSet, \code{keepComp} must refer to the index of the components in the original icaSet object.
## @title subset an \code{\link{IcaSet}} object
## @param icaSet An IcaSet object
## @param keepSamples The sample IDs to be selected, must be included in \code{sampleNames(icaSet)}. If missing, all samples are kept.
## @param keepComp The component indices to be selected, must  included in \code{indComp(icaSet)}. If missing, all components are kept.
## @return This function returns an \code{IcaSet} object.
## @author Anne Biton
## @keywords internal
## @seealso \code{\link{IcaSet-class}}
## @examples 
## ## load an example of IcaSet
## data(icaSetCarbayo)
## 
## #keep a subset of samples (samples of stage T2, T3, T4)
## invSamples <- rownames(subset(pData(icaSetCarbayo),stage %in% c("T2","T3","T4")))
## icaSetCarbayo_subSamples <- subIcaSet(icaSetCarbayo, keepSamples=invSamples)
## 
## #keep a subset of components (components of index 5:9 and 11,12)
## icaSetCarbayo_subComp <- subIcaSet(icaSetCarbayo, keepComp=c(5:9,11,12))
## 
subIcaSet <- function(
                      icaSet
                      ,keepSamples
                      ,keepComp
                      ) {

    if (!missing(keepSamples)) {
        diffSamples <- setdiff(keepSamples,sampleNames(icaSet))
        if (length(diffSamples) == length(keepSamples))
            stop("The sample ids provided in 'keepSamples' are not available in 'icaSet'.")
        else if (length(diffSamples) > 0)
            warning(paste("The samples:",paste(diffSamples,collapse=", "),"are not available in icaSet."))
        keepSamples <- intersect(keepSamples,sampleNames(icaSet)) 
             
        icaSet@A <- A(icaSet)[keepSamples,]
        varLabelsor <- varLabels(icaSet)
        pData(icaSet) <- pData(icaSet)[keepSamples,,drop=FALSE]
        varLabels(icaSet) <-  varLabelsor
        pData(protocolData(icaSet))<- pData(protocolData(icaSet))[keepSamples,,drop=FALSE]
        assayDataElement(icaSet,"dat") <- assayDataElement(icaSet,"dat")[,keepSamples,drop=FALSE]
        icaSet@datByGene <- datByGene(icaSet)[,keepSamples,drop=FALSE]
        refSamples(icaSet) <- intersect(refSamples(icaSet),keepSamples)
                    
    }
    if (!missing(keepComp)) {
        if (!is.numeric(keepComp))
            stop ("keepComp must be numeric")
        nbcomp <- nbComp(icaSet)
        if (sum(sort(keepComp) %in% indComp(icaSet)) != nbcomp) {
            diffcomp <- setdiff(keepComp,indComp(icaSet))
            if (length(diffcomp) == length(keepComp))
                stop("The component numbers provided in 'keepComp' are not included in 'indComp(icaSet)'.")
            else if (length(diffcomp) > 0)
                warning(paste("The component(s):",paste(diffcomp,collapse=", ",sep=""),"are not available in 'icaSet'."))

            keepComp <- intersect(keepComp,indComp(icaSet))
            keepCompor <- keepComp
            keepComp <- match(keepComp, indComp(icaSet))

            icaSet@compNames <- compNames(icaSet)[keepComp]

            if (length(witGenes(icaSet))>0)
                icaSet@witGenes <- witGenes(icaSet)[keepComp]


            icaSet@A <- A(icaSet)[,keepComp, drop=FALSE]
            icaSet@S <- S(icaSet)[,keepComp, drop=FALSE]
            icaSet@SByGene <- SByGene(icaSet)[,keepComp, drop=FALSE]
            
            indComp(icaSet) <- keepCompor
        }
        
    }
    validObject(icaSet)
    return(icaSet)
}

##' Extract projection values of a given set of IDs on a subset of components.
##'
##' 
##' @title Extract projection values
##' @param icaSet An object of class \code{\link{IcaSet}}
##' @param ids feature or gene IDs
##' @param keepComp Index of the components to be conserved, must be in \code{indComp(icaSet)}
##' @param level The level of projections to be extracted, either \code{"features"} or \code{"genes"}
##' @return A vector or a list of projection values
##' @author Anne Biton
##' @export 
##' @examples 
##'
##' ## load an example of IcaSet
##' data(icaSetCarbayo)
##' 
##' ##get the projection of your favorite proliferation genes
##' #on all components
##' getProj(icaSetCarbayo, ids=c("TOP2A","CDK1","CDC20"), level="genes")
##' 
##' #on some components 
##' getProj(icaSetCarbayo, ids=c("TOP2A","CDK1","CDC20"),
##' keepComp=c(1,6,9,12),level="genes")
##' 
##' ##get the gene projection values on the sixth component
##' getProj(icaSetCarbayo, keepComp=6,level="genes")
##' 
getProj <- function(icaSet,ids,keepComp, level = c("features","genes")) {

    level <- match.arg(tolower(level), choices=c("features","genes"))

    switch(level,
           features={ data=S(icaSet)},
           genes={data=SByGene(icaSet)}
           )
    
    if (!missing(ids)) {
        interG <- intersect(rownames(data),ids)
    
        if (length(interG) == 0)
            stop("The given ids are not available in icaSet")
        else if (length(interG) < length(ids)) {
            diff <- setdiff(ids,rownames(data))
            warning(paste("The IDs",paste(diff,collapse = ", "),"are not avaible in icaSet"))
        }
    }
    else interG <- rownames(data)
    
    if (!missing(keepComp)) {
        diffcomp <- setdiff(keepComp,indComp(icaSet))
        if (length(diffcomp) == length(keepComp))
            stop("The component numbers provided in 'keepComp' are not included in 'indComp(icaSet)'.")
        else if (length(diffcomp) > 0)
            warning(paste("The component(s):",paste(diffcomp,collapse=", ",sep=""),"are not available in 'icaSet'."))
        keepComp <- intersect(keepComp,indComp(icaSet))
    }
    else
        keepComp <- indComp(icaSet)
    
        out <- data[interG,keepComp]
    
    if (length(keepComp)==1) 
            names(out) <- interG
    
    return(out)
}



##' This function annotates the features of an \code{\link{IcaSet}} object and fills its attributes \code{SByGene} and \code{datByGene}.
##'
##' When attribute \code{annotation} of \code{icaSet} is not specified (of length \code{0}), \code{biomaRt} is used to annotate the features through function \code{\link{annotFeaturesWithBiomaRt}}.
##' 
##' When specified, attribute \code{annotation} of argument \code{icaSet} must be an annotation package and will be used to annotate the \code{featureNames} of \code{icaSet}. In addition, the attribute \code{typeID} (a vector) of argument \code{icaSet} must contain a valid element  \code{geneID_annotation} that determines the object of the package to be used for the annotation, see \code{\link{IcaSet}}.
##' 
##' When argument \code{annot} is TRUE, this function fills the attributes \code{SByGene} and \code{datByGene} of \code{icaSet}. When several feature IDs are available for a same gene ID, the median value of the corresponding features IDs is attributed to the gene (the median of the projection values is used for attribute \code{SByGene}, and the median of the expression values is used for attribute \code{datByGene}).
##' 
##' When attribute \code{chipManu} of the argument \code{icaSet} is "illumina", the features are first converted into nuID using the package 'lumi*Mapping' and then annotated into genes.
##' In that case, features can only be annotated in ENTREZID or SYMBOL. It means that \code{typeID(icaSet)['geneID_annotation']} must be either  'ENTREZID' or 'SYMBOL'. You will need to annotate yourself the IcaSet object if you want to use different IDs.
##' 
##' @title Features annotation of an object of class IcaSet. 
##' @param icaSet An object of class \code{\link{IcaSet}} to be annotated, must contain a valid \code{annotation} attribute.
##' @param params An object of class \code{\link{MineICAParams}} containing the parameters of the analysis.
##' @param annot TRUE (default) if the IcaSet object must indeed be annotated
##' @return The modified argument \code{icaSet}, with filled attributes \code{SByGene} and \code{datByGene}.
##' @seealso \code{\link{annotFeaturesComp}}
##' @author Anne Biton
##' @export
##' @examples
##' #load data
##' data(icaSetCarbayo)
##' require(hgu133a.db)
##' 
##' # run annotation of the features into gene Symbols as specified in 'typeID(icaSetCarbayo)["geneID_annotation"]',
##' # using package hgu133a.db as defined in 'annotation(icaSetMainz)' 
##' icaSetCarbayo <- annotInGene(icaSet=icaSetCarbayo, params=buildMineICAParams())
##' 
##' \dontrun{
##' #load data
##' library(breastCancerMAINZ)
##' data(mainz)
##' #run ICA
##' resJade <- runICA(X=exprs(mainz), nbComp=5, method = "JADE", maxit=10000) 
##' 
##' #build params
##' params <- buildMineICAParams(resPath="mainz/")
##' 
##' #build a new IcaSet object, omitting annotation of the features (runAnnot=FALSE)
##' #but specifying the element "geneID_annotation" of argument 'typeID'
##' icaSetMainz <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S),
##'                              dat=exprs(mainz), pData=pData(mainz),
##'                              annotation="hgu133a.db", typeID= c(geneID_annotation = "SYMBOL",
##'                              geneID_biomart = "hgnc_symbol", featureID_biomart = "affy_hg_u133a"),
##'                              chipManu = "affymetrix", runAnnot=FALSE,
##'                              mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl"))
##'
##' #Attributes SByGene is empty and attribute datByGene refers to assayData 
##' SByGene(icaSetMainz)
##' head(datByGene(icaSetMainz))
##' 
##' # run annotation of the features into gene Symbols as specified in 'typeID(icaSetMainz)["geneID_annotation"]',
##' # using package hgu133a.db as defined in 'annotation(icaSetMainz)' 
##' icaSetMainz <- annotInGene(icaSet=icaSetMainz, params=params)
##' }
annotInGene <- function(icaSet,
                        params,
                        annot = TRUE
                        ) { 
   
        chip <- annotation(icaSet)

        if (annot) {
            if (length(chip)==0 || (chip == ""))
                message("No annotation package is available, biomaRt is used for the annotation.")
            else {
                if (!("geneID_annotation" %in% names(typeID(icaSet))))
                    stop(paste("Since you want to annotate the features using ", chip,", typeID attribute of IcaSet object must contain a 'geneID_annotation' element which determines the object from the package you want to use.", sep = ""))
                else {
                    listIDs <- gsub(x=ls(paste("package:",chip,sep="")), pattern = gsub(chip,pattern=".db",replacement = ""), replacement = "")
                    if (!(typeID(icaSet)["geneID_annotation"] %in% listIDs))
                        stop(paste("The element 'geneID_annotation' of attribute 'typeID' of object IcaSet is not available in annotation package,",chip))

                    library(chip, character.only = TRUE)
                }
            }
            icaSet <- annotFeaturesComp(icaSet = icaSet, params=params, type = toupper(typeID(icaSet)["geneID_annotation"]))

        }
        

        return(icaSet)
                        
}

##' This function annotates a set of features
##'
##' @title Annotation of features using an annotation package
##' @param features Feature IDs to be annotated
##' @param type The object from the package used to annotate the features, must be available in \code{ls("package:package_name")}
##' @param annotation An annotation package
##' @return A vector of gene/object IDs indexed by the feature IDs.
##' @author Anne Biton
##' @export
##' @examples
##' library(hgu133a.db)
##' annotFeatures(features = c("1007_s_at", "1053_at", "117_at", "121_at", "1255_g_at"),
##'               type="SYMBOL", annotation="hgu133a.db")
##' 
annotFeatures <- function(features,
                          type,
                          annotation
                          ) {
    
        pack.annot <- eval(as.name(paste(gsub(".db", "", annotation), type, sep = "")))
        probeToGene = unlist(AnnotationDbi::mget(features, pack.annot, ifnotfound = NA))
        names(probeToGene) = features
        probeToGene = na.omit(probeToGene)
        message(paste("The number of probe sets not associated with a",type,"ID is",length(features)-length(probeToGene)))
	message(paste("The remaining number of probe sets is",length(probeToGene)))

	
	return(probeToGene)
}


## ## annotation of an IcaSet object (icaSetKim) containing data from an illumina microarray
## # chipManu must be available and equal to "illumina"
## chipManu(icaSetKim)
## # annotation must be available and equal to "lumi*All.db" where * is the organism, in this case 'Human'.
## annotation(icaSetKim)
## icaSetKim_annot <- annotFeaturesComp(icaSet=icaSetKim, params=params, type="SYMBOL")

##' ##' This function annotates the features of an object of class \code{\link{IcaSet}}, and fills its attributes \code{SByGene} and \code{datByGene}.
##'
##' This function is called by function \code{\link{annotInGene}} which will check the validity of the attributes \code{annotation, typeID, chipManu} and eventually \code{chipVersion} of \code{icaSet}.
##' If available, the attribute \code{annotation} of argument \code{icaSet} must be an annotation package and will be used to annotate the \code{featureNames} of \code{icaSet}.
##' If attribute \code{annotation} of argument \code{icaSet} is not available (of length 0), \code{biomaRt} is used to annotate the features. 
##' 
##' This function fills the attributes \code{SByGene} and \code{datByGene} of the argument \code{icaSet}. When several feature IDs are available for a same gene ID, the median value of the corresponding features IDs is attributed to the gene (the median of projection values is used for attribute \code{SByGene}, and the median of expression values is used for attribute \code{datByGene}).
##' 
##' When attribute \code{chipManu} of the argument \code{icaSet} is "illumina", the features are first converted into nuID using the package 'lumi*Mapping' and then annotated into genes.
##' In that case, features can only be annotated in ENTREZID or SYMBOL. It means that \code{typeID(icaSet)['geneID_annotation']} must be either  'ENTREZID' or 'SYMBOL'. You will need to annotate yourself the \code{\link{IcaSet}} object if you want to use different IDs.
##'
##' 
##' @title Features annotation
##' @param icaSet An object of class \code{\link{IcaSet}} whose features have to be annotated. The attribute \code{annotation} of this object contains the annotation package to be used.
##' @param params An object of class \code{\link{MineICAParams}} containing the parameters of the analysis.
##' @param type The ID of the object of the annotation package to be used for the annotation, must be available in \code{ls("package:package_name")}
##' @param featureId The type of the feature IDs, in the \code{biomaRt} way (type \code{listFilters(mart)} to choose one). Used when \code{annotation(icaSet)} is of length 0.
##' @param geneId The type of the gene IDs, in the \code{biomaRt} way (type \code{listAttributes(mart)} to choose one). Used when \code{annotation(icaSet)} is of length 0.
##' @return This function returns the argument \code{icaSet} with attributes
##' \code{SByGene} and \code{datByGene} filled.
##' @author Anne Biton
##' @export
##' @seealso \code{\link{annotFeatures}}, \code{\link{annotFeaturesWithBiomaRt}}, \code{\link{annotInGene}}
##' @examples 
##'
##' ## load an example of IcaSet
##' data(icaSetCarbayo)
##' params <- buildMineICAParams()
##' require(hgu133a.db)
##' ####===================================================
##' ## Use of annotation package contained in annotation(icaSet)
##' ####====================================================
##' ## annotation in SYMBOL 
##' icaSetCarbayo_annot <- annotFeaturesComp(icaSet=icaSetCarbayo, params=params, type="SYMBOL")
##' # arg 'type' is optional since the function uses contents of typeID(icaSet) as the defaults,
##' # it is specified in these examples for pedagogy views 
##'
##' ## annotation in Entrez Gene 
##' icaSetCarbayo_annot <- annotFeaturesComp(icaSet=icaSetCarbayo, params=params, type="ENTREZID") 
##'
##' \dontrun{
##' ####===================================================
##' ## Use of biomaRt, when annotation(icaSet) is of length 0
##' ####====================================================
##' ## empty attribute 'annotation' of the IcaSet object
##' # when this attribute is not specified, biomaRt is used for annotation
##' annotation(icaSetCarbayo) <- character()
##' 
##' # make sure the mart attribute is correctly defined
##' mart(icaSetCarbayo) <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
##' 
##' ## make sure elements "featureID_biomaRt" and "geneID_biomaRt" of typeID(icaSet) are correctly filled 
##' # they will be used by function 'annotFeaturesComp' through biomaRt to query the database
##' typeID(icaSetCarbayo)
##' 
##' ## run annotation of HG-U133A probe set IDs into Gene Symbols using biomaRt
##' icaSetCarbayo_annot <- annotFeaturesComp(icaSet=icaSetCarbayo, params=params)
##' 
##' }
annotFeaturesComp <- function(icaSet,
                              params,
                              type = toupper(typeID(icaSet)["geneID_annotation"]),
                              featureId = typeID(icaSet)["featureID_biomart"],
                              geneId = typeID(icaSet)["geneID_biomart"]
                              ) {
	
    if (nrow(S(icaSet))==0)
        stop("Attribute 'S' is missing in icaSet")

    if (length(chipManu(icaSet))>0)
        manu <- match.arg(tolower(chipManu(icaSet)), c("illumina","affymetrix","agilent"))
    else
        manu <- ""
    
    ## if annotation package not provided, annotation with biomaRt
    if (length(annotation(icaSet))==0 || annotation(icaSet)=="") {
        message("..Features annotation with biomaRt..","\n")
        
        if (mart(icaSet)@dataset == "")
            stop("Please fill the mart attribute of the 'icaSet' object using the function 'useMart'.")
 
        message("...Annotation of ",featureId," into ",geneId,"\n")
 	probe2gene <- annotFeaturesWithBiomaRt(features=featureNames(icaSet), featureId=featureId, geneId=geneId, mart=mart(icaSet))
    }
    else {
        message("..Features annotation with package ",annotation(icaSet),"..","\n")
        message("...Annotation into ",typeID(icaSet)["geneID_annotation"],"..","\n")

      if (manu == "illumina") {
        if ((!(toupper(type) %in% c("ENTREZID", "SYMBOL") )))
            stop("For illumina, features can only be annotated in ENTREZID or SYMBOL. You will need to annotate yourself the IcaSet object if you want to use different IDs.")

	ilmn_features = unique(featureNames(icaSet))
        
	
	message("....Mapping features<->nuID.... \n")
        if (length(agrep(organism(icaSet),c("Homo Sapiens", "HomoSapiens")))>0)
            spec <- "Human"
        else
            spec <- c("Human","Rat","Mouse")[agrep(organism(icaSet), c("Human","Rat","Mouse"))]
        if (length(spec) == 0)
            spec <- "Unknown"

        message("....Loading illumina mapping package needed for annotation....\n")

        eval(parse(text=c(Human = "library('lumiHumanIDMapping');library('lumiHumanAll.db')", Rat =  "library('lumiRatIDMapping');library('lumiRatAll.db')", Mouse =  "library('lumiMouseIDMapping');library('lumiMouseAll.db')")[spec]))

        ### Mapping Illumina features <-> nuIDs
	chipInfo <- lumi:::getChipInfo(ilmn_features, species = spec, lib.mapping = c(Human = "lumiHumanIDMapping", Rat =  "lumiRatIDMapping", Mouse =  "lumiMouseIDMapping")[spec], chipVersion = if (length(chipVersion(icaSet))>0 | chipVersion(icaSet) != "") chipVersion(icaSet) else NULL, idMapping = TRUE, verbose = FALSE) 
	# nuID to illumina probe 
	mapping <- chipInfo$idMapping
	newId <- mapping[, "nuID"]
	names(newId) = rownames(mapping)

	## restrict to the common genes
	newId = newId[intersect(names(newId),featureNames(icaSet))]

	message("......mapping nuID<->genes......\n")
        
	## Mapping nuIDs <-> genes
	## Now I need to retrieve the annotations of the illumina features thanks to their nuIDs.
	## Using lumiHumanAll.db package
	if (type == "ENTREZID") {
		nuid2genes = getEG(newId, "lumiHumanAll.db")
	}
	else if (type == "SYMBOL") {
		nuid2genes = getSYMBOL(newId, "lumiHumanAll.db")	
	}

	naInd <- which(is.na(nuid2genes))
	nuid2genes <- nuid2genes[which(!is.na(nuid2genes))]
	message(paste("Number of features = ",length(ilmn_features), "\n"))
	message(paste("Number of features with no associated gene id = ",length(naInd), "\n"))
	message(paste("Number of features with an associated gene id = ",length(nuid2genes), "\n"))
	message(paste("Number of gene ids = ",length(unique(nuid2genes)), "\n"))

	message(".......mapping illumina features<->genes....... \n")
	### gene <-> illumina probe 
	# nuid to probe
	probe2nuid = newId
	#names(probe2nuid) = names(newId)
	probe2gene = nuid2genes[probe2nuid] #long
	names(probe2gene) = names(probe2nuid)
    }
    else {
	# all features
	allfeatures = featureNames(icaSet)

	# probe -> gene
	probe2gene = annotFeatures(allfeatures, type = type, annotation = annotation(icaSet))
	noms = names(probe2gene)
	probe2gene = as.character(probe2gene) #seul moyen de le transformer en vecteur
	names(probe2gene) = noms					
	probe2gene = probe2gene[intersect(names(probe2gene),featureNames(icaSet))]
    }
   }

    probe2gene <- probe2gene[!is.na(probe2gene)]
    allGenes = unique(unlist(probe2gene))

    ## List : Gene -> c(features)
    gene2probe <- lapply(allGenes,
                           function(geneId, probe2gene) {
                                selecProbeGenes = probe2gene[probe2gene %in% geneId]	
                                return(names(selecProbeGenes))
                           }
                           , probe2gene = probe2gene
                    )
    names(gene2probe) <- allGenes



    ### Is there only one probe associated to each gene (use of BrainArray for example)?
    gene2nbfeatures <- unlist(llply(gene2probe,length))
    if (all(gene2nbfeatures == 1)) {
        oneProbeByGene = TRUE
        message("Each gene is annotated by only one probe set. \n")
        gene2probe = unlist(gene2probe)
    }
    else
        oneProbeByGene = FALSE

    names(gene2nbfeatures) <- names(gene2probe) <- allGenes
    genes2oneprobe = names(gene2nbfeatures[gene2nbfeatures == 1])

    
    dat <- assayDataElement(icaSet,"dat")
    dat_genes = dat[unlist(gene2probe[genes2oneprobe]), ]
    rownames(dat_genes) = as.character(genes2oneprobe)

    #### Gene -> median features
    if (!oneProbeByGene) {

       genesmanyfeatures = names(gene2nbfeatures[gene2nbfeatures > 1])
       message(paste ("The number of genes annotated by several features is", length(genesmanyfeatures[!is.na(genesmanyfeatures)]), "\n"))
       featuresmanygenes = unlist(gene2probe[genesmanyfeatures])
       gene2probe_manyfeatures = gene2probe[genesmanyfeatures]
       
       ### Gene -> Median of the probe sets projection
       genesComp = llply(Slist(icaSet), 
                             function (probe2proj, probe2gene, genesmanyfeatures, featuresmanygenes, gene2probe_manyfeatures) {
                                 
                                 orOrd <- unique(probe2gene)
                                 # restrict to features annotated by a gene (since probe2gene was restricted to non NA values)
                                 probe2proj <- probe2proj[names(probe2gene)]

                                 # take median projection for genes annotated by several features
                                 gene2proj <- 
                                     sapply(genesmanyfeatures,
                                            function(gene, probe2proj, probe2gene) {
                                                median(probe2proj[names(probe2gene[probe2gene %in% gene])])
                                            }
                                            , probe2proj = probe2proj, probe2gene = probe2gene
                                            )
                                                
                                 # merge with probe projection for genes annotated by only one probe
                                 gene2proj2 <- probe2proj[!(names(probe2proj) %in% featuresmanygenes)]
                                 names(gene2proj2) <- probe2gene[!(names(probe2gene) %in% featuresmanygenes)]
                                 gene2proj <- c(gene2proj,gene2proj2)[orOrd]
                                 gene2proj <- gene2proj[!is.na(gene2proj)]

                                 
                                 return(gene2proj)

                             }
                             , probe2gene = probe2gene[!is.na(probe2gene)]
                             , genesmanyfeatures = genesmanyfeatures
                             , featuresmanygenes = featuresmanygenes
                             , gene2probe_manyfeatures = gene2probe_manyfeatures
                  )



       gc()
       
       dat_medianGenes <- lapply(gene2probe[genesmanyfeatures], 
                                       function (features, dat) {
                                           # restriction to the probe sets annotating the genes
                                           dataSub =  dat[features,]
                                           return(apply(dataSub, 2, median))
                                       }
                                       , dat = dat
                                       )

       dat_medianGenes = as.data.frame(t(as.data.frame(dat_medianGenes)))
       rownames(dat_medianGenes) =  as.character(genesmanyfeatures)
       dat_genes = rbind(dat_genes, dat_medianGenes)

       
    }
    else {
        genesComp <- llply(Slist(icaSet),
                            function(probe2proj, gene2probe) {
                                probe2proj <- probe2proj[gene2probe]
                                names(probe2proj) <- names(gene2probe)
                                return(probe2proj)
                            }
                            , gene2probe = gene2probe)
    }

     gc()

    dat_genes <- dat_genes[names(genesComp[[1]]),]
    
    icaSet@datByGene <- as.data.frame(dat_genes, stringsAsFactors = FALSE, check.names = FALSE)
    icaSet@SByGene <- as.data.frame(genesComp, stringsAsFactors = FALSE, check.names = FALSE)

    # redefine witness genes if they were selected based on a different annotation
    if (length(witGenes(icaSet))>0 && length(intersect(witGenes(icaSet),geneNames(icaSet)))<nbComp(icaSet)) {
        message(".....Redefine gene witnesses.....","\n")
        witGenes(icaSet) <- selectWitnessGenes(icaSet=icaSet, params=params, level="gene", maxNbOcc=1, selectionByComp=NULL)
    }
    
    colnames(SByGene(icaSet)) <- colnames(S(icaSet))
    gc()

    return(icaSet)
}


##' This function annotates a set of features using \code{biomaRt}
##'
##' @title Annotation of features using \code{biomaRt}
##' @param features Feature IDs to be annotated
##' @param featureId The type of the feature IDs, in the \code{biomaRt} way (type \code{listFilters(mart)} to choose one)
##' @param geneId The type of the gene IDs, in the \code{biomaRt} way (type \code{listAttributes(mart)} to choose one)
##' @param mart The mart object (database and dataset) used for annotation, see function \code{useMart} of package \code{biomaRt} 
##' @return A vector of gene IDs indexed by the feature IDs.
##' @author Anne Biton
##' @export
##' @examples if (interactive()) {
##' # define the database to be queried by biomaRt
##' mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
##'
##' # annotate a set of HG-U133a probe sets IDs into Gene Symbols
##' annotFeaturesWithBiomaRt(features = c("1007_s_at", "1053_at", "117_at", "121_at", "1255_g_at"),
##' featureId="affy_hg_u133a", geneId="hgnc_symbol", mart=mart)
##' 
##' # annotate a set of Ensembl Gene IDs into Gene Symbols
##' annotFeaturesWithBiomaRt(features = c("ENSG00000101412", "ENSG00000112242",
##'                                       "ENSG00000148773", "ENSG00000131747", "ENSG00000170312",
##'                                       "ENSG00000117399"), featureId="ensembl_gene_id", geneId="hgnc_symbol", mart=mart)
##' }
annotFeaturesWithBiomaRt <- function(features,
                                     featureId, 
                                     geneId, 
                                     mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
                                     ) {
	
    #input = features
    featureId <- match.arg(featureId,c(NULL,listFilters(mart)[,1]))
    #output = genes
    geneId <- match.arg(geneId,c(NULL,listAttributes(mart)[,1]))

    if (length(geneId)>1)
        warning("Only the first element of 'geneID' will be used.")
    if (length(featureId)>1)
        warning("Only the first element of 'featureID' will be used.")
    
    geneId <- geneId[1]
    featureId <- featureId[1]
    
    # probe -> gene
    d <- getBM(attributes = c(featureId,geneId),
               filters = featureId,
               values = features,
               mart = mart)

    d <- subset(d, d[[geneId]] != "")
    ## features without available annotation which are not returned by getBM
    noAnnot <- features[!(features %in% d[[featureId]])]
    if (length(noAnnot)>0)
        message("Number of features without associated annotation is ", length(noAnnot), "\n")

    probe2gene <- as.character(d[[geneId]])
    names(probe2gene) <- as.character(d[[featureId]])
    
    probe2gene <- probe2gene[!is.na(probe2gene)]
    allGenes <- unique(unlist(probe2gene))

    return(probe2gene)
}






##' This function annotates IDs (typically gene IDs) provided by the user and returns an html file with their description.
##' 
##'
##' \code{"hgnc_symbol", "ensembl_gene_id", "description", "chromosome_name", "start_position", "end_position", "band"}, and \code{"strand"},
##' are automatically added to the list of fields available in argument \code{typeRetrieved} queried on biomaRt.
##' The web-links to www.genecards.org and www.proteinatlas.org are automatically added in the columns of the output respectively
##' corresponding to \code{hgnc_symbol} and \code{ensembl_gene_id}.
##' @title Description of features using package \code{biomaRt}.
##' @param data Either a data.frame whose rownames or one of its columns contain the IDs to be annotated, or a vector of IDs.
##' @param filename The name of the HTML file where gene annotations are written.
##' @param mart Output of function \code{useMart} from package \code{biomaRt}.
##' @param typeId The type of IDs available in \code{data}, in the biomaRt way (type \code{listFilters(mart)} to choose one).
##' @param typeRetrieved The descriptors uses to annotate the features of \code{data} (type \code{listAttributes(mart)} to choose one or several).
##' @param sortBy Name of a column of \code{data} used to order the output.
##' @param sortAbs If TRUE absolute value of column \code{sortBy} is used to order the output.
##' @param colAnnot The column containing the IDs to be annotated, if NULL or missing and argument \code{data} is a data.frame, then rownames of \code{data} must contain the IDs.
##' @param decreasing If TRUE, the output is sorted by decreasing values of the \code{sortBy} column
##' @param highlight IDs to be displayed in colour red in the returned table 
##' @param caption A title for the HTML table
##' @return This function returns a data.frame which contains annotations of the input data.
##' @examples
##' if (interactive()) {
##' ## define the database to be used 
##' mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
##'
##' ### Describe:
##' ## a set of hgnc symbols with default descriptions (typeRetrieved=NULL)
##' genes <- c("TOP2A","E2F3","E2F1","CDK1","CDC20","MKI67")
##' writeGenes(data=genes, filename="foo", mart=mart, typeId = "hgnc_symbol")
##' 
##' ## a data.frame indexed by hngc symbols, sort output according to column "values", add a title to the HTML output 
##' datagenes <- data.frame(values=rnorm(6),row.names = genes)
##' writeGenes(data=datagenes, filename="foo", sortBy = "values", caption = "Description of some proliferation genes.")
##' 
##' ## a set of Entrez Gene IDs with default descriptions 
##' genes <- c("7153","1871","1869","983","991","4288")
##' writeGenes(data=genes, filename="foo", mart=mart, typeId = "entrezgene")
##' }
##' \dontrun{
##' ## add the GO category the genes belong to
##' ## search in listAttributes(mart)[,1] which filter correspond to the Gene Ontology -> "go_id"
##' writeGenes(data=genes, filename="foo", mart=mart, typeId = "entrezgene", typeRetrieved = "go_id")
##' }
##' @author Anne Biton
##' @export
##' @seealso \code{\link[biomaRt]{getBM}}, \code{\link[biomaRt]{listFilters}}, \code{\link[biomaRt]{listAttributes}}, \code{\link[biomaRt]{useMart}}
writeGenes <- function (
                        data
                        ,filename = NULL
                        ,mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
                        ,typeId = "hgnc_symbol"
                        ,typeRetrieved=NULL  #c("hgnc_symbol", "entrezgene")
                        ,sortBy = NULL
                        ,sortAbs = TRUE
                        ,colAnnot = NULL
                        ,decreasing = TRUE
                        ,highlight = NULL
                        ,caption = ""
                        ) {

        typeId <- match.arg(typeId,c(NULL,listFilters(mart)[,1]), several.ok = TRUE)

        if (!is.null(typeRetrieved))
            typeRetrieved <- match.arg(typeRetrieved,listAttributes(mart)[,1], several.ok = TRUE)
        else
            typeRetrieved <- c()
        
        addDefaultDescr <- setdiff(c("description",  "chromosome_name", "start_position", "end_position", "band","strand", "ensembl_gene_id"),c(typeRetrieved,typeId))
        
        if (length(addDefaultDescr)>0)
            typeRetrieved <- c(typeRetrieved,addDefaultDescr)
        

        
	if (is.vector(data)) {	
            namesGenes <- data 
	    data <- data.frame(row.names = namesGenes)
	    data[[typeId]] <- namesGenes
            colData <- NULL
	}
	else {
            namesGenes <- if (is.null(colAnnot)) rownames(data) else data[[colAnnot]] 
            data[[typeId]] <- namesGenes
                    colData <- colnames(data)
        
            if (!is.null(sortBy))
                colData <- colData[colData != sortBy]

        }

        ## add symbol id by default
        if (!("hgnc_symbol" %in% c(typeId, typeRetrieved))) {
            typeRetrieved <- c("hgnc_symbol",typeRetrieved)
            addSymbol <- TRUE
        }
        if (!("ensembl_gene_id" %in% c(typeId, typeRetrieved))) {
            typeRetrieved <- c(typeRetrieved, "ensembl_gene_id")
            addSymbol <- TRUE
        }
        
        d <- getBM(attributes = c(typeId, typeRetrieved),
                   filters = typeId,
                   values = namesGenes,
                   mart = mart)
        
        if (!is.null(colData)) {
            colData <- setdiff(colData,c(typeId,typeRetrieved))
            if (length(colData)==0)
                colData <- NULL
        }

        ## add genes without available annotation which are not returned by getBM
	noAnnot <- namesGenes[!(namesGenes %in% d[[typeId]])]
	if (length(noAnnot) > 0) {
		nr = nrow(d)
		d[(nr+1):(nr+length(noAnnot)),typeId] <- noAnnot
                
	}

        ## merge to original data
        d <- merge(data, d, by = typeId, all.x = TRUE)

	## Rank data frame according to the values in the sortBy column, or by start position
	if (!is.null(sortBy)) {
            if (!is.numeric(d[[sortBy]])) stop("The column used to sort the data.frame is not numeric")
            if (sortAbs) 
                d <- d[order(abs(as.numeric(d[[sortBy]])), decreasing = decreasing),]
            else
                d <- d[order(as.numeric(d[[sortBy]]), decreasing = decreasing),]
            
          d[[sortBy]] <- signif(d[[sortBy]],4) 
          d[[sortBy]] <- as.character(d[[sortBy]]) 
          
        }
        else d <- d[order(as.numeric(d[["start_position"]]),decreasing = FALSE),]
	
	if (!is.null(filename)) {
          dhtml = d

          if (!is.null(highlight)) {
              classGene <- rep("gene",times=length(dhtml$hgnc_symbol))
              classGene[which(dhtml$hgnc_symbol %in% highlight)] <- "geneh"
              dhtml$hgnc_symbol <- paste("<a class =", classGene, " href=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=", dhtml$hgnc_symbol, "\">", dhtml$hgnc_symbol,"</a>",sep = "")
              
          }
          else dhtml$hgnc_symbol <- paste("<a class = 'gene' href=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                                          dhtml$hgnc_symbol, "\">",
                                          dhtml$hgnc_symbol,"</a>",
                                          sep = "")
          dhtml$ensembl_gene_id <- paste("<a class = 'normal' href='http://www.proteinatlas.org/gene_info.php?ensembl_gene_id=",
                                         dhtml$ensembl_gene_id, "'>",
                                         dhtml$ensembl_gene_id,"</a>",
                                         sep = "")

          dhtml = dhtml[,c(typeId, if(!is.null(colData)) c(sortBy, colData) else sortBy, if (!is.null(typeRetrieved)) typeRetrieved)]

          x <- xtable(dhtml,
                      caption = caption,
                      align = c("l","l","l","c","c","l",rep("c",ncol(dhtml)-5)),
		      digits = 4)
          
          x <- capture.output(print(x,
                                    type = "html",
                                    sanitize.text.function = force,
                                    include.rownames=FALSE,
                                    caption.placement = "top")
                              )
 

	style <- 
	"
	<head>
	 <style type='text/css'>

	     a.comp:link{text-decoration:none;color:white;}
	     a.comp:visited{text-decoration:none;color:white;}
	     a.comp:hover{text-decoration:none;font-weight:bold;color:black;}

	     a.gene:link{text-decoration:none;color:black;font-weight:bold;}
	     a.gene:visited{text-decoration:none;color:black;font-weight:bold;}
	     a.gene:hover{text-decoration:none;font-weight:bold;color:#B22222;}

	     a.geneh:link{text-decoration:none;color:red;font-weight:bold;}
	     a.geneh:visited{text-decoration:none;color:red;font-weight:bold;}
	     a.geneh:hover{text-decoration:none;font-weight:bold;color:red;}

	     a.normal:link{text-decoration:none;color:black;}
	     a.normal:visited{text-decoration:none;color:black;}
	     a.normal:hover{text-decoration:none;font-weight:bold;color:black;}

	     TH.geneTable {
	        padding: 2px;
	            border:2px solid;
	            color: white;
	            border-color: #8B1A1A;
	            text-align: center;
	            font-weight: bold;
	            background-color: #B22222
	     }

	 </style>
	</head>"

        x <- gsub(x,pattern="<TABLE",replacement = "<TABLE style='font-family:Helvetica; font-size:12px; border-top:solid thin black; border=2px;' ", ignore.case = TRUE)
        x <- gsub(x,pattern="<CAPTION",replacement="<CAPTION  style='color:black; text-align:left; font-size:14px;' ", ignore.case = TRUE)
        x <- gsub(x,pattern="<TH>",replacement = "<TH  class='geneTable'>", ignore.case = TRUE)
        x <- paste(style,x,sep="")                
	x <- write(x, file = paste(gsub(filename,pattern=".htm",replacement="",fixed=TRUE),".htm",sep=""))


	}


	return(d)
}




##### =====================================
####   Write contributing genes for each
####   component in separate files
##### =====================================

writeProjByComp <- function(icaSet
                             ### The IcaSet-object 
                             ,params
                             ### The MineICAParams-object which contains the parameters of the analysis, the files are written in the path provided in the attribute \code{genesPath} of \code{params}. The attribute \code{selCutoff} is used to select the features/genes that are annotated and writte.
                             ,mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
                             ### The mart object used for annotation, see function \code{useMart} of package biomaRt 
                             ,typeRetrieved=NULL #= listFilters(mart)[,1] #c("hgnc_symbol", "entrezgene")
                             ### The type of the annotation you want to use to annotate the IDs of argument \code{data}
                             ,addNbOcc = TRUE
                             ### If TRUE, the number of components the features/genes contribute to is added to the output. A gene/feature is considered as a contributor of a component if its scaled projection value is higher than attribute \code{selCutoff} of \code{icaSet}.
                             ,selectionByComp = NULL
                             ### A list which contains the feature/gene projections on each component, already restricted to the ones considered as contributors.
                             ,level = c("features","genes")
                             ### The attribute of \code{icaSet} that will be annotated: either the feature projections, or the gene projections
                             ,typeId
                             ### The type of the IDs the features or the genes of \code{icaSet} correspond to. It must be provided in the biomaRt way (type listFilters(mart) to choose the typeId).
                             ,selCutoffWrite=2.5
                             )
{
    ##details<< One file is created by component, each file is named by the index of the components (\code{indComp(icaSet)}) and located in the path \code{genePath(params)}.
    ## The genes are ranked according to their absolute projection values.
    ##
    ##
    ## In addition, this function write an html file named "genes2comp" which provides, for each feature/gene, the number of components in which the gene is contributor (according to the threshold \code{cutoffSel(params)}), and its projection value on all the components.  
## The projection values are scaled.
##
    sel <- i <- comp <- cutoff <- NULL
    
    level <- match.arg(tolower(level), choices = c("genes","features"))
    
    switch(level,
           features={SlistByGene=Slist(icaSet); object="S"},
           genes={SlistByGene=SlistByGene(icaSet); object="SByGene"}
           )

    if (missing(typeId))
        typeId <- typeID(icaSet)[c(features="featureID_biomart",genes="geneID_biomart")[level]]
    
    
    if (length(SlistByGene)==0)
        stop(paste("The attribute",object,"is empty."))
    
    system(paste("mkdir",genesPath(params)), ignore.stderr=TRUE)
    
    selCutoff <- selCutoff(params)

    if ((length(selCutoff)>1) & (length(selCutoff)<nbComp(icaSet))) {
        warning("The length of 'selCutoff' must be either 1 or equal to the number of components. Only the first value is used.")
        selCutoff <- selCutoff[1]
        params["selCutoff"] <- selCutoff
    }
        
        
    if (is.null(selectionByComp))
        selectionByComp <- selectContrib(SlistByGene, cutoff = selCutoff)

    selectionByComp_write <- selectContrib(SlistByGene, cutoff = selCutoffWrite)

        

    #### ******* Write one table by component with the gene projections and annotations
    gene2comp = NULL
    if (addNbOcc) {
        ## If selectionByComp is NULL, it will be recomputed by nbOccByGeneInComp
        ## sometimes we recompute the selection or we use a different selection from data because in the global script the done selection is made with a minimal cutoff in order to write a lot of genes
        gene2comp <- nbOccByGeneInComp(sel = selectionByComp, cutoff = selCutoff)
        gene2comp_nb <- unlist(lapply(gene2comp, length))
        names(gene2comp_nb) <- names(gene2comp)
        gene2comp_char <- as.character(lapply(gene2comp, paste, collapse = ",")); names(gene2comp_char) = names(gene2comp)
    }
    else {
        gene2comp_nb <- NULL
        gene2comp_char <- NULL
    }

    if (length(selCutoff)==1)
        selCutoff <- rep(selCutoff,nbComp(icaSet))

    annotD <- 
    foreach(sel=selectionByComp_write,i=1:length(SlistByGene)) %do% {
             d <- data.frame(scaled_proj = sel);
             rownames(d) <- names(sel);
             if (addNbOcc) {
               gene2comp_nb[is.na(gene2comp_nb)] <- 0

               d[[paste("nbOcc_forThreshold:",selCutoff[i],sep="")]] <- "0"
               d[[paste("comp_forThreshold:",selCutoff[i],sep="")]] <- ""
               
               dbis <- d
               d[intersect(rownames(d),names(gene2comp_nb)),paste("nbOcc_forThreshold:",selCutoff[i],sep="")] <-
                   paste("<a  class='normal' href='genes2comp_threshold",
                         if (length(selCutoff)>1) paste(selCutoff,collapse="_") else selCutoff,".htm#",
                         intersect(rownames(d),names(gene2comp_nb)),"'>",
                         gene2comp_nb[intersect(rownames(d),names(gene2comp_nb))],"</a>",sep="")
               dbis[intersect(rownames(d),names(gene2comp_nb)),paste("nbOcc_forThreshold:",selCutoff[i],sep="")] <- gene2comp_nb[intersect(rownames(d),names(gene2comp_nb))]

               #d[[paste("nbOcc_forThreshold:",selCutoff,sep="")]][is.na(d[[paste("nbOcc_forThreshold:",selCutoff,sep="")]])] <- "0"
               #dbis[[paste("nbOcc_forThreshold:",selCutoff,sep="")]][is.na(dbis[[paste("nbOcc_forThreshold:",selCutoff,sep="")]])] <- "0"
               
               d[intersect(rownames(d),names(gene2comp_nb)),paste("comp_forThreshold:",selCutoff[i],sep="")] <-
                   sapply(gene2comp_char[intersect(rownames(d),names(gene2comp_nb))], 
                          function(x) {
                             ssplit <- strsplit(x,split=",")[[1]]
                             #if (is.na(ssplit[1])) ssplit <- "0"
                             paste(paste("<a  class='normal' href='",ssplit,".htm'>",ssplit,"</a>",sep=""),collapse = ",")
                           }
                          )
               dbis[intersect(rownames(d),names(gene2comp_nb)),paste("comp_forThreshold:",selCutoff[i],sep="")] <- gene2comp_char[intersect(rownames(d),names(gene2comp_nb))]


             }

             if (length(sel) != 0) {
                 resW <- writeGenes(data = d,
                                    filename = paste(genesPath(params),i,sep=""),
                                    mart = mart,
                                    typeId = typeId,
                                    typeRetrieved = typeRetrieved,
                                    sortBy = "scaled_proj",
                                    caption = paste("<b>Component ",comp," : this table provides the scaled genes projections exceeding ",selCutoffWrite,". <br> The genes are ranked according to their scaled projection on the component.</b> <br>Click on the gene symbol and ensembl gene id to respectively access to the corresponding GeneCards and Protein Atlas pages. <br> Click on the components index to open the corresponding genes tables. Click on the number of occurrences of a gene to see its projections values on every component.", sep = ""))
                 if (addNbOcc) {
                     
                     resW[[paste("comp_forThreshold:",selCutoff[i],sep="")]] <- dbis[[paste("comp_forThreshold:",selCutoff[i],sep="")]][match(resW[[typeId]], rownames(dbis))]
                     resW[[paste("nbOcc_forThreshold:",selCutoff[i],sep="")]] <- dbis[[paste("nbOcc_forThreshold:",selCutoff[i],sep="")]][match(resW[[typeId]], rownames(dbis))]
                 }
             }
             else
                 resW <- NULL

             return(resW)
                
      }

    nbOcc.inallComp <- nbOccInComp (icaSet = icaSet,
                                    params = params,
                                    selectionByComp = NULL, 
                                    file = paste(genesPath(params),"/genes2comp_threshold",if (length(selCutoff)>1) paste(selCutoff,collapse="_") else selCutoff,".htm",sep = ""),
                                    level = level)
    nbOcc.inallComp <- as.data.frame(apply(apply(nbOcc.inallComp,2,gsub,pattern="<b>",replacement=""),2,gsub,pattern="</b>",replacement=""), stringsAsFactors=FALSE)
    colnames(nbOcc.inallComp)[(ncol(nbOcc.inallComp)-nbComp(icaSet)+1):ncol(nbOcc.inallComp)] <- indComp(icaSet)
    

    
    return(list(listAnnotComp=annotD, nbOccInComp=nbOcc.inallComp))
### This function returns a list of two elements: the first element is named 'listAnnotComp' and contain a list with the output of the function 'writeGenes' for each component, and the second element is named 'nbOccInComp' and provides a data.frame which gives for each feature/gene (row) its projection values across all the components (columns).
    ##seealso<<writeGenes, getBM, listFilters, listAttributes, useMart, selectContrib, nbOcc.inallComp
}





nbOccByGeneInComp <- function(
### For each feature/gene, this function returns the indices of the components they contribute to. 
                              Slist,
                              ### The list of components, each element contains the projection values
                              cutoff,
                              ### The threshold used to define a gene as contributor
                              sel
                              ### The list of components already restricted to the contributing genes
                              ) {

  if (missing(sel) || is.null(sel))
      sel <- selectContrib(Slist, cutoff = cutoff)
  names(sel) = NULL
  ## all genes (or features) that have a zvalue of projection higher than cutoffOnZval on components
  allGenes = unique(unlist(lapply(sel,names)))
  ## for each gene, return the ids of the components it belongs to
  genes2comp = lapply(allGenes,
    function(gene, sel) {
      isPresent = lapply(sel,
             function(comp, gene) {
               if (gene %in% names(comp)) return (TRUE)
               else return(FALSE)
             }
             , gene = gene)
      presentInComp = which(isPresent == TRUE)
      return(presentInComp)
             
    }
    , sel = sel)
  
  names(genes2comp) <- allGenes
  
  return(genes2comp)
}

##' For each feature/gene, this function returns the indices of the components they contribute to.
##'
##' @title nbOccInComp_simple
##' @param icaSet An object of class \code{\link{IcaSet}}. 
##' @param params An object of class \code{\link{MineICAParams}} containing the parameters of the analysis. \code{cutoffSel(params)} is used as a threshold on the absolute projections to select the contributing features/genes. 
##' @param selectionByComp The list of components already restricted to the contributing features/genes (each element is a vector of projection values indexed by features or genes).
##' @param level The attribute of \code{icaSet} to be used, the occurences of either the \code{"features"} (using \code{S(icaSet)}) or the \code{"genes"} (using \code{SByGene(icaSet)}) will be reported.
##' @return Returns a data.frame whose columns are: \code{gene} the feature or gene IDs, \code{nbOcc} the number of components the gene contributes to, \code{components} the indices of those components.
##' @author Anne Biton
##' @keywords internal 
##' @examples 
##'  data(icaSetCarbayo)
##' params <- buildMineICAParams(resPath="carbayo/")
##' nbOcc <- MineICA:::nbOccInComp_simple(icaSet=icaSetCarbayo, params=params, level="genes")
##' 
nbOccInComp_simple <- function(
                               icaSet,
                               params,
                               selectionByComp = NULL,
                               level = c("features","genes")
                               )
{
    level <- match.arg(tolower(level), choices = c("features","genes"))
    
    switch(level,
           features={Slist=Slist(icaSet)},
           genes={Slist=SlistByGene(icaSet)}
           )
    
    cutoff <- selCutoff(params)
    
    if (is.null(selectionByComp))
        selectionByComp <- selectContrib(Slist, cutoff = cutoff)
    
    genesToComps_index <- nbOccByGeneInComp(Slist = Slist, sel = selectionByComp, cutoff = cutoff)
    genesToComps_nb <- unlist(lapply(genesToComps_index, length));
    names(genesToComps_nb) <- names(genesToComps_index)
    genesToComps <- unlist(lapply(genesToComps_index, paste, collapse = ","));
    names(genesToComps) <- names(genesToComps_index)
    d.genesToComps <- data.frame(gene = names(genesToComps),  nbOcc = genesToComps_nb[names(genesToComps)], components = as.character(genesToComps), stringsAsFactors = FALSE)
    return(d.genesToComps)
    ### Returns a data.frame whose columns are: 'gene' the feature or gene ID, 'nbOcc' the number of components on which the gene contributes according to the threshold, 'components' the indices of these components.
}

##' selectWitnessGenes
##' 
##' This function selects a gene per component. 
##' 
##' Selects as feature/gene witness, for each component, the first gene whose
##' absolute projection is greater than a given threshold in at the most
##' \code{maxNbOcc} components. 
##' These witnesses can then be used as representatives of the expression behavior of the contributing genes of the components.
##'
##' When a feature/gene respecting the given constraints is not found, \code{maxNbOcc} is incremented of one until a gene is found.
##' 
##' 
##' @param icaSet An object of class \code{\link{IcaSet}}
##' @param params An object of class \code{\link{MineICAParams}} containing the parameters of the analysis, the attribute \code{cutoffSel} is used as the threshold.
##' @param level The attribute of \code{icaSet} to be used, the witness
##' elements will be either selected within the \code{"features"} or the \code{"genes"}
##' @param maxNbOcc The maximum number of components where the genes can have an
##' absolute projection value higher than \code{cutoffSel(params)} in order to
##' be selected.
##' @param selectionByComp The list of components already restricted to the
##' contributing genes
##' @return This function returns a vector of IDs.
##' @export
##' @author Anne Biton
##' @examples 
##' ## load an example of IcaSet
##' data(icaSetCarbayo)
##'
##' ## define parameters: features or genes are considered to be contributor
##' # when their absolute projection value exceeds a threshold of 4. 
##' params <- buildMineICAParams(resPath="carbayo/", selCutoff=4)
##' 
##' ## selection, as gene witnesses, of the genes whose absolute projection is greater than 4
##' # in at the most one component. I.e, a gene is selected as a gene witness of a component
##' # if he has a large projection on this component only.
##' selectWitnessGenes(icaSet=icaSetCarbayo, params=params, level="genes", maxNbOcc=1)
##' 
##' ## selection, as gene witnesses, of the genes whose absolute projection is greater than 4
##' # in at the most two components.
##' # I.e, a gene is selected as a gene witness of a given component if he has a large projection
##' # in this component and at the most another. 
##' selectWitnessGenes(icaSet=icaSetCarbayo, params=params, level="genes", maxNbOcc=2)
##' 
##' 
selectWitnessGenes <- function(icaSet,
                               params,
                               level = c("genes","features"),
                               maxNbOcc = 1,
                               selectionByComp = NULL
                               ) {

    level <- match.arg(level)
    
    if ((maxNbOcc) < 1)
        stop("maxNbOcc must at least equal 1.")
    
     switch(level,
           features={Slist=Slist(icaSet)},
           genes={Slist=SlistByGene(icaSet)}
           )

    
    cutoff <- selCutoff(params)

     if (is.null(selectionByComp) | missing(selectionByComp))
        selectionByComp <- selectContrib(Slist, cutoff = cutoff)

    ## For components with no genes selected using the cutoff,  
    indEmptySel <- which(unlist(lapply(selectionByComp,length))==0)
    selectionByComp[indEmptySel] <- Slist[indEmptySel]
        
    docc <- nbOccInComp_simple(icaSet = icaSet, params = params, selectionByComp = selectionByComp)
    nbocc <- docc$nbOcc
    names(nbocc) <- docc[,1]
    ## Select as witness genes the first gene whose contribution is greater than cutoff in at the most maxNbOcc components.
    witGenes <- 
        lapply(selectionByComp,
           function(sel, nbocc, maxNbOcc) {
                   sel <- sort(abs(sel), decreasing = TRUE)
                   ind <- which(nbocc[names(sel)] <= maxNbOcc)
                   while (is.na(ind[1])) {
                       maxNbOcc <- maxNbOcc+1
                       ind <- which(nbocc[names(sel)] <= maxNbOcc)
                   }
                   return(names(sel[ind[1]]))
           }, nbocc = nbocc
            , maxNbOcc = maxNbOcc
           )
                   
    return(unlist(witGenes))
    ### This function returns a vector of IDs.
}

##' For each feature/gene, this function returns the components they contribute
##' to and their projection values across all the components.
##' 
##' A feature/gene is considered as a contributor when its scaled projection value exceeds the threshold \code{selCutoff(icaSet)}.
##' 
##' This function plots the number of times the feature/gene is a contributor as a function of the standard deviation of its expression profile.
##' 
##' The created files are located in \code{genePath(params)}. An extensiom '.htm' and '.pdf' is respectively added to the \code{file} name for the data.frame and the plot outputs.
##' @title Select components the features contribute to
##' @param icaSet An object of class \code{\link{IcaSet}}
##' @param params An object of class \code{\link{MineICAParams}} containing the parameters of the
##' analysis, the attribute \code{cutoffSel} is used as a threshold on the
##' absolute projections to determine which genes contribute to the components.
##' @param selectionByComp The list of components already restricted to the
##' contributing genes
##' @param level The attribute of \code{icaSet} to be used, are reported the
##' occurences of either the \code{"features"} or the \code{"genes"}.
##' @param file The file where the output data.frame and plots are written. 
##' @return Returns a data.frame whose columns are: 'gene' the feature or gene
##' ID, 'nbOcc' the number of components on which the gene contributes according
##' to the threshold, 'components' the indices of these components, and then the
##' component indices which contain its projection values.
##' @export
##' @author Anne Biton
##' @examples 
##' data(icaSetCarbayo)
##' params <- buildMineICAParams(resPath="carbayo/")
##' nbOcc <- nbOccInComp(icaSet=icaSetCarbayo, params=params, level="genes", file="gene2MixingMatrix")
##' 
nbOccInComp <- function(
                        icaSet
                        ,params
                        ,selectionByComp = NULL
                        ,level = c("features","genes")
                        ,file  = NULL
                        )
{
    
    level <- match.arg(tolower(level), choices = c("features","genes"))
     switch(level,
           features={Slist=Slist(icaSet)},
           genes={Slist=SlistByGene(icaSet)}
           )

    ##Slist <- icaSet[object]
    cutoff <- selCutoff(params)
    if (is.null(selectionByComp))
        selectionByComp <- selectContrib(Slist, cutoff = cutoff)

    d.genesToComps <- nbOccInComp_simple(icaSet = icaSet, params = params, selectionByComp = selectionByComp, level = level)

    

            if (is.null(datByGene(icaSet)) | nrow(datByGene(icaSet))==0)
                dat <- read.table (datfile(params), header = TRUE)
            else {
                if (level == "genes")
                    dat <- datByGene(icaSet)
                else
                    dat <- assayDataElement(icaSet,"dat")
            }
            dat.sd = getSdExpr(d.genesToComps$gene, dat)
            d.genesToComps$sd_expr = signif(dat.sd[d.genesToComps$gene],5)
            d.genesToComps = d.genesToComps[order(d.genesToComps$sd_expr, decreasing = TRUE),]

	    ## plot sd expr vs projection
            d.aggregate <- aggregate(d.genesToComps[,c("nbOcc","sd_expr")],by=list(occ=as.factor(d.genesToComps$nbOcc)),FUN=mean)
    
	    pdf(file=gsub(file,pattern="htm",replacement="pdf"),width = 6, height = 6, title=gsub(file,pattern="htm",replacement="pdf"))
            if (anyDuplicated(d.aggregate$nbOcc) != 0)
                d.aggregate <- d.aggregate[-which(duplicated(d.aggregate$nbOcc)),]

	    plot(d.aggregate$sd_expr,
                 d.aggregate$nbOcc,
                 pch = 19,
                 col = "black",
                 ylab = "number of occurence",
                 xlab = "mean of sd(expr)" ,
                 type = "b",
                 main = paste("Number of gene contributions to a component \n as a function of the average standard deviation of their profiles","\n cutoff ",if (length(unique(cutoff))>1) paste(cutoff,collapse=",") else cutoff, sep=""),
                 cex.main = 1) 
	    dev.off()
    

   if (!is.null(file)) {

       ## merge with a table providing the genes scaled projections on the components
       compsel <- proj <- NULL
       genesProj <- foreach(proj = Slist, compsel = selectionByComp, .combine = cbind) %dopar% {
	   # proj of the contributing genes	
	   p <- signif(proj[d.genesToComps$gene],4)
	   names(p) <- d.genesToComps$gene
	   # mark the contributing genes in an other color(i.e. genes exceeding the threshold)
	   p[names(compsel)] <- paste("<b>",p[names(compsel)],"</b>", sep = "")	
	   return(p)
       }

       genesProj <- as.data.frame(genesProj, row.names = d.genesToComps$gene)

       d.genesToCompsbis <-cbind(d.genesToComps, genesProj[d.genesToComps$gene,]) 
           
       names(genesProj) <- paste("<a  class='comp' href='",c(1:ncol(genesProj)),".htm'>",c(1:length(Slist)),"</a>",sep="") 
       d.genesToComps <- cbind(d.genesToComps, genesProj[d.genesToComps$gene,]) 
       d.genesToComps$gene <- paste("<a class='gene' name='",d.genesToComps$gene,"' href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", d.genesToComps$gene, "'>", d.genesToComps$gene,"</a>",sep = "")
       
       d.genesToComps$components <-
	sapply(as.character(d.genesToComps$components), 
	       function(x) {
			ssplit <- strsplit(x,split=",")[[1]]
			paste(paste("<a  class='normal' href='",ssplit,".htm'>",ssplit,"</a>",sep=""),collapse = ",") 
		})

       

	style <- 
	"
	<head>
	 <style type='text/css'>

	     a.comp:link{text-decoration:none;color:white;}
	     a.comp:visited{text-decoration:none;color:white;}
	     a.comp:hover{text-decoration:none;font-weight:bold;color:black;}

	     a.gene:link{text-decoration:none;color:black;font-weight:bold;}
	     a.gene:visited{text-decoration:none;color:black;font-weight:bold;}
	     a.gene:hover{text-decoration:none;font-weight:bold;color:#B22222;}

	     a.normal:link{text-decoration:none;color:black;}
	     a.normal:visited{text-decoration:none;color:black;}
	     a.normal:hover{text-decoration:none;font-weight:bold;color:black;}

	     TH.geneTable {
	        padding: 2px;
	            border:2px solid;
	            color: white;
	            border-color: #8B1A1A;
	            text-align: center;
	            font-weight: bold;
	            background-color: #B22222
	     }

	 </style>
	</head>"

      x <- xtable(d.genesToComps, caption = paste("<b>This table describes the genes contributing to at least one component, given the threshold(s) ", if (length(cutoff>1)) paste(cutoff,collapse=",") else cutoff, ". <br> The genes are ranked according to the standard deviation of their expression profiles.</b> <br>Click on the gene id to access to its GeneCards page. <br> Click on the components index to access the scaled projections of the complete list of genes. ", sep = ""),align = c("l","l","c","l",rep("l",ncol(d.genesToComps)-3)))

      x <- capture.output(print( x,
                                type = "html",
                                sanitize.text.function = force,
                                include.rownames=FALSE,
                                caption.placement = "top")
                          )
      x <- gsub(x,pattern="<TABLE",replacement = "<TABLE style='font-family:Helvetica; font-size:12px; border-top:solid thin black; border=2px;' ", ignore.case = TRUE)
      x <- gsub(x,pattern="<CAPTION",replacement="<CAPTION  style='color:black; text-align:left; font-size:14px;' ", ignore.case = TRUE)
      x <- gsub(x,pattern="<TH>",replacement = "<TH  class='geneTable'>", ignore.case = TRUE)
      x <- paste(style,x,sep="")                
      write(x, file = file)

   }
    else {
        d.genesToComps <- cbind(d.genesToComps, as.data.frame(Slist)[d.genesToComps$gene,])
        colnames(d.genesToComps)[ncol(d.genesToComps):(ncol(d.genesToComps)-length(compNames(icaSet)))] <- rev(compNames(icaSet))
    }

    return(d.genesToCompsbis)

}



##' getSdExpr
##' 
##' Compute standard deviation of the gene expression
##' 
##' 
##' @param features IDs
##' @param dat Expression data indexed by IDs
##' @return Returns a vector
##' @author Anne Biton
##' @export
##' @keywords internal
##' @examples 
##' dat <- matrix(rnorm(1000),ncol=10,nrow=100)
##' rownames(dat) <- 1:100
##' MineICA:::getSdExpr(features = 2:20, dat = dat)
##' 
getSdExpr <- function (
                       ### Compute standard deviation of the gene expression
                       features
                       ### IDs 
                       ,dat
                       ### Expression data indexed by IDs
                       ) {
		genes.sdExpr <- as.numeric( sapply(features, function (gene, dat) { return(sd(as.numeric(dat[gene,]))) }, dat))
		names(genes.sdExpr) <- features
		return(genes.sdExpr)
                ### Returns a vector 
}




##' This function builds an object of class \code{\link[MineICA:MineICAParams-class]{MineICAParams}}.
##' It contains the parameters that will be used by function \code{\link{runAn}} to analyze the ICA decomposition contained in an object of class \code{\link[MineICA:IcaSet-class]{IcaSet}}. 
##'
##' 
##' @title Creates an object of class MineICAParams
##' @param Sfile A txt file containing the Source matrix S. 
##' @param Afile A txt file containing the Mixing matrix A. 
##' @param datfile A txt file containing the data (e.g expression data) on which the decomposition was calculated.
##' @param annotfile Either a "rda" or "txt" file containing the annotation data for the samples (must be of dimensions samples x annotations). 
##' @param resPath The path where the outputs of the analysis will be written, default is the current directory. 
##' @param genesPath The path _within_ the resPath where the gene projections will be written. If missing, will be automatically attributed as \code{resPath}/ProjByComp/.
##' @param annot2col A vector of colors indexed by annotation levels. If missing, will be automatically attributed using function \code{annot2Color}.
##' @param pvalCutoff The cutoff used to consider a p-value significant, default is 0.05.
##' @param selCutoff The cutoff applied to the absolute feature/gene projection values to consider them as contributors, default is 3. Must be either of length 1 and the same treshold is applied to all components, or of length equal to the number of components in order to a specific threshold is for each component.
##' @return An object of class \code{\link{MineICAParams}}
##' @export
##' @author Anne Biton
##' @seealso \code{\link[MineICA:MineICAParams-class]{MineICAParams}}, \code{\link[MineICA]{runAn}}
##' @examples 
##'
##' ## define default parameters and fill resPath
##' params <- buildMineICAParams(resPath="resMineICACarbayo/")
##' 
##' ## change the default cutoff for selection of contribugint genes/features 
##' params <- buildMineICAParams(resPath="resMineICACarbayo/", selCutoff=4)
##' 
buildMineICAParams <- function (Sfile = new("character"), Afile=new("character"), datfile= new("character"), annotfile=new("character"), resPath="", genesPath, annot2col=new("character"), pvalCutoff= 0.05, selCutoff=3) {

 if (resPath != "") {
     home <- system("echo $HOME",intern=TRUE);
     resPath <- gsub(resPath,pattern="~",replacement=home);
 }

 if (missing(genesPath))
     genesPath <- paste(resPath,"ProjByComp/",sep="")
 else
     genesPath <- paste(resPath,genesPath,sep="/")
 

  
 params <- new("MineICAParams",
           Sfile = Sfile,
           Afile= Afile,
           datfile= datfile,
           annotfile=annotfile,
           resPath= resPath ,
           genesPath=genesPath,
           annot2col=annot2col,
           pvalCutoff=pvalCutoff,
           selCutoff=selCutoff)

 if (resPath != "")
     system(paste("mkdir",resPath(params)), ignore.stderr=TRUE)
 system(paste("mkdir",genesPath(params)), ignore.stderr=TRUE)

 return(params)

       }


##' This function builds an object of class \code{\link{IcaSet}}.
##' 
##' 
##' @param params An object of class \code{\link{MineICAParams}} containing the parameters of the analysis
##' @param A The mixing matrix of the ICA decomposition (of dimension samples x components).
##' @param S The source matrix of the ICA decomposition (of dimension features x components).
##' @param dat The data matrix the ICA was applied to  (of dimension features x samples).
##' @param pData Phenotype data, a data.frame which contains the sample
##' informations of dimension samples x annotations.
##' @param fData Feature data, a data.frame which contrains the feature
##' descriptions of dimensions features x annotations.
##' @param witGenes A vector of witness genes. They are representative of the
##' expression behavior of the contributing genes of each component. 
##' If missing or NULL, they will be automatically attributed using function \code{\link{selectWitnessGenes}}.
##' @param compNames A vector of component labels.
##' @param refSamples A vector of reference sample IDs (e.g the "normal" samples).
##' @param annotation An annotation package (e.g a ".db" package specific to the
##' microarray used to generate \code{dat})
##' @param chipManu If microarray data, the manufacturer: either 'affymetrix' or 'illumina'.
##' @param chipVersion For illumina microarrays: the version of the microarray.
##' @param alreadyAnnot TRUE if the feature IDs contained in the row names of \code{dat} and \code{S} already correspond to the final level of annotation (e.g if they are already gene IDs). In that case, no annotation is performed.
##' @param typeID A character vector specifying the annotation IDs, it includes
##' three elements : \describe{ 
##' \item{geneID_annotation}{the IDs from the
##' package to be used to annotate the features into genes. It will be used to
##' fill the attributes \code{datByGene} and \code{SByGene} of the \code{icaSet}.
##' It must match one of the objects the corresponding package supports
##' (you can access the list of objects by typing ls("package:packagename")). If
##' no annotation package is provided, this element is not useful.}
##' \item{geneID_biomart}{the type of gene IDs, as available in
##' \code{listFilters(mart)}; where mart is specified as described in \code{\link[biomaRt]{useMart}}.
##' If you have directly built the IcaSet at the
##' gene level (i.e if no annotation package is used), \code{featureID_biomart} and
##' \code{geneID_biomart} will be identical.} 
##' \item{featureID_biomart}{the
##' type of feature IDs, as available in \code{listFilters(mart)}; where
##' \code{mart} is specified as described in function \code{\link[biomaRt]{useMart}}.
##' Not useful if you work at the gene level.} }
##' @param runAnnot If TRUE, \code{icaSet} is annotated with function \code{annotInGene}.
##' @param organism The organism the data correspond to.
##' @param mart The mart object (database and dataset) used for annotation, see function \code{useMart} of package \code{biomaRt} 
##' @seealso \code{\link{selectWitnessGenes}}, \code{\link{annotInGene}}
##' @return An object of class IcaSet 
##' @export
##' @author Anne Biton
##' @examples
##'
##' dat <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
##' rownames(dat) <- paste("g", 1:1000, sep="")
##' colnames(dat) <- paste("s", 1:10, sep="")
##' 
##' ## build a data.frame containing sample annotations
##' annot <- data.frame(type=c(rep("a",5),rep("b",5)))
##' rownames(annot) <- colnames(dat)
##' 
##' ## run ICA
##' resJade <- runICA(X=dat, nbComp=3, method = "JADE")
##' 
##' ## build params
##' params <- buildMineICAParams(resPath="toy/")
##'
##' ## build IcaSet object
##' icaSettoy <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S),
##'                          dat=dat, pData=annot, alreadyAnnot=TRUE)
##' params <- icaSettoy$params
##' icaSettoy <- icaSettoy$icaSet
##' 
##' \dontrun{
##' ## load data
##' library(breastCancerMAINZ)
##' data(mainz)
##' 
##' ## run ICA
##' resJade <- runICA(X=dataMainz, nbComp=10, method = "JADE", maxit=10000) 
##' 
##' ## build params
##' params <- buildMineICAParams(resPath="mainz/")
##' 
##' ## build IcaSet object
##' 
##' # fill typeID, Mainz data originate from affymetrix HG-U133a  microarray and are indexed by probe sets
##' # we want to annotate the probe sets into Gene Symbols
##' typeIDmainz <-  c(geneID_annotation="SYMBOL", geneID_biomart="hgnc_symbol", featureID_biomart="affy_hg_u133a")
##' 
##' icaSetMainz <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S),
##'                              dat=exprs(mainz), pData=pData(mainz),
##'                              annotation="hgu133a.db", typeID= c(geneID_annotation = "SYMBOL",
##'                              geneID_biomart = "hgnc_symbol", featureID_biomart = "affy_hg_u133a"),
##'                              chipManu = "affymetrix", runAnnot=TRUE,
##'                              mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl"))
##' }
buildIcaSet <- function(params,
                        A,
                        S,
                        dat,
                        pData = new("data.frame"),
                        fData = new("data.frame"),
                        witGenes= new("character"),
                        compNames= new("character"),
                        refSamples= new("character"),
                        annotation = new("character"),
                        chipManu = new("character"),
                        chipVersion = new("character"),
                        alreadyAnnot = FALSE,
                        typeID = c(geneID_annotation = "SYMBOL", geneID_biomart = "hgnc_symbol", featureID_biomart = ""),
                        runAnnot = TRUE,
                        organism = "Human",
                        mart=new("Mart")
                        ) { 
    
    if (missing(dat)) {
        if (is.null(datfile(params)) | datfile(params) == "" | length(datfile(params))==0)
            stop("Either the expression matrix or the corresponding file must be provided.")
        dat <- read.table(datfile(params), header = TRUE,  stringsAsFactors= FALSE, check.names = FALSE)
    }

   

    if (missing(pData) & (length(annotfile(params))>0)) {
        if (length(grep(pattern = ".rda", x = tolower(annotfile(params))))>0)
            pData <- get(load(annotfile(params)))
        else
            pData <- read.table(annotfile(params), header = TRUE,  stringsAsFactors= FALSE, check.names = FALSE)
        
    }
    if (missing(A)) 
        A <- readA(Afile = Afile(params), dat = dat) #, sep = sep)
    

    if (missing(S))
        S <- readS(Sfile = Sfile(params), dat = dat) #, sep = sep)
    

    if(length(compNames) <= ncol(A))
        compNames <- as.character(1:ncol(A))
    
    if (length(annotation)==0 || (annotation == "")) {
        if (is.null(mart))
            stop("A 'mart' object must be defined using function 'useMart'.")
    }
    
    icaSet <- new("IcaSet",
                  A = A,
                  S = S,
                  ##witGenes = witGenes, ## to be initialized after annotation in gene
                  annotation=annotation,
                  dat=as.matrix(dat),
                  compNames=compNames,
                  chipManu = chipManu,
                  chipVersion=chipVersion,
                  refSamples=refSamples,
                  indComp = 1:ncol(A),
                  typeID = typeID,
                  mart=mart)
    if (!missing(pData) || ncol(pData)>0)
        pData(icaSet) <- pData
    
    organism(icaSet) <- organism

    
    annot2col(params) <- annot2Color(pData)   


    if (alreadyAnnot) {
        SByGene(icaSet) <- S
        datByGene(icaSet) <- dat
        runAnnot <- FALSE
    }
    else {
        SByGene(icaSet) <- datByGene(icaSet) <- data.frame()
    }
    
    
    if (runAnnot) {
        if (length(annotation(icaSet))==0 || (annotation(icaSet) == "")) {
            if (missing(mart) || is.null(mart))
                stop("A 'mart' object must be defined using function 'useMart'.")
            #warning("Since no annotation package is provided, annotation will be performed using biomaRt.")
            #stop(paste("You must provide an 'annotation' (i.e package) attribute for the annotation."))
        }
        else {
        if (!("geneID_annotation" %in% names(typeID(icaSet))))
            stop(paste("Since you want to annotate the features using ", annotation(icaSet),", typeID attribute of IcaSet object must contain an element 'geneID_annotation' which determines the object from the package you want to use."))

        if ((length(chipManu(icaSet))>1 && chipManu(icaSet) == "illumina") & (!(toupper(typeID(icaSet)["geneID_annotation"]) %in% c("ENTREZID", "SYMBOL") )))
            stop("For illumina, features can only be annotated in ENTREZID or SYMBOL. You will need to annotate yourself the IcaSet object if you want to use different IDs.")
        }

    }


    ## if no annotation package is provided, still the biomart ID of the features or genes must be provided
    if (length(annotation(icaSet))==0 || (annotation(icaSet) == "")) {
        if (length(typeID(icaSet)["geneID_biomart"]) == 0 || typeID(icaSet)["geneID_biomart"] == "" || is.na(typeID(icaSet)["geneID_biomart"]))
            if (length(typeID(icaSet)["featureID_biomart"]) == 0 || typeID(icaSet)["featureID_biomart"] == "" || is.na(typeID(icaSet)["featureID_biomart"]) )
                stop("Since no annotation package is available, you must fill at least 'featureID_biomart' or 'geneID_biomart' in typeID attribute.")
            else
                typeID(icaSet)["geneID_biomart"] <- typeID(icaSet)["featureID_biomart"]
        else
            if (length(typeID(icaSet)["featureID_biomart"]) == 0 || typeID(icaSet)["featureID_biomart"] == "" || is.na(typeID(icaSet)["featureID_biomart"]))
                typeID(icaSet)["featureID_biomart"] <- typeID(icaSet)["geneID_biomart"]
    }


    if (length(chipManu(icaSet))>1 && chipManu(icaSet) == "illumina" & runAnnot) {
        message("Since you are using illumina microarrays, additional packages are needed to complete the annotation: lumi, lumiHumanAll.db and either lumiRatIDMapping or  lumiMouseIDMapping depending on the organism. \n")
        message("...You must load packages 'lumi' and 'lumiHumanAll.db if not done yet'...")

    }
 
    ## launched even if runAnnot = FALSE in order to write rnk files if needed
    icaSet <- annotInGene(icaSet = icaSet, params = params, annot = runAnnot)
    witGenes(icaSet) <- witGenes
    
    if (length(witGenes(icaSet))==0)
        witGenes(icaSet) <- selectWitnessGenes(icaSet=icaSet, params=params, level = if (nrow(SByGene(icaSet))>0) "genes" else "features", maxNbOcc = 1, selectionByComp = NULL)


    return(list(icaSet=icaSet,params = params))    
}


##' This function selects the features having the largest Inter Quartile Range (IQR).
##'
##' @title Selection of features based on their IQR
##' @param data Measured data of dimension features x samples (e.g, gene expression data)
##' @param nb The number of features to be selected  
##' @return A subset of \code{data} restricted to the features having the \code{nb} highest IQR value
##' @author Pierre Gestraud
##' @export
##' @examples 
##' dat <- matrix(rnorm(10000),ncol=10,nrow=1000)
##' rownames(dat) <- 1:1000
##' selectFeatures_IQR(data=dat, nb=500)
##' 
selectFeatures_IQR <- 
function(data, nb){
   iqrValues <- apply(data, 1, IQR, na.rm=TRUE)
   sortedIqr <- sort(iqrValues, decreasing=TRUE)
   maxIqr <- sortedIqr[nb]
   index <- which(iqrValues >= maxIqr)
   nbGenes <- length(index)
   message("Number of selected genes is ", nbGenes, "\n")
   message("Max IQR is ", signif(maxIqr,4), "\n")   
   return(data[index,])
}
 
##' Writes the gene projection values of each component in a '.rnk' file for GSEA.
##' 
##' The .rnk format requires two columns, the first containing the gene IDs, the second containing the projection values.
##' The genes are ordered by projection values. 
##' The files are named "index-of-component_abs.rnk" if \code{abs=TRUE}, or  "index-of-component.rnk" if \code{abs=FALSE}. 
##' @title Write rnk files containing gene projections
##' @param icaSet An object of class \code{IcaSet}
##' @param abs If TRUE (default) the absolute projection values are used.
##' @param path The path that will contain the rnk files.
##' @return NULL
##' @author Anne
##' @export
writeRnkFiles <- function(icaSet, abs=TRUE, path) {
            
            dirGsea <- paste(path,"/GSEA/",sep="")
            system(paste("mkdir",dirGsea), ignore.stderr=TRUE)

            Slist <- SlistByGene(icaSet)
            
               if (abs) {
                   system(paste("mkdir ",dirGsea, "rnk_abs",sep=""), ignore.stderr=TRUE)
                   rnk_path= paste(dirGsea,"rnk_abs/",sep="")
               }
               else {
                   system(paste("mkdir ",dirGsea, "rnk",sep=""), ignore.stderr=TRUE)
                   rnk_path= paste(dirGsea,"rnk/",sep="")
               }

            if (abs)
                Slist <- lapply(Slist, function(x) {abs(x)})
            
		sapply (1:length(Slist),
                        function (indcomp, path, Slist, icaSet) {
                            comp.hugo <- Slist[[indcomp]]
			    dat = as.data.frame(comp.hugo)
			    write.table(dat,file = paste(path,"/",indComp(icaSet)[indcomp],"_",if(abs) "abs",".rnk",sep=""), sep = "\t", quote = FALSE, col.names = FALSE)
			 }
			 , path = rnk_path
                         , Slist =  Slist
                         , icaSet = icaSet
		)
	}


