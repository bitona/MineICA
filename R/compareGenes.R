##' Compute and annotate the intersection or union between contributiong genes of components originating from different IcaSet objects.
##'
##' @title Union and intersection of contributing genes
##' @param icaSets List of \code{IcaSet} objects, the \code{geneNames} of the \code{IcaSet} objects must be from the same type (e.g, gene Symbols).
##' @param keepCompByIcaSet Indices of the components to be considered in each \code{IcaSet}.
##' @param lab The names of the icaSets (e.g the names of the datasets they originate from).
##' @param cutoff The cutoff (on the absolute centered and scaled projections) above which the genes have to be considered.
##' @param type \code{"intersection"} to restrict the list of genes to the ones that are common between all datasets, or \code{"union"} to consider all the union of genes available across the datasets.
##' @param annotate If TRUE (default) the genes are annotated using function \code{writeGenes}.
##' @param file The HTML file name where the genes and their annotations are written, default is \code{type}Genes_lab1-i_lab2-j_... where i and j are the component indices contained in \code{keepCompByIcaSet}.
##' @param mart The mart object (database and dataset) used for annotation, see function \code{useMart} of package \code{biomaRt}.
##' @return A data.frame containing \describe{  
##' \item{\code{typeID(icaSets[[1]])['geneID_biomart']}:}{the gene IDs,}
##' \item{median_rank}{the median of the ranks of each gene across the \code{IcaSet} objects,} 
##' \item{analyses}{the labels of the \code{IcaSet} objects in which each gene is above the given \code{cutoff}} 
##' \item{min_rank}{the minimum of the ranks of each gene across the \code{IcaSet} objects,}
##' \item{ranks}{the ranks of each gene in each \code{IcaSet} where it is available,}
##' \item{scaled_proj}{the centered and reduced projection of each gene in each \code{IcaSet} where it is available.}}
##' @author Anne Biton
##' @export
##' @seealso \code{\link[MineICA]{writeGenes}}
##' @examples \dontrun{
##' data(icaSetCarbayo)
##' mart <- useMart("ensembl", "hsapiens_gene_ensembl")
##'
##' ## comparison of two components
##' ## here the components come from the same IcaSet for convenience
##' ## but they must come from different IcaSet in practice.
##' compareGenes(keepCompByIcaSet = c(9,4), icaSets = list(icaSetCarbayo, icaSetCarbayo),
##'              lab=c("Carbayo", "Carbayo2"), cutoff=3, type="union",  mart=mart)
##' 
##' }
##' 
compareGenes <- function(keepCompByIcaSet, icaSets, lab,  cutoff=0, type = c("union","intersection"), annotate=TRUE, file, mart=useMart("ensembl", "hsapiens_gene_ensembl")) {

    type <- match.arg(type)
    
    if (missing(file) || is.null(file))        
        file <- paste(type,"Genes_", paste(lab,"-",unlist(keepCompByIcaSet),collapse="_",sep=""),sep="",collapse="")

    
    paths <- lab

    icaSet <- p <- NULL
    comps <- 
        foreach(icaSet=icaSets, comp=keepCompByIcaSet) %dopar% {
            comp <- SByGene(icaSet)[,comp]
            names(comp) <- geneNames(icaSet)
            return(comp)
        }

    names(comps) <- paths

    comps <- llply(comps,function(x) (x-mean(x))/sd(x))

    #### intersection of the available genes
    inter <-names(comps[[1]])
    for (comp in comps)
        inter <- intersect(inter,names(comp))


    compsinter <- llply(comps,function(x,inter) x[inter], inter = inter)
    #plot(as.data.frame(compsinter))
            
    if (!missing(cutoff) && !is.null(cutoff)) {
        interCut <- names(comps[[1]][abs(comps[[1]])>cutoff])
        for (comp in comps) interCut <- intersect(interCut,names(comp)[abs(comp)>cutoff])
    }
    ### use the signs of the genes in the intersection to determine the signs of the genes in the union
    ## I use the signs in one analysis and assume that the signs are the same in the other analysis...

    if (FALSE) {
    intersigns <- sign(comps[[1]][inter])
    if (type == "union") {
        #### union of the available genes
        union <-c()
        union2signs <- c()
        for (comp in comps) {
          if (sum(sign(comp[names(intersigns[intersigns==1])]))<0) comp <- -comp
          union <- unique(c(union,names(comp)))
          union2signs <- c(union2signs,sign(comp))
        }
        print(length(union))
    }
    else if (type == "intersection") {
        if (is.null(cutoff))
            stop("When intersection is used, cutoff should be specified")
        union2signs <- c()
        for (comp in comps) {
          if (sum(sign(comp[names(intersigns[intersigns==1])]))<0) comp <- -comp
          union2signs <- c(union2signs,sign(comp[interCut]))
        }
        union <- interCut
              
    }
}

        switch(type,
              "intersection"={union <- interCut},
              "union"={union <- c();for (comp in comps) {union <- unique(c(union,names(comp)))}}
              )
    compRank <- foreach(comp=comps) %dopar% {
        r <- rank(-abs(comp))
    }
    names(compRank) <- paths #if(!is.null(path2name)) path2name[paths] else basename(paths)



    ## add analysis where the gene is present
    union2an <- 
    foreach (p=paths, comp=comps, .combine = cbind) %dopar% {
            union2an <- match(union,names(comp))
            union2an[!is.na(union2an)] <- p
            union2an[is.na(union2an)] <- ""
            names(union2an) <- union
            return(union2an)

    }
    union2nban <- apply(union2an,1,function(x) {x=as.character(x);x=x[x != ""]; length(x)})
    
    union2an <- (apply(union2an,1,paste,collapse=",",sep=""))


    ## add rank details where the gene is present
    union2proj <- 
    foreach (p=paths, comp=comps, .combine = cbind) %dopar% {
            union2p <- signif(comp[match(union,names(comp))],2)
            names(union2p) <- union
            #union2p[!is.na(union2p)] <- comp[names(union2p[!is.na(union2p)])]
            union2p[is.na(union2p)] <- ""
            #names(union2p) <- union
            return(union2p)

    }
    union2proj <- (apply(union2proj,1,paste,collapse=",",sep=""))

    ## add rank details where the gene is present
    union2rank <- as.data.frame(llply(compRank, function(comp,union) comp[union], union = union))
    union2rankdetails <-apply(union2rank,1,paste,collapse=",",sep="")
    union2rankdetails <-gsub(union2rankdetails,pattern="NA",replacement="")
    union2rankmin <- apply(union2rank, 1,min, na.rm = TRUE)
    union2rankmedian <- apply(union2rank, 1,median, na.rm = TRUE)
    names(union2rankmin) <- union
    names(union2rankmedian) <- union

    if (type=="intersection") 
        d <- data.frame(min_rank = as.character(union2rankmin), median_rank = union2rankmedian, ranks = union2rankdetails, scaled_proj = union2proj, row.names = union)
    else {
        d <- data.frame(analyses = union2an, min_rank = as.character(union2rankmin), median_rank = union2rankmedian, ranks = union2rankdetails, scaled_proj = union2proj, row.names = union)
            d$analyses <- gsub(d$analyses,pattern="^,,",replacement="", perl = TRUE)
        d$analyses <- gsub(d$analyses,pattern="^,",replacement="", perl = TRUE)
        d$analyses <- gsub(d$analyses,pattern=",,,$",replacement="", perl = TRUE)
        d$analyses <- gsub(d$analyses,pattern=",,$",replacement="", perl = TRUE)
        d$analyses <- gsub(d$analyses,pattern=",,",replacement=",")
        d$nbAn <- union2nban
    }

    d <- d[order(abs(d$median_rank), decreasing = FALSE),]
    if (nrow(d)>1000) d <- d[1:1000,]

    d$ranks <- gsub(d$ranks,pattern="^,,",replacement="", perl = TRUE)
    d$ranks <- gsub(d$ranks,pattern="^,",replacement="", perl = TRUE)
    d$ranks <- gsub(d$ranks,pattern=",,,$",replacement="", perl = TRUE)
    d$ranks <- gsub(d$ranks,pattern=",,$",replacement="", perl = TRUE)
    d$ranks <- gsub(d$ranks,pattern=",,",replacement=",")
    d$scaled_proj <- gsub(d$scaled_proj,pattern="^,",replacement="", perl = TRUE)
    d$scaled_proj <- gsub(d$scaled_proj,pattern=",,,$",replacement="", perl = TRUE)
    d$scaled_proj <- gsub(d$scaled_proj,pattern=",,$",replacement="", perl = TRUE)
    d$scaled_proj <- gsub(d$scaled_proj,pattern=",,",replacement=",")

    
    
    ### annotate genes and order by mean rank value
    #if (FALSE) {

    if (annotate) {
        writeGenes(data = d,
               filename = file,
               mart = mart, 
               typeId = typeID(icaSets[[1]])["geneID_biomart"],
               typeRetrieved = NULL,
               sortBy = "median_rank",
               colAnnot = NULL,
               decreasing = FALSE
               )
    }
    
    ## attribute meanRank to genes in the union
    return(d)


}

