# A list of contributing gene projection values per component.  Each element of the list corresponds to a component and is restricted to the genes exceeding a given threshold.
# laisser writeGostatsHtmltable

##' doEnrichment
##' This internal function is called by \code{hypergeoAn} and runs hypergeometric tests through function \code{hyperGTest} to associate the contributing genes of a component to gene sets.
##'
##' 
##' @title Runs enrichment analysis of contributing genes 
##' @param compSel A list containing three elements
##' \describe{
##' \item{compSel}{the projection values of contributing genes that were selected based on their absolute projection}
##' \item{compSel_neg}{the projection values of contributing genes that have negative projections}
##' \item{compSel_pos}{the projection values of contributing genes that have positive projections}}
##' @param chip The annotation package 
##' @param onto A string specifying the GO ontology to use. Must be one of 'BP', 'CC', or 'MF', see \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}. Only used when argument \code{db} is 'GO'.
##' @param hgCutoff The p-value threshold
##' @param cond A logical indicating whether the calculation should conditioned on the GO structure, see \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}. 
##' @param universe The universe for the hypergeometric tests, see \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}.
##' @param path The path where results will be saved
##' @param db The used database to use ('GO' or 'KEGG')
##' @param pack.annot.EID The name of the environment of the annotation package containing the annotation for Entrez Gene. 
##' @param Slist The list of gene projections across all components
##' @param it The index of the component
##' @param cutoff The threshold applied on the gene projections, used to select the contributing genes
##' @param entrez2symbol A vector of all gene Symbols involved in the analysis indexed by their Entrez Gene IDs. It is only used when \code{annotation(params)} is empty, and allows to associate gene sets to Symbols. 
##' @return Object of class \code{GOHyperGResult-class}
##' @keywords internal
##' @author Anne Biton
doEnrichment <- function  (compSel,
                           chip,
                           onto,
                           hgCutoff,
                           cond,
                           universe,
                           path,
                           db,
                           pack.annot.EID,
                           Slist,
                           it,
                           cutoff,
                           entrez2symbol
                           ) { 	

    db <- tolower(db)
    ## if the component has previously been annotated and decomposed into 3 vectors, compSel would be a list
    if (!is.list(compSel)) {
        ## Annote and Select negative and positive projections in two different vectors. 
        compSel_neg =  na.omit(unique(unlist(AnnotationDbi::mget(names(compSel[which(compSel < 0)]), pack.annot.EID, ifnotfound = NA))))
        compSel_pos =  na.omit(unique(unlist(AnnotationDbi::mget(names(compSel[which(compSel > 0)]), pack.annot.EID, ifnotfound = NA))))
        compSel = na.omit(unique(unlist(AnnotationDbi::mget(names(compSel), pack.annot.EID, ifnotfound = NA))))
    }
    else {
        compSel_neg = compSel$compSel_neg
        compSel_pos = compSel$compSel_pos
        compSel = compSel$compSel
    }
    
    system(paste("mkdir",path), ignore.stderr=TRUE)
    file_resComp <- paste(path,it,"_",hgCutoff,"cond",as.character(cond), sep = "", collapse = "")

    #First we treat genes with negative and positive projections at the same time and then we treat separately the both groups
    compSel_list = list(both = compSel, left = compSel_neg, right = compSel_pos)

    hg.res = list(); length(hg.res) = length(compSel_list); names(hg.res) = names(compSel_list)
    for (j in names(compSel_list)) {
        geneIds =  as.character(compSel_list[[j]])
        if (length(geneIds)>1) {

          if (db == "go") {
                   params <- new("GOHyperGParams", geneIds = geneIds,
                   universeGeneIds = universe, annotation = chip,
                   ontology = onto, pvalueCutoff = hgCutoff, conditional = cond,
                   testDirection = "over")
           }
           else if (db == "kegg") {
                   params <- new ("KEGGHyperGParams", geneIds = geneIds,
                   universeGeneIds = universe, annotation = chip,
                   pvalueCutoff = hgCutoff, 
                   testDirection = "over")

           }
           else if (db == "pfam") {
                   params <- new ("PFAMHyperGParams", geneIds = geneIds,
                   universeGeneIds = universe, annotation = chip,
                   pvalueCutoff = hgCutoff, 
                   testDirection = "over")
          }

          hgOver <- hyperGTest(params)

          # Add Genes implicated in selection
          hg.res[[j]] <- addGenesToGoReport(hgOver = hgOver, universe = if (missing(entrez2symbol) || is.null(entrez2symbol)) names(Slist[[1]]) else universe, db = db, onto = onto, annotation = chip, entrez2symbol = if (missing(entrez2symbol)) NULL else entrez2symbol) ## because when biomart was used for annotation (entrez2symbol provided) universe corresponds to entrez ID indexed by the feature ID of the IcaSet and these Entrez IDs will be used in report, while Slist  names corresponds to feature ID like probe sets that will be annotated using the annotation package "chip"
         }
         else hg.res[[j]] <- data.frame()

         if (!is.null(file_resComp)) {
              file = paste(file_resComp,"_",j,".htm",sep="")
              writeGostatsHtmltable(d = hg.res[[j]], label = paste("component",it), side = j, db = db, file= file, cutoff = cutoff) 
         }

        }

        return(hg.res)
}

##' Runs an enrichment analysis of the contributing genes associated with each component, using the function \code{hyperGTest} of package \code{\link[GOstats]{GOstats}}.
##' The easiest way to run enrichment analysis is to use function \code{\link{runEnrich}}. 
##'
##' An annotation package must be available in \code{annotation(icaSet)} to provide the contents of the gene sets. If none corresponds to the technology you deal with, please choose the org.*.eg.db package according to the organism (for example org.Hs.eg.db for Homo sapiens).
##' Save results of the enrichment tests in a '.rda' file located in \code{path}/\code{db}/\code{onto}/\code{zvalCutoff(params)}.
##' @title Runs an enrichment analysis per component using package \code{\link{GOstats}}.
##' @param icaSet An object of class \code{IcaSet}
##' @param params An object of class \code{\link[MineICA:MineICAParams-class]{MineICAParams}} containing the parameters of the analysis
##' @param path The path where results will be saved
##' @param SlistSel A list of contributing gene projection values per component.  Each element of the list corresponds to a component and is restricted to the features or genes exceeding a given threshold. If missing, is computed by the function.
##' @param hgCutoff The p-value threshold
##' @param db The database to be used (\code{"GO"} or \code{"KEGG"})
##' @param onto A character specifying the GO ontology to use. Must be one of \code{"BP"}, \code{"CC"}, or \code{"MF"}, see \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}. Only used when argument \code{db} is \code{"GO"}.
##' @param cond A logical indicating whether the calculation should conditioned on the GO structure, see \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}. 
##' @param universe The universe for the hypergeometric tests, see \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}.
##' @param entrez2symbol A vector of all gene Symbols involved in the analysis indexed by their Entrez Gene IDs. It is only used when \code{annotation(params)} is empty, and allows to associate gene sets to Symbols.
##' @author Anne Biton
##' @seealso \code{\link{runEnrich}}, \code{\link[xtable]{xtable}}, \code{\link[biomaRt]{useMart}}, \code{\link[Category]{hyperGTest}}, \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}, \code{\link{mergeGostatsResults}}
##' @examples \dontrun{
##' ## load an example of IcaSet 
##' data(icaSetCarbayo)
##'
##' ## define params
##' # Use threshold 3 to select contributing genes.
##' # Results of enrichment analysis will be written in path 'resPath(params)/GOstatsEnrichAnalysis'
##' params <- buildMineICAParams(resPath="~/resMineICACarbayo/", selCutoff=3)
##'
##' ## Annotation package for IcaSetCarbayo is hgu133a.db.
##' # check annotation package
##' annotation(icaSetCarbayo)
##' 
##' ## Define universe, i.e the set of EntrezGene IDs mapping to the feature IDs of the IcaSet object.
##' universe <- as.character(na.omit(unique(unlist(AnnotationDbi::mget(featureNames(icaSetCarbayo),
##'                          hgu133aENTREZID, ifnotfound = NA)))))
##' 
##' ## Apply enrichement analysis (of the contributing genes) to the first components using gene sets from KEGG.
##' # Since an annotation package is available, we don't need to fill arg 'entrez2symbol'.
##' # run the actual enrichment analysis
##' hypergeoAn(icaSet=icaSetCarbayo[,,1], params=params, db="GO",onto="BP", universe=universe)
##' }
##' @export
hypergeoAn <- function ( icaSet,
                        params,
                        path = paste(resPath(params),"GOstatsEnrichAnalysis/",sep="/"),
                        SlistSel,
                        hgCutoff = 0.01,
                        db = "go",
                        onto = "BP", 
                        cond = TRUE,
                        universe,
                        entrez2symbol
                        
	) {

    message("..Add Hugo Symbols to GOstats report..") 
    db <- tolower(db)

    pack.annot <- annotation(icaSet)
    selCutoff <- params["selCutoff"]
        
    ## 1. Selection of the probesets on the component
    if (missing(SlistSel)) 
        SlistSel = selectContrib(Slist(icaSet), cutoff = selCutoff)
    

    ## 2. Annotation of the selected probesets using Entrez ids
    if (length(pack.annot) > 0 && pack.annot != "" &&  substr(pack.annot, start = 1, stop = 3) != "org")  {
        pack.annot.EID <- eval(as.name(paste(gsub(".db", "", pack.annot), "ENTREZID", sep = "")))
        library(pack.annot, character.only = TRUE)
    }

    ## 3. Universe = all probe sets available on the component
    if (missing(universe))
        universe <- as.character(na.omit(unique(unlist(AnnotationDbi::mget(featureNames(icaSet), pack.annot.EID, ifnotfound = NA)))))
    
    ## Always tests for over-representation of a category
    hg.res = list()
    testDir = "over"
    selCutoff <- params["selCutoff"]
    if (length(selCutoff)==1)
        selCutoff <- rep(selCutoff,length(SlistSel))
    
    ## For each component, make the enrichment tests
    hg.res <- llply (c(1:length(SlistSel)),
                         function(it, sel, chip, onto, hgCutoff, cond, universe, path, db, pack.annot.EID, Slist,  selCutoff, entrez2symbol, indComp) {      
                             #print(it)
                             compSel <- sel[[it]]
                             doEnrichment(compSel = compSel
                                          , chip = chip
                                          , onto = onto
                                          , hgCutoff = hgCutoff
                                          , cond = cond
                                          , universe = universe
                                          , path = path
                                          , db = db
                                          , pack.annot.EID = pack.annot.EID
                                          , Slist = Slist
                                          , it = indComp[it]
                                          , cutoff = selCutoff[it]
                                          , entrez2symbol = entrez2symbol
                                          )
                         }
                         , sel = SlistSel
                         , chip = annotation(icaSet)
                         , onto = onto
                         , hgCutoff = hgCutoff
                         , cond = cond
                         , universe = universe
                         , path = path
                         , db = db
                         , pack.annot.EID = pack.annot.EID
                         , Slist = Slist(icaSet)
                         , selCutoff = selCutoff
                         , entrez2symbol = if (missing(entrez2symbol)) NULL else entrez2symbol
                         , indComp = indComp(icaSet) 
                         )
        names(hg.res) <- indComp(icaSet) 
	save(hg.res, file = paste(path,"hgres.rda", sep = "", collapse = ""))
	gc()
	return(hg.res)

}


## @examples \dontrun{
## ## load an example of output of function hyperGTest
## data(hgOver)
## ## look at the summary of this object
## summary(hgOver)
## ## add Gene Symbols to this summary
## addGenesToGoReport(hgOver=hgOver,db="KEGG",annotation="hgu133a.db",
##                    universe=unique(unlist(geneIdUniverse(hgOver))))
## }

##'  Add gene Symbols contained in gene sets selected as significant by \code{\link[Category]{hyperGTest}} function
##'
##' This function takes as inputs the outputs of \code{\link[Category]{hyperGTest}} which takes Entrez Gene IDs as inputs to perform the enrichment analysis.
##' The goal of this function is to select the Entrez Gene IDs responsible for the significant enrichment of a given gene set and annotate them in to gene Symbol IDs.
##' When the annotation package  \code{annotation} was used to map feature IDs to Entrez Gene ID, it can also be used here to map Entrez and Symbol IDs. 
##' If the annotation package was not used, but the Entrez Gene IDs were directly provided to the hyperGtest function, \code{annotation} is expected to be NULL and \code{entrez2symbol} must be specified.
##' 
##' This function returns the outputs of function \code{\link[Category]{hyperGTest}} which contain: 
##' \describe{
##' \item{DB, ID, Term}{The database, the gene set ID, and the gene Set name,}
##' \item{P-value}{probability of observing the number of genes annotated for the gene set among the selected gene list, knowing the total number of annotated genes among the universe}, 
##' \item{Expected counts}{expected number of genes in the selected gene list to be found at each tested category term/gene set,}
##' \item{Odds ratio}{odds ratio for each category term tested which is an indicator of the level of enrichment of genes within the list as against the universe,} 
##' \item{Counts}{number of genes in the selected gene list which are annotated for the gene set,}
##' \item{Size}{number of genes from the universe annotated for the gene set.}}
##' 
##' @title Add Symbol IDs to hyperGTest results
##' @param hgOver Output of function \code{\link[Category]{hyperGTest}}
##' @param universe A vector including all IDs on which enrichment analysis was applied
##' @param db The database to use, default is c("GO","KEGG")
##' @param onto A string specifying the GO ontology to use. Must be one of "BP", "CC", or "MF", see \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}. Only used when argument \code{db} is "GO".
##' @param annotation An annotation package
##' @param entrez2symbol A vector indexed by Entrez Gene ID and filled with the corresponding Gene Symbols
##' @return A data.frame containing the summary of the output of function hyperGTest (\code{summary(hgOver)}) with an additional column providing the gene Symbols included in the significant gene sets.
##' @seealso \code{\link[Category]{hyperGTest}}, \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}
##' @author Anne Biton
##' @keywords internal
`addGenesToGoReport` <-
function(hgOver, universe, db = c("GO","KEGG"), onto = c("CC", "MF", "BP"), annotation = NULL, entrez2symbol = NULL) {


  db <- match.arg(tolower(db), choices = c("go","kegg"))
  onto <- match.arg(toupper(onto), choices = c("CC", "MF", "BP"))

  a <- GOstats::geneIdsByCategory(hgOver)
  b <- GOstats::geneIdUniverse(hgOver, cond=conditional(hgOver))

  a <- a[sigCategories(hgOver)]
  b <- b[sigCategories(hgOver)]

  pack.annot <- annotation
  if (!is.null(pack.annot) & pack.annot != "" &  substr(pack.annot, start = 1, stop = 3) != "org") {
      pack.annot.EID <- eval(as.name(paste(gsub(".db", "", pack.annot), "ENTREZID", sep="")))
      pack.annot.SYMBOL <- eval(as.name(paste(gsub(".db", "", pack.annot), "SYMBOL", sep = "")))
  
      if (tolower(db) == "go") {
	  pack.annot.ALLPROBES <- eval(as.name(paste(gsub(".db", "", pack.annot), "GO2ALLPROBES", sep="")))
	  if (onto=="CC") {
	    probes <- AnnotationDbi::mget('GO:0005575', pack.annot.ALLPROBES)
	  } else if (onto=="BP") {
	    probes <- AnnotationDbi::mget('GO:0008150', pack.annot.ALLPROBES)
	  } else if (onto=="MF") {
	    probes <- AnnotationDbi::mget('GO:0003674', pack.annot.ALLPROBES)
	  }
	  entrez <- AnnotationDbi::mget(unique(unlist(probes)), pack.annot.EID)

      }
      else if (tolower(db) == "kegg") {
  	  entrez <- unlist(AnnotationDbi::as.list(pack.annot.EID))
          # If Kegg gene sets are used, summary function of "HyperGResult-class" doesn't work because of the ".db"
	  hgOver@annotation = gsub(".db", "", pack.annot)
      }


      pbset <- lapply( a, function (genes_a, entrez) { return(unique(names(unlist(entrez[entrez%in%genes_a])))) }, entrez = entrez)
      pbset <- lapply( pbset, intersect, universe)

      gs <- lapply(pbset, function(pbset, pack.annot.SYMBOL) { return(unique(unlist(AnnotationDbi::mget(pbset, pack.annot.SYMBOL)))) }, pack.annot.SYMBOL = pack.annot.SYMBOL)

  }
  else {
      if (is.null(entrez2symbol))
          stop("Mapping between Hugo Symbols and Entrez Gene ID is missing and must be provided when argument 'annotation' is empty.")
      
      
      pbset <- lapply( a, function (genes_a, entrez) { return(unique(unlist(entrez[entrez %in% genes_a]))) }, entrez = names(entrez2symbol))
      pbset <- lapply( pbset, intersect, universe)

      gs <- lapply(pbset, function(pbset, entrez2symbol) { return(unique(unlist(entrez2symbol[pbset]))) }, entrez2symbol = entrez2symbol)
  }

  hg_tmp = AnnotationDbi::summary(hgOver)

  if (db == "go") {
      hg_tmp$In_geneSymbols <- unlist(lapply(gs[hg_tmp[[paste("GO",onto,"ID",sep="")]]], paste, collapse = ","))
  }
  else if (db == "kegg") {
      hg_tmp$In_geneSymbols = unlist(lapply(gs[hg_tmp$KEGGID], paste, collapse = ","))
  }

  return(hg_tmp)
    

}

##' This function takes as input in argument \code{d} the output of function \code{\link{addGenesToGoReport}} whose goal is to add genes included in gene sets detected as significantly enriched by \code{\link[Category]{hyperGTest}} function.
##' It writes the enrichment results in an HTML file which redirects each gene set ID to its web-description and each gene to its Gene Card web-page.   
##'
##' @title Writes enrichment results in a HTML file
##' @param d A data.frame describing enrichment results, output of function \code{\link[Category]{hyperGTest}}
##' @param label The label of the data the results originate from
##' @param side The side of the component used for enrichment analysis
##' @param db The database used ("GO" or "KEGG")
##' @param file File name for output
##' @param cutoff The threshold used to select the genes used to run the enrichment analysis
##' @return NULL
##' @seealso \code{\link[xtable]{xtable}}, \code{\link{addGenesToGoReport}}, \code{\link[Category]{hyperGTest}}
##' @author Anne Biton
##' @keywords internal
##' @examples
##'
##' hgOver <- structure(list(GOBPID = c("GO:0003012", "GO:0030049"),
##'                     Pvalue = c(1.70848789161935e-10, 6.62508415367712e-05),
##'                     OddsRatio = c(22.1043956043956, 26.4190476190476),
##'                     ExpCount = c(1.19549929676512, 0.246132208157525),
##'                     Count = c(12L, 4L), Size = c(68L, 14L),
##'                     Term = c("muscle system process", "muscle filament sliding"),
##'                     In_geneSymbols = c("ACTA2,ACTC1,ACTG2,CASQ2,CNN1,DES,MYH3,MYLK,PTGS1,TPM2,MYL9,LMOD1","ACTC1,DES,MYH3,TPM2")), 
##'                     .Names = c("GOBPID", "Pvalue", "OddsRatio", "ExpCount", "Count", "Size", "Term", "In_geneSymbols"),
##'                     class = "data.frame", row.names=1:2)
##' 
##' MineICA:::writeGostatsHtmltable(d=hgOver, label="Example of enrichment analysis", db="KEGG",
##'                                 file="outputHyper_example.htm")
##' 
writeGostatsHtmltable <- function(d, label, side="both", db, file, cutoff=3) {

    if (side == "both")
        signgenes <- "both positive and negative"
    else
        if (side == "right")
            signgenes <- "positive"
        else
            if (side == "left")
                signgenes <- "negative"

    style <- 
    "
    <head>
     <style type='text/css'>

         a.normal:link{text-decoration:none;color:black;}
         a.normal:visited{text-decoration:none;color:purple;}
         a.normal:hover{text-decoration:none;font-weight:bold;color:black;}

         a.sign:link{text-decoration:none;color:blue;font-weight:bold;}
         a.sign:visited{text-decoration:none;color:#0000CD;font-weight:bold;}
         a.sign:hover{text-decoration:none;font-weight:bold;color:#00EE00;}

         TH.richTable {
            padding: 2px;
                border:2px solid;
                color: white;
                border-color: #6495ED;
                text-align: center;
                font-weight: bold;
                background-color: #C6E2FF
         }

     </style>
    </head>"

    if (nrow(d)>0) {

     #col 1 = KEGGID GOBPID ....
	    if (tolower(db) == "go") {
	          d$Term <- paste("<a  class='normal' href='http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=",d[,1],"'>", d$Term,"</a>",sep="")
	          d[,1] <- paste("<a  class='normal' href='http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=",d[,1],"'>", d[,1],"</a>",sep="")

	    }
	    else if (tolower(db) == "kegg") {
	          d$Term <-  paste("<a  class='normal' href='http://www.genome.jp/dbget-bin/www_bget?map",d[,1],"'>",d$Term,"</a>",sep="")
	          d[,1] <- paste("<a  class='normal' href='http://www.genome.jp/dbget-bin/www_bget?map",d[,1],"'>",d[,1],"</a>",sep="")
	    }

	    ## add gene cards link
	    d$In_geneSymbols <- sapply(d$In_geneSymbols, 
	                            function(x) {
	                                 ssplit <- strsplit(x,split=",")[[1]]
	                                 paste(paste("<a  class='normal' href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=",ssplit,"'>",ssplit,"</a>",sep=""),collapse = ",")
	                                 })


	    dnames <- colnames(d)
	    dnames <- c(dnames[1],"Term","Pvalue",dnames[!(dnames %in% c(dnames[1],"Term","Pvalue"))])
	    d <- d[,dnames]
	    d$Pvalue <- as.character(signif(d$Pvalue,4))

            
	    x <- xtable(d, caption = paste("<b>Significant <a  class='normal' href='http://bioinformatics.oxfordjournals.org/cgi/content/short/23/2/257'>GOstats</a> results for ", label, ".</b><br> The listed results were computed using <b>",signgenes," genes projections</b> whose scaled values exceed the threshold ",cutoff,".",sep=""),digits = 3)
	    x <- capture.output(print(x,
                                      type = "html",
                                      sanitize.text.function = force,
                                      include.rownames=FALSE,
                                      caption.placement = "top")
                                )

	    x <- paste(style,x,sep="")          
	    x <- gsub(x,pattern="<TABLE",replacement = "<TABLE style='font-family:Helvetica; font-size:12px; border-top:solid thin black; border=2px;' ", ignore.case = TRUE)
	    x <- gsub(x,pattern="<CAPTION",replacement="<CAPTION  style='color:black; text-align:left; font-size:14px;' ")
	    x <- gsub(x,pattern="<TH>",replacement = "<TH  class='richTable'>", ignore.case = TRUE)

    }
    else x = paste("<b>No Significant <a  class='normal' href='http://bioinformatics.oxfordjournals.org/cgi/content/short/23/2/257'>GOstats</a> results for ", label, ".</b><br> The listed results were computed using <b>",signgenes," genes projections</b> whose scaled values exceed the threshold ",cutoff,".",sep="")
    write(x,file=file)
    
}




####### WARNING : 1. I assume that the same cutoffs are available for all DB
#######           2. I assume that there are only GO and KEGG tested, if not change the part "Construction of the matrix gene sets -> components"
######            3. We won't have pvalue for gene sets which are not significant on a component, because the threshold is directly given to GOstats
######            4. hgCutoff and cond have to be provided because they will be used in the file names of the results summary files
######               Because if one new analysis on the same data is run, the results will be put in the same place
######               I chose not to delete the previous file and to precise with parameters have been used to compute each summary file.
##' This function is internal and called by function \code{runEnrich}. It merges enrichment results obtained with either KEGG, GO, or both databases into one file.
##' 
##' This function writes an HTML file per component, containing the outputs of the enrichment tests computed through the function \code{\link[Category]{hyperGTest}}.
##' The results of the enrichment tests are loaded from .rda files located in \code{resPath(icaSet)}/GOstatsEnrichAnalysis/byDb/'db-name'/('ontology-name'/).
##' The results obtained for the different databases/ontologies are then merged into an array for each component, this array is written as an HTML file in the directory \code{resPath(icaSet)}/\code{zvalCutoff(params)}.
##' The arguments \code{hgCutoff} and \code{cond} have to be provided because they will be used in the file names of the resulting files. 
##' 
##' This function makes several important assumptions: only databases GO and KEGG have been tested, p-values are not available for gene sets that have not been selected as significant.
##' 
##' The outputs of \code{\link[Category]{hyperGTest}} that are given in each table are: 
##' \describe{
##' \item{DB, ID, Term}{The database, the gene set ID, and the gene set name,}
##' \item{P-value}{probability of observing the number of genes annotated for the gene set among the selected gene list, knowing the total number of annotated genes among the universe}, 
##' \item{Expected counts}{expected number of genes in the selected gene list to be found at each tested category term/gene set,}
##' \item{Odds ratio}{odds ratio for each category term tested which is an indicator of the level of enrichment of genes within the list as against the universe,} 
##' \item{Counts}{number of genes in the selected gene list which are annotated for the gene set,}
##' \item{Size}{number of genes from the universe annotated for the gene set.}}
##' 
##' @title Merge enrichment results obtained for different databases into one file per component. 
##' @param resPath The global path where results of ICA analysis are written 
##' @param GOstatsPath The path within argument \code{resPath} where files will be written
##' @param rdata The name of the rdata file containing the enrichment analysis of all components
##' @param cutoffs The threshold(s) used to select genes used in enrichment analysis
##' @param hgCutoff The p-value threshold 
##' @param cond A logical indicating whether the calculation has been conditioned on the GO structure, see \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}. 
##' @param pathGenes The path where HTML files containing gene projections for each component are located
##' @seealso  \code{\link[xtable]{xtable}}, \code{\link[biomaRt]{useMart}}, \code{\link[Category]{hyperGTest}}, \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}, \code{\link{hypergeoAn}}, \code{\link{mergeGostatsResults}}
##' @return NULL
##' @keywords internal
##' @author Anne Biton
mergeGostatsResults <- function (resPath, GOstatsPath, rdata="hgres", cutoffs = NULL, hgCutoff = 0.01, cond = TRUE, pathGenes) {

        d <- NULL
        GOstatsPath <- gsub(x=GOstatsPath,pattern = resPath,replacement = "")
        pathRacine <- paste(resPath,GOstatsPath,"/",sep="")
        
        path <- paste(pathRacine,"/byDB/",sep="")

        pathLink <- "byDB/"
        
        ### ******* ------------------------------------------------------
        ### Extraction of the DB available in the directory and the cutoffs available for each DB
        db <- list.files(path = path,  recursive = FALSE, full.names = TRUE)
        db <- gsub(db, pattern = pathRacine, replacement = "")
        dblink <- list.files(path = path, recursive = FALSE, full.names = FALSE)
        dblink <- sapply(dblink, function(x) paste(x,pathLink,dblink,sep=""))
        
        if (length(db) == 0) print(paste("No functional analysis is available for", path))
        ## if go is one of the available db, one more level has to be examined
        indgo <- grep(db, pattern = "/go$", ignore.case = TRUE)
        if (length(indgo) > 0) {
          # if go in one of the tested databases
          # GO categories could have been tested (BP, CC, MF)
          # and would correspond to directories in the GO/ path
          onto <- system(paste("ls ", pathRacine, db[indgo], sep = ""), intern = TRUE)
          dbgo <- paste(db[indgo], "/",onto, "/", sep = "")
          dbgobis <- paste(dblink[indgo],"/", onto, "/", sep = "")
          db <- c(db[-indgo],dbgo)
          dblink <- c(dblink[-indgo],dbgo)
          
        }

        
style <- 
"
<head>
 <style type='text/css'>

     a.normal:link{text-decoration:none;color:black;}
     a.normal:visited{text-decoration:none;color:purple;}
     a.normal:hover{text-decoration:none;font-weight:bold;color:black;}

     a.sign:link{text-decoration:none;color:blue;font-weight:bold;}
     a.sign:visited{text-decoration:none;color:#0000CD;font-weight:bold;}
     a.sign:hover{text-decoration:none;font-weight:bold;color:#00EE00;}

     a.comp:link{text-decoration:none;color:white;}
     a.comp:visited{text-decoration:none;color:white;}
     a.comp:hover{text-decoration:none;font-weight:bold;color:black;}


     TH.richTable {
        padding: 2px;
            border:2px solid;
            color: white;
            border-color: #6495ED;
            text-align: center;
            font-weight: bold;
            background-color: #C6E2FF
     }

 </style>
</head>"

        rdataName <- load(paste(pathRacine,db[1],"/",rdata,".rda",sep=""))

        indComp <- c(1:length(eval(as.name(rdataName))))
        rm(rdataName)
        


        path2db <-
            sapply(db, function(path) {
                if (length(grep(pattern="/go/",x=path,ignore.case = TRUE))>0) {
                    onto <- strsplit(path,split="/")[[1]]
                    db <- onto[length(onto)-1]
                    onto <- onto[length(onto)]
                    return(paste(db,"/",onto,sep=""))
                }
                else {
                    db <- basename(path)
                    return(db)
                }
            })
        names(path2db) <- db


        db2comppath <- lapply(db,
                              function(dir, indComp, cutoffs, hgCutoff, cond, resPath, GOstatsPath, pathRacine) {
                                 dir <- paste(pathRacine,"/",dir,"/", sep="")
                                 listFiles <- list.files(path = dir, pattern = paste("*_",hgCutoff,"cond",as.character(cond),"*",sep=""), full.names = TRUE, ignore.case = TRUE)
                                 listFiles <- sapply(listFiles,gsub, pattern = "///", replacement = "/")
                                 listFiles <- sapply(listFiles, gsub, pattern = "////", replacement = "/")
                                 listFiles <- sapply(listFiles, gsub, pattern = "//", replacement = "/")
                                 listFiles <- sapply(listFiles, gsub, pattern = pathRacine, replacement = "")
                                 names(listFiles) <- NULL

                                 trueindComp <- unique(sapply(listFiles, function(x) strsplit(basename(x),split="_")[[1]][1]))
                                 ## order components by increasing index
                                 ordcomp <- order(as.numeric(trueindComp))
                                 trueindComp <- trueindComp[ordcomp]
                                 
                                 filesbycomp <- 
                                  lapply(trueindComp,function(comp,dir,listFiles, hgCutoff, cond) {
                                      listFilesComp <- listFiles[grep(listFiles,pattern=paste(comp,"_",hgCutoff,"cond",as.character(cond),sep=""))]
                                      listFilesComp <- gsub(listFilesComp, pattern = GOstatsPath, replacement = "")
                                      listFilesComp <- gsub(listFilesComp, pattern = resPath, replacement = "") 
                                      listFilesComp <- gsub(listFilesComp, pattern = pathRacine, replacement = "") 

                                      names(listFilesComp) <- NULL
                                      
                                       return(list(both=listFilesComp[grep(listFilesComp,pattern="_both")],
                                                   left=listFilesComp[grep(listFilesComp,pattern="_left")],
                                                   right=listFilesComp[grep(listFilesComp,pattern="_right")]))
                                      
                                  }
                                         ,dir=dir
                                         ,listFiles = listFiles
                                         ,hgCutoff = hgCutoff
                                         ,cond = cond
                                    )
                                   names(filesbycomp) <- trueindComp

                                 return(filesbycomp)
                               }
                              , indComp = indComp
                              , cutoffs = cutoffs
                              , hgCutoff = hgCutoff
                              , cond = cond
                              , resPath = resPath
                              , GOstatsPath = GOstatsPath
                              , pathRacine = pathRacine
                              )
        names(db2comppath) <- db
        
        trueindComp <- names(db2comppath[[1]])

        db2res <- foreach(path = db) %dopar% {
		dd <- paste(pathRacine,path,"/",sep="")
                rdataName <- load(paste(dd,"/",rdata,".rda",sep=""))
                res <- eval(as.name(rdataName))
	  }
	  names(db2res) <- db
       
        resbyComp <-
            lapply(indComp, function(comp, db, dir, cutoff, listFilesComp, path2db, hgCutoff, cond, db2res, trueindComp) {
                comp2gs <- foreach(path = db,  .combine = rbind) %dopar% {

                                   dbid <- paste(gsub(path2db[path],pattern="/",replacement=""),"ID",sep="")
                                   dbb <- path2db[path]
                                   listFiles <- listFilesComp[[path]][[comp]]
				   res <- db2res[[path]][[comp]]$both

                                   ### ******* ------------------------------------------------------
                                   ### Construction of the matrix sgnificant gene sets -> res for the component
                                   ### ******* ------------------------------------------------------                                   
                                   if (nrow(res)>0) {
                                      res$DB <- res$DBi <- rep("",nrow(res))
                                      res$DBi[1] <- dbb

                                      ## link to the results for the db on the component using the 3 cases (enrich analysis on negative, positive, and absolute gene projections)
                                      res$DB[1] <- paste("<b>",dbb,"</b>", "<br>(", "<a  class='normal' href='",listFiles$left,"'>","neg","</a>","/",
                                                         "<a  class='normal' href='",listFiles$right,"'>","pos","</a>",")",sep="")
                                      
                                      res <- res[c("DB", dbid, "Term","Pvalue","OddsRatio","ExpCount","Count","Size","In_geneSymbols","DBi")]
                                      colnames(res)[2] <- "ID"
                                      res$Termi <- res$Term
                                      res$IDi <- res$ID
                                      
                                      if (length(grep(db,pattern="go",ignore.case=TRUE)) != 0) {
                                          res$Term <- paste("<a  class='normal' href='http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=",res$ID,"'>", res$Term,"</a>",sep="")
                                          res$ID <- paste("<a  class='normal' href='http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=",res$ID,"'>", res$ID,"</a>",sep="")
                                          
                                      }
                                      else if (tolower(db) == "kegg") {
                                          res$Term <-  paste("<a  class='normal' href='http://www.genome.jp/dbget-bin/www_bget?map",res$ID,"'>",res$Term,"</a>",sep="")
                                          res$ID <- paste("<a  class='normal' href='http://www.genome.jp/dbget-bin/www_bget?map",res$ID,"'>",res$ID,"</a>",sep="")
                                      }
                                      res$Pvalue <- as.character(signif(res$Pvalue,4))
                                      gc()
                                      return(res)
                                   }
                                   else return(NULL)
                                      
                               }

            if (length(comp2gs) !=0) {
                ##### ------------------------
                ##### Write all significant gene sets related to the component in an html file
                ##### ------------------------
                dcomp <- comp2gs[,colnames(comp2gs)[!(colnames(comp2gs) %in% c("DBi","Termi","IDi"))]]           
                ## add gene cards link
                dcomp$In_geneSymbols <- sapply(dcomp$In_geneSymbols, 
                                               function(x) {
                                                   ssplit <- strsplit(x,split=",")[[1]]
                                                   paste(paste("<a  class='normal' href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=",ssplit,"'>",ssplit,"</a>",sep=""),collapse = ",")
                                               })
                
                x <- xtable(dcomp, caption = paste("<b>Significant <a  class='normal' href='http://bioinformatics.oxfordjournals.org/cgi/content/short/23/2/257'>GOstats</a> results for component ", trueindComp[comp], "</b><br> The listed results were computed using genes projections whose absolute scaled values exceed the threshold(s) ",if (length(unique(cutoff))==1) cutoff[1] else cutoff[comp],". <br> Click on 'neg' and 'pos' below each DB name to see the results of the enrichment analysis made using  genes with negative and resp. positive projections only. <br>Click on a gene id to open the corresponding GeneCards page.",sep=""))

                x <- capture.output(print(x,
                                      type = "html",
                                      sanitize.text.function = force,
                                      include.rownames=FALSE,
                                      caption.placement = "top",
                                      digits = 6)
                                    )

                x <- paste(style,x,sep="")          
                x <- gsub(x,pattern="<TABLE",replacement = "<TABLE style='font-family:Helvetica; font-size:12px; border-top:solid thin black; border=2px;' ", ignore.case = TRUE)
                x <- gsub(x,pattern="<CAPTION",replacement="<CAPTION  style='color:black; text-align:left; font-size:14px;' ")
                x <- gsub(x,pattern="<TH>",replacement = "<TH  class='richTable'>", ignore.case = TRUE)
       
                ##### ------------------------
                ##### Return 2 vectors for each component :
                ##### 1. Pvalue indexed by gene set id (kegg or go )
                ##### 2. term indexed by gene set id 
                ##### ------------------------
                gs2pval <- comp2gs$Pvalue
                names(gs2pval) <- comp2gs$IDi
                gs2term <- comp2gs$Termi
                names(gs2term) <- comp2gs$IDi
                gs2db <- foreach(d=comp2gs$DBi, .combine=c) %do% { if (d != "") state = d; return(state)}
                names(gs2db) <- comp2gs$IDi
            
                
            }
            else {
                x = paste("<b>No results.</b>", sep="")
                gs2pval=NULL
                gs2term=NULL
                gs2db = NULL

            }
            write(x,file= paste(pathRacine,"/",trueindComp[comp],"_comp2gs_",hgCutoff,"cond",as.character(cond),"_cutoff",if (length(unique(cutoff))==1) cutoff[1] else cutoff[comp],".htm", sep = ""))

            return(list(gs2pval=gs2pval, gs2term=gs2term, gs2db = gs2db))

            
          }
          , db = db
          , dir = path
          , cutoff = cutoffs
          , listFilesComp = db2comppath
          , path2db = path2db
          , hgCutoff = hgCutoff
          , cond = cond
	  , db2res = db2res
          , trueindComp = trueindComp)
        
        ### ******* ------------------------------------------------------
        ### Construction of the matrix gene sets -> components
        ### ******* ------------------------------------------------------

        gs2pval_bycomp <- lapply(resbyComp,function(x) return(x$gs2pval))
        gs2term_bycomp <- lapply(resbyComp,function(x) return(x$gs2term))
        gs2db <- unlist(lapply(resbyComp,function(x) return(x$gs2db)))
        #gs2db <- gs2db[-duplicated(gs2db)]
        
        ## components without significant results
        noSignComp <- which(unlist(lapply(gs2pval_bycomp,is.null)))
        if (length(noSignComp)>0) {
            gs2pval_bycomp <- gs2pval_bycomp[-noSignComp]
            gs2term_bycomp <- gs2term_bycomp[-noSignComp]
            trueindComp <- trueindComp[-noSignComp]
            indComp <- indComp[-noSignComp]
        }
        
             
        if (length(gs2db)>0) {

            ## extraction of of gene sets significant at least once
            gsSign <- unique(unlist(lapply(gs2pval_bycomp, names)))
            names(gs2term_bycomp) <- NULL # to keep only geen set ids when applying unlist
            gsSign2term <- unlist(gs2term_bycomp);
            if (sum(duplicated(gsSign2term))>0)
                gsSign2term <- gsSign2term[-which(duplicated(gsSign2term))]
            
            ## in order to have the same number of elements in each component
            ## NA will be attributed if the gene set is not significant on the component
            gs2pval_bycomp <- lapply(gs2pval_bycomp, function(x,n) x[n], n = gsSign)
            names(gs2pval_bycomp) <- trueindComp ##m indComp
            ## transformation in a data.frame
            gs2comp <- as.data.frame(gs2pval_bycomp, row.names = gsSign)
            names(gs2comp) <- trueindComp #m as.character(indComp)
            gs2comp$ID <- rownames(gs2comp)
            gs2comp$DB <- gs2db[gs2comp$ID]
            gs2comp$Term <- gsSign2term[gs2comp$ID]
                                        #colnames(gs2comp)[4:(length(indComp)+3)] <- paste("<a  class='comp' href='",pathGenes,indComp,".htm'>",indComp,"</a>",sep="")

 
            ## order by database
            gs2comp <- gs2comp[order(gs2comp$DB),]
            ## transform pvalue into a link toward the complete results of the component on the database correspondign to the line
            gs2compHtml <- foreach (comp=trueindComp, .combine = cbind) %dopar% {

                comp <- as.character(comp)
                colHtml <- as.character(gs2comp[,comp])
                indNoNA <- which(!is.na(colHtml))
                                        # to avoid the covering element by element : do it by column
                db2path <- sapply(gs2comp[indNoNA,]$DB, function(db, path2db) {names(path2db)[match(db,path2db)]},path2db=path2db)
                
                filesCompB <- (sapply(paste("db2comppath[['",db2path,"']][['",comp,"']][['both']]",sep=""),function(x) eval(parse(text=x))))
                filesCompL <- (sapply(paste("db2comppath[['",db2path,"']][['",comp,"']][['left']]",sep=""),function(x) eval(parse(text=x))))
                filesCompR <- (sapply(paste("db2comppath[['",db2path,"']][['",comp,"']][['right']]",sep=""),function(x) eval(parse(text=x))))
                colHtml[indNoNA] <- paste(" <a  class='normal' href='",filesCompB,"'>",colHtml[indNoNA], "</a><br>",
                                          "(<a  class='normal' href='",filesCompL,"'>","neg", "</a>",
                                          "/<a  class='normal' href='",filesCompR,"'>","pos","</a>)",sep="")
                return(colHtml)
            }
                    
            gs2comp[,trueindComp] <- gs2compHtml
          
            gs2comp$DB[which(duplicated(gs2comp$DB))] <- ""
            gs2comp$DB[gs2comp$DB != ""] <- paste("<b>",gs2comp$DB[gs2comp$DB != ""],"</b>",sep="")
            gs2comp <- gs2comp[,c("DB","ID","Term",trueindComp)]
            names(gs2comp)[4:(length(indComp)+3)] <- paste("<a  class='comp' href='",relativePath(pathGenes, pathRacine),basename(pathGenes),"/",trueindComp,".htm'>",trueindComp,"</a>",sep="")                    



        ### **** ----------------
        ### Write gs2comp in html format
        ### Pvalue is computed using genes with both negative and positive projections, access to the results on genes with
        ### only negative and positive projections exceeding the threshold are provided by a link below each pvalue
        ### The link on pvalue results provide all the results of the component on the provided gene sets
        ### **** ----------------

        x <- xtable(gs2comp, caption = paste("<b>Pvalue of the (at least once) significant gene sets on the different components. The enrichment analysis were performed with <a  class='normal' href='http://bioinformatics.oxfordjournals.org/cgi/content/short/23/2/257'>GOstats</a>.", "</b><br> The listed results were computed using genes projections whose absolute scaled values exceed the threshold(s) ",if (length(unique(cutoffs))==1) cutoffs[1] else paste(cutoffs,collapse=","),". <br> Click on the pvalue to access to the description of the GOstats results of the component on the corresponding database. <br> Click on 'neg' and 'pos' below each pvalue to access to the results of the enrichment analysis of the component on the corresponding database made using only genes with negative and resp. positive projections only.",sep=""))
        x <- capture.output(print(x,
                                  type = "html",
                                        #file = paste(dir,"/",comp,"_comp2gs.htm", sep = ""),
                                  sanitize.text.function = force,
                                  include.rownames=FALSE,
                                  caption.placement = "top",
                                  digits = 6)
                            )

            x <- paste(style,x,sep="")          
            x <- gsub(x,pattern="<TABLE",replacement = "<TABLE style='font-family:Helvetica; font-size:12px; border-top:solid thin black; border=2px;' ", ignore.case = TRUE)
            x <- gsub(x,pattern="<CAPTION",replacement="<CAPTION  style='color:black; text-align:left; font-size:14px;' ")
            x <- gsub(x,pattern="<TH>",replacement = "<TH  class='richTable'>", ignore.case = TRUE)
            write(x,file= paste(pathRacine,"gs2comp_pvalue_",hgCutoff,"cond",as.character(cond),"_cutoff",if (length(unique(cutoffs))==1) unique(cutoffs) else paste(cutoffs,collapse="_"),".htm", sep = ""))
        }

    }
      

