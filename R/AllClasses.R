setClass(Class = "IcaSet",
         contains  = "eSet",
         representation = representation(
           A="data.frame",
           S="data.frame", 
           SByGene="data.frame", 
           compNames="character", 
           indComp="numeric",
           witGenes="character", 
           datByGene="data.frame",
           chipManu="character",
           chipVersion="character",
           refSamples="character",
           typeID="character", 
           organism = "character",
           mart="Mart"
         ),
         prototype   = prototype(new("VersionedBiobase", versions = c(classVersion("eSet"), IcaSet="0.1.0")))
         )


setClass(Class = "MineICAParams",
         representation = representation(
           Sfile="character",
           Afile="character",
           datfile="character",
           annotfile="character",
           resPath="character",
           genesPath="character",
           annot2col="character",
           pvalCutoff="numeric",
           selCutoff="numeric"
         )
         )

     
