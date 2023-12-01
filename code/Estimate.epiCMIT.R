##
## epiCMIT version v.2. Illumina 450k, EPIC and NGS supported.
##

# last update 31/05/2021


options(stringsAsFactors = F,error=NULL)
##needed packages
library(GenomicRanges)
library(data.table)



## Although not mandatory, minimum knowledge of R and GenomicRanges package is expected to run these lines. Please, contact me if you need help! ;-)

##load necessary info
download.file("https://github.com/Duran-FerrerM/Pan-B-cell-methylome/raw/master/data/Estimate.epiCMIT.RData", destfile = "Estimate.epiCMIT.RData", method="libcurl")
load("Estimate.epiCMIT.RData")

##
## Working functions to calculate the epiCMIT 
##


# ------------------------------------------------
# DNAm.to.epiCMIT()
# this function converts a DNA methylation matrix (DNAm) into a suitable GRange object to calculate the epiCMIT with epiCMIT() function.
# 
# ARGUMENTS
# DNAm --> DNA methylation matrix. Should be a matrix, data.frame or GRanges object.
# DNAm.genome.assembly --> The genome assembly of the DNAm matrix. Supported versions are "hg19" and "hg38".
# map.DNAm.to --> Map the data to Illumina-450k epiCMIT-CpGs or WGBS epiCMIT-CpGs. Should be "Illumina.450K.epiCMIT" or "WGBS.epiCMIT"
# CpGs.as.DNAm.rownames --> logical. Are rownames CpG names? If so they will be used for mapping. If FALSE, coordinates should be provided.
# If DNAm is of class matrix and CpGs.as.DNAm.rownames=FALSE, CpG coordinates and DNAm info columns from DNAm should be provided:
# DNAm.chr=NULL
# DNAm.start=NULL
# DNAm.end=NULL
# DNAm.meth=NULL
# min.epiCMIT.CpGs=800, integer. Recommended a minimum of 800 epiCMIT-CpGs, less CpGs have not been thoroughly tested. 
# ------------------------------------------------

# ------------------------------------------------
# epiCMIT()
# This function calculates the epiCMIT after running DNAm.to.epiCMIT() function.
# ARGUMENTS:
# DNAm.epiCMIT = The GRange object after running DNAm.to.epiCMIT() function.
# return.epiCMIT.annot = logical indicating whether to export the epiCMIT-CpGs used and metadata info.
# export.results = logical indicating wether to export the results.
# export.results.dir = character indicating the path to export the results.
# export.results.name = character string with the name for the exported file.
# ------------------------------------------------


##
## Calculate epiCMIT in example Illumina 450K data from Duran-Ferrer 2020, Nature Cancer
## 

head(Illumina.450k.hg19.example)
DNAm.epiCMIT <- DNAm.to.epiCMIT(DNAm = Illumina.450k.hg19.example,
                         DNAm.genome.assembly = "hg19",
                         map.DNAm.to = "Illumina.450K.epiCMIT",
                         min.epiCMIT.CpGs = 800 # minimum recommended
                         )
DNAm.epiCMIT

##calculate epiCMIT
epiCMIT.Illumina <- epiCMIT(DNAm.epiCMIT = DNAm.epiCMIT,
                             return.epiCMIT.annot = FALSE,
                             export.results = TRUE,
                             export.results.dir = ".",
                             export.results.name = "Illumina.450k.example_"
                             )
head(epiCMIT.Illumina$epiCMIT.scores)
epiCMIT.Illumina$epiCMIT.run.info


##
## Example with RRBS-SE hg19 data
##

RRBS.SE.hg19.Examples

epiCMIT.RRBS <- lapply(names(RRBS.SE.hg19.Examples),function(sample.i){
  
  ##transform data to run epiCMIT() function.
  DNAm.epiCMIT <- DNAm.to.epiCMIT(DNAm = RRBS.SE.hg19.Examples[[sample.i]],
                                  DNAm.genome.assembly = "hg19",
                                  map.DNAm.to = "WGBS.epiCMIT",
                                  min.epiCMIT.CpGs = 800 # minimum recommended
  )
  
  ##calculate epiCMIT
  epiCMIT.results <- epiCMIT(DNAm.epiCMIT = DNAm.epiCMIT,
                             return.epiCMIT.annot = TRUE,
                             export.results = TRUE,
                             export.results.dir = ".",
                             export.results.name = paste0("RRBS-SE.",sample.i,"_example_")
  )
})

epiCMIT.RRBS.scores <- do.call(rbind,lapply(epiCMIT.RRBS,function(x){x[["epiCMIT.scores"]]}))
epiCMIT.RRBS.scores$epiCMIT.CpGs <- as.numeric(do.call(rbind,lapply(epiCMIT.RRBS,function(x){x[["epiCMIT.run.info"]][["epiCMIT.CpGs"]]})))
epiCMIT.RRBS.scores$epiCMIT.hyper.CpGs <- as.numeric(do.call(rbind,lapply(epiCMIT.RRBS,function(x){x[["epiCMIT.run.info"]][["epiCMIT.hyper.CpGs"]]})))
epiCMIT.RRBS.scores$epiCMIT.hypo.CpGs <- as.numeric(do.call(rbind,lapply(epiCMIT.RRBS,function(x){x[["epiCMIT.run.info"]][["epiCMIT.hypo.CpGs"]]})))
epiCMIT.RRBS.scores[,-1]

fwrite(epiCMIT.RRBS.scores,"epiCMIT.RRBS.scores.tsv",sep="\t")



