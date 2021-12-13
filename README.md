# The proliferative history shapes the DNA methylome of B-cell tumors and predicts clinical outcome

## Abstract
We report a systematic analysis of the DNA methylation variability in 1,595 samples of normal cell subpopulations and 14 tumor subtypes spanning the entire human B-cell lineage. Differential methylation among tumor entities relates to differences in cellular origin and to de novo epigenetic alterations, which allowed us to build an accurate machine learning-based diagnostic algorithm. We identify extensive patient-specific methylation variability in silenced chromatin associated with the proliferative history of normal and neoplastic B cells. Mitotic activity generally leaves both hyper- and hypomethylation imprints, but some B-cell neoplasms preferentially gain or lose DNA methylation. Subsequently, we construct a DNA methylation-based mitotic clock called epiCMIT, whose lapse magnitude represents a strong independent prognostic variable in B-cell tumors and is associated with particular driver genetic alterations. Our findings reveal DNA methylation as a holistic tracer of B-cell tumor developmental history, with implications in the differential diagnosis and prediction of clinical outcome. <br />
Here, we provide the diagnostic algorithm and the epiCMIT mitotic clock calculator.

## Graphical summary
![alt text](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/Graphical.abstract.png)

## Citation
If you use any data or code derived from this study, please cite:<br />
Duran-Ferrer, M. et al. The proliferative history shapes the DNA methylome of B-cell tumors and predicts clinical outcome. Nat Cancer (2020). https://doi.org/10.1038/s43018-020-00131-2. <br />
The pdf of the article can be found [here](https://www.nature.com/articles/s43018-020-00131-2.epdf?sharing_token=XRuBq8qwGeJcZf6SuE08pdRgN0jAjWel9jnR3ZoTv0Mh5o7ypQ2fGNatxzZC0VSATkPfrfQL1kKKlFISzIABgfdmJeGPofsBj1UFSVxn5ru5tgRQoqXwF63VqsH5u33nJ-Zp1gzOEZNuXb-F6VcxSAiniSABihzhc5dJ9z5PP1M%3D). <br />
Please, also check the [comment](https://www.nature.com/articles/s43018-020-00132-1) of our manuscipt by Paolo Strati and Michael R. Green.

## LICENSE
LICENSE terms can be found [here](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/LICENSE)

## Code availability

### epiCMIT mitotic clock calculator
### NOTE: Please, make sure you have the latest R version and R packages (R≥4.1.1 and GenomicRanges≥1.44.0) to run the following code. 
#### Illumina arrays (450k and EPIC) and NGS (RRBS, ERRBS and WGBS) are currently supported.
The epiCMIT (epigenetically-determined Cumulative MIToses) mitotic clock represents a relative measure of the total proliferative history of normal and neoplastic B cells. It is built considering the highest score from their two underlying hyper- and hypomethylation-based mitotic clocks, called epiCMIT-hyper and the epiCMIT-hypo, respectively. The code we provide here calculates the three mitotic clocks for each sample. All of them range from 0 to 1, depending on low or high realtive proliferative history. Based on the data analyzed in the manuscipt, considering hyper- or hypomethyaltion separately may not be sufficient to capture the entire mitotic history of cancer cells. We performed a comprehensive selection of CpGs to buid the epiCMIT and showed it as an accurate mitotic clock for normal and neoplastic B cells. Nonetheless, given our careful CpG filtering we strongly belive that epiCMIT represent a pan-cancer mitotic clock. As a final note, the proliferative history of B-cell tumors include the proliferative history of normal B-cell development and the proliferative history of malignant transformation and progression. The epiCMIT should be compared then among B-cell tumors with the same cellular origin.
    
Current version: v.2.0.  
The source R file can be found [here](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/Estimate.epiCMIT.R)

#### Code:
```
##
## epiCMIT version v.2. Illumina 450k, EPIC and NGS supported.
##

# last update 19/04/2021


options(stringsAsFactors = F,error=NULL)
##needed packages
library(GenomicRanges)
library(data.table)


## Although not mandatory, minimum knowledge of R and GenomicRanges package is expected to run these lines. Please, contact me if you need help! ;-)

##load necessary info
download.file("https://github.com/Duran-FerrerM/Pan-B-cell-methylome/raw/master/Estimate.epiCMIT.RData", destfile = "Estimate.epiCMIT.RData", method="libcurl")
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

```

### Pan B-cell tumor classifier algorithm
#### The CpGs of the classifier are present in both Illumina 450k and EPIC arrays, and thus both platforms are supported. A first step to predict B-cell tumor content in the samples is available (strongly recommended), as the predictor assumes a minimum of 60% tumor cell content. Please, note that the tumor cell content prediction based on DNA methylation may not be accurate for some DLBCL and for MM cases, as we reported in Duran-Ferrer 2020.
B-cell tumors comprise a variety of neoplasias derived from normal B cells of the heamopoietic system. Collectively, ALL, MCL, CLL, DLBCL and MM represent the majority of diagnosed B-cell neoplasias. They are further classified in subtypes with different clinicobiological features. The code we provide here represents a Pan B-cell tumor classifier algorithm, which contains two steps: in the firt step, an unknown sample is classified into one of the five major B-cell tumor entities previously mentioned; in the second step, these major entities are further classified into their subtypes: ALL subtypes, HeH; 11q23/MLL; t(12;21); t(1;19); t(9;22); dic(9;20). MCL subtypes, C1/cMCL and C2/nnMCL. CLL subtypes, n-CLL/low programmed CLL, i-CLL/intermediate programmed CLL and m-CLL/high programmed CLL and  DLBCL subtypes, ABC and GCB. Although an unbiased prediction can serve as an internal control, if trusty knowledge of the entity is avaialble, we recommend to specify it to the classifier algorithm.

Current version: v.2.0.  
The source R file can be found [here](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/B.cell.tumor.classifier.R)  
Additional details on how the classifier was built can be found [here](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/predictor.code.R)  

#### Code:

```
##
## Pan B cell tumor classifier. v.2
##

# last update, 19/04/2021

options(max.print = 10000,stringsAsFactors = F,error=NULL)


##
## Load necessary info
##

# download a file with the classifier, load it in R, and delete the file
download.file("https://github.com/Duran-FerrerM/Pan-B-cell-methylome/raw/master/B.cell.tumor.classifier.RData", destfile = "B.cell.tumor.classifier.RData", method="libcurl")
load("B.cell.tumor.classifier.RData")


# Description of the function Deconvolute.Bcell.tumor() (which substitutes the previous function methy.classifier()).

# This function predicts the entity/subtype of a tumor sample based on reduced CpG set published in Duran-Ferrer, M, 2020.

# Output
# One table with the cellular compositions of B-cell tumor samples depending on the assumed microenvironmental cells, 
# and the others the estimated probabilities of B-cell tumor entities and subtypes, and
# a final table with the predicted entity/subtype for each sample. The output includes a raw svm prediction and a
# suggested prediction that takes into account confusion among entities/subtypes. If only one predictor is applied (see arguments) the ouptut
# only includes a table with the predictions of that specific predictor

# Function arguments
# data: data.frame or matrix with the methylation beta values and named rows(CpGs) and columns(Samples)
# predict.tumor.cell.content: logical indicating whether predict B-cell tumor content based on DNA methylation values.
# microenvironment.CpGs: Type of microenvironmental cells used for deconvolution. Bcells.2CpGs only allow the prediction of B-cell proportion.
# predict.Bcell.tumors: logical indicating whether predict B-cell tumor entities and subtypes.
# which.predictor: character string that specifies if the full prediction (entity+subtype) or which one of the five predictors will be applied.
# Possible values are "entity.subtype","entity","ALL","CLL","DLBCL" or "MCL". Defaults to "entity.subtype" (full prediction)
# impute.missings: impute missing beta values? Default to FALSE. Use with caution, the effect of missings to the predictions
# have not been extensively tested
# export: export results to a .xlsx file? The file will be named "methy.classifier.results.xlsx" and stored in the working directory.

# example
# example.betas contain differnt B-cell tumors and microenvironmental cells for illustrative purposes
Deconvolute.Bcell.tumor(data=example.betas)
Deconvolute.Bcell.tumor(data=example.betas,impute.missings = T)
Results <- Deconvolute.Bcell.tumor(data=example.betas,
                        predict.tumor.cell.content = T,
                        microenvironment.CpGs = "Pan.Bcell.microenvironment",
                        predict.Bcell.tumors = T,
                        which.predictor="entity.subtype",
                        impute.missings = T,
                        export=F
                        )

Results

cbind(Results$Cellular.proportions$Bcells,Results$combined.prediction)
```

## Data availability
DNA methylation and gene expression data that support the findings of this study have been deposited at the European Genome-phenome Archive (EGA) under accession number [EGAS00001004640](https://ega-archive.org/studies/EGAS00001004640). <br />

Previously published DNA methylation data re-analyzed in this study can be found under accession codes: B cells, [EGAS00001001196](https://www.ebi.ac.uk/ega/studies/EGAS00001001196); ALL, [GSE56602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56602), [GSE49032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49032), [GSE76585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76585), [GSE69229](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69229); MCL, [EGAS00001001637](https://www.ebi.ac.uk/ega/studies/EGAS00001001637), [EGAS00001004165](https://ega-archive.org/studies/EGAS00001004165); CLL, [EGAD00010000871](https://www.ebi.ac.uk/ega/datasets/EGAD00010000871), [EGAD00010000948](https://www.ebi.ac.uk/ega/datasets/EGAD00010000948); MM, [EGAS00001000841](https://ega-archive.org/studies/EGAS00001000841); In vitro B-cell differentiation model of naïve B cells from human primary samples, [GSE72498](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72498). <br />

Normalized DNA methylation matrices used for all the analyses in this study are available [here](http://resources.idibaps.org/paper/the-proliferative-history-shapes-the-DNA-methylome-of-B-cell-tumors-and-predicts-clinical-outcome). <br />

Published gene expression datasets can be found under the accession codes: B cells, [EGAS00001001197](https://www.ebi.ac.uk/ega/studies/EGAS00001001197); ALL, [GSE49032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49032); MCL, [GSE36000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36000); CLL, [EGAS00000000092](https://www.ebi.ac.uk/ega/studies/EGAS00000000092), [EGAD00010000252](https://www.ebi.ac.uk/ega/datasets/EGAD00010000252); MM, [GSE19784](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19784); In vitro B-cell differentiation model of naïve B cells from human primary samples, [GSE72498](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72498).<br />
ChIP-seq datasets that were re-analyzed here can be found under the accession codes: [GSE109377](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109377) (NALM6 ALL cell line, n=1) and [EGAS00001000326](https://www.ebi.ac.uk/ega/studies/EGAS00001000326) (15 normal B cells donors, and 5 MCL, 7 CLL and 4 MM patients) available from [Blueprint](https://www.blueprint-epigenome.eu/).

## Contact
If you have any question, comment or suggestions please contact me at: *maduran@clinic.cat* :-)
