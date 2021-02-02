# The proliferative history shapes the DNA methylome of B-cell tumors and predicts clinical outcome

## Abstract
We report a systematic analysis of the DNA methylation variability in 1,595 samples of normal cell subpopulations and 14 tumor subtypes spanning the entire human B-cell lineage. Differential methylation among tumor entities relates to differences in cellular origin and to de novo epigenetic alterations, which allowed us to build an accurate machine learning-based diagnostic algorithm. We identify extensive patient-specific methylation variability in silenced chromatin associated with the proliferative history of normal and neoplastic B cells. Mitotic activity generally leaves both hyper- and hypomethylation imprints, but some B-cell neoplasms preferentially gain or lose DNA methylation. Subsequently, we construct a DNA methylation-based mitotic clock called epiCMIT, whose lapse magnitude represents a strong independent prognostic variable in B-cell tumors and is associated with particular driver genetic alterations. Our findings reveal DNA methylation as a holistic tracer of B-cell tumor developmental history, with implications in the differential diagnosis and prediction of clinical outcome.  

## Graphical summary
![alt text](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/Graphical.abstract.png)

## Citation
If you use any data or code derived from this study, please cite:<br />
Duran-Ferrer, M., Clot, G., Nadeu, F. et al. The proliferative history shapes the DNA methylome of B-cell tumors and predicts clinical outcome. Nat Cancer (2020). https://doi.org/10.1038/s43018-020-00131-2. The pdf of the article can be found [here](https://www.nature.com/articles/s43018-020-00131-2.epdf?sharing_token=XRuBq8qwGeJcZf6SuE08pdRgN0jAjWel9jnR3ZoTv0Mh5o7ypQ2fGNatxzZC0VSATkPfrfQL1kKKlFISzIABgfdmJeGPofsBj1UFSVxn5ru5tgRQoqXwF63VqsH5u33nJ-Zp1gzOEZNuXb-F6VcxSAiniSABihzhc5dJ9z5PP1M%3D). <br />
Please, also check the [comment](https://www.nature.com/articles/s43018-020-00132-1) of our manuscipt by Paolo Strati and Michael R. Geen.

## LICENSE
LICENSE terms can be found [here](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/LICENSE)

## Code availability

### epiCMIT mitotic clock calculator 
#### Illumina 450k and EPIC arrays currently supported. An epiCMIT based on Next Generation Sequencing will be introduced very soon. Please, feel free to contact me if you have NGS-based data.
The epiCMIT (epigenetically-determined Cumulative MIToses) mitotic clock represents a relative measure of the total proliferative history of normal and neoplastic B cells. It is build considering the highest score from their two underlying hyper- and hypomethylation-based mitotic clocks, called epiCMIT-hyper and the epiCMIT-hypo, respectively. The code we provide here calculates the three mitotic clocks for each sample. All of them range from 0 to 1, depending on low or high realtive proliferative history. Based on our data, considering hyper- or hypomethyaltion separately may not be sufficient to capture the entire mitotic history of cancer cells. We performed a comprehensive selection of CpGs to buid the epiCMIT and showed it as an accurate mitotic clock for normal and neoplastic B cells. Nonetheless, given our careful CpG filtering we strongly belive that epiCMIT represent a pan-cancer mitotic clock.  
  
As a final note, the proliferative history of B-cell tumors include the proliferative history of normal B-cell development and the proliferative history of malignant transformation and progression. The epiCMIT should be compared then among B-cell tumors with the same cellular origin.
    
Current version: v.1.0.  
The source R file can be found [here](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/Estimate.epiCMIT.R)

#### Code:
```
##
## epiCMIT calculation example
##

# Version v.1.0, last update 07/08/2020

# download a file necessary data load it in R, and delete the file
download.file("https://github.com/Duran-FerrerM/Pan-B-cell-methylome/raw/master/Estimate.epiCMIT.RData", destfile = "Estimate.epiCMIT.RData", method="libcurl")
load("Estimate.epiCMIT.RData")
file.remove("Estimate.epiCMIT.RData")

# Description of the function Estimate.epiCMIT()
# This function calculates the epiCMIT, epiCMIT-hyper and epiCMIT-hypo in B-cell tumors as published in Duran-Ferrer, M, 2020.

#INPUT
# Function arguments
#betas: DNA methylation matrix of samples to calculate the epiCMIT. CpGs on rows, samples on columns. Rownames of matrix should be names of CpGs. See betas.example.
#epiCMIT.annot: Annotation necessary for calculation of epiCMIT. epiCMIT have been built with 450K aray data.
#export: Whether to export or not results in the current directory. Default to FALSE


# OUTPUT
# A data.frame with epiCMIT,epiCMIT-hyper and epiCMIT-hypo in all the samples. If export=T, export an .xlsx file in yout current directory with epiCMIT results.

# Example execution:
Estimate.epiCMIT(betas = betas.example,epiCMIT.annot = epiCMIT.annot,export = T)
```

### Pan B-cell tumor classifier algorithm
B-cell tumors comprise a variety of neoplasias derived from normal B cells of the heamopoietic system. Collectively, ALL, MCL, CLL, DLBCL and MM represent the majority of diagnosed B-cell neoplasias. They are further classified in subtypes with different clinicobiological features. The code we provide here represents a Pan B-cell tumor classifier algorithm, which contains two steps: in the firt step, an unknown sample is classified into one of the five major B-cell tumor entities; in the second step, these major entities are further classified into their subtypes. If prior knowledge of the entity is avaialble, we recommend to specify it to the classifier algorithm.

Current version: v.1.0.  
The source R file can be found [here](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/B.cell.tumor.classifier.R)  
Additional details on how the classifier was built can be found [here](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/predictor.code.R)  

#### Code:

```
##
## B-cell tumor classifier example
##

# Version v.1.0 last update 07/08/2020.

# download a file with the classifier, load it in R, and delete the file
download.file("https://github.com/Duran-FerrerM/Pan-B-cell-methylome/raw/master/B.cell.tumor.classifier.RData", destfile = "B.cell.tumor.classifier.RData", method="libcurl")
load("B.cell.tumor.classifier.RData")
file.remove("B.cell.tumor.classifier.RData")

# Description of the function methy.classifier()
# This function predicts the entity/subtype of a tumor sample based on reduced CpG set published in Duran-Ferrer, M, 2020.

# Output
# List of tables (one for each predictor) with the estimated probabilities of the classes, and
# a final table with the predicted entity/subtype for each sample. The output includes a raw svm prediction and a
# suggested prediction that takes into account confusion among entities/subtypes. If only one predictor is applied (see arguments) the ouptut
# only includes a table with the predictions of that specific predictor

# Function arguments
# data: data.frame or matrix with the methylation beta values and named rows(samples) and columns(cpgs)
# which.predictor: character string that specifies if the full prediction (entity+subtype) or which one of the five predictors will be applied.
# Possible values are "entity.subtype","entity","ALL","CLL","DLBCL" or "MCL". Defaults to "entity.subtype" (full prediction)
# impute.missings: impute missing beta values? Default to FALSE. Use with caution, the effect of missings to the predictions
# have not been extensively tested
# export: export results to a .xlsx file? The file will be named "methy.classifier.results.xlsx" and stored in the working directory.

# example

methy.classifier(data=example)
methy.classifier(data=example,which.predictor="MCL",export=TRUE)
methy.classifier(data=example,which.predictor="entity.subtype",export=TRUE)
```

## Data availability
DNA methylation and gene expression data that support the findings of this study have been deposited at the European Genome-phenome Archive (EGA) under accession number [EGAS00001004640](https://ega-archive.org/studies/EGAS00001004640). <br />

Previously published DNA methylation data re-analyzed in this study can be found under accession codes: B cells, [EGAS00001001196](https://www.ebi.ac.uk/ega/studies/EGAS00001001196); ALL, [GSE56602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56602), [GSE49032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49032), [GSE76585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76585), [GSE69229](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69229); MCL, [EGAS00001001637](https://www.ebi.ac.uk/ega/studies/EGAS00001001637), [EGAS00001004165](https://ega-archive.org/studies/EGAS00001004165); CLL, [EGAD00010000871](https://www.ebi.ac.uk/ega/datasets/EGAD00010000871), [EGAD00010000948](https://www.ebi.ac.uk/ega/datasets/EGAD00010000948); MM, [EGAS00001000841](https://ega-archive.org/studies/EGAS00001000841); In vitro B-cell differentiation model of naïve B cells from human primary samples, [GSE72498](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72498). <br />

Normalized DNA methylation matrices used for all the analyses in this study are available [here](http://resources.idibaps.org/paper/the-proliferative-history-shapes-the-DNA-methylome-of-B-cell-tumors-and-predicts-clinical-outcome). <br />

Published gene expression datasets can be found under the accession codes: B cells, [EGAS00001001197](https://www.ebi.ac.uk/ega/studies/EGAS00001001197); ALL, [GSE49032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49032); MCL, [GSE36000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36000); CLL, [EGAS00000000092](https://www.ebi.ac.uk/ega/studies/EGAS00000000092), [EGAD00010000252](https://www.ebi.ac.uk/ega/datasets/EGAD00010000252); MM, [GSE19784](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19784); In vitro B-cell differentiation model of naïve B cells from human primary samples, [GSE72498](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72498).<br />
ChIP-seq datasets that were re-analyzed here can be found under the accession codes: [GSE109377](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109377) (NALM6 ALL cell line, n=1) and [EGAS00001000326](https://www.ebi.ac.uk/ega/studies/EGAS00001000326) (15 normal B cells donors, and 5 MCL, 7 CLL and 4 MM patients) available from [Blueprint](https://www.blueprint-epigenome.eu/).

## Contact
If you have any question, comment or suggestions please contact me at: *maduran@clinic.cat* :-)
