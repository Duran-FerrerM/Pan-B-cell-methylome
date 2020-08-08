# The proliferative history shapes the DNA methylome of B-cell tumors and predicts clinical outcome

## Abstract
We report a systematic analysis of the DNA methylation variability in 1,595 samples of normal cell subpopulations and 14 tumor subtypes spanning the entire human B-cell lineage. Differential methylation among tumor entities relates to differences in cellular origin and to de novo epigenetic alterations, which allowed us to build an accurate machine learning-based diagnostic algorithm. We identify extensive patient-specific methylation variability in silenced chromatin associated with the proliferative history of normal and neoplastic B cells. Mitotic activity generally leaves both hyper- and hypomethylation imprints, but some B-cell neoplasms preferentially gain or lose DNA methylation. Subsequently, we construct a DNA methylation-based mitotic clock called epiCMIT, whose lapse magnitude represents a strong independent prognostic variable in B-cell tumors and is associated with particular driver genetic alterations. Our findings reveal DNA methylation as a holistic tracer of B-cell tumor developmental history, with implications in the differential diagnosis and prediction of clinical outcome.  
For further details, please check the original manuscript: https://www.biorxiv.org/content/10.1101/2020.02.06.937383v3

## Citation
If you use the Pan B-cell tumor classifier or the epiCMIT mitotic clock calculator, please cite : https://www.biorxiv.org/content/10.1101/2020.02.06.937383v3

## Graphical summary
![alt text](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/Graphical.abstract.png)


## Code availability
### Pan B-cell tumor classifier algorithm:
Current version: v.1.0.  
The source R file can be found here: https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/B.cell.tumor.classifier.R  
Here you can find additional details on how the classifier was built: https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/predictor.code.R  

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

### epiCMIT mitotic clock calculator: 
Current version: v.1.0.  
The R source file can be found here: https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/Estimate.epiCMIT.R

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
## Data availability
All the normalized DNA methylation matrices used for the sutdy can be found here: http://resources.idibaps.org/paper/the-proliferative-history-shapes-the-dna-methylome-of-b-cell-tumors-and-predicts-clinical-outcome

## Contact
If you have any question, comment or suggestions please contact me at: *maduran@clinic.cat* :-)
