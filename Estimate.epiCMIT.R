##
## epiCMIT calculation example
##


# download a file necessary data load it in R, and delete the file
download.file("https://github.com/Duran-FerrerM/Pan-B-cell-methylome/raw/master/epiCMIT.example.RData", destfile = "methy.classifier.RData", method="libcurl")
load("methy.classifier.RData")
file.remove("methy.classifier.RData")

# Description of the function Estimate.epiCMIT()
# This function calculates the epiCMIT, epiCMIT-hyper and epiCMIT-hypo in B-cell tumors as published in Duran-Ferrer, M, 2020.

# Output
# A data.frame with epiCMIT,epiCMIT-hyper and epiCMIT-hypo in all the samples


# Function arguments
#betas: DNA methylation matrix of samples to calculate the epiCMIT. CpGs on rows, samples on columns. Rownames of matrix should be names of CpGs. See betas.example.
#epiCMIT.annot: Annotation necessary for calculation of epiCMIT. epiCMIT have been built with 450K aray data.
#export: Whether to export or not results in the current directory. Default to FALSE

Estimate.epiCMIT(betas = betas.example,epiCMIT.annot = epiCMIT.annot,export = T)
