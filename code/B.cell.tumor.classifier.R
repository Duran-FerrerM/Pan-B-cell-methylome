##
## Pan B cell tumor classifier. v.3
##

# last update, 31/04/2022

options(max.print = 10000,stringsAsFactors = F,error=NULL)


##
## Load necessary info
##

# download a file with the classifier, load it in R, and delete the file
download.file("https://github.com/Duran-FerrerM/Pan-B-cell-methylome/raw/master/data/B.cell.tumor.classifier.RData", destfile = "B.cell.tumor.classifier.RData", method="libcurl")
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


