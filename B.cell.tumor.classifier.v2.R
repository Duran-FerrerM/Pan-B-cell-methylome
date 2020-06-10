##
## B-cell tumor classifier example
##

# B-cell tumor classifier, last updated 10/06/2020.

# download a file with the classifier, load it in R, and delete the file
download.file("https://github.com/Duran-FerrerM/Pan-B-cell-methylome/raw/master/methy.classifier.RData", destfile = "methy.classifier.RData", method="libcurl")
load("methy.classifier.RData")
file.remove("methy.classifier.RData")

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

