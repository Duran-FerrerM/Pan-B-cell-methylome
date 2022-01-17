# The proliferative history shapes the DNA methylome of B-cell tumors and predicts clinical outcome

## Abstract
We report a systematic analysis of the DNA methylation variability in 1,595 samples of normal cell subpopulations and 14 tumor subtypes spanning the entire human B-cell lineage. Differential methylation among tumor entities relates to differences in cellular origin and to de novo epigenetic alterations, which allowed us to build an accurate machine learning-based diagnostic algorithm. We identify extensive patient-specific methylation variability in silenced chromatin associated with the proliferative history of normal and neoplastic B cells. Mitotic activity generally leaves both hyper- and hypomethylation imprints, but some B-cell neoplasms preferentially gain or lose DNA methylation. Subsequently, we construct a DNA methylation-based mitotic clock called epiCMIT, whose lapse magnitude represents a strong independent prognostic variable in B-cell tumors and is associated with particular driver genetic alterations. Our findings reveal DNA methylation as a holistic tracer of B-cell tumor developmental history, with implications in the differential diagnosis and prediction of clinical outcome. <br />
Here, we provide the diagnostic algorithm and the epiCMIT mitotic clock calculator.

## Graphical summary
![alt text](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/Images/Graphical.abstract.png)

## Citation
If you use any data or code derived from this study, please cite:<br />
Duran-Ferrer, M. et al. The proliferative history shapes the DNA methylome of B-cell tumors and predicts clinical outcome. Nat Cancer (2020). https://doi.org/10.1038/s43018-020-00131-2. <br />
The pdf of the article can be found [here](https://www.nature.com/articles/s43018-020-00131-2.epdf?sharing_token=XRuBq8qwGeJcZf6SuE08pdRgN0jAjWel9jnR3ZoTv0Mh5o7ypQ2fGNatxzZC0VSATkPfrfQL1kKKlFISzIABgfdmJeGPofsBj1UFSVxn5ru5tgRQoqXwF63VqsH5u33nJ-Zp1gzOEZNuXb-F6VcxSAiniSABihzhc5dJ9z5PP1M%3D). <br />
Please, also check the [comment](https://www.nature.com/articles/s43018-020-00132-1) of our manuscipt by Paolo Strati and Michael R. Green.

## LICENSE
LICENSE terms can be found [here](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/blob/master/LICENSE)

##  epiCMIT mitotic clock calculator

A complete tutorial to estimate the epiCMIT mitotic clock score in DNA methylation data can be found [here](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/tree/master/Code/Estimate.epiCMIT.html)

## Pan B-cell tumor classifier algorithm
A complete tutorial to estimate cellular composition of your DNA methylation data in Illumina arrays as well as to predict the main B-cell tumor entities ans subtypes can be found [here](https://github.com/Duran-FerrerM/Pan-B-cell-methylome/tree/master/Code/B.cell.tumor.classifier.html).

## Data availability
DNA methylation and gene expression data that support the findings of this study have been deposited at the European Genome-phenome Archive (EGA) under accession number [EGAS00001004640](https://ega-archive.org/studies/EGAS00001004640). <br />

Previously published DNA methylation data re-analyzed in this study can be found under accession codes: B cells, [EGAS00001001196](https://www.ebi.ac.uk/ega/studies/EGAS00001001196); ALL, [GSE56602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56602), [GSE49032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49032), [GSE76585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76585), [GSE69229](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69229); MCL, [EGAS00001001637](https://www.ebi.ac.uk/ega/studies/EGAS00001001637), [EGAS00001004165](https://ega-archive.org/studies/EGAS00001004165); CLL, [EGAD00010000871](https://www.ebi.ac.uk/ega/datasets/EGAD00010000871), [EGAD00010000948](https://www.ebi.ac.uk/ega/datasets/EGAD00010000948); MM, [EGAS00001000841](https://ega-archive.org/studies/EGAS00001000841); In vitro B-cell differentiation model of naïve B cells from human primary samples, [GSE72498](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72498). <br />

Normalized DNA methylation matrices used for all the analyses in this study are available [here](http://resources.idibaps.org/paper/the-proliferative-history-shapes-the-DNA-methylome-of-B-cell-tumors-and-predicts-clinical-outcome). <br />

Published gene expression datasets can be found under the accession codes: B cells, [EGAS00001001197](https://www.ebi.ac.uk/ega/studies/EGAS00001001197); ALL, [GSE49032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49032); MCL, [GSE36000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36000); CLL, [EGAS00000000092](https://www.ebi.ac.uk/ega/studies/EGAS00000000092), [EGAD00010000252](https://www.ebi.ac.uk/ega/datasets/EGAD00010000252); MM, [GSE19784](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19784); In vitro B-cell differentiation model of naïve B cells from human primary samples, [GSE72498](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72498).<br />
ChIP-seq datasets that were re-analyzed here can be found under the accession codes: [GSE109377](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109377) (NALM6 ALL cell line, n=1) and [EGAS00001000326](https://www.ebi.ac.uk/ega/studies/EGAS00001000326) (15 normal B cells donors, and 5 MCL, 7 CLL and 4 MM patients) available from [Blueprint](https://www.blueprint-epigenome.eu/).

## Contact
If you have any question, comment or suggestions please contact me at: *maduran@clinic.cat* :-)
