Functions useful to calculate some performance measures of a survival model and the corresponding clinical utility

The R/SAS functions are:

+ internal_cv.R calculates the bootstrap optimism-corrected internal validation;

+ stdca.r calculates the net benefit and the corresponding decision curves as proposed by [Vickers et al](https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/decision-curve-analysis). The corresponding SAS function is stdca.sas.

+ Rsquared.R: calculation of R<sup>2</sup> and Royston's D as additional overall and discrimation measures to evaluate a model over the entire follow-up and not at a fixed time horizon;

+ UnoC.R: calculation of Uno's C-index as a discrimination measure for all follow-up time as proposed by [Uno et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3079915/);

+ File RCSPLINE macro.sas is a SAS macro that calculates the restricted cubic spline.







