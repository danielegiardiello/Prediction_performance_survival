# Assessing prediction performance of survival models using SAS

The corresponding SAS code created by [David McLernon](https://twitter.com/davemclernon?lang=en) are available here:

+ [01_model_development_internal_validation](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/SAS_code/01_model_development_internal_validation.sas) develops the prediction model and provides the corresponding internal bootstrap optimism-corrected validation;

+ [02_model_development_internal_validation_withPGR](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/SAS_code/02_model_development_internal_validation_withPGR.sas) develops the prediction model including the PGR marker and provides the corresponding internal bootstrap optimism-corrected validation;

+ [03_external_validation_when_original_model_dataset_available](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/SAS_code/03_external_validation_when_original_model_dataset_available.sas) provides the SAS code to externally validate the prediction model in a new dataset when the individual patient data of the development set are available;

+ [04_external_validation_when_model_estimates_and_baseline_survival_available](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/SAS_code/04_external_validation_when_model_estimates_and_baseline_survival_available.sas) provides the SAS code to externally validate the prediction model in a new dataset when only the model equation of the prediction model is provided;

Data useful to run SAS codes are also available: [development data](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/Data/rotterdam.csv) and [validation data](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/Data/gbsg.csv). 

NOTE -  SAS code requires some external macros:
+ [RCSPLINE macro](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/Functions/RCSPLINE%20macro.sas) to calculate restricted cubic splines;

+ [stdca](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/Functions/stdca.sas) to calculate net benefit and decision curves.
