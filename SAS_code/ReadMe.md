The corresponding SAS code are provided here created by [David McLernon](https://twitter.com/davemclernon?lang=en):

+ [01_model_development_internal_validation](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/SAS_code/01_model_development_internal_validation.sas) develops the prediction model and provides the corresponding internal bootstrap optimism-corrected validation;

+ [02_model_development_internal_validation_withPGRmarker](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/SAS_code/02_model_development_internal_validation_withPGRmarker.sas) develops the prediction model and provides the corresponding internal bootstrap optimism-corrected validation  including the PGR marker;

+ [03_external_validation_when_original_model_dataset_available](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/SAS_code/03_external_validation_when_original_model_dataset_available.sas) provides the SAS code to externally validate the prediction model in a new data when individual patient data are of the development set are available.

+ [04_external_validation_when_model_estimates_and_baseline_survival_available](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/SAS_code/04_external_validation_when_model_estimates_and_baseline_survival_available.sas) provided the SAS code to externally validated the prediction model in a new data when only the model equation of the prediction model is provided.

NOTE: SAS code requires some external macros:
+ [RCSPLINE macro](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/Functions/RCSPLINE%20macro.sas) to calculate restricted cubic splines;

+ [stdca](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/Functions/stdca.sas) to calculate net benefit and decision curves.
