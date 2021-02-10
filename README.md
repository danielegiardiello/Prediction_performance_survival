# Assessing prediction performance of survival models
The following repository illustrates how to develop and validate a risk prediction models with survival outcome.
Different guidelines are available according to the different level of statistical education, experience and interests of the user according to the [STRATOS](https://www.stratos-initiative.org/) initiative.

The collaboration with STRATOS topic [groups](https://www.stratos-initiative.org/groups) 6 and 8 was essential to provide this work.

The repository (will) contain the following code:  

+ [01_predsurv_simplified](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_simplified.md) illustrates how to develop and validate a risk prediction model with survival outcomes when both development and validation data are available. Basically, the prediction performance of survival models are evaluated at a fixed time horizon (for example at 5 years). The RMarkdown source code (.Rmd) is [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_simplified.Rmd).People with basic or low statistical knowledge are encouraged to use these files. 

+ [02_presurv](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/02_predsurv.md) illustrates how to validate a risk prediction model with survival outcomes when model equation of a developed risk prediction model (using Cox regression) is available and validation data are available. The RMarkdown source code (.Rmd) is [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/02_predsurv.Rmd).

+ [03_predsurv_extended](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/03_predsurv_extended.md) is a more extensive illustration of [01_predsurv_simplified](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_simplified.md) for users with more advanced statistical knowledge. The prediction performance of survival models is not only evaluated at a fixed time horizon but also further assessments of the entire follow-up are provided.
The RMarkdown source code (.Rmd) is [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/03_predsurv_extended.Rmd).

R/Rstudio and packages version information used to generate the documentation and code are available [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/sessionInfo.md).

Data, external functions and figures are available in the corresponding subfolders.
A general SAS code will be also available.  

