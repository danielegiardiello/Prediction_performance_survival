# Assessing prediction performance of survival models
The following repository illustrates how to develop and validate a risk prediction models with survival outcome.
Different guidelines are available according to the different level of statistical education, experience and interests of the user according to the [STRATOS](https://www.stratos-initiative.org/) initiative.

The collaboration with STRATOS topic [groups](https://www.stratos-initiative.org/groups) 6 and 8 was essential to provide this work.

The repository (will) contain the following code:  

+ [01_predsurv_simplified](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_simplified.md) illustrates how to develop and validate a risk prediction model with survival outcomes when both development and validation data are available. Basically, the prediction performance of survival models are evaluated at a fixed time horizon (for example at 5 years). The RMarkdown source code (.Rmd) is [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_simplified.Rmd). People with basic or low statistical knowledge are encouraged to use these files. Another nice R code was motivated by the paper of [Royston & Altman](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-13-33) and is provided by [Martijn Heymans](https://github.com/mwheymans) [here](https://missingdatasolutions.rbind.io/2021/02/cox-external-validation/).

+ [02_presurv](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/02_predsurv.md) illustrates how to validate a risk prediction model with survival outcomes when model equation of a developed risk prediction model (using Cox regression) is available and validation data are available. The RMarkdown source code (.Rmd) is [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/02_predsurv.Rmd). 

+ [03_predsurv_extended](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/03_predsurv_extended.md) is a more extensive illustration of [01_predsurv_simplified](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_simplified.md) for users with more advanced statistical knowledge. The prediction performance of survival models is not only evaluated at a fixed time horizon but also further assessments of the entire follow-up are provided.
The RMarkdown source code (.Rmd) is [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/03_predsurv_extended.Rmd).

R/Rstudio and packages version information used to generate the documentation and code are available [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/sessionInfo.md).

Data, external functions and figures are available in the corresponding subfolders.  

SAS code is available [here](https://github.com/danielegiardiello/Prediction_performance_survival/tree/main/SAS_code).

## Usage

You can either download a zip file containing the directory, or you can clone it by using

```bash
git clone https://github.com/danielegiardiello/Prediction_performance_survival.git
```

In either case, you can then use the `Prediction_performance_survival.Rproj` file to open
and Rstudio session in the directory you have just downloaded. You may then knit
both rmarkdown files, or run them line-by-line.

## Contributions

| Name                                                         | Affiliation                           | Role                  |
| ------------------------------------------------------------ | ------------------------------------- | ----------------------|
| [Daniele Giardiello](https://github.com/danielegiardiello/)  | The Netherlands Cancer Institute (NL) | Author/maintainer     |
| [David McLernorn](https://twitter.com/davemclernon?lang=en) | University of Aberdeen (UK) | Author of SAS code              |
| [Edouard Bonneville](https://www.lumc.nl/org/bds/medewerkers/1968807) | Leiden University Medical Center (NL) | Contributor  |
| [Terry Therneau](https://www.mayo.edu/research/faculty/therneau-terry-m-ph-d/bio-00025991) | Mayo Clinic (US)| Contributor  |


