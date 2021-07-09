# Assessing prediction performance of survival models
The following repository illustrates how to develop and validate a risk prediction models with survival outcome.
Different guidelines are available according to the different level of statistical education, experience and interests of the user according to the [STRATOS](https://www.stratos-initiative.org/) initiative.

The collaboration with STRATOS topic [groups](https://www.stratos-initiative.org/groups) 6 and 8 was essential to provide this work.

The repository contains the following code:  

+ Minimal and essential [code](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_minimal.R) to develop and validate a risk prediction model with survival outcomes when both development and validation data are available. More elaborated output can be found [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_minimal_RMD.md). The corresponding .Rmd source code is [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_minimal_RMD.Rmd). People with basic or low statistical knowledge are encouraged to use these files. Another nice R code was motivated by the paper of [Royston & Altman](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-13-33) and is provided by [Martijn Heymans](https://github.com/mwheymans) [here](https://missingdatasolutions.rbind.io/2021/02/cox-external-validation/).

+ Minimal and essential [code](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/02_predsurv_minimal.R) to validate a risk prediction model in a external data when model equation of a developed risk prediction model is available.More elaborated output can be found [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/02_predsurv.md). The corresponding .Rmd source code is [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/02_predsurv.Rmd).

+ Extensive output and [code](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/03_predsurv_extended.md) to develop and validate a risk prediction model with a survival outcome. The .Rmd source code is [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/03_predsurv_extended.Rmd). People with advanced knowledge in statistics are encouraged to use these files.

External functions and figures are available in the corresponding subfolders.  

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
| [David McLernon](https://twitter.com/davemclernon?lang=en) | University of Aberdeen (UK) | Author of SAS code              |
| [Edouard Bonneville](https://www.lumc.nl/org/bds/medewerkers/1968807) | Leiden University Medical Center (NL) | Contributor  |
| [Terry Therneau](https://www.mayo.edu/research/faculty/therneau-terry-m-ph-d/bio-00025991) | Mayo Clinic (US)| Contributor  |


