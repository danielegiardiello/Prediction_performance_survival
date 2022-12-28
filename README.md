# Assessing performance and clinical usefulness in prediction models with survival outcomes: practical guidance for Cox proportional hazards models

R and SAS Code repository for the manuscript 'Assessing performance and clinical usefulness in prediction models with survival outcomes: practical guidance for Cox proportional hazards models' published in Annals of Internal Medicine.

**Journal version of the manuscript is [here](https://www.acpjournals.org/doi/pdf/10.7326/M22-0844).**  
**A preprint version of the manuscript is [here](https://www.medrxiv.org/content/10.1101/2022.03.17.22272411v1).**


The repository contains the following code:  

+ Minimal and essential [code](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_minimal.R) to develop and validate a risk prediction model with survival outcomes when both development and validation data are available. People with basic or low statistical knowledge and basic R programming knowledge are encouraged to use these files. **To reproduce the main results of the manuscript, this script is sufficient**.  
More elaborated output can be found [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_minimal_RMD.md). The corresponding .Rmd source code is [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_minimal_RMD.Rmd).  
Another nice R code was motivated by the paper of [Royston & Altman](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-13-33) is provided by [Martijn Heymans](https://github.com/mwheymans) [here](https://missingdatasolutions.rbind.io/2021/02/cox-external-validation/).
Minimal and essential code is also available in Python [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_minimal.py). More elaborated output in RMarkdown based on Python code is available [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_minimal_RMD_py.md) with the corresponding .Rmd source code [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_minimal_RMD_py.Rmd).

+ Minimal and essential [code](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/02_predsurv_minimal.R) to validate a risk prediction model in a external data when model equation of a developed risk prediction model is available. More elaborated output can be found [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/02_predsurv.md). The corresponding .Rmd source code is [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/02_predsurv.Rmd).

+ Extensive output and [code](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/03_predsurv_extended.md) to develop and validate a risk prediction model with a survival outcome. The .Rmd source code is [here](https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/03_predsurv_extended.Rmd). People with advanced knowledge in statistics are encouraged to use these files.

External [functions](https://github.com/danielegiardiello/Prediction_performance_survival/tree/main/Functions) and [figures](https://github.com/danielegiardiello/Prediction_performance_survival/tree/main/imgs) are available in the corresponding subfolders.  

SAS code is available [here](https://github.com/danielegiardiello/Prediction_performance_survival/tree/main/SAS_code).

## Usage

You can either download a zip file containing the directory, or you can clone it by using

```bash
git clone https://github.com/danielegiardiello/Prediction_performance_survival.git
```

In either case, you can then use the `Prediction_performance_survival.Rproj` file to open
and Rstudio session in the directory you have just downloaded. You may then knit
both rmarkdown files, or run them line-by-line.

The collaboration with STRATOS topic [groups](https://www.stratos-initiative.org/groups) 6 and 8 was essential to provide this work.

## Contributions

| Name                                                         | Affiliation                           | Role                  |
| ------------------------------------------------------------ | ------------------------------------- | ----------------------|
| [Daniele Giardiello](https://github.com/danielegiardiello/)  | The Netherlands Cancer Institute (NL) <br /> Leiden University Medical Center (NL) <br /> EURAC research (IT) | Author/maintainer and creator of R and Python code     |
| [David McLernon](https://twitter.com/davemclernon?lang=en) | University of Aberdeen (UK) | Author of SAS code              |
| [Laure Wynants](https://www.maastrichtuniversity.nl/laure.wynants) | Maastrict University (NL)  <br /> KU Leuven (BE) | Review of minimal .R/.Rmd codes |
| [Nan van Geloven](https://www.lumc.nl/org/bds/medewerkers/1216536) | Leiden University Medical Center (NL) | Review of 02_.Rmd code |
| [Maarten van Smeden](https://www.umcutrecht.nl/en/research/researchers/van-smeden-maarten-m) | University Medical Centre Utrecht (NL) |Review of 03_.Rmd code     |
| [Edouard Bonneville](https://www.lumc.nl/org/bds/medewerkers/1968807) | Leiden University Medical Center (NL) | Contributor  |
| [Terry Therneau](https://www.mayo.edu/research/faculty/therneau-terry-m-ph-d/bio-00025991) | Mayo Clinic (US)| Contributor  |


