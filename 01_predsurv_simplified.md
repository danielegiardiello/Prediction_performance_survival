Performance assessment of survival prediction models - simplified code
================

-   [Goals](#goals)
    -   [Install/load packages and import
        data](#installload-packages-and-import-data)
    -   [Descriptive statistics](#descriptive-statistics)
-   [Goal 1: develop a risk prediction model with a time-to-event
    outcome](#goal-1-develop-a-risk-prediction-model-with-a-time-to-event-outcome)
    -   [1.1 Preliminary investigation: survival and censoring curves in
        the development and validation
        data](#preliminary-investigation-survival-and-censoring-curves-in-the-development-and-validation-data)
    -   [1.2 Secondary investigation: check non-linearity of continuous
        predictors](#secondary-investigation-check-non-linearity-of-continuous-predictors)
    -   [1.3 Model development: first check - the proportional hazard
        (PH)
        assumption](#model-development-first-check---the-proportional-hazard-ph-assumption)
    -   [1.4 Model development: fit the risk prediction
        models](#model-development-fit-the-risk-prediction-models)
-   [Goal 2: Assessing performance in survival prediction
    models](#goal-2-assessing-performance-in-survival-prediction-models)
    -   [2.1 Overall performance
        measures](#overall-performance-measures)
    -   [2.2 Discrimination measures](#discrimination-measures)
    -   [2.3 Calibration](#calibration)
    -   [2.3.1 Observed/Expected ratio](#observedexpected-ratio)
    -   [2.3.2 Calibration plot using restricted cubic
        splines](#calibration-plot-using-restricted-cubic-splines)
-   [Goal 3: Clinical utility](#goal-3-clinical-utility)
-   [References](#references)
-   [Reproducibility ticket](#reproducibility-ticket)

## Goals

In this document, we assume that individual data of the development and
validation set are both available. This file illustrates in a simplified
way how to develop a survival prediction model and how to assess the
corresponding prediction performance using internal and external
validation.

The goals are:  
1. To develop a risk prediction model with a time-to-event outcome;  
2. To assess the prediction performance of a model with a time-to-event
outcome;  
3. To assess the potential clinical utility of a risk prediction model
with time-to-event outcome;

### Install/load packages and import data

We following libraries are needed to achieve the following goals, if you
have not them installed, please use install.packages(’‘)
(e.g. install.packages(’survival’)) or use the user-friendly approach if
you are using RStudio.

``` r
# Use pacman to check whether packages are installed, if not load
if (!require("pacman")) install.packages("pacman")
library(pacman)

pacman::p_load(
  rio,
  survival,
  rms,
  mstate,
  sqldf,
  pec,
  riskRegression,
  survAUC,
  survivalROC,
  timeROC,
  plotrix,
  splines,
  knitr,
  table1,
  kableExtra,
  gtsummary,
  boot,
  tidyverse,
  rsample,
  gridExtra,
  webshot
)
# webshot::install_phantomjs()


rdata <- readRDS(here::here("Data/rdata.rds"))
vdata <- readRDS(here::here("Data/vdata.rds"))
```

We loaded the development (rdata) and the validation data (vdata) from
.rds format. The Rotterdam breast cancer data was used to predict the
risk of recurrence or death using size, stage and tumor size as
predictors. These three predictors were used in the Nottingham
Prognostic Index, one of the most popular index to determine prognosis
following surgery of breast cancer.  
The Germany Breast Cancer Study Group data was used as an external
validation of the model developed in the Rotterdam breast cancer data.
The prediction model will be then extended using the progesterone (PGR)
marker measured at primary surgery.  
The improvement in prediction performance will be evaluated internally
in the Rotterdam data (development data) and in German Breast Cancer
Study data (validation data).

### Descriptive statistics

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Characteristic
</th>
<th style="text-align:left;">
Development dataset, N = 2,982
</th>
<th style="text-align:left;">
Validation dataset, N = 686
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Size (cm)
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left; padding-left:  2em;" indentlevel="1">
&lt;=20
</td>
<td style="text-align:left;">
1,387 (47%)
</td>
<td style="text-align:left;">
180 (26%)
</td>
</tr>
<tr>
<td style="text-align:left; padding-left:  2em;" indentlevel="1">
21-50
</td>
<td style="text-align:left;">
1,291 (43%)
</td>
<td style="text-align:left;">
453 (66%)
</td>
</tr>
<tr>
<td style="text-align:left; padding-left:  2em;" indentlevel="1">
&gt;50
</td>
<td style="text-align:left;">
304 (10%)
</td>
<td style="text-align:left;">
53 (7.7%)
</td>
</tr>
<tr>
<td style="text-align:left;">
Number of nodes
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left; padding-left:  2em;" indentlevel="1">
0
</td>
<td style="text-align:left;">
1,436 (48%)
</td>
<td style="text-align:left;">
0 (0%)
</td>
</tr>
<tr>
<td style="text-align:left; padding-left:  2em;" indentlevel="1">
1-3
</td>
<td style="text-align:left;">
764 (26%)
</td>
<td style="text-align:left;">
376 (55%)
</td>
</tr>
<tr>
<td style="text-align:left; padding-left:  2em;" indentlevel="1">
4+
</td>
<td style="text-align:left;">
782 (26%)
</td>
<td style="text-align:left;">
310 (45%)
</td>
</tr>
<tr>
<td style="text-align:left;">
Grade of tumor
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left; padding-left:  2em;" indentlevel="1">
1-2
</td>
<td style="text-align:left;">
794 (27%)
</td>
<td style="text-align:left;">
525 (77%)
</td>
</tr>
<tr>
<td style="text-align:left; padding-left:  2em;" indentlevel="1">
3+
</td>
<td style="text-align:left;">
2,188 (73%)
</td>
<td style="text-align:left;">
161 (23%)
</td>
</tr>
<tr>
<td style="text-align:left;">
Age (years)
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left; padding-left:  2em;" indentlevel="1">
Mean (SD)
</td>
<td style="text-align:left;">
55 (13)
</td>
<td style="text-align:left;">
53 (10)
</td>
</tr>
<tr>
<td style="text-align:left; padding-left:  2em;" indentlevel="1">
Median (Range)
</td>
<td style="text-align:left;">
54 (24, 90)
</td>
<td style="text-align:left;">
53 (21, 80)
</td>
</tr>
<tr>
<td style="text-align:left;">
PGR (ng/mL)
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left; padding-left:  2em;" indentlevel="1">
Mean (SD)
</td>
<td style="text-align:left;">
162 (291)
</td>
<td style="text-align:left;">
110 (202)
</td>
</tr>
<tr>
<td style="text-align:left; padding-left:  2em;" indentlevel="1">
Median (Range)
</td>
<td style="text-align:left;">
41 (0, 5,004)
</td>
<td style="text-align:left;">
32 (0, 2,380)
</td>
</tr>
</tbody>
</table>

## Goal 1: develop a risk prediction model with a time-to-event outcome

Prediction models are useful to provide the estimated probability of a
specific outcome using personal information. In many studies, especially
in medicine, the main outcome under assessment is the time to an event
of interest defined generally as survival time. Prognostic models for
survival end points, such as recurrence or progression of disease, need
to account for drop out during follow-up. Patients who have not
experienced the event of interest are censored observations. Cox
regression analysis is the most popular statistical model to deal with
such data in oncology and other medical research.

### 1.1 Preliminary investigation: survival and censoring curves in the development and validation data

First, we draw the survival and the censoring curves of the development
and validation data

``` r
# Development set
sfit1 <- survfit(Surv(ryear, status == 1) ~ 1, data = rdata) # survival
sfit2 <- survfit(Surv(ryear, status == 0) ~ 1, rdata) # censoring

par(xaxs = "i", yaxs = "i", las = 1)
plot(sfit1, conf.int = FALSE, lwd = 2, xlab = "Years", bty = "n")
lines(sfit2, conf.int = FALSE, col = 2, lwd = 2)
legend(11, .9, c("Death", "Censoring"), col = 1:2, lwd = 2, bty = "n")
title("Development set")
```

<img src="imgs/01_predsurv_simplified/surv-1.png" width="672" style="display: block; margin: auto;" />

``` r
# Validation set
sfit3 <- survfit(Surv(ryear, status == 1) ~ 1, data = vdata) # survival
sfit4 <- survfit(Surv(ryear, status == 0) ~ 1, vdata) # censoring

par(xaxs = "i", yaxs = "i", las = 1)
plot(sfit3, conf.int = FALSE, lwd = 2, xlab = "Years", bty = "n", xlim = c(0, 8))
lines(sfit4, conf.int = FALSE, col = 2, lwd = 2)
legend("bottomleft", c("Death", "Censoring"), col = 1:2, lwd = 2, bty = "n")
title("Validation set")
```

<img src="imgs/01_predsurv_simplified/surv-2.png" width="672" style="display: block; margin: auto;" />
A number of 2982 patients were included to develop the risk prediction
model for survival with a median follow-up of 9 years. The median
survival in the development data was 8 years with the corresponding 95%
confidence intervals (CIs) of 7 and 9 years. The 5-year survival was 59%
(95% CI: 58-61%). A number of 686 patients were selected to externally
validate the risk prediction model.The median survival in the validation
data was 4 years. The median survival was 5 years while the 5-year
survival was 49% (95% CI: 45-54%).

### 1.2 Secondary investigation: check non-linearity of continuous predictors

The potential non-linear relation between continuous predictors
(i.e. progesterone level) and the outcome should be investigated before
developing a risk prediction model. Non-linearity of continuous
predictors can be checked using splines.  
Physically, a spline is a flexible wood or metal strip, which is passed
through a set of fixed points (knots) in order to approximate a curve.
The most common computational approximation to this is a cubic smoothing
spline which is cubic between the knot points, and constrained to be
linear beyond the two end knots. For the restricted cubic spline using
rms::rcs() R package::function(), the position of the knots are defined
at 10<sup>th</sup>,50<sup>th</sup> and 90<sup>th</sup> quantile of the
continuous predictor distribution. For more details see Frank Harrell’s
book ‘Regression Model Strategies’ on page 27 (second edition).  
The user can specify the positions of the knots instead of using the
default calculation of the knots proposed in the book of Frank Harrell.
To deal with very large influential value, we winzorize progesterone
level to the 90<sup>th</sup> percentile.

``` r
# Winzorise PGR to the 99th percentile to deal with very large influential values;
p99 <- quantile(rdata$pgr, probs = .99)
rdata$pgr2 <- pmin(rdata$pgr, p99)
vdata$pgr2 <- pmin(vdata$pgr, p99)

dd <- datadist(rdata)
options(datadist = "dd")
fit_pgr <- cph(Surv(ryear, status) ~ rcs(pgr2), data = rdata, x = T, y = T, surv = T)
plot(Predict(fit_pgr))
```

<img src="imgs/01_predsurv_simplified/ff-1.png" width="672" style="display: block; margin: auto;" />

``` r
options(datadist = NULL)
```

We should model the progesterone level using a three-knot restricted
cubic spline. We save the spline in the development and validation data.

``` r
# Save in the date the restricted cubic spline term using rms::rcspline.eval() package
# Development set
rcs3_pgr <- rcspline.eval(rdata$pgr2, knots = c(0, 41, 486))
attr(rcs3_pgr, "dim") <- NULL
attr(rcs3_pgr, "knots") <- NULL
rdata$pgr3 <- rcs3_pgr
rm(rcs3_pgr)

# Validation set
rcs3_pgr <- rcspline.eval(vdata$pgr2, knots = c(0, 41, 486))
attr(rcs3_pgr, "dim") <- NULL
attr(rcs3_pgr, "knots") <- NULL
vdata$pgr3 <- rcs3_pgr
rm(rcs3_pgr)
```

### 1.3 Model development: first check - the proportional hazard (PH) assumption

We now examine the fits in a more careful way by checking the
proportionality of the hazards of the Cox regression model. Firstly, we
fit the first prediction model in the development data using size, node,
grade. Then, we check the PH assumption.

``` r
dd <- datadist(rdata)
options(datadist = "dd")
options(prType = "html")
fit1_cph <- cph(Surv(ryear, status) ~ csize + cnode + cgrade, rdata, x = T, y = T, surv = T)


zp1 <- cox.zph(fit1_cph, transform = "identity")
kable(round(zp1$table, 3)) %>% kable_styling("striped", position = "center")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
chisq
</th>
<th style="text-align:right;">
df
</th>
<th style="text-align:right;">
p
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
csize
</td>
<td style="text-align:right;">
28.320
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.000
</td>
</tr>
<tr>
<td style="text-align:left;">
cnode
</td>
<td style="text-align:right;">
16.798
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.000
</td>
</tr>
<tr>
<td style="text-align:left;">
cgrade
</td>
<td style="text-align:right;">
3.779
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.052
</td>
</tr>
<tr>
<td style="text-align:left;">
GLOBAL
</td>
<td style="text-align:right;">
35.323
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.000
</td>
</tr>
</tbody>
</table>

``` r
oldpar <- par(mfrow = c(2, 2), mar = c(5, 5, 1, 1))
for (i in 1:3) {
  plot(zp1[i], resid = F)
  abline(0, 0, lty = 3)
}
par(oldpar)
```

<img src="imgs/01_predsurv_simplified/ph-1.png" width="672" style="display: block; margin: auto;" />

``` r
options(datadist = NULL)
```

The statistical tests show strong evidence of non-proportionality. Since
the number of death is large the formal tests are quite sensitive,
however, and it is important to also examine the graphs.  
These show an estimated coefficient as a function of time. As a further
follow-up we will divide the data into 3 epochs of 0-5, 5-10, and 10+
years, fitting a separate model to each.

``` r
# Development
edata <- survSplit(Surv(ryear, status) ~ .,
  data = rdata, cut = c(5, 10),
  episode = "epoch"
)
efit1 <- cph(Surv(ryear, status) ~ csize + cnode + cgrade,
  data = edata[edata$epoch == 1, ], x = T, y = T, surv = T
)
efit2 <- cph(Surv(ryear, status) ~ csize + cnode + cgrade,
  data = edata[edata$epoch == 2, ], x = T, y = T, surv = T
)
efit3 <- cph(Surv(ryear, status) ~ csize + cnode + cgrade,
  data = edata[edata$epoch == 3, ], x = T, y = T, surv = T
)
res_efit <- round(rbind(e1 = coef(efit1), e2 = coef(efit2), e3 = coef(efit3)), 2)
rownames(res_efit) <- c("Epoch 1: 0-5 years", "Epoch 2: 5-10 years", "Epoch 3: >10 years")
kable(res_efit) %>% kable_styling("striped", position = "center")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
csize=21-50
</th>
<th style="text-align:right;">
csize=&gt;50
</th>
<th style="text-align:right;">
cnode=1-3
</th>
<th style="text-align:right;">
cnode=4+
</th>
<th style="text-align:right;">
cgrade=3+
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Epoch 1: 0-5 years
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
1.09
</td>
<td style="text-align:right;">
0.41
</td>
</tr>
<tr>
<td style="text-align:left;">
Epoch 2: 5-10 years
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.31
</td>
</tr>
<tr>
<td style="text-align:left;">
Epoch 3: &gt;10 years
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.41
</td>
</tr>
</tbody>
</table>

``` r
options(prType = "html")
tt1 <- table(edata$epoch, edata$status, dnn = c("Epoch", "Status"))
rownames(tt1) <- c("Epoch 1: 0-5 years", "Epoch 2: 5-10 years", "Epoch 3: >10 years")
kable(tt1,
  col.names = c("Censored", "Event"),
  row.names = TRUE
) %>% kable_styling("striped", position = "center")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Censored
</th>
<th style="text-align:right;">
Event
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Epoch 1: 0-5 years
</td>
<td style="text-align:right;">
1801
</td>
<td style="text-align:right;">
1181
</td>
</tr>
<tr>
<td style="text-align:left;">
Epoch 2: 5-10 years
</td>
<td style="text-align:right;">
1275
</td>
<td style="text-align:right;">
296
</td>
</tr>
<tr>
<td style="text-align:left;">
Epoch 3: &gt;10 years
</td>
<td style="text-align:right;">
440
</td>
<td style="text-align:right;">
41
</td>
</tr>
</tbody>
</table>

A drastic change in the size coefficients across all epochs is apparent,
along with a major reduction in the nodes coefficient in epoch 3. As an
ameleoration of this we will refit the model using only the first epoch,
which includes most of the recurrences and deaths.  
We applied the administrative censoring at 5 years in the development
data and we assessed the prediction performance of the prognostic model
at 5 years. The 5-year horizon was chosen because it is also a common
prediction time horizon in clinical practice. The hazards in the
development data seem not totally proportional within 5 years but minor
deviation of proportionality were considered acceptable. In R the
summary function creates a more detailed report of a fit.  
Possible alternatives to address the violation of proportional hazards
are:

-   a stratified Cox proportional hazard models for the predictor(s)
    that violates the proportional hazard. The offending predictor will
    be absorbed into the baseline (as a strata), and then use the full
    follow-up. In this scenario, it seems that the size predictor
    violates the proportional hazard assumption. So, three baseline
    hazards will be non-parametrically estimated. For prediction, the
    absolute risks will be calculated considering the different baseline
    hazards. The hazard ratio for size will not be estimated since it is
    incorporated in the baseline;

-   fitting the model with interaction between the offending predictor
    (in this case size) and the logarithm of the time. This complicates
    the interpretability of the model and risk prediction might be not
    fully estimated for new patients at the time of breast cancer
    diagnosis since size interacts with the time which is unknown for
    new patients.

If we ignored the non-proportional hazards entirely, the prediction
performance of the model may lead to increase bias in estimating the
risk especially in terms of calibration performances.

### 1.4 Model development: fit the risk prediction models

We develop the risk prediction model in the development data considering
the first 5-year follow-up to minimize the violation of proportional
hazard including size, nodel and grade. The second model also includes
the progesterone level modelled using a 3-knot restricted cubic
spline.  
We also administratively censored the validation data at 5 years.

``` r
# Consider the first 5-year epoch in the development set
edata1 <- edata[edata$epoch == 1, ]
# Refit the model
efit1 <- coxph(Surv(ryear, status) ~ csize + cnode + cgrade,
  data = edata1, x = T, y = T
)
# Additional marker
efit1b <- coxph(Surv(ryear, status) ~ csize + cnode + cgrade + pgr2 + pgr3,
  data = edata1, x = T, y = T
)


# Validation
evdata <- survSplit(Surv(ryear, status) ~ ., data = vdata, cut = 5, episode = "epoch")
evdata1 <- evdata[evdata$epoch == 1, ]
```

Below the results of the models:

``` r
dd <- datadist(edata1)
options(datadist = "dd")
options(prType = "html")
s5 <- cph(Surv(ryear, status) ~ csize + cnode + cgrade,
  data = edata1, x = T, y = T, surv = T
)
print(s5)
```

 <strong>Cox Proportional Hazards Model</strong>
 
 <pre>
 cph(formula = Surv(ryear, status) ~ csize + cnode + cgrade, data = edata1, 
     x = T, y = T, surv = T)
 </pre>
 
 <table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'></th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; border-right: 1px solid black; text-align: center;'>Model Tests</th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; border-right: 1px solid black; text-align: center;'>Discrimination<br>Indexes</th>
</tr>
</thead>
<tbody>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'>Obs 2982</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>LR χ<sup>2</sup> 465.78</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>R</i><sup>2</sup> 0.145</td>
</tr>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'>Events 1181</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>d.f. 5</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>D</i><sub>xy</sub> 0.356</td>
</tr>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'>Center 0.9167</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>Pr(>χ<sup>2</sup>) 0.0000</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>g</i> 0.702</td>
</tr>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'></td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>Score χ<sup>2</sup> 533.49</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>g</i><sub>r</sub> 2.017</td>
</tr>
<tr>
<td style='min-width: 9em; border-bottom: 2px solid grey; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'></td>
<td style='min-width: 9em; border-bottom: 2px solid grey; border-right: 1px solid black; text-align: center;'>Pr(>χ<sup>2</sup>) 0.0000</td>
<td style='min-width: 9em; border-bottom: 2px solid grey; border-right: 1px solid black; text-align: center;'></td>
</tr>
</tbody>
</table>

 
 <table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr><th style='border-bottom: 1px solid grey; font-weight: 900; border-top: 2px solid grey; min-width: 7em; text-align: center;'></th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>β</th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>S.E.</th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>Wald <i>Z</i></th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>Pr(>|<i>Z</i>|)</th>
</tr>
</thead>
<tbody>
<tr>
<td style='min-width: 7em; text-align: left;'>csize=21-50</td>
<td style='min-width: 7em; text-align: right;'> 0.3943</td>
<td style='min-width: 7em; text-align: right;'> 0.0676</td>
<td style='min-width: 7em; text-align: right;'> 5.83</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; text-align: left;'>csize=>50</td>
<td style='min-width: 7em; text-align: right;'> 0.6230</td>
<td style='min-width: 7em; text-align: right;'> 0.0956</td>
<td style='min-width: 7em; text-align: right;'> 6.52</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; text-align: left;'>cnode=1-3</td>
<td style='min-width: 7em; text-align: right;'> 0.3613</td>
<td style='min-width: 7em; text-align: right;'> 0.0788</td>
<td style='min-width: 7em; text-align: right;'> 4.59</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; text-align: left;'>cnode=4+</td>
<td style='min-width: 7em; text-align: right;'> 1.0897</td>
<td style='min-width: 7em; text-align: right;'> 0.0730</td>
<td style='min-width: 7em; text-align: right;'>14.92</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: left;'>cgrade=3+</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'> 0.4145</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'> 0.0750</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'> 5.53</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'><0.0001</td>
</tr>
</tbody>
</table>

``` r
options(prType = "html")
s6 <- cph(Surv(ryear, status) ~ csize + cnode + cgrade + rcs(pgr2, c(0, 41, 486)),
  data = edata1, x = T, y = T, surv = T
)
print(s6)
```

 <strong>Cox Proportional Hazards Model</strong>
 
 <pre>
 cph(formula = Surv(ryear, status) ~ csize + cnode + cgrade + 
     rcs(pgr2, c(0, 41, 486)), data = edata1, x = T, y = T, surv = T)
 </pre>
 
 <table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'></th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; border-right: 1px solid black; text-align: center;'>Model Tests</th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; border-right: 1px solid black; text-align: center;'>Discrimination<br>Indexes</th>
</tr>
</thead>
<tbody>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'>Obs 2982</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>LR χ<sup>2</sup> 504.93</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>R</i><sup>2</sup> 0.156</td>
</tr>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'>Events 1181</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>d.f. 7</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>D</i><sub>xy</sub> 0.374</td>
</tr>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'>Center 0.6584</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>Pr(>χ<sup>2</sup>) 0.0000</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>g</i> 0.751</td>
</tr>
<tr>
<td style='min-width: 9em; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'></td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'>Score χ<sup>2</sup> 569.09</td>
<td style='min-width: 9em; border-right: 1px solid black; text-align: center;'><i>g</i><sub>r</sub> 2.119</td>
</tr>
<tr>
<td style='min-width: 9em; border-bottom: 2px solid grey; border-left: 1px solid black; border-right: 1px solid black; text-align: center;'></td>
<td style='min-width: 9em; border-bottom: 2px solid grey; border-right: 1px solid black; text-align: center;'>Pr(>χ<sup>2</sup>) 0.0000</td>
<td style='min-width: 9em; border-bottom: 2px solid grey; border-right: 1px solid black; text-align: center;'></td>
</tr>
</tbody>
</table>

 
 <table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr><th style='border-bottom: 1px solid grey; font-weight: 900; border-top: 2px solid grey; min-width: 7em; text-align: center;'></th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>β</th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>S.E.</th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>Wald <i>Z</i></th>
<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>Pr(>|<i>Z</i>|)</th>
</tr>
</thead>
<tbody>
<tr>
<td style='min-width: 7em; text-align: left;'>csize=21-50</td>
<td style='min-width: 7em; text-align: right;'>  0.3695</td>
<td style='min-width: 7em; text-align: right;'> 0.0677</td>
<td style='min-width: 7em; text-align: right;'> 5.46</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; text-align: left;'>csize=>50</td>
<td style='min-width: 7em; text-align: right;'>  0.5983</td>
<td style='min-width: 7em; text-align: right;'> 0.0955</td>
<td style='min-width: 7em; text-align: right;'> 6.26</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; text-align: left;'>cnode=1-3</td>
<td style='min-width: 7em; text-align: right;'>  0.3856</td>
<td style='min-width: 7em; text-align: right;'> 0.0788</td>
<td style='min-width: 7em; text-align: right;'> 4.89</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; text-align: left;'>cnode=4+</td>
<td style='min-width: 7em; text-align: right;'>  1.0856</td>
<td style='min-width: 7em; text-align: right;'> 0.0729</td>
<td style='min-width: 7em; text-align: right;'>14.88</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; text-align: left;'>cgrade=3+</td>
<td style='min-width: 7em; text-align: right;'>  0.3491</td>
<td style='min-width: 7em; text-align: right;'> 0.0758</td>
<td style='min-width: 7em; text-align: right;'> 4.60</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; text-align: left;'>pgr2</td>
<td style='min-width: 7em; text-align: right;'> -0.0033</td>
<td style='min-width: 7em; text-align: right;'> 0.0006</td>
<td style='min-width: 7em; text-align: right;'>-5.47</td>
<td style='min-width: 7em; text-align: right;'><0.0001</td>
</tr>
<tr>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: left;'>pgr2'</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'>  0.0143</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'> 0.0030</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'> 4.78</td>
<td style='min-width: 7em; border-bottom: 2px solid grey; text-align: right;'><0.0001</td>
</tr>
</tbody>
</table>

``` r
options(datadist = NULL)
```

The coefficients of the models indicated that higher size, higher number
of positive lymph nodes and higher grade is more associate with poorer
prognosis. The association of the progesterone biomarker and the outcome
is non-linear as investigated previously.

## Goal 2: Assessing performance in survival prediction models

The performance of a risk prediction models may be evaluated through:

-   discrimination: the ability of the model to identify patients with
    and without the outcome and it requires the coefficients (or the log
    of the hazard ratios) of the developed risk prediction model to be
    evaluated.

-   calibration: the agreement between observed and predicted
    probabilities. It requires the baseline (cumulative) hazard or
    survival.

-   overall performance measures: as a combination of discrimination and
    calibration and/or as a measure of the explained variation;

Unfortunately, a few publications report the complete baseline
(cumulative) hazard or survival or even the baseline (cumulative) hazard
or survival at fixed time horizon *t*. If we had both individual data of
the development and validation, a complete assessment of discrimination
and calibration would be possible. We could evaluate the prediction
performance of a risk prediction model at a fixed time horizon(s) *t*
and for the complete follow-up time. In risk prediction, physicians
typically focus on one or more clinically relevant time horizons to
inform subjects about their risk. For this reasons, according to
information available, different levels of validation assessment are
possible. Here we aim to assess the prediction performance of a risk
prediction model with time-to-event outcome in case all individual data
are available and in case of only the model equation of a fixed time
horizon (i.e. at 5 years) is provided including the baseline survival.

### 2.1 Overall performance measures

We calculate the Brier Score and the Index of Prediction Accuracy (IPA,
the scaled Brier) as a overall performance measure.

First we save elements need to calculate the performance measures as the
linear predictor and the predicted survival at 5 years in the
development and validation data. Secondly, we create B bootstrap data to
calculate percentile bootstrap confidence intervals needed for some
performance measures.

``` r
edata1$lp <- predict(efit1) # linear predictor
edata1$predsurv5 <- predictSurvProb(efit1, newdata = edata1, times = 5) # predicted survival
edata1$lp_1b <- predict(efit1b)
edata1$predsurv5_1b <- predictSurvProb(efit1b, newdata = edata1, times = 5)


evdata1$lp_1b <- predict(efit1b, newdata = evdata1)
evdata1$predsurv5 <- predictSurvProb(efit1, newdata = evdata1, times = 5)
evdata1$lp_1b <- predict(efit1b, newdata = evdata1)
evdata1$predsurv5_1b <- predictSurvProb(efit1b, newdata = evdata1, times = 5)

set.seed(20200416)
eboot <- bootstraps(edata1, times = 10)
vboot <- bootstraps(evdata1, times = 10)
# NOTE: B=10 otherwise the computation time will be too long
```

We calculate the overall performance measures: Brier score, IPA and the
corresponding confidence intervals.

``` r
# Development set (apparent Brier and IPA) without pgr
score_rdata1 <-
  Score(list("Cox development" = efit1),
    formula = Surv(ryear, status) ~ 1, data = edata1, conf.int = TRUE, times = 4.95,
    cens.model = "km", metrics = "brier",
    summary = "ipa"
  )

# Validation set without pgr
score_vdata1 <-
  Score(list("Cox development" = efit1),
    formula = Surv(ryear, status) ~ 1, data = evdata1, conf.int = TRUE, times = 4.95,
    cens.model = "km", metrics = "brier",
    summary = "ipa"
  )

# Development set (apparent Brier and IPA) with pgr
score_rdata1b <-
  Score(list("Cox development" = efit1b),
    formula = Surv(ryear, status) ~ 1, data = edata1, conf.int = TRUE, times = 4.95,
    cens.model = "km", metrics = "brier",
    summary = "ipa"
  )

# Validation set with pgr
score_vdata1b <-
  Score(list("Cox development" = efit1b),
    formula = Surv(ryear, status) ~ 1, data = evdata1, conf.int = TRUE, times = 4.95,
    cens.model = "km", metrics = "brier",
    summary = "ipa"
  )


# IPA confidence intervals
# NOTE: computational time is too long to compute them when the number of bootstrap replications is high
score_boot_1 <- function(split) {
  Score(list("Cox development" = efit1),
    formula = Surv(ryear, status) ~ 1, data = analysis(split), conf.int = FALSE, times = 4.95,
    cens.model = "km", metrics = "brier",
    summary = "ipa"
  )$Brier[[1]]$IPA[2]
}

score_boot_1b <- function(split) {
  Score(list("Cox development" = efit1b),
    formula = Surv(ryear, status) ~ 1, data = analysis(split), conf.int = FALSE, times = 4.95,
    cens.model = "km", metrics = "brier",
    summary = "ipa"
  )$Brier[[1]]$IPA[2]
}
eboot <- eboot %>% mutate(
  IPA1 = map_dbl(splits, score_boot_1),
  IPA1b = map_dbl(splits, score_boot_1b)
) # Development dataset

vboot <- vboot %>% mutate(
  IPA1 = map_dbl(splits, score_boot_1),
  IPA1b = map_dbl(splits, score_boot_1b)
) # Validation dataset
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Apparent

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Apparent + PGR

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

External

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

External + PGR

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Brier
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.23
</td>
</tr>
<tr>
<td style="text-align:left;">
IPA
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.19
</td>
</tr>
</tbody>
</table>

As expected the overall performance measures were lower in the external
validation. Including information about PGR slightly improved the
overall performance.

### 2.2 Discrimination measures

Discrimination is the ability to differentiate between subjects who have
the outcome and subjects who do not. In prognostic modelling,
discrimination reflects separation between survival curves for
individuals or groups. We propose to calculate Tte Uno’s time-dependent
AUC measure of discrimination at a specific time horizon *t*.  
More details are in the paper and in the references below.  
The time horizon to calculate the time-dependent measures was set to 5
years. Values close to 1 indicate good discrimination ability, while
values close to 0.5 indicated poor discrimination ability.  
We used the time horizon at 4.95 and not 5 years since controls are
considered patients at risk after the time horizon and we
administratively censored at 5 years to minimize the violation of PH
assumption (see paragraph 1.3).

``` r
# Time-dependent AUC (in Table 3 called Uno's TD AUC at 5 years) ###
# Uno's time-dependent Area Under the Curve
# Apparent
Uno_rdata1 <-
  timeROC(
    T = edata1$ryear, delta = edata1$status,
    marker = predict(efit1, newdata = edata1),
    cause = 1, weighting = "marginal", times = 4.95,
    iid = TRUE
  )

# External validation
Uno_vdata1 <-
  timeROC(
    T = evdata1$ryear, delta = evdata1$status,
    marker = predict(efit1, newdata = evdata1),
    cause = 1, weighting = "marginal", times = 4.95,
    iid = TRUE
  )

# Apparent with pgr
Uno_rdata1b <-
  timeROC(
    T = edata1$ryear, delta = edata1$status,
    marker = predict(efit1b, newdata = edata1),
    cause = 1, weighting = "marginal", times = 4.95,
    iid = TRUE
  )

# External validation with pgr
Uno_vdata1b <-
  timeROC(
    T = evdata1$ryear, delta = evdata1$status,
    marker = predict(efit1b, newdata = evdata1),
    cause = 1, weighting = "marginal", times = 4.95,
    iid = TRUE
  )
# NOTE: if you have a lot of data n > 2000, standard error computation may be really long.
# In that case, please use bootstrap percentile to calculate confidence intervals.
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Apparent

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Apparent + PGR

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

External

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

External + PGR

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Uno AUC
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.7
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.78
</td>
</tr>
</tbody>
</table>

The time-dependent AUCs at 5 years were in the external validation were
between 0.68 and 0.72 showing moderate discrimination. The additional
information of the PGR marker increased the discrimination by 0.04.

### 2.3 Calibration

Calibration is the agreement between observed outcomes and predicted
probabilities. For example, in survival models, a predicted survival
probability at a fixed time horizon *t* of 80% is considered reliable if
it can be expected that 80 out of 100 will survive among patients
received a predicted survival probability of 80%.

Calibration is measured by:

-   Observed and Expected ratio at time horizon (*t*):

    -   the number of observed event is calculated as the one minus the
        Kaplan-Meier curve at time *t*;

    -   the number of expected event is calculated as the mean of the
        predicted risk at time *t*;

    -   Confidence intervals are calculated using the Normal
        approximation of the Poisson distribution.

-   Calibration plot: it is a graphical representation of calibration
    in-the-large and calibration. It shows:

    -   on the *x-axis* the predicted survival (or risk) probabilities
        at a fixed time horizon (e.g. at 5 years);

    -   on the *y-axis* the observed survival (or risk) probabilities at
        a fixed time horizon (e.g. at 5 years);

    -   The 45-degree line indicates the good overall calibration.
        Points below the 45-degree line indicates that the model
        overestimate the observed risk. If points are above the
        45-degree line, the model underestimate the observed risk; The
        observed probabilities estimated by the Kaplan-Meier curves (in
        case of survival) or by the complementary of the Kaplan-Meier
        curves (in case of risk in absence of competing risks) are
        represented in terms of percentiles of the predicted survival
        (risk) probabilities.

Other calibration measures are proposed in the literature. More details
are provided in the references at the end of the document.

### 2.3.1 Observed/Expected ratio

We calculate the observed/ expected ratio (OE) at 5 years in the
development and validation data. In the development data the OE should
be (close to) 1.

``` r
#' Calculate the O/E ratio
#'
#' @param tfup follow-up time
#' @param event variable indication event (0=no, 1=yes)
#' @param pred.risk predicted risk at time horizon
#' @param thorizon time horizon to calculate the observed/expected ratio at time t
#' @param alpha type I error (to calculate confidence intervals, default: 0.05)
#'
#' @return 
#'
#' @examples

OE.ratio.t <- function(tfup, event, pred.risk, thorizon, alpha = .05) {
  obs_t <- 1 - summary(survfit(Surv(tfup, event) ~ 1), times = thorizon)$surv
  exp_t <- mean(pred.risk)

  OE <- obs_t / exp_t
  OE.lower <- OE * exp(-qnorm(1 - alpha / 2) * sqrt(1 / sum(event)))
  OE.upper <- OE * exp(+qnorm(1 - alpha / 2) * sqrt(1 / sum(event)))
  res <- c(OE, OE.lower, OE.upper)
  names(res) <- c("Estimate", "Lower .95", "Upper .95")
  return(res)
}

OE.rdata <- round(OE.ratio.t(rdata$ryear, event = rdata$status, pred.risk = 1 - edata1$predsurv5, thorizon = 5), 2)
OE.rdata1b <- round(OE.ratio.t(rdata$ryear, event = rdata$status, pred.risk = 1 - edata1$predsurv5_1b, thorizon = 5), 2)

OE.vdata <- round(OE.ratio.t(vdata$ryear, event = vdata$status, pred.risk = 1 - evdata1$predsurv5, thorizon = 5), 2)
OE.vdata1b <- round(OE.ratio.t(vdata$ryear, event = vdata$status, pred.risk = 1 - evdata1$predsurv5_1b, thorizon = 5), 2)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Apparent

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Apparent + PGR

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

External

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

External + PGR

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
<th style="text-align:right;">
Estimate
</th>
<th style="text-align:right;">
Lower .95
</th>
<th style="text-align:right;">
Upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
O/E ratio
</td>
<td style="text-align:right;">
0.99
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:right;">
0.99
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:right;">
1.11
</td>
<td style="text-align:right;">
0.99
</td>
<td style="text-align:right;">
1.25
</td>
<td style="text-align:right;">
1.08
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
1.21
</td>
</tr>
</tbody>
</table>

The models PGR tend to slightly underpredict the risk of mortality in
the validation data (OE: 1.11, 95%CI: 0.99-1.25 without PGR; OE: 1.08
95%CI: 0.96-1.21).

### 2.3.2 Calibration plot using restricted cubic splines

Calibration plots of the external validation data with and without PGR
are calculated and shown using restricted cubic splines.  
The interpretation of the calibration plot were provided in the section
2.3 of this document, in the corresponding paper and in the literature
provided in the paper and at the end of this document. More details
about the method are given in the references below.

``` r
## External data without PGR
evdata1$predmort5 <- 1 - survest(s5, newdata = evdata1, times = 5)$surv
evdata1$predmort5.cll <- log(-log(1 - evdata1$predmort5))

# Estimate
vcalibrate.efit1 <- cph(Surv(ryear, status) ~ rcs(predmort5.cll, 3),
  x = T, y = T,
  data = evdata1, surv = T
) # better rms::cph?
predict.grid <- seq(quantile(evdata1$predmort5, prob = 0.01), quantile(evdata1$predmort5, prob = 0.99), length = 100)
predict.grid.cll <- log(-log(1 - predict.grid))
predict.grid.df <- data.frame(predict.grid)
predict.grid.cll.df <- data.frame(predict.grid.cll)
names(predict.grid.df) <- "predmort5"
names(predict.grid.cll.df) <- "predmort5.cll"

# Plot
pred.vcalibrate.efit1 <- 1 - survest(vcalibrate.efit1, newdata = predict.grid.cll.df, times = 5)$surv
pred.vcalibrate.efit1.upper <- 1 - survest(vcalibrate.efit1, newdata = predict.grid.cll.df, times = 5)$lower
pred.vcalibrate.efit1.lower <- 1 - survest(vcalibrate.efit1, newdata = predict.grid.cll.df, times = 5)$upper
par(xaxs = "i", yaxs = "i", las = 1)
plot(predict.grid, pred.vcalibrate.efit1,
  type = "l", lty = 1, xlim = c(0, 1),
  ylim = c(0, 1), lwd = 2,
  xlab = "Predicted probability",
  ylab = "Observed probability", bty = "n"
)
lines(predict.grid, pred.vcalibrate.efit1.lower, type = "l", lty = 2, lwd = 2)
lines(predict.grid, pred.vcalibrate.efit1.upper, type = "l", lty = 2, lwd = 2)
abline(0, 1, lwd = 2, lty = 2, col = "red")
title("A External data without PGR", adj = 0)
```

<img src="imgs/01_predsurv_simplified/cal_rcs-1.png" width="672" style="display: block; margin: auto;" />

``` r
# par(new=T)
# plot(density(evdata1$predmort5),axes=F,xlab=NA,ylab=NA,
#     main="")

## External data with PGR
evdata1$predmort5 <- 1 - survest(s6, newdata = evdata1, times = 5)$surv
evdata1$predmort5.cll <- log(-log(1 - evdata1$predmort5))

# Estimate
vcalibrate.efit1b <- cph(Surv(ryear, status) ~ rcs(predmort5.cll, 3),
  x = T, y = T, data = evdata1, surv = T
)
predict.grid <- seq(quantile(evdata1$predmort5, prob = 0.01), quantile(evdata1$predmort5, prob = 0.99), length = 100)
predict.grid.cll <- log(-log(1 - predict.grid))
predict.grid.df <- data.frame(predict.grid)
predict.grid.cll.df <- data.frame(predict.grid.cll)
names(predict.grid.df) <- "predmort5"
names(predict.grid.cll.df) <- "predmort5.cll"

# Plot
pred.vcalibrate.efit1b <- 1 - survest(vcalibrate.efit1b, newdata = predict.grid.cll.df, times = 5)$surv
pred.vcalibrate.efit1b.lower <- 1 - survest(vcalibrate.efit1b, newdata = predict.grid.cll.df, times = 5)$upper
pred.vcalibrate.efit1b.upper <- 1 - survest(vcalibrate.efit1b, newdata = predict.grid.cll.df, times = 5)$lower
par(xaxs = "i", yaxs = "i", las = 1)
plot(predict.grid, pred.vcalibrate.efit1b,
  type = "l", lty = 1, xlim = c(0, 1),
  ylim = c(0, 1), lwd = 2,
  xlab = "Predicted probability of recurrence at 5 years",
  ylab = "Observed probability of recurrence at 5 years", bty = "n"
)
lines(predict.grid, pred.vcalibrate.efit1b.lower, type = "l", lty = 2, lwd = 2)
lines(predict.grid, pred.vcalibrate.efit1b.upper, type = "l", lty = 2, lwd = 2)
abline(0, 1, lwd = 2, lty = 2, col = "red")
title("B External data with PGR", adj = 0)
```

<img src="imgs/01_predsurv_simplified/cal_rcs-2.png" width="672" style="display: block; margin: auto;" />

``` r
# par(new=T)
# plot(density(evdata1$predmort5),axes=F,xlab=NA,ylab=NA,
#     main="")
```

Both plots identified good calibration although probabilities of
recurrence were slightly underestimated especially for the lowest and
the highest values of the observed probabilities of recurrence.  
The additional information of PGR improved the overall calibration,
especially for the highest values, as shown in the two calibration plots
above.

## Goal 3: Clinical utility

Discrimination and calibration measures are essential to assess the
prediction performance but insufficient to evaluate the potential
clinical utility of a risk prediction model for decision making. When
new markers are available, clinical utility assessment evaluates whether
the extended model helps to improve decision making.  
Clinical utility is measured by the net benefit that includes the number
of true positives and the number of false positives. For example, in
time-to-event models, the true positives reflect the benefit of being
event free for a given time horizon using additional interventions such
as additional treatments, personalized follow-up or additional
surgeries. The false positives represent the harms of unnecessary
interventions.  
Generally, in medicine, clinicians accepts to treat a certain number of
patients for which interventions are unnecessary to be event free for a
given time horizon. So, false negatives (the harm of not being event
free for a given time horizon) are more important than false positives
(the harm of unnecessary interventions). Thus, net benefit is the number
of true positives classifications minus the false positives
classifications weighted by a factor related to the harm of not
preventing the event versus unnecessary interventions. The weighting is
derived from the threshold probability to death (one minus survival
probability) using a defined time horizon (for example 5 years since
diagnosis). For example, a threshold of 10% implies that additional
interventions for 10 patients of whom one would have experience the
event in 5 years if untreated is acceptable (thus treating 9 unnecessary
patients). This strategy is compared with the strategies of treat all
and treat none patients. If overtreatment is harmful, a higher threshold
should be used.  
The net benefit is calculated as:

<img src="https://render.githubusercontent.com/render/math?math=%5Chuge%7B%5Cfrac%7BTP%7D%7Bn%7D-%5Cfrac%7BFP%7D%7Bn%7D*%5Cfrac%7Bp_t%7D%7B1-p_t%7D%7D">

*TP*=true positive patients  
*FP*=false positive patients  
*n*=number of patients and *p*<sub>t</sub> is the risk threshold.

For survival data *TP* and *FP* is calculated as follows:  
<img src="https://render.githubusercontent.com/render/math?math=%5CLarge%7BTP%20%3D%20%5B1-S(t)%7C%20X%3D1%5D*P(X%3D1)*n%7D">

<img src="https://render.githubusercontent.com/render/math?math=%5CLarge%7BFP%20%3D%20%5BS(t)%7C%20X%3D1%5D*P(X%3D1)*n%7D">

where  
*S(t)* survival at time *t*  
*X=1* where the predicted probability at time *t* is *p*<sub>t</sub>

And the the decision curve is calculated as follows:

1.  Choose a time horizon (in this case 5 years);
2.  Specify a risk threshold which reflects the ratio between harms and
    benefit of an additional intervention;
3.  Calculate the number of true positive and false positive given the
    threshold specified in (2);
4.  Calculate the net benefit of the survival model;
5.  Plot net benefit on the *y-axis* against the risk threshold on the
    *x-axis*;
6.  Repeat steps 2-4 for each model consideration;
7.  Repeat steps 2-4 for the strategy of assuming all patients are
    treated;
8.  Draw a straight line parallel to the *x-axis* at y=0 representing
    the net benefit associated with the strategy of assuming that all
    patients are not treated.

Given some thresholds, the model/strategy with higher net benefit
represents the one that potentially improves clinical decision making.
However, poor discrimination and calibration lead to lower net benefit.

``` r
# We could use the rms::survest() or pec::predictEventProb() to get the mortality at observed time
# Development data
# Predicted probability calculation
edata1$mort5_efit1 <- 1 - predictSurvProb(efit1, newdata = edata1, times = 5)

# Extended model with PGR
# Predicted probability calculation
edata1$mort5_efit1b <- 1 - predictSurvProb(efit1b, newdata = edata1, times = 5)

# Run decision curve analysis

# Development data
# Model without PGR
edata1 <- as.data.frame(edata1)
dca_rdata_1 <- stdca(
  data = edata1, outcome = "status", ttoutcome = "ryear",
  timepoint = 5, predictors = "mort5_efit1", xstop = 1.0,
  ymin = -0.01, graph = FALSE
)
```

    ## [1] "mort5_efit1: No observations with risk greater than 81%, and therefore net benefit not calculable in this range."

``` r
# Model with PGR
dca_rdata_1b <- stdca(
  data = edata1, outcome = "status", ttoutcome = "ryear",
  timepoint = 5, predictors = "mort5_efit1b", xstop = 1.0,
  ymin = -0.01, graph = FALSE
)
```

    ## [1] "mort5_efit1b: No observations with risk greater than 86%, and therefore net benefit not calculable in this range."

``` r
# Decision curves plot
par(xaxs = "i", yaxs = "i", las = 1)
plot(dca_rdata_1$net.benefit$threshold,
  dca_rdata_1$net.benefit$mort5_efit1,
  type = "l", lwd = 2, lty = 1,
  xlab = "Threshold probability in %", ylab = "Net Benefit",
  xlim = c(0, 1), ylim = c(-0.10, 0.45), bty = "n",
  cex.lab = 1.2, cex.axis = 1
)
# legend('topright',c('Treat all','Treat none','Prediction model'),
#        lwd=c(2,2,2),lty=c(1,1,2),col=c('darkgray','black','black'),bty='n')
lines(dca_rdata_1$net.benefit$threshold, dca_rdata_1$net.benefit$none, type = "l", lwd = 2, lty = 4)
lines(dca_rdata_1$net.benefit$threshold, dca_rdata_1$net.benefit$all, type = "l", lwd = 2, col = "darkgray")
lines(dca_rdata_1b$net.benefit$threshold, dca_rdata_1b$net.benefit$mort5_efit1b, type = "l", lwd = 2, lty = 5)
legend("topright",
  c(
    "Treat All",
    "Original model",
    "Original model + PGR",
    "Treat None"
  ),
  lty = c(1, 1, 5, 4), lwd = 2, col = c("darkgray", "black", "black", "black"),
  bty = "n"
)
title("A Development data", adj = 0, cex = 1.5)
```

<img src="imgs/01_predsurv_simplified/dca-1.png" width="672" style="display: block; margin: auto;" />

``` r
# External data
# Validation data
# Predicted probability calculation
evdata1$mort5_model1 <- 1 - predictSurvProb(efit1, newdata = evdata1, times = 5)

# Extended model with PGR
# Predicted probability calculation
evdata1$mort5_model1b <- 1 - predictSurvProb(efit1b, newdata = evdata1, times = 5)

# Run decision curve analysis

# Development data
# Model without PGR
evdata1 <- as.data.frame(evdata1)
dca_vdata_model1 <- stdca(
  data = evdata1, outcome = "status", ttoutcome = "ryear",
  timepoint = 5, predictors = "mort5_model1", xstop = 1.0,
  ymin = -0.01, graph = FALSE
)
```

    ## [1] "mort5_model1: No observations with risk greater than 81%, and therefore net benefit not calculable in this range."

``` r
# Model with PGR
dca_vdata_model1b <- stdca(
  data = evdata1, outcome = "status", ttoutcome = "ryear",
  timepoint = 5, predictors = "mort5_model1b", xstop = 1,
  ymin = -0.01, graph = FALSE
)
```

    ## [1] "mort5_model1b: No observations with risk greater than 86%, and therefore net benefit not calculable in this range."

``` r
# Decision curves plot
par(xaxs = "i", yaxs = "i", las = 1)
plot(dca_vdata_model1$net.benefit$threshold,
  dca_vdata_model1$net.benefit$mort5_model1,
  type = "l", lwd = 2, lty = 1,
  xlab = "Threshold probability in %", ylab = "Net Benefit",
  xlim = c(0, 1), ylim = c(-0.10, 0.60), bty = "n",
  cex.lab = 1.2, cex.axis = 1
)
# legend('topright',c('Treat all','Treat none','Prediction model'),
#        lwd=c(2,2,2),lty=c(1,1,2),col=c('darkgray','black','black'),bty='n')
lines(dca_vdata_model1$net.benefit$threshold, dca_vdata_model1$net.benefit$none, type = "l", lwd = 2, lty = 4)
lines(dca_vdata_model1$net.benefit$threshold, dca_vdata_model1$net.benefit$all, type = "l", lwd = 2, col = "darkgray")
lines(dca_vdata_model1b$net.benefit$threshold, dca_vdata_model1b$net.benefit$mort5_model1b, type = "l", lwd = 2, lty = 5)
legend("topright",
  c(
    "Treat All",
    "Original model",
    "Original model + PGR",
    "Treat None"
  ),
  lty = c(1, 1, 5, 4), lwd = 2, col = c("darkgray", "black", "black", "black"),
  bty = "n"
)
title("B External data", adj = 0, cex = 1.5)
```

<img src="imgs/01_predsurv_simplified/dca-2.png" width="672" style="display: block; margin: auto;" />

Based on previous research we used a range of thresholds from 14% to 23%
for adjuvant chemotherapy. If we choose a threshold of 20% the model had
a net benefit of 0.262 in the development data using the basic model.
This means that the model would identify 26 patients per 100 who will
have recurrent breast cancer or die within 5 years since diagnosis and
thus adjuvant chemotherapy is really needed.  
The decision curve shows that the net benefit would be much larger for
higher threshold values, i.e., patients accepting higher risks of
recurrence compared to the treat all strategy. The same interpretation
may be used in the validation data: choosing a threshold of 20% the
basic model had a net benefit of 0.385 for the basic and the extended
model.  
Moreover, net benefit can be defined in terms of reduction of avoidable
interventions (e.g adjuvant chemotherapy per 100 patients) by:

<img src="https://render.githubusercontent.com/render/math?math=%5Chuge%7B%5Cfrac%7BNB_%7Bmodel%7D%20-%20NB_%7Ball%7D%7D%7B(p_t%2F%20(1-p_t))%7D*100%7D%0A">

where *NB*<sub>model</sub> is the net benefit of the prediction model,
*NB*<sub>all</sub> is the net benefit of the strategy treat all and
*p*<sub>*t*</sub> is the risk threshold.

## References

-   Overall measures  
    Reference:
    <https://diagnprognres.biomedcentral.com/articles/10.1186/s41512-018-0029-2*/>  
    R Vignette:
    <https://cran.r-project.org/web/packages/riskRegression/vignettes/IPA.html#fn.1>  

-   Discrimination measures  
    References:  
    <https://www.jstor.org/stable/27639883>  
    <https://onlinelibrary.wiley.com/doi/10.1002/sim.5958>  

-   Calibration  
    <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3933449/pdf/nihms542648.pdf>  
    <https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8281>  
    <https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8570>  

-   Clinical utility (decision curves)  
    R/SAS/STATA code and references:
    <https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/decision-curve-analysis>  
    More guidelines about net benefit assessment and interpretation  
    <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6261531/>  
    <https://diagnprognres.biomedcentral.com/articles/10.1186/s41512-019-0064-7>  

-   Other useful references  
    <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6728752/>  
    <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7100774/>  

## Reproducibility ticket

``` r
sessioninfo::session_info()
```

    ## - Session info ---------------------------------------------------------------
    ##  setting  value                       
    ##  version  R version 4.0.3 (2020-10-10)
    ##  os       Windows 10 x64              
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_Netherlands.1252    
    ##  ctype    English_Netherlands.1252    
    ##  tz       Europe/Berlin               
    ##  date     2021-04-03                  
    ## 
    ## - Packages -------------------------------------------------------------------
    ##  package        * version    date       lib source        
    ##  assertthat       0.2.1      2019-03-21 [1] CRAN (R 4.0.3)
    ##  backports        1.2.0      2020-11-02 [1] CRAN (R 4.0.3)
    ##  base64enc        0.1-3      2015-07-28 [1] CRAN (R 4.0.3)
    ##  bit              4.0.4      2020-08-04 [1] CRAN (R 4.0.3)
    ##  bit64            4.0.5      2020-08-30 [1] CRAN (R 4.0.3)
    ##  blob             1.2.1      2020-01-20 [1] CRAN (R 4.0.3)
    ##  boot           * 1.3-25     2020-04-26 [2] CRAN (R 4.0.3)
    ##  broom            0.7.4      2021-01-29 [1] CRAN (R 4.0.3)
    ##  broom.helpers    1.2.1      2021-02-26 [1] CRAN (R 4.0.4)
    ##  cachem           1.0.1      2021-01-21 [1] CRAN (R 4.0.3)
    ##  cellranger       1.1.0      2016-07-27 [1] CRAN (R 4.0.3)
    ##  checkmate        2.0.0      2020-02-06 [1] CRAN (R 4.0.3)
    ##  chron            2.3-56     2020-08-18 [1] CRAN (R 4.0.3)
    ##  cli              2.3.0      2021-01-31 [1] CRAN (R 4.0.3)
    ##  cluster          2.1.0      2019-06-19 [2] CRAN (R 4.0.3)
    ##  cmprsk           2.2-10     2020-06-09 [1] CRAN (R 4.0.3)
    ##  codetools        0.2-16     2018-12-24 [2] CRAN (R 4.0.3)
    ##  colorspace       2.0-0      2020-11-11 [1] CRAN (R 4.0.3)
    ##  conquer          1.0.2      2020-08-27 [1] CRAN (R 4.0.3)
    ##  crayon           1.4.0      2021-01-30 [1] CRAN (R 4.0.3)
    ##  curl             4.3        2019-12-02 [1] CRAN (R 4.0.3)
    ##  data.table       1.13.6     2020-12-30 [1] CRAN (R 4.0.3)
    ##  DBI              1.1.1      2021-01-15 [1] CRAN (R 4.0.3)
    ##  dbplyr           2.1.0      2021-02-03 [1] CRAN (R 4.0.3)
    ##  digest           0.6.27     2020-10-24 [1] CRAN (R 4.0.3)
    ##  dplyr          * 1.0.3      2021-01-15 [1] CRAN (R 4.0.3)
    ##  ellipsis         0.3.1      2020-05-15 [1] CRAN (R 4.0.3)
    ##  evaluate         0.14       2019-05-28 [1] CRAN (R 4.0.3)
    ##  fastmap          1.1.0      2021-01-25 [1] CRAN (R 4.0.3)
    ##  forcats        * 0.5.1      2021-01-27 [1] CRAN (R 4.0.3)
    ##  foreach          1.5.1      2020-10-15 [1] CRAN (R 4.0.3)
    ##  foreign          0.8-80     2020-05-24 [2] CRAN (R 4.0.3)
    ##  Formula        * 1.2-4      2020-10-16 [1] CRAN (R 4.0.3)
    ##  fs               1.5.0      2020-07-31 [1] CRAN (R 4.0.3)
    ##  furrr            0.2.2      2021-01-29 [1] CRAN (R 4.0.3)
    ##  future           1.21.0     2020-12-10 [1] CRAN (R 4.0.3)
    ##  generics         0.1.0      2020-10-31 [1] CRAN (R 4.0.3)
    ##  ggplot2        * 3.3.3      2020-12-30 [1] CRAN (R 4.0.3)
    ##  globals          0.14.0     2020-11-22 [1] CRAN (R 4.0.3)
    ##  glue             1.4.2      2020-08-27 [1] CRAN (R 4.0.3)
    ##  gridExtra      * 2.3        2017-09-09 [1] CRAN (R 4.0.3)
    ##  gsubfn         * 0.7        2018-03-16 [1] CRAN (R 4.0.3)
    ##  gt               0.2.2      2020-08-05 [1] CRAN (R 4.0.4)
    ##  gtable           0.3.0      2019-03-25 [1] CRAN (R 4.0.3)
    ##  gtsummary      * 1.3.7      2021-02-26 [1] CRAN (R 4.0.4)
    ##  haven            2.3.1      2020-06-01 [1] CRAN (R 4.0.3)
    ##  here             1.0.1      2020-12-13 [1] CRAN (R 4.0.4)
    ##  highr            0.8        2019-03-20 [1] CRAN (R 4.0.3)
    ##  Hmisc          * 4.4-2      2020-11-29 [1] CRAN (R 4.0.3)
    ##  hms              1.0.0      2021-01-13 [1] CRAN (R 4.0.3)
    ##  htmlTable        2.1.0      2020-09-16 [1] CRAN (R 4.0.3)
    ##  htmltools        0.5.1.1    2021-01-22 [1] CRAN (R 4.0.3)
    ##  htmlwidgets      1.5.3      2020-12-10 [1] CRAN (R 4.0.3)
    ##  httr             1.4.2      2020-07-20 [1] CRAN (R 4.0.3)
    ##  iterators        1.0.13     2020-10-15 [1] CRAN (R 4.0.3)
    ##  jpeg             0.1-8.1    2019-10-24 [1] CRAN (R 4.0.3)
    ##  jsonlite         1.7.2      2020-12-09 [1] CRAN (R 4.0.3)
    ##  kableExtra     * 1.3.1      2020-10-22 [1] CRAN (R 4.0.3)
    ##  knitr          * 1.31       2021-01-27 [1] CRAN (R 4.0.3)
    ##  lattice        * 0.20-41    2020-04-02 [2] CRAN (R 4.0.3)
    ##  latticeExtra     0.6-29     2019-12-19 [1] CRAN (R 4.0.3)
    ##  lava             1.6.8.1    2020-11-04 [1] CRAN (R 4.0.3)
    ##  lifecycle        0.2.0      2020-03-06 [1] CRAN (R 4.0.3)
    ##  listenv          0.8.0      2019-12-05 [1] CRAN (R 4.0.3)
    ##  lubridate        1.7.9.2    2020-11-13 [1] CRAN (R 4.0.3)
    ##  magrittr         2.0.1      2020-11-17 [1] CRAN (R 4.0.3)
    ##  MASS             7.3-53     2020-09-09 [2] CRAN (R 4.0.3)
    ##  Matrix           1.2-18     2019-11-27 [2] CRAN (R 4.0.3)
    ##  MatrixModels     0.4-1      2015-08-22 [1] CRAN (R 4.0.3)
    ##  matrixStats      0.58.0     2021-01-29 [1] CRAN (R 4.0.3)
    ##  memoise          2.0.0      2021-01-26 [1] CRAN (R 4.0.3)
    ##  mets             1.2.8.1    2020-09-28 [1] CRAN (R 4.0.3)
    ##  modelr           0.1.8      2020-05-19 [1] CRAN (R 4.0.3)
    ##  mstate         * 0.3.1      2020-12-17 [1] CRAN (R 4.0.3)
    ##  multcomp         1.4-15     2020-11-14 [1] CRAN (R 4.0.3)
    ##  munsell          0.5.0      2018-06-12 [1] CRAN (R 4.0.3)
    ##  mvtnorm          1.1-1      2020-06-09 [1] CRAN (R 4.0.3)
    ##  nlme             3.1-149    2020-08-23 [2] CRAN (R 4.0.3)
    ##  nnet             7.3-14     2020-04-26 [2] CRAN (R 4.0.3)
    ##  numDeriv         2016.8-1.1 2019-06-06 [1] CRAN (R 4.0.3)
    ##  openxlsx         4.2.3      2020-10-27 [1] CRAN (R 4.0.4)
    ##  pacman         * 0.5.1      2019-03-11 [1] CRAN (R 4.0.4)
    ##  parallelly       1.23.0     2021-01-04 [1] CRAN (R 4.0.3)
    ##  pec            * 2020.11.17 2020-11-16 [1] CRAN (R 4.0.3)
    ##  pillar           1.4.7      2020-11-20 [1] CRAN (R 4.0.3)
    ##  pkgconfig        2.0.3      2019-09-22 [1] CRAN (R 4.0.3)
    ##  plotrix        * 3.8-1      2021-01-21 [1] CRAN (R 4.0.3)
    ##  png              0.1-7      2013-12-03 [1] CRAN (R 4.0.3)
    ##  polspline        1.1.19     2020-05-15 [1] CRAN (R 4.0.3)
    ##  prodlim        * 2019.11.13 2019-11-17 [1] CRAN (R 4.0.3)
    ##  proto          * 1.0.0      2016-10-29 [1] CRAN (R 4.0.3)
    ##  purrr          * 0.3.4      2020-04-17 [1] CRAN (R 4.0.3)
    ##  quantreg         5.83       2021-01-22 [1] CRAN (R 4.0.3)
    ##  R6               2.5.0      2020-10-28 [1] CRAN (R 4.0.3)
    ##  RColorBrewer     1.1-2      2014-12-07 [1] CRAN (R 4.0.3)
    ##  Rcpp             1.0.6      2021-01-15 [1] CRAN (R 4.0.3)
    ##  readr          * 1.4.0      2020-10-05 [1] CRAN (R 4.0.3)
    ##  readxl           1.3.1      2019-03-13 [1] CRAN (R 4.0.3)
    ##  reprex           1.0.0      2021-01-27 [1] CRAN (R 4.0.3)
    ##  rio            * 0.5.26     2021-03-01 [1] CRAN (R 4.0.4)
    ##  riskRegression * 2020.12.08 2020-12-09 [1] CRAN (R 4.0.3)
    ##  rlang            0.4.10     2020-12-30 [1] CRAN (R 4.0.3)
    ##  rmarkdown        2.6        2020-12-14 [1] CRAN (R 4.0.3)
    ##  rms            * 6.1-0      2020-11-29 [1] CRAN (R 4.0.3)
    ##  rpart            4.1-15     2019-04-12 [2] CRAN (R 4.0.3)
    ##  rprojroot        2.0.2      2020-11-15 [1] CRAN (R 4.0.3)
    ##  rsample        * 0.0.8      2020-09-23 [1] CRAN (R 4.0.3)
    ##  RSQLite        * 2.2.3      2021-01-24 [1] CRAN (R 4.0.3)
    ##  rstudioapi       0.13       2020-11-12 [1] CRAN (R 4.0.3)
    ##  rvest            0.3.6      2020-07-25 [1] CRAN (R 4.0.3)
    ##  sandwich         3.0-0      2020-10-02 [1] CRAN (R 4.0.3)
    ##  scales           1.1.1      2020-05-11 [1] CRAN (R 4.0.3)
    ##  sessioninfo      1.1.1      2018-11-05 [1] CRAN (R 4.0.4)
    ##  SparseM        * 1.78       2019-12-13 [1] CRAN (R 4.0.3)
    ##  sqldf          * 0.4-11     2017-06-28 [1] CRAN (R 4.0.3)
    ##  stringi          1.5.3      2020-09-09 [1] CRAN (R 4.0.3)
    ##  stringr        * 1.4.0      2019-02-10 [1] CRAN (R 4.0.3)
    ##  survAUC        * 1.0-5      2012-09-04 [1] CRAN (R 4.0.3)
    ##  survival       * 3.2-7      2020-09-28 [1] CRAN (R 4.0.3)
    ##  survivalROC    * 1.0.3      2013-01-13 [1] CRAN (R 4.0.3)
    ##  table1         * 1.2.1      2020-11-26 [1] CRAN (R 4.0.3)
    ##  TH.data          1.0-10     2019-01-21 [1] CRAN (R 4.0.3)
    ##  tibble         * 3.0.6      2021-01-29 [1] CRAN (R 4.0.3)
    ##  tidyr          * 1.1.2      2020-08-27 [1] CRAN (R 4.0.3)
    ##  tidyselect       1.1.0      2020-05-11 [1] CRAN (R 4.0.3)
    ##  tidyverse      * 1.3.0      2019-11-21 [1] CRAN (R 4.0.3)
    ##  timereg          1.9.8      2020-10-05 [1] CRAN (R 4.0.3)
    ##  timeROC        * 0.4        2019-12-18 [1] CRAN (R 4.0.3)
    ##  usethis          2.0.1      2021-02-10 [1] CRAN (R 4.0.4)
    ##  vctrs            0.3.6      2020-12-17 [1] CRAN (R 4.0.3)
    ##  viridisLite      0.3.0      2018-02-01 [1] CRAN (R 4.0.3)
    ##  webshot        * 0.5.2      2019-11-22 [1] CRAN (R 4.0.3)
    ##  withr            2.4.1      2021-01-26 [1] CRAN (R 4.0.3)
    ##  xfun             0.20       2021-01-06 [1] CRAN (R 4.0.3)
    ##  xml2             1.3.2      2020-04-23 [1] CRAN (R 4.0.3)
    ##  yaml             2.2.1      2020-02-01 [1] CRAN (R 4.0.3)
    ##  zip              2.1.1      2020-08-27 [1] CRAN (R 4.0.4)
    ##  zoo              1.8-8      2020-05-02 [1] CRAN (R 4.0.3)
    ## 
    ## [1] C:/Users/danie/Documents/R/win-library/4.0
    ## [2] C:/Program Files/R/R-4.0.3/library
