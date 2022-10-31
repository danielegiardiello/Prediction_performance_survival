Performance assessment of survival prediction models - only validation
individual data available
================

- <a href="#goals" id="toc-goals">Goals</a>
- <a href="#installload-packages-and-import-data"
  id="toc-installload-packages-and-import-data">Install/load packages and
  import data</a>
  - <a href="#data-preparation" id="toc-data-preparation">Data
    preparation</a>
- <a
  href="#goal-1-assessing-performance-of-a-developed-survival-model-in-a-new-data"
  id="toc-goal-1-assessing-performance-of-a-developed-survival-model-in-a-new-data">Goal
  1: Assessing performance of a developed survival model in a new data</a>
  - <a
    href="#11-calculate-the-absolute-risk-prediction-at-5-years-in-the-validation-data"
    id="toc-11-calculate-the-absolute-risk-prediction-at-5-years-in-the-validation-data">1.1
    Calculate the absolute risk prediction at 5 years in the validation
    data</a>
  - <a
    href="#12-histograms-of-predictions-with-and-without-the-additional-marker"
    id="toc-12-histograms-of-predictions-with-and-without-the-additional-marker">1.2
    Histograms of predictions with and without the additional marker</a>
  - <a href="#13-discrimination-measures"
    id="toc-13-discrimination-measures">1.3 Discrimination measures</a>
  - <a href="#14-calibration" id="toc-14-calibration">1.4 Calibration</a>
    - <a
      href="#141-mean-calibration---observedexpected-ratio-for-fixed-time-point"
      id="toc-141-mean-calibration---observedexpected-ratio-for-fixed-time-point">1.4.1
      Mean calibration - observed/expected ratio for fixed time point</a>
    - <a href="#142-weak-calibration---calibration-slope-for-fixed-time-point"
      id="toc-142-weak-calibration---calibration-slope-for-fixed-time-point">1.4.2
      Weak calibration - calibration slope for fixed time point</a>
    - <a href="#143-moderate-calibration---fixed-time-point"
      id="toc-143-moderate-calibration---fixed-time-point">1.4.3 Moderate
      calibration - fixed time point</a>
  - <a href="#15-overall-performance-measures"
    id="toc-15-overall-performance-measures">1.5 Overall performance
    measures</a>
- <a href="#goal-2-clinical-utility" id="toc-goal-2-clinical-utility">Goal
  2. Clinical utility</a>
- <a href="#reproducibility-ticket"
  id="toc-reproducibility-ticket">Reproducibility ticket</a>

## Goals

When a risk prediction model has been developed and published in the
literature, individual data that was used during model development are
not always available. In this document, we assume the scenario that a
risk prediction model was already developed and is available in the
literature. We assume that the author(s) developed a risk prediction
model using a Cox proportional hazard regression and provided the model
equation in terms of coefficients and the baseline survival at a fixed
time horizon *t* (e.g. five years).

The goals are:  
1. to assess the prediction performance of a published prediction model
with a time-to-event outcome in a new independent (external) data;  
2. to assess the potential clinical utility of a prediction model with
time-to-event outcome in the new data;

## Install/load packages and import data

First of all, install the R packages essential for the analyses. We
following libraries are needed to achieve the following goals. If you
don’t have them installed, please use install.packages(’‘)
(e.g. install.packages(’survival’)) or use the user-friendly approach if
you are using RStudio.

``` r
# Use pacman to check whether packages are installed, if not load them
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(survival,
               Hmisc,
               pec,
               timeROC,
               riskRegression,
               rms,
               knitr,
               kableExtra,
               tidyverse)

options(show.signif.stars = FALSE)  # display statistical intelligence
palette("Okabe-Ito")  # color-blind friendly  (needs R 4.0)
```

### Data preparation

Outcome and predictors in the new data must be coded as provided in the
model equation of the developed model. In our case study, the
time-to-event outcome should be in years and the predictors should be
categorized exactly as in the developed model.

In the prediction model developed using the Rotterdam data, the data was
administratively censored at 5 years. For this reason, we also
administratively censor the data from patients in the new (validation)
data at 5 years.

<details>
<summary>
Click to expand code
</summary>

``` r
# Validation data
gbsg$ryear <- gbsg$rfstime/365.25
gbsg$rfs   <- gbsg$status           # the GBSG data contains RFS
gbsg$cnode <- cut(gbsg$nodes, 
                  c(-1,0, 3, 51),
                  c("0", "1-3", ">3"))   # categorized node
gbsg$csize <- cut(gbsg$size,  
                  c(-1, 20, 50, 500), #categorized size
                  c("<=20", "20-50", ">50"))
pgr99 <- 1360 #p99 from development data 
gbsg$pgr2 <- pmin(gbsg$pgr, pgr99) # Winsorized value of PGR
nodes99 <- 19 #p99 from development data 
gbsg$nodes2 <- pmin(gbsg$nodes, nodes99) # Winsorized value of continuous nodes
gbsg$grade3 <- as.factor(gbsg$grade)
levels(gbsg$grade3) <- c("1-2", "1-2", "3")

# Restricted cubic spline 
# Continuous nodes
rcs3_nodes <- Hmisc::rcspline.eval(gbsg$nodes2, knots = c(0, 1, 9))
attr(rcs3_nodes, "dim") <- NULL
attr(rcs3_nodes, "knots") <- NULL
gbsg$nodes3 <- rcs3_nodes

# PGR
rcs3_pgr <- Hmisc::rcspline.eval(gbsg$pgr2, knots = c(0, 41, 486))
attr(rcs3_pgr, "dim") <- NULL
attr(rcs3_pgr, "knots") <- NULL
gbsg$pgr3 <- rcs3_pgr


# Much of the analysis will focus on the first 5 years: create
#  data sets that are censored at 5
temp <- survSplit(Surv(ryear, rfs) ~ ., data = gbsg, cut = 5,
                  episode ="epoch")
gbsg5 <- subset(temp, epoch == 1)

# Relevel
gbsg5$cnode <- relevel(gbsg$cnode, "0")
```

</details>

## Goal 1: Assessing performance of a developed survival model in a new data

The performance of a risk prediction models may be evaluated through:

- discrimination: the ability of the model to correctly rank patients
  with and without the outcome by a certain time point. This requires
  the coefficients (or the log of the hazard ratios) of the developed
  Cox prediction model to be evaluated;

- calibration: the agreement between observed and predicted
  probabilities. It additionally requires the baseline (cumulative)
  hazard or survival;

- overall performance measures: a combination of discrimination and
  calibration.

Unfortunately, few publications report the complete baseline
(cumulative) hazard or survival or even the baseline (cumulative) hazard
or survival at fixed time horizon *t*.  
It is common that physicians focus on one or more clinically relevant
time horizons to inform subjects about their risk. We aim to assess the
prediction performance of a risk prediction model with time-to-event
outcome in a new data when information at a fixed time horizon(s) (here
at 5 years) of a developed prediction model were provided.

When the baseline is not available (unfortunately not uncommon in the
literature), only a graphical representation of calibration is possible.
We assume here to know the coefficients *and the baseline survival at 5
years *S*<sub>0</sub>*(t = 5)\* of the developed prediction model.

If the model equation is provided including the coefficients and the
baseline at fixed time point *t* (e.g. 5 years), we could validate the
risk prediction model in our external data. Typically, the model
equation is provided in terms of predicted survival at a fixed time
point *t*.

<img src="https://render.githubusercontent.com/render/math?math=%5Clarge%7B%5Csf%7BS(t)%20%3D%20S_0(t)%5E%7Bexp(%5Cbeta_1X_1%2B%5Cbeta_2X_2%2B%5Cbeta_3X_3%2B%5Ccdots%2B%5Cbeta_pX_p)%7D%7D%3DS_0(t)%5E%7Bexp(PI)%7D%7D">

where:  
*S(t)* is the probability of surviving by time *t*.  
*S*<sub>0</sub>*(t)* is the baseline survival probability by time *t*.  
<img src="https://render.githubusercontent.com/render/math?math=%5Csf%7BPI%20%3D%20%5Cbeta_1X_1%2B%5Ccdots%2B%5Cbeta_pX_p%7D">
is the prognostic index: the combination of the model coefficients and
the value of the predictors.  

In some software, the baseline survival that is reported might relate to
the baseline survival when covariate values are all at the mean value.
Beware of that. See for example, the function `rms::cph()` and
`rms::cph()$center` in the `rms` package and in `survival` package
`help(basehaz)`, especially the argument `centered`. If the centercept
is mentioned in the model equation, this can be used to rescaled the
baseline using some algebraic steps.

<img src="https://render.githubusercontent.com/render/math?math=%5Clarge%7BS(t)%20%3D%20%7BS_%7B0%7D(t)%7D%5E%7Bexp(PI-c)%7D%20%3D%20%5B%7BS_%7B0%7D(t)%7D%5E%7Bexp(-c)%7D%5D%5E%7Bexp(PI)%7D%20%3D%7BS_%7B0%7D(t)_%7Bresc%7D%7D%5E%7Bexp(PI)%7D%7D">

### 1.1 Calculate the absolute risk prediction at 5 years in the validation data

This part must be run for all following parts of the code. After running
this, the user can also focus on one particular performance measure only
(e.g. discrimination).

<details>
<summary>
Click to expand code
</summary>

``` r
# Absolute risk calculation --------------
# Baseline survival at 5 years - basic model
S0.5yr <- 0.804 
# Design matrix of predictors
des_matr <- as.data.frame(model.matrix(~ csize + nodes2 + nodes3 + grade3,                                    data = gbsg5))
des_matr$`(Intercept)` <- NULL
# Coefficients
coef <- c(0.342, 0.574, 0.304, -0.811,  0.362)
# Prognostic index (PI)
gbsg5$PI <- as.vector(as.matrix(des_matr) %*% cbind(coef))
# Estimated absolute risk at 5 years (1 - S(t), 1 - survival at time t)
gbsg5$pred5 <-  as.vector(1 - S0.5yr**exp(gbsg5$PI))

#Absolute risk  at 5 yrs - Extended model with PGR ---------------
# Baseline survival at 5 years - extended model
S0.5yr_pgr <- 0.761 
# Design matrix of predictors
des_matr <- as.data.frame(model.matrix(~ csize + 
                                         nodes2 + nodes3 + grade3 + 
                                         I(pgr2) + I(pgr3), data = gbsg5))
des_matr$`(Intercept)` <- NULL
# Coefficients
coef <- c(0.320, 0.554, 0.305, 
          -0.820, 0.305, -0.003, 0.013)
# Prognostic index (PI)
gbsg5$PI_pgr <- as.vector(as.matrix(des_matr) %*% cbind(coef))
# Absolute risk at 5 years (1 - S(t), 1 - survival at time t)
gbsg5$pred5_pgr <-  as.vector(1 - S0.5yr_pgr**exp(gbsg5$PI_pgr))
```

</details>

### 1.2 Histograms of predictions with and without the additional marker

<details>
<summary>
Click to expand code
</summary>

``` r
# Validation data

par(las = 1)
xlab <- c(paste0('Basic model\nvariance = ', 
                 round(var(gbsg5$pred5), 3)),
          paste0('Extended model with PGR\nvariance = ',
                round(var(gbsg5$pred5_pgr), 3)))
Hmisc::histbackback(gbsg5$pred5,
                    gbsg5$pred5_pgr, 
                    brks = seq(0.01, 0.99, by = 0.02), 
                    xlab = xlab, 
                    ylab = 'Predicted probability')
title("Validation data")
```

</details>

<img src="imgs/02_predsurv/hist_pred-1.png" width="672" style="display: block; margin: auto;" />

### 1.3 Discrimination measures

Discrimination is the ability to differentiate between subjects who have
the outcome by a certain time point and subjects who do not. Concordance
can be assessed over several different time intervals:

- the entire range of the data. Two concordance measures are suggested:

  - Harrell’s C quantifies the degree of concordance as the proportion
    of evaluable pairs where the patient with a longer survival time has
    better predicted survival;

  - Uno’s C uses a time dependent weighting that more fully adjusts for
    censoring;

- a 5 year window corresponding to our target assessment point. Uno’s
  cumulative/dynamic time-dependent Area Under the Curve (AUC) is
  suggested. Uno’s time-dependent AUC summarizes discrimination at
  specific fixed time points. At any time point of interest, *t*, a
  patient is classified as having an event if the patient experienced
  the event between baseline and *t* (5 years in our case study), and as
  a non-event if the patient remained event-free at *t*. The
  time-dependent AUC evaluates whether predicted probabilities were
  higher for cases than for non-cases.

There is some uncertainty in the literature about the original Harrell
formulation versus Uno’s suggestion to re-weight the time scale by the
factor $1/G^2(t)$ where $G$ is the censoring distribution. There is more
detailed information in the concordance vignette found in the survival
package.

For all three measures, values close to 1 indicate good discrimination
ability, while values close to 0.5 indicated poor discrimination
ability.

<details>
<summary>
Click to expand code
</summary>

``` r
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(survival,
               Hmisc,
               pec,
               timeROC)

harrell_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ PI, 
                               gbsg5, 
                               reverse = TRUE)

harrell_C_gbsg5_pgr <- concordance(Surv(ryear, rfs) ~ PI_pgr, 
                               gbsg5, 
                               reverse = TRUE)

# Uno's C
Uno_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ PI, 
                           gbsg5, 
                           reverse = TRUE,
                           timewt = "n/G2")

Uno_C_gbsg5_pgr <- concordance(Surv(ryear, rfs) ~ PI_pgr, 
                           gbsg5, 
                           reverse = TRUE,
                           timewt = "n/G2")
```

</details>
<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
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
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Harrell C - Validation data
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.71
</td>
</tr>
<tr>
<td style="text-align:left;">
Uno C - Validation data
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.69
</td>
</tr>
</tbody>
</table>

Concordance was between 0.64 and 0.68. The extended model slightly
improved discrimination ability compared to the basic model.

<details>
<summary>
Click to expand code
</summary>

``` r
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(survival,
               Hmisc,
               pec,
               timeROC)

# External validation
Uno_gbsg5 <-
  timeROC::timeROC(
    T = gbsg5$ryear, delta = gbsg5$rfs,
    marker = gbsg5$PI,
    cause = 1, weighting = "marginal", times = 4.99,
    iid = TRUE
  )

# External validation with pgr
Uno_gbsg5_pgr <-
  timeROC::timeROC(
    T = gbsg5$ryear, delta = gbsg5$rfs,
    marker = gbsg5$PI_pgr,
    cause = 1, weighting = "marginal", times = 4.99,
    iid = TRUE
  )
# NOTE: if you have a lot of data n > 2000, standard error computation may be really long. In that case, please use bootstrap percentile to calculate confidence intervals.
```

</details>
<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
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
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Uno AUC
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.7
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.76
</td>
</tr>
</tbody>
</table>

The time-dependent AUCs at 5 years in the external validation were 0.68
and 0.70 for the basic and extended model, respectively.

### 1.4 Calibration

Calibration is the agreement between observed outcomes and predicted
probabilities. For example, in survival models, a predicted survival
probability at a fixed time horizon *t* of 80% is considered reliable if
it can be expected that 80 out of 100 will survive among patients who
received a predicted survival probability of 80%. Calibration can be
assessed at a fixed time point (e.g. at 5 years), and globally
(considering the entire range of the data). In addition, different level
of calibration assessment can be estimated according to the level of
information available in the data. When individual data of development
and validation set are available, full assessment of calibration is
possible. Calibration at fixed time point is possible when baseline
hazard at fixed time point and coefficient are available. When only
coefficients are available, assessment of calibration is limited.

In the scenario we consider here, we can evaluate calibration only at
fixed time point *t* (i.e. 5 years) since we may have baseline survival
at time *t* (5 years) and coefficients of the model.

- Mean calibration at a fixed time point can be estimated using the
  Observed versus Expected ratio at time t;

- Weak calibration can be estimated by additionally calculating
  calibration slope.

- Moderate calibration can estimated at a fixed time point using a
  flexible calibration curve, complemented with ICI, E50, E90.

More detailed explanations are available in the paper.

#### 1.4.1 Mean calibration - observed/expected ratio for fixed time point

The mean calibration at fixed time point (e.g. at 5 years) can be
estimated using the Observed versus Expected ratio. The observed is
estimated using the complementary of the Kaplan-Meier curve at the fixed
time point. The expected is estimated using the average predicted risk
of the event at the fixed time point.

<details>
<summary>
Click to expand code
</summary>

``` r
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(survival,
               Hmisc,
               pec,
               timeROC,
               rms)
##  Observed / Expected ratio at time t ------------
# Observed: 1-Kaplan Meier at time (t)
horizon <- 5
obj <- summary(survfit(Surv(ryear, rfs) ~ 1, 
                       data = gbsg5), 
               times = horizon)

OE <- (1 - obj$surv) / mean(gbsg5$pred5)
OE_pgr <- (1 - obj$surv) / mean(gbsg5$pred5_pgr)
```

</details>
<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
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
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
OE ratio
</td>
<td style="text-align:right;">
1.02
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.14
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.12
</td>
</tr>
</tbody>
</table>

Observed and Expected ratio is 1.07 (95% CI: 0.93 - 1.17) for the basic
model and 1.02 (95% CI: 0.91 - 1.14) for the extended model.

#### 1.4.2 Weak calibration - calibration slope for fixed time point

<details>
<summary>
Click to expand code
</summary>

``` r
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(survival,
               Hmisc,
               pec,
               timeROC,
               rms)

# cloglog and center for the basic and extended model
lp.val <- log(-log(1 - gbsg5$pred5))   # lp = cloglog
lp.val_pgr <- log(-log(1 - gbsg5$pred5_pgr)) 
center <- mean(lp.val)  # center
center_pgr <- mean(lp.val_pgr)  # center


### Model with a slope and an intercept
horizon <- 5
f.val <- coxph(Surv(gbsg5$ryear, gbsg5$rfs) ~ lp.val)  
slope <- f.val$coefficients[1]
slope.se <- sqrt(vcov(f.val)[[1, 1]])

f.val_pgr <- coxph(Surv(gbsg5$ryear, gbsg5$rfs) ~ lp.val_pgr)  
slope_pgr <- f.val_pgr$coefficients[1]
slope.se_pgr <- sqrt(vcov(f.val_pgr)[[1, 1]])
```

</details>
<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
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
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Calibration slope
</td>
<td style="text-align:right;">
1.06
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
1.14
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
1.37
</td>
</tr>
</tbody>
</table>

Calibration slope was 1.07 and 1.20 for the basic and extended model,
respectively.

#### 1.4.3 Moderate calibration - fixed time point

Moderate calibration at fixed time point can be assessed using flexible
calibration curve, complemented with ICI, E50, E90 as suggested by
Austin et al.

- Calibration curve is a graphical representation of moderate
  calibration. It shows:

  - on the *x-axis* the predicted survival (or risk) probabilities at a
    fixed time horizon (e.g. at 5 years);

  - on the *y-axis* the observed survival (or risk) probabilities at a
    fixed time horizon (e.g. at 5 years);

  - The 45-degree line indicates perfect calibration. Points below the
    45-degree line indicate that the model overestimates the observed
    risk. If points are above the 45-degree line, the model
    underestimate the observed risk; The observed probabilities
    estimated by the Kaplan-Meier curves (in case of survival) or by the
    complementary of the Kaplan-Meier curves (in case of risk) are
    represented in terms of percentiles of the predicted survival (risk)
    probabilities.

- Integrated Calibration Index (ICI) is the weighted mean of absolute
  difference between smoothed observed proportions and predicted
  probabilities in which observations are weighted by the empirical
  density function of the predicted probabilities;

- E50 and E90 denote the median and the 90th percentile of the absolute
  differences between observed and predicted probabilities of the
  outcome at time *t*;

<details>
<summary>
Click to expand code
</summary>

``` r
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(survival,
               Hmisc,
               pec,
               timeROC,
               rms)

# Calibration plot --------
# Basic model
gbsg5 <- data.frame(gbsg5)
gbsg5$pred.cll <- log(-log(1 - gbsg5$pred5))

# Extended model
gbsg5$pred.cll_pgr <- log(-log(1 - gbsg5$pred5_pgr))


# Estimate actual risk - basic model
vcal <- rms::cph(Surv(ryear, rfs) ~ rcs(pred.cll, 3),
                 x = T,
                 y = T,
                 surv = T,
                 data = gbsg5
) 

# Estimate actual risk - extended model
vcal_pgr <- rms::cph(Surv(ryear, rfs) ~ rcs(pred.cll_pgr, 3),
                     x = T,
                     y = T,
                     surv = T,
                     data = gbsg5
) 


dat_cal <- cbind.data.frame(
  "obs" = 1 - rms::survest(vcal,
                           times = 5,
                           newdata = gbsg5)$surv,
  
  "lower" = 1 - rms::survest(vcal,
                             times = 5,
                             newdata = gbsg5)$upper,
  
  "upper" = 1 - rms::survest(vcal,
                             times = 5,
                             newdata = gbsg5)$lower,
  
  "pred" = as.vector(gbsg5$pred5),
  
  
   "obs_pgr" = 1 - rms::survest(vcal_pgr,
                                times = 5,
                                newdata = gbsg5)$surv,
  
  "lower_pgr" = 1 - rms::survest(vcal_pgr,
                                 times = 5,
                                 newdata = gbsg5)$upper,
  
  "upper_pgr" = 1 - rms::survest(vcal_pgr,
                                 times = 5,
                                 newdata = gbsg5)$lower,
  
  "pred_pgr" = as.vector(gbsg5$pred5_pgr)
  
)


# Flexible calibration curve - basic model
dat_cal <- dat_cal[order(dat_cal$pred), ]

par(xaxs = "i", yaxs = "i", las = 1)
plot(
  dat_cal$pred, 
  dat_cal$obs,
  type = "l", 
  lty = 1, 
  xlim = c(0, 1),
  ylim = c(0, 1), 
  lwd = 2,
  xlab = "Predicted risk from developed model",
  ylab = "Predicted risk from refitted model", bty = "n"
)
lines(dat_cal$pred, 
      dat_cal$lower, 
      type = "l", 
      lty = 2, 
      lwd = 2)
lines(dat_cal$pred, 
      dat_cal$upper,
      type = "l", 
      lty = 2, 
      lwd = 2)
abline(0, 1, lwd = 2, lty = 2, col = 2)
legend("bottomright",
        c("Ideal calibration",
          "Calibration curve based on secondary Cox model",
          "95% confidence interval"),
        col = c(2, 1, 1),
        lty = c(2, 1, 2),
        lwd = c(2, 2, 2),
        bty = "n",
        cex = 0.85)
title("Basic model - validation data ")


# Flexible calibration curve - extended model
dat_cal <- dat_cal[order(dat_cal$pred_pgr), ]
par(xaxs = "i", yaxs = "i", las = 1)
plot(
  dat_cal$pred_pgr, 
  dat_cal$obs_pgr,
  type = "l", 
  lty = 1, 
  xlim = c(0, 1),
  ylim = c(0, 1), 
  lwd = 2,
  xlab = "Predicted risk from developed model",
  ylab = "Predicted risk from refitted model", bty = "n"
)
lines(dat_cal$pred_pgr, 
      dat_cal$lower_pgr, 
      type = "l", 
      lty = 2, 
      lwd = 2)
lines(dat_cal$pred_pgr, 
      dat_cal$upper_pgr,
      type = "l", 
      lty = 2, 
      lwd = 2)
abline(0, 1, lwd = 2, lty = 2, col = 2)
legend("bottomright",
        c("Ideal calibration",
          "Calibration curve based on secondary Cox model",
          "95% confidence interval"),
        col = c(2, 1, 1),
        lty = c(2, 1, 2),
        lwd = c(2, 2, 2),
        bty = "n",
        cex = 0.85)
title("Extended model - validation data ")

# Numerical measures ---------------
# Basic model
absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)

numsum_cph <- c(
  "ICI" = mean(absdiff_cph),
  setNames(quantile(absdiff_cph, c(0.5, 0.9)), c("E50", "E90"))
)

# Extended model ------
absdiff_cph_pgr <- abs(dat_cal$pred_pgr - dat_cal$obs_pgr)

numsum_cph_pgr <- c(
  "ICI" = mean(absdiff_cph_pgr),
  setNames(quantile(absdiff_cph_pgr, c(0.5, 0.9)), c("E50", "E90"))
)
```

</details>

<img src="imgs/02_predsurv/cal_rcs_metrics-1.png" width="576" style="display: block; margin: auto;" /><img src="imgs/02_predsurv/cal_rcs_metrics-2.png" width="576" style="display: block; margin: auto;" />

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
ICI
</th>
<th style="text-align:right;">
E50
</th>
<th style="text-align:right;">
E90
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
External data
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
External data + PGR
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
</tbody>
</table>

In the validation, ICI at 5 years was 0.03 and 0.02 for the basic and
extended model, respectively.

### 1.5 Overall performance measures

Two overall performance measures are proposed for prediction models with
a survival outcome:

- Brier score: it is the mean squared difference between observed event
  indicators and predicted risks at a fixed time point (e.g. at 5
  years), lower is better;

- Scaled Brier score, also known as Index of Prediction Accuracy (IPA):
  it improves interpretability by scaling the Brier Score. It is the
  decrease in Brier compared to a null model, expressed as a percentage,
  higher is better.

<details>
<summary>
Click to expand code
</summary>

``` r
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(riskRegression)

# For basic model
score_gbsg5 <-
  riskRegression::Score(list("Validation" = gbsg5$pred5),
                        formula = Surv(ryear, rfs) ~ 1, 
                        data = gbsg5, 
                        conf.int = TRUE, 
                        times = 4.99,
                        cens.model = "km", 
                        metrics = "brier",
                        summary = "ipa")

# For extended code including PGR
score_gbsg5_pgr <-
  riskRegression::Score(list("Validation with PGR" = gbsg5$pred5_pgr),
                        formula = Surv(ryear, rfs) ~ 1, 
                        data = gbsg5, 
                        conf.int = TRUE, 
                        times = 4.99,
                        cens.model = "km", 
                        metrics = "brier",
                        summary = "ipa")

# Extra: bootstrap confidence intervals for IPA ------
B <- 100
horizon <- 4.99
data_boot <- list()
b_score <- list()
b_score_pgr <- list()
ipa_boot <- c()
ipa_boot_pgr <- c()

for (j in 1:B) {
  data_boot[[j]] <- gbsg5[sample(nrow(gbsg5), replace = TRUE), ]
  
  b_score[[j]] <- riskRegression::Score(
    list("Validation" =  data_boot[[j]]$pred5),
    formula = Surv(ryear, rfs) ~ 1,
    cens.model = "km", 
    data = data_boot[[j]], 
    conf.int = FALSE, 
    times = horizon,
    metrics = c("brier"),
    summary = c("ipa")
  )

   b_score_pgr[[j]] <- riskRegression::Score(
    list("Validation with PGR" =  data_boot[[j]]$pred5_pgr),
    formula = Surv(ryear, rfs) ~ 1,
    cens.model = "km", 
    data = data_boot[[j]], 
    conf.int = FALSE, 
    times = horizon,
    metrics = c("brier"),
    summary = c("ipa")
  )
   
   ipa_boot[j] <- b_score[[j]]$Brier$score[model == "Validation"][["IPA"]]
   ipa_boot_pgr[j] <- b_score_pgr[[j]]$Brier$score[model == "Validation with PGR"][["IPA"]]
}
```

</details>
<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
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
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Brier
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.23
</td>
</tr>
<tr>
<td style="text-align:left;">
Scaled Brier
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.19
</td>
</tr>
</tbody>
</table>

Overall performances for the basic and extended model are quite similar
in the validation data.

## Goal 2. Clinical utility

Discrimination and calibration measures are essential to assess the
prediction performance but insufficient to evaluate the potential
clinical utility of a risk prediction model for decision making.
Clinical utility assessment evaluates whether a prediction model helps
to improve decision making.  
Clinical utility is measured by the net benefit that includes the number
of true positives and the number of false positives. The true positives
reflect the benefit of correctly predicting patients who will experience
an event during a given time horizon, giving the opportunity to use
additional interventions such as additional treatments, personalized
follow-up or additional surgeries. The false positives represent the
harms of incorrectly predicting events possibly leading to unnecessary
additional interventions. Net benefit is the number of true positives
classifications minus the false positives classifications weighted by a
factor related to the harm of failing to predict an event versus and the
harm of falsely predicting an event. The weighting is derived from a
predefined threshold probability using a defined time horizon (for
example 5 years since diagnosis). For example, a threshold of 10%
implies that additional interventions for 10 patients of whom one would
have experienced the event in the next 5 years if untreated is
acceptable (thus treating 9 patients unnecessarily). This strategy is
compared with the strategies of treating all and treating none of the
patients. If overtreatment is harmful, a higher threshold should be
used.  
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
*X=1* when the predicted probability at time *t* exceeds *p*<sub>t</sub>

And the decision curve is calculated as follows:

1.  Choose a time horizon (in this case 5 years);
2.  Specify a risk threshold which reflects the ratio between harms and
    benefit of an additional intervention;
3.  Calculate the number of true positives and false positives given the
    threshold specified in (2);
4.  Calculate the net benefit of the survival model;
5.  Plot net benefit on the *y-axis* against the risk threshold on the
    *x-axis*;
6.  Repeat steps 2-4 for each model consideration;
7.  Repeat steps 2-4 for the strategy assuming all patients are treated;
8.  Draw a straight line parallel to the *x-axis* at y=0 representing
    the net benefit associated with the strategy assuming that none of
    the patients are treated.

We smoothed the decision curves based on the risk prediction models to
reduce the visual impact of random noise using `stats::smooth()`
function.

<details>
<summary>
Click to expand code
</summary>

``` r
# External data
# Run decision curve analysis

# Model without PGR
gbsg5 <- as.data.frame(gbsg5)
dca_gbsg5 <- stdca(
  data = gbsg5, outcome = "rfs", ttoutcome = "ryear",
  timepoint = 5, predictors = "pred5", xstop = 1.0,
  ymin = -0.01, graph = FALSE
)
# Model with PGR
dca_gbsg5_pgr <- stdca(
  data = gbsg5, outcome = "rfs", ttoutcome = "ryear",
  timepoint = 5, predictors = "pred5_pgr", xstop = 1,
  ymin = -0.01, graph = FALSE
)

# Smoothing DCA without PGR
dca_gbsg5_smooth <- smooth(dca_gbsg5$net.benefit$pred5
                           [!is.na(dca_gbsg5$net.benefit$pred5)],
                           twiceit = TRUE)
dca_gbsg5_smooth <- c(dca_gbsg5_smooth, 
                      rep(NA, sum(is.na(dca_gbsg5$net.benefit$pred5))))

# Smoothing DCA with PGR
dca_gbsg5_pgr_smooth <- smooth(dca_gbsg5_pgr$net.benefit$pred5_pgr
                            [!is.na(dca_gbsg5_pgr$net.benefit$pred5_pgr)],
                               twiceit = TRUE)
dca_gbsg5_pgr_smooth <- c(dca_gbsg5_pgr_smooth,
                  rep(NA,sum(is.na(dca_gbsg5_pgr$net.benefit$pred5_pgr))))

# Decision curves plot
par(xaxs = "i", yaxs = "i", las = 1)
plot(dca_gbsg5$net.benefit$threshold,
  dca_gbsg5_smooth,
  type = "l", 
  lwd = 3, 
  lty = 2,
  xlab = "Threshold probability in %", 
  ylab = "Net Benefit",
  xlim = c(0, 1), 
  ylim = c(-0.10, 0.60), 
  bty = "n",
  cex.lab = 1.2, 
  cex.axis = 1,
  col = 4
)
lines(dca_gbsg5$net.benefit$threshold, 
      dca_gbsg5$net.benefit$none, 
      type = "l", 
      lwd = 3, 
      lty = 4,
      col = 8)
lines(dca_gbsg5$net.benefit$threshold, 
      dca_gbsg5$net.benefit$all, 
      type = "l", 
      lwd = 3, 
      col = 2)
lines(dca_gbsg5_pgr$net.benefit$threshold,
      dca_gbsg5_pgr_smooth, 
      type = "l", 
      lwd = 3, 
      lty = 5,
      col = 7)
legend("topright",
  c(
    "Treat All",
    "Original model",
    "Original model + PGR",
    "Treat None"
  ),
  lty = c(1, 2, 5, 4), lwd = 3, 
  col = c(2, 4, 7, 8),
  bty = "n"
)
title("Validation set", cex = 1.5)
```

</details>

    ## [1] "pred5: No observations with risk greater than 93%, and therefore net benefit not calculable in this range."

    ## [1] "pred5_pgr: No observations with risk greater than 95%, and therefore net benefit not calculable in this range."

<img src="imgs/02_predsurv/dca-1.png" width="672" style="display: block; margin: auto;" />

Moreover, net benefit can be defined in terms of reduction of avoidable
interventions (e.g adjuvant chemotherapy per 100 patients) by:

<img src="https://render.githubusercontent.com/render/math?math=%5Chuge%7B%5Cfrac%7BNB_%7Bmodel%7D%20-%20NB_%7Ball%7D%7D%7B(p_t%2F%20(1-p_t))%7D*100%7D%0A">

where *NB*<sub>model</sub> is the net benefit of the prediction model,
*NB*<sub>all</sub> is the net benefit of the strategy treat all and
*p*<sub>*t*</sub> is the risk threshold.

## Reproducibility ticket

``` r
sessioninfo::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.2.1 (2022-06-23 ucrt)
    ##  os       Windows 10 x64 (build 22000)
    ##  system   x86_64, mingw32
    ##  ui       RTerm
    ##  language (EN)
    ##  collate  English_Netherlands.utf8
    ##  ctype    English_Netherlands.utf8
    ##  tz       Europe/Berlin
    ##  date     2022-10-31
    ##  pandoc   2.19.2 @ C:/Program Files/RStudio/bin/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package        * version    date (UTC) lib source
    ##  assertthat       0.2.1      2019-03-21 [1] CRAN (R 4.2.1)
    ##  backports        1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
    ##  base64enc        0.1-3      2015-07-28 [1] CRAN (R 4.2.0)
    ##  broom            1.0.1      2022-08-29 [1] CRAN (R 4.2.1)
    ##  cellranger       1.1.0      2016-07-27 [1] CRAN (R 4.2.1)
    ##  checkmate        2.1.0      2022-04-21 [1] CRAN (R 4.2.1)
    ##  cli              3.4.1      2022-09-23 [1] CRAN (R 4.2.1)
    ##  cluster          2.1.3      2022-03-28 [2] CRAN (R 4.2.1)
    ##  cmprsk           2.2-11     2022-01-06 [1] CRAN (R 4.2.1)
    ##  codetools        0.2-18     2020-11-04 [2] CRAN (R 4.2.1)
    ##  colorspace       2.0-3      2022-02-21 [1] CRAN (R 4.2.1)
    ##  crayon           1.5.2      2022-09-29 [1] CRAN (R 4.2.1)
    ##  data.table       1.14.2     2021-09-27 [1] CRAN (R 4.2.1)
    ##  DBI              1.1.3      2022-06-18 [1] CRAN (R 4.2.1)
    ##  dbplyr           2.2.1      2022-06-27 [1] CRAN (R 4.2.1)
    ##  deldir           1.0-6      2021-10-23 [1] CRAN (R 4.2.0)
    ##  digest           0.6.29     2021-12-01 [1] CRAN (R 4.2.1)
    ##  dplyr          * 1.0.10     2022-09-01 [1] CRAN (R 4.2.1)
    ##  ellipsis         0.3.2      2021-04-29 [1] CRAN (R 4.2.1)
    ##  evaluate         0.17       2022-10-07 [1] CRAN (R 4.2.1)
    ##  fansi            1.0.3      2022-03-24 [1] CRAN (R 4.2.1)
    ##  fastmap          1.1.0      2021-01-25 [1] CRAN (R 4.2.1)
    ##  forcats        * 0.5.2      2022-08-19 [1] CRAN (R 4.2.1)
    ##  foreach          1.5.2      2022-02-02 [1] CRAN (R 4.2.1)
    ##  foreign          0.8-82     2022-01-16 [2] CRAN (R 4.2.1)
    ##  Formula        * 1.2-4      2020-10-16 [1] CRAN (R 4.2.0)
    ##  fs               1.5.2      2021-12-08 [1] CRAN (R 4.2.1)
    ##  future           1.28.0     2022-09-02 [1] CRAN (R 4.2.1)
    ##  future.apply     1.9.1      2022-09-07 [1] CRAN (R 4.2.1)
    ##  gargle           1.2.1      2022-09-08 [1] CRAN (R 4.2.1)
    ##  generics         0.1.3      2022-07-05 [1] CRAN (R 4.2.1)
    ##  ggplot2        * 3.3.6      2022-05-03 [1] CRAN (R 4.2.1)
    ##  globals          0.16.1     2022-08-28 [1] CRAN (R 4.2.1)
    ##  glue             1.6.2      2022-02-24 [1] CRAN (R 4.2.1)
    ##  googledrive      2.0.0      2021-07-08 [1] CRAN (R 4.2.1)
    ##  googlesheets4    1.0.1      2022-08-13 [1] CRAN (R 4.2.1)
    ##  gridExtra        2.3        2017-09-09 [1] CRAN (R 4.2.1)
    ##  gtable           0.3.1      2022-09-01 [1] CRAN (R 4.2.1)
    ##  haven            2.5.1      2022-08-22 [1] CRAN (R 4.2.1)
    ##  here             1.0.1      2020-12-13 [1] CRAN (R 4.2.1)
    ##  highr            0.9        2021-04-16 [1] CRAN (R 4.2.1)
    ##  Hmisc          * 4.7-1      2022-08-15 [1] CRAN (R 4.2.1)
    ##  hms              1.1.2      2022-08-19 [1] CRAN (R 4.2.1)
    ##  htmlTable        2.4.1      2022-07-07 [1] CRAN (R 4.2.1)
    ##  htmltools        0.5.3      2022-07-18 [1] CRAN (R 4.2.1)
    ##  htmlwidgets      1.5.4      2021-09-08 [1] CRAN (R 4.2.1)
    ##  httr             1.4.4      2022-08-17 [1] CRAN (R 4.2.1)
    ##  interp           1.1-3      2022-07-13 [1] CRAN (R 4.2.1)
    ##  iterators        1.0.14     2022-02-05 [1] CRAN (R 4.2.1)
    ##  jpeg             0.1-9      2021-07-24 [1] CRAN (R 4.2.0)
    ##  jsonlite         1.8.2      2022-10-02 [1] CRAN (R 4.2.1)
    ##  kableExtra     * 1.3.4      2021-02-20 [1] CRAN (R 4.2.1)
    ##  knitr          * 1.40       2022-08-24 [1] CRAN (R 4.2.1)
    ##  lattice        * 0.20-45    2021-09-22 [2] CRAN (R 4.2.1)
    ##  latticeExtra     0.6-30     2022-07-04 [1] CRAN (R 4.2.1)
    ##  lava             1.6.10     2021-09-02 [1] CRAN (R 4.2.1)
    ##  lifecycle        1.0.3      2022-10-07 [1] CRAN (R 4.2.1)
    ##  listenv          0.8.0      2019-12-05 [1] CRAN (R 4.2.1)
    ##  lubridate        1.8.0      2021-10-07 [1] CRAN (R 4.2.1)
    ##  magrittr         2.0.3      2022-03-30 [1] CRAN (R 4.2.1)
    ##  MASS             7.3-57     2022-04-22 [2] CRAN (R 4.2.1)
    ##  Matrix           1.5-1      2022-09-13 [1] CRAN (R 4.2.1)
    ##  MatrixModels     0.5-1      2022-09-11 [1] CRAN (R 4.2.1)
    ##  mets             1.3.1      2022-10-02 [1] CRAN (R 4.2.1)
    ##  modelr           0.1.9      2022-08-19 [1] CRAN (R 4.2.1)
    ##  multcomp         1.4-20     2022-08-07 [1] CRAN (R 4.2.1)
    ##  munsell          0.5.0      2018-06-12 [1] CRAN (R 4.2.1)
    ##  mvtnorm          1.1-3      2021-10-08 [1] CRAN (R 4.2.0)
    ##  nlme             3.1-157    2022-03-25 [2] CRAN (R 4.2.1)
    ##  nnet             7.3-17     2022-01-16 [2] CRAN (R 4.2.1)
    ##  numDeriv         2016.8-1.1 2019-06-06 [1] CRAN (R 4.2.0)
    ##  pacman         * 0.5.1      2019-03-11 [1] CRAN (R 4.2.1)
    ##  parallelly       1.32.1     2022-07-21 [1] CRAN (R 4.2.1)
    ##  pec            * 2022.05.04 2022-05-04 [1] CRAN (R 4.2.1)
    ##  pillar           1.8.1      2022-08-19 [1] CRAN (R 4.2.1)
    ##  pkgconfig        2.0.3      2019-09-22 [1] CRAN (R 4.2.1)
    ##  png              0.1-7      2013-12-03 [1] CRAN (R 4.2.0)
    ##  polspline        1.1.20     2022-04-25 [1] CRAN (R 4.2.0)
    ##  prodlim        * 2019.11.13 2019-11-17 [1] CRAN (R 4.2.1)
    ##  purrr          * 0.3.5      2022-10-06 [1] CRAN (R 4.2.1)
    ##  quantreg         5.94       2022-07-20 [1] CRAN (R 4.2.1)
    ##  R6               2.5.1      2021-08-19 [1] CRAN (R 4.2.1)
    ##  RColorBrewer     1.1-3      2022-04-03 [1] CRAN (R 4.2.0)
    ##  Rcpp             1.0.9      2022-07-08 [1] CRAN (R 4.2.1)
    ##  readr          * 2.1.3      2022-10-01 [1] CRAN (R 4.2.1)
    ##  readxl           1.4.1      2022-08-17 [1] CRAN (R 4.2.1)
    ##  reprex           2.0.2      2022-08-17 [1] CRAN (R 4.2.1)
    ##  riskRegression * 2022.09.23 2022-09-26 [1] CRAN (R 4.2.1)
    ##  rlang            1.0.6      2022-09-24 [1] CRAN (R 4.2.1)
    ##  rmarkdown        2.17       2022-10-07 [1] CRAN (R 4.2.1)
    ##  rms            * 6.3-0      2022-04-22 [1] CRAN (R 4.2.1)
    ##  rpart            4.1.16     2022-01-24 [2] CRAN (R 4.2.1)
    ##  rprojroot        2.0.3      2022-04-02 [1] CRAN (R 4.2.1)
    ##  rstudioapi       0.14       2022-08-22 [1] CRAN (R 4.2.1)
    ##  rvest            1.0.3      2022-08-19 [1] CRAN (R 4.2.1)
    ##  sandwich         3.0-2      2022-06-15 [1] CRAN (R 4.2.1)
    ##  scales           1.2.1      2022-08-20 [1] CRAN (R 4.2.1)
    ##  sessioninfo      1.2.2      2021-12-06 [1] CRAN (R 4.2.1)
    ##  SparseM        * 1.81       2021-02-18 [1] CRAN (R 4.2.0)
    ##  stringi          1.7.8      2022-07-11 [1] CRAN (R 4.2.1)
    ##  stringr        * 1.4.1      2022-08-20 [1] CRAN (R 4.2.1)
    ##  survival       * 3.3-1      2022-03-03 [2] CRAN (R 4.2.1)
    ##  svglite          2.1.0      2022-02-03 [1] CRAN (R 4.2.1)
    ##  systemfonts      1.0.4      2022-02-11 [1] CRAN (R 4.2.1)
    ##  TH.data          1.1-1      2022-04-26 [1] CRAN (R 4.2.1)
    ##  tibble         * 3.1.8      2022-07-22 [1] CRAN (R 4.2.1)
    ##  tidyr          * 1.2.1      2022-09-08 [1] CRAN (R 4.2.1)
    ##  tidyselect       1.2.0      2022-10-10 [1] CRAN (R 4.2.1)
    ##  tidyverse      * 1.3.2      2022-07-18 [1] CRAN (R 4.2.1)
    ##  timereg          2.0.2      2022-04-11 [1] CRAN (R 4.2.1)
    ##  timeROC        * 0.4        2019-12-18 [1] CRAN (R 4.2.1)
    ##  tzdb             0.3.0      2022-03-28 [1] CRAN (R 4.2.1)
    ##  utf8             1.2.2      2021-07-24 [1] CRAN (R 4.2.1)
    ##  vctrs            0.4.2      2022-09-29 [1] CRAN (R 4.2.1)
    ##  viridisLite      0.4.1      2022-08-22 [1] CRAN (R 4.2.1)
    ##  webshot          0.5.4      2022-09-26 [1] CRAN (R 4.2.1)
    ##  withr            2.5.0      2022-03-03 [1] CRAN (R 4.2.1)
    ##  xfun             0.33       2022-09-12 [1] CRAN (R 4.2.1)
    ##  xml2             1.3.3      2021-11-30 [1] CRAN (R 4.2.1)
    ##  yaml             2.3.5      2022-02-21 [1] CRAN (R 4.2.1)
    ##  zoo              1.8-11     2022-09-17 [1] CRAN (R 4.2.1)
    ## 
    ##  [1] C:/Users/danie/AppData/Local/R/win-library/4.2
    ##  [2] C:/Program Files/R/R-4.2.1/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
