Performance assessment of survival prediction models - simplified code
================

-   [Goals](#goals)
    -   [Set up - load packages and import
        data](#set-up---load-packages-and-import-data)
-   [Goal 1 - Develop a risk prediction model with a time to event
    outcome](#goal-1---develop-a-risk-prediction-model-with-a-time-to-event-outcome)
    -   [1.1 Model development - fit the risk prediction
        models](#11-model-development---fit-the-risk-prediction-models)
-   [Goal 2 - Assessing performance in survival prediction
    models](#goal-2---assessing-performance-in-survival-prediction-models)
    -   [2.1 Discrimination measures](#21-discrimination-measures)
    -   [2.2 Calibration](#22-calibration)
    -   [2.2.1 Observed Expected ratio](#221-observed-expected-ratio)
    -   [2.2.2 Calibration plot using restricted cubic
        splines](#222-calibration-plot-using-restricted-cubic-splines)
    -   [2.3 Overall performance
        measures](#23-overall-performance-measures)
-   [Goal 3 - Clinical utility](#goal-3---clinical-utility)
-   [Additional notes](#additional-notes)
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

### Set up - load packages and import data

Please run the following code to set up the data used in the following
document. The following libraries are needed to achieve the following
goals, if you have not them installed, please use install.packages(’‘)
(e.g. install.packages(’survival’)) or use the user-friendly approach if
you are using RStudio.

We loaded the development (rotterdam) and the validation data (gbsg)
from survival package. The Rotterdam breast cancer data was used to
predict the risk of recurrence or death using size, stage and tumor size
as predictors. These three predictors were used in the Nottingham
Prognostic Index, one of the most popular indeces to determine prognosis
following surgery of breast cancer.  
The Germany Breast Cancer Study Group data was used as an external
validation of the model developed in the Rotterdam breast cancer data.
The prediction model will be then extended using the progesterone (PGR)
marker measured at primary surgery.  
The improvement in prediction performance will be evaluated internally
in the Rotterdam data (development data) and in German Breast Cancer
Study data (validation data).

## Goal 1 - Develop a risk prediction model with a time to event outcome

Prediction models are useful to provide the estimated probability of a
specific outcome using personal information. In many studies, especially
in medicine, the main outcome under assessment is the time to an event
of interest defined generally as survival time. Prognostic models for
survival end points, such as recurrence or progression of disease, need
to account for drop out during follow-up. Patients who have not
experienced the event of interest are censored observations. Cox
regression analysis is the most popular statistical model to deal with
such data in oncology and other medical research.

### 1.1 Model development - fit the risk prediction models

We develop the risk prediction model in the development data considering
the first 5-year follow-up to minimize the violation of proportional
hazard including size, node and grade. We also administratively censored
the validation data at 5 years.

<details>
<summary>
Click to expand code
</summary>

``` r
# Libraries needed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(survival,
              Hmisc,
              pec)

# Fit the model without PGR
efit1 <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade3,
  data = rott5, 
  x = T, 
  y = T)
efit1

# Baseline at 5 years
bh <- basehaz(efit1, centered = FALSE) # uncentered
bh$surv <- exp(-bh$hazard) # baseline survival
S0_t5 <- bh$surv[bh$time == 5] 
# NOTE: this can be used to calculate S(t = 5) = S0(t = 5)**exp(X*beta)

# The model with additional PGR marker
efit1_pgr  <- update(efit1, . ~ . + pgr2 + pgr3)
```

</details>

    ## Call:
    ## coxph(formula = Surv(ryear, rfs) ~ csize + cnode + grade3, data = rott5, 
    ##     x = T, y = T)
    ## 
    ##               coef exp(coef) se(coef)      z        p
    ## csize20-50 0.38342   1.46729  0.06504  5.895 3.74e-09
    ## csize>50   0.66355   1.94167  0.09126  7.271 3.57e-13
    ## cnode1-3   0.35998   1.43330  0.07534  4.778 1.77e-06
    ## cnode>3    1.06278   2.89440  0.07035 15.108  < 2e-16
    ## grade33    0.37477   1.45466  0.07130  5.256 1.47e-07
    ## 
    ## Likelihood ratio test=483.6  on 5 df, p=< 2.2e-16
    ## n= 2982, number of events= 1275

The coefficients of the models indicated that higher size, higher number
of positive lymph nodes and higher grade is more associate with poorer
prognosis. The association of the progesterone marker and the outcome is
non-linear as investigated previously.

## Goal 2 - Assessing performance in survival prediction models

The performance of a risk prediction models may be evaluated through:

-   discrimination: the ability of the model to identify patients with
    and without the outcome. It requires the coefficients (or the log of
    the hazard ratios) of the developed risk prediction model to be
    evaluated.

-   calibration: the agreement between observed and predicted
    probabilities. It requires the baseline (cumulative) hazard or
    survival.

-   overall performance measures: as a combination of discrimination and
    calibration and/or as a measure of the explained variation;

Unfortunately, only few publications report the complete baseline
(cumulative) hazard or survival or even the baseline (cumulative) hazard
or survival at fixed time horizon *t*. If we had both individual data of
the development and validation, a complete assessment of discrimination
and calibration would be possible. We could evaluate the prediction
performance of a risk prediction model at a fixed time horizon(s) *t*
and for the complete follow-up time. In risk prediction, physicians
typically focus on one or more clinically relevant time horizons to
inform subjects about their risk. For this reason, according to
information available, different levels of validation assessment are
possible. Here we aim to assess the prediction performance of a risk
prediction model with time-to-event outcome in case all individual data
are available and in case of only the model equation of a fixed time
horizon (i.e. at 5 years) is provided including the baseline survival.

### 2.1 Discrimination measures

Discrimination is the ability to differentiate between subjects who have
the outcome and subjects who do not. Concordance can be assessed over
several different time intervals:

-   the entire range of the data. Two concordance measures are
    suggested:

    -   Harrell’s C quantifies the degree of concordance as the
        proportion of such pairs where the patient with a longer
        survival time has better predicted survival;

    -   Uno’s C uses a time dependent weighting that more fully adjusts
        for censoring;

-   a 5 year window corresponding to our target assessment point. Uno’s
    time-dependent Area Under the Curve (AUC) is suggested. Uno’s
    time-dependent AUC summarizes discrimination at specific fixed time
    points. At any time point of interest, *t*, a patient is classified
    as having an event if the patient experienced the event between
    baseline and *t* (5 years in our case study), and as a non-event if
    the patient remained event-free at *t*. The time-dependent AUC
    evaluates whether predicted probabilities were higher for cases than
    for non-cases.

Clearly the last of these is most relevant.

This is easy to compute using the concordance function in the survival
package. There is some uncertainty in the literature about the original
Harrell formulation versus Uno’s suggestion to re-weigh the time scale
by the factor 1/*G*<sup>2</sup>(*t*) where *G* is the censoring
distribution. There is more detailed information in the concordance
vignette found in the survival package.

We also propose to calculate Uno’s time-dependent AUC at a specific time
horizon *t*.  
More explanations and details are in the paper.

The time horizon to calculate the time-dependent measures was set to 5
years. Values close to 1 indicate good discrimination ability, while
values close to 0.5 indicated poor discrimination ability.  
We used the time horizon at 4.99 and not 5 years since non-cases are
considered patients at risk after the time horizon and we
administratively censored at 5 years to minimize the violation of PH
assumption.

<details>
<summary>
Click to expand code
</summary>

``` r
# Libraries needed
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(survival,
               Hmisc,
               pec,
               timeROC)

# Add linear predictor in the validation set
gbsg5$lp <- predict(efit1, newdata = gbsg5)

### Harrell and Uno's concordance index 
# Harrell's C


## Validation data
# Harrell's C
harrell_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ lp, 
                               gbsg5, 
                               reverse = TRUE)
# Uno's C
Uno_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ lp, 
                           gbsg5, 
                           reverse = TRUE,
                           timewt = "n/G2")
```

</details>

    ##                              Estimate Lower .95 Upper .95
    ## Harrell C - Validation data      0.65      0.62      0.68
    ## Uno C - Validation data          0.64      0.61      0.67

Harrell C and Uno C were 0.65 and 0.64, respectively.

    ##   Uno AUC Lower .95 Upper .95 
    ##      0.69      0.63      0.74

The time-dependent AUCs at 5 years were in the external validation was
0.69.

### 2.2 Calibration

Calibration is the agreement between observed outcomes and predicted
probabilities. For example, in survival models, a predicted survival
probability at a fixed time horizon *t* of 80% is considered reliable if
it can be expected that 80 out of 100 will survive among patients who
received a predicted survival probability of 80%.

Calibration is measured by:

-   Observed and Expected ratio at time horizon (*t*):

    -   the number of observed events (per 100) is calculated as one
        minus the Kaplan-Meier curve at time *t*;

    -   the number of expected events (per 100) is calculated as the
        mean of the predicted risk at time *t*;

    -   Confidence intervals are calculated using the Normal
        approximation of the Poisson distribution.

-   Calibration plot: it is a graphical representation of calibration.
    It shows:

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

### 2.2.1 Observed Expected ratio

We calculate the observed/ expected ratio (OE) at 5 years in the
development and validation data. In the development data the OE should
be (close to) 1.

<details>
<summary>
Click to expand code
</summary>

``` r
# Libraries needed
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(survival,
               Hmisc)

# Observed / Expected ratio
t_horizon <- 5

# Observed
obj <- summary(survfit(
  Surv(ryear, rfs) ~ 1, 
  data = gbsg5),
  times = t_horizon)

obs_t <- 1 - obj$surv

# Predicted risk 
gbsg5$pred <- 1 - predictSurvProb(efit1, 
                                  newdata = gbsg5,
                                  times = t_horizon)
# Expected
exp_t <- mean(gbsg5$pred)

# Observed / Expected ratio
OE_t <- obs_t / exp_t

alpha <- .05
OE_summary <- c(
  "OE" = OE_t,
  "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
  "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
)

OE_summary
```

</details>

    ##        OE     2.5 %    97.5 % 
    ## 1.0444489 0.9299645 1.1730270

Observed and expected ratio was 1.04.

### 2.2.2 Calibration plot using restricted cubic splines

Calibration plots of the external validation data with and without PGR
are calculated and shown using restricted cubic splines.  
The interpretation of the calibration plot was provided in the section
2.2 reported above, in the corresponding paper and in the literature
provided in the paper and at the end of this document.

<details>
<summary>
Click to expand code
</summary>

``` r
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(survival,
               Hmisc,
               rms)

gbsg5$pred <- 1 - predictSurvProb(efit1, 
                                  newdata = gbsg5, 
                                  times = 5)
gbsg5$pred.cll <- log(-log(1 - gbsg5$pred))


# Estimate actual risk
vcal <- cph(Surv(ryear, rfs) ~ rcs(pred.cll, 3),
            x = T,
            y = T,
            surv = T,
            data = gbsg5
) 

dat_cal <- cbind.data.frame(
  "obs" = 1 - survest(vcal, 
                      times = 5, 
                      newdata = gbsg5)$surv,
  
  "lower" = 1 - survest(vcal, 
                        times = 5, 
                        newdata = gbsg5)$upper,
  
  "upper" = 1 - survest(vcal, 
                        times = 5, 
                        newdata = gbsg5)$lower,
  "pred" = gbsg5$pred
)

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
  xlab = "Predicted probability",
  ylab = "Observed probability", bty = "n"
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
abline(0, 1, lwd = 2, lty = 2, col = "red")

# Numerical measures
absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)

numsum_cph <- c(
  "ICI" = mean(absdiff_cph),
  setNames(quantile(absdiff_cph, c(0.5, 0.9)), c("E50", "E90")),
  "Emax" = max(absdiff_cph)
)
numsum_cph
```

</details>

<img src="imgs/01_predsurv_simplified/cal_rcs-1.png" width="672" style="display: block; margin: auto;" />

    ##        ICI        E50        E90       Emax 
    ## 0.02724664 0.02973399 0.06101062 0.06912009

Good calibration was estimated using calibration plot and calibration
measures.

### 2.3 Overall performance measures

We calculate the Brier Score and the Index of Prediction Accuracy (IPA,
the scaled Brier) as a overall performance measure.

We calculate the overall performance measures: Brier score, Scaled Brier
(IPA) and the corresponding confidence intervals.

<details>
<summary>
Click to expand code
</summary>

``` r
# Libraries needed
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(survival,
               Hmisc,
               pec)

# Fit the model without PGR
efit1 <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade3,
  data = rott5, 
  x = T, 
  y = T)

# The model with additional PGR marker
efit1_pgr  <- update(efit1, . ~ . + pgr2 + pgr3)

# Brier Score and IPA in the validation set (model without PGR)
score_gbsg5 <-
  Score(list("cox_validation" = efit1),
    formula = Surv(ryear, rfs) ~ 1, 
    data = gbsg5, 
    conf.int = TRUE, 
    times = 4.99,
    cens.model = "km", 
    metrics = "brier",
    summary = "ipa"
)

# Extra: bootstrap confidence intervals for IPA ------
B <- 100
horizon <- 4.99
boots_ls <- lapply(seq_len(B), function(b) {
  
  # Resample validation data
  data_boot <- gbsg5[sample(nrow(gbsg5), replace = TRUE), ]

  
  # Get IPA on boot validation data
  score_boot <- Score(
    list("cox_validation" = efit1),
    formula = Surv(ryear, rfs) ~ 1,
    cens.model = "km", 
    data = data_boot, 
    conf.int = FALSE, 
    times = horizon,
    metrics = c("brier"),
    summary = c("ipa")
  )
  
  #.. can add other measure heres, eg. concordance
  
  ipa_boot <- score_boot$Brier$score[model == "cox_validation"][["IPA"]]
  cbind.data.frame("ipa" = ipa_boot)
})

df_boots <- do.call(rbind.data.frame, boots_ls)
```

</details>

    ##                                Estimate Lower .95  Upper .95
    ## Brier - Validation data            0.22       0.21      0.24
    ## Scaled Brier - Validation data     0.10       0.04      0.16

Brier and scaled Brier score were 0.22 and 0.11, respectively.

## Goal 3 - Clinical utility

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
and treat none of the patients. If overtreatment is harmful, a higher
threshold should be used.  
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

<details>
<summary>
Click to expand code
</summary>

``` r
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(survival,
               Hmisc)

# Run decision curve analysis

# Development data
# Model without PGR
gbsg5 <- as.data.frame(gbsg5)
dca_gbsg5 <- stdca(
  data = gbsg5, 
  outcome = "rfs", 
  ttoutcome = "ryear",
  timepoint = 5, 
  predictors = "pred", 
  xstop = 1.0,
  ymin = -0.01, 
  graph = FALSE
)

# Decision curves plot
par(xaxs = "i", yaxs = "i", las = 1)
plot(dca_gbsg5$net.benefit$threshold,
  dca_gbsg5$net.benefit$pred,
  type = "l", 
  lwd = 2,
  lty = 1,
  xlab = "Threshold probability in %", 
  ylab = "Net Benefit",
  xlim = c(0, 1), 
  ylim = c(-0.10, 0.60), 
  bty = "n",
  cex.lab = 1.2, 
  cex.axis = 1
)

lines(dca_gbsg5$net.benefit$threshold, 
      dca_gbsg5$net.benefit$none, 
      type = "l", 
      lwd = 2, 
      lty = 4)
lines(dca_gbsg5$net.benefit$threshold, 
      dca_gbsg5$net.benefit$all, 
      type = "l", 
      lwd = 2, 
      col = "darkgray")
lines(dca_gbsg5$net.benefit$threshold, 
      dca_gbsg5$net.benefit$pred, 
      type = "l", lwd = 2, lty = 5)
legend("topright",
  c(
    "Treat All",
    "Original model",
    "Treat None"
  ),
  lty = c(1, 1, 4), 
  lwd = 2, 
  col = c("darkgray", "black", "black"),
  bty = "n"
)
```

</details>

    ## [1] "pred: No observations with risk greater than 84%, and therefore net benefit not calculable in this range."

<img src="imgs/01_predsurv_simplified/dca-1.png" width="672" style="display: block; margin: auto;" />

The potential benefit at 23% threshold of the prediction model is 0.36.
This means that the model might identify a net 36 patients out of 100
who will have recurrent breast cancer or die within 5 years of surgery
and thus require adjuvant chemotherapy.

Potential benefit can be also defined in terms of net reduction of
avoidable interventions (e.g adjuvant chemotherapy per 100 patients) by:

<img src="https://render.githubusercontent.com/render/math?math=%5Chuge%7B%5Cfrac%7BNB_%7Bmodel%7D%20-%20NB_%7Ball%7D%7D%7B(p_t%2F%20(1-p_t))%7D*100%7D%0A">

where *NB*<sub>model</sub> is the net benefit of the prediction model,
*NB*<sub>all</sub> is the net benefit of the strategy treat all and
*p*<sub>*t*</sub> is the risk threshold.

## Additional notes

1.  To run the apparent validation find in any performance measure
    calculation find “gbsg5” and replace with “rott5” except for model
    development part;

2.  To run validation of the extended model in any performance find
    “efit1” and replace with “efit1_pgr”.

## Reproducibility ticket

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19044)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] here_1.0.1                timeROC_0.4              
    ##  [3] riskRegression_2021.10.10 pec_2021.10.11           
    ##  [5] prodlim_2019.11.13        rms_6.2-0                
    ##  [7] SparseM_1.81              Hmisc_4.6-0              
    ##  [9] ggplot2_3.3.5             Formula_1.2-4            
    ## [11] lattice_0.20-45           survival_3.2-13          
    ## [13] pacman_0.5.1             
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-153         cmprsk_2.2-10        matrixStats_0.61.0  
    ##  [4] lubridate_1.8.0      RColorBrewer_1.1-2   rprojroot_2.0.2     
    ##  [7] numDeriv_2016.8-1.1  tools_4.1.2          backports_1.3.0     
    ## [10] utf8_1.2.2           R6_2.5.1             rpart_4.1-15        
    ## [13] DBI_1.1.1            colorspace_2.0-2     nnet_7.3-16         
    ## [16] withr_2.4.3          tidyselect_1.1.1     gridExtra_2.3       
    ## [19] compiler_4.1.2       quantreg_5.86        htmlTable_2.3.0     
    ## [22] sandwich_3.0-1       scales_1.1.1         checkmate_2.0.0     
    ## [25] polspline_1.1.19     mvtnorm_1.1-3        stringr_1.4.0       
    ## [28] digest_0.6.29        foreign_0.8-81       rmarkdown_2.11      
    ## [31] base64enc_0.1-3      jpeg_0.1-9           pkgconfig_2.0.3     
    ## [34] htmltools_0.5.2      parallelly_1.29.0    highr_0.9           
    ## [37] fastmap_1.1.0        htmlwidgets_1.5.4    rlang_0.4.12        
    ## [40] rstudioapi_0.13      generics_0.1.1       zoo_1.8-9           
    ## [43] dplyr_1.0.7          ModelMetrics_1.2.2.2 magrittr_2.0.1      
    ## [46] Matrix_1.3-4         Rcpp_1.0.7           munsell_0.5.0       
    ## [49] fansi_0.5.0          lifecycle_1.0.1      stringi_1.7.6       
    ## [52] multcomp_1.4-17      pROC_1.18.0          yaml_2.2.1          
    ## [55] MASS_7.3-54          plyr_1.8.6           recipes_0.1.17      
    ## [58] grid_4.1.2           parallel_4.1.2       listenv_0.8.0       
    ## [61] mets_1.2.9           crayon_1.4.2         splines_4.1.2       
    ## [64] timereg_2.0.1        knitr_1.36           pillar_1.6.4        
    ## [67] future.apply_1.8.1   reshape2_1.4.4       codetools_0.2-18    
    ## [70] stats4_4.1.2         glue_1.5.1           evaluate_0.14       
    ## [73] latticeExtra_0.6-29  data.table_1.14.2    png_0.1-7           
    ## [76] vctrs_0.3.8          foreach_1.5.1        MatrixModels_0.5-0  
    ## [79] gtable_0.3.0         purrr_0.3.4          future_1.23.0       
    ## [82] assertthat_0.2.1     xfun_0.28            gower_0.2.2         
    ## [85] class_7.3-19         timeDate_3043.102    tibble_3.1.6        
    ## [88] conquer_1.2.1        iterators_1.0.13     cluster_2.1.2       
    ## [91] lava_1.6.10          globals_0.14.0       TH.data_1.1-0       
    ## [94] ellipsis_0.3.2       caret_6.0-90         ipred_0.9-12
