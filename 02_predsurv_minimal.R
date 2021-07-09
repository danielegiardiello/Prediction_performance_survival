


# Libraries and options ----------------------------------

library(survival)
library(rms)
library(pec)
library(riskRegression)
library(timeROC)

options(show.signif.stars = FALSE)  # display statistical intelligence
palette("Okabe-Ito")  # color-blind friendly  (needs R 4.0)



# Data and recoding ----------------------------------
# Validation data
gbsg$ryear <- gbsg$rfstime/365.25
gbsg$rfs   <- gbsg$status           # the GBSG data contains RFS
gbsg$cnode <- cut(gbsg$nodes, 
                  c(-1,0, 3, 51),
                  c("0", "1-3", ">3"))   # categorized node
gbsg$csize <- cut(gbsg$size,  
                  c(-1, 20, 50, 500), #categorized size
                  c("<=20", "20-50", ">50"))
pgr99 <- 1347.85 
gbsg$pgr2 <- pmin(gbsg$pgr, pgr99) # Winsorized value


# Restricted cubic spline for PGR
rcs3_pgr <- rcspline.eval(gbsg$pgr2, 
                          knots = c(0, 41, 486))
attr(rcs3_pgr, "dim") <- NULL
attr(rcs3_pgr, "knots") <- NULL
gbsg$pgr3 <- rcs3_pgr


# The analysis will focus on the first 5 years: create
# data sets that are censored at 5 years
temp <- survSplit(Surv(ryear, rfs) ~ ., 
                  data = rotterdam, 
                  cut=5,
                  episode="epoch")
rott5 <- subset(temp, epoch==1)  # only the first 5 years
temp <- survSplit(Surv(ryear, rfs) ~ ., 
                  data = gbsg, 
                  cut=5,
                  episode ="epoch")
gbsg5 <- subset(temp, epoch==1)

# Relevel
rott5$cnode <- relevel(rotterdam$cnode, "1-3")
gbsg5$cnode <- relevel(gbsg$cnode, "1-3")


# Absolute risk calculation -------------
# Baseline survival at 5 years - basic model
S0.5yr <- .8604385 
# Design matrix of predictors
des_matr <- as.data.frame(model.matrix(~ csize + cnode + grade, data = gbsg5))
des_matr$`(Intercept)` <- NULL
# Coefficients
coef <- c(0.3922098, 0.6656456, -0.3538533, 0.6936283,  0.3766110)
# Prognostic index (PI)
gbsg5$PI <- as.vector(as.matrix(des_matr) %*% cbind(coef))
# Absolute risk at 5 years (1 - S(t), 1 - survival at time t)
gbsg5$pred5 <-  as.vector(1 - S0.5yr**exp(gbsg5$PI))


#Absolute risk  at 5 yrs - Extended model with PGR ---------------
# Baseline survival at 5 years - basic model
S0.5yr_pgr <- .8084657  
# Design matrix of predictors
des_matr <- as.data.frame(model.matrix(~ csize + cnode + grade + 
                                         I(pgr2) + I(pgr3), data = gbsg5))
des_matr$`(Intercept)` <- NULL
# Coefficients
coef <- c(0.37096136, 0.64380754, -0.37445717, 
          0.66988919, 0.32169223, -0.00293231, 0.01281538)
# Prognostic index (PI)
gbsg5$PI_pgr <- as.vector(as.matrix(des_matr) %*% cbind(coef))
# Absolute risk at 5 years (1 - S(t), 1 - survival at time t)
gbsg5$pred5_pgr <-  as.vector(1 - S0.5yr_pgr**exp(gbsg5$PI_pgr))

# Discrimination ---------------------------------------

## Validation data
# Harrell's C
harrell_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ PI, 
                               gbsg5, 
                               reverse = TRUE)
# Uno's C
Uno_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ PI_pgr, 
                           gbsg5, 
                           reverse = TRUE,
                           timewt = "n/G2")
alpha <- .05
res_C <- matrix(
  c(
    harrell_C_gbsg5$concordance,
    harrell_C_gbsg5$concordance - 
      qnorm(1 - alpha/2) * sqrt(harrell_C_gbsg5$var),
    harrell_C_gbsg5$concordance + 
      qnorm(1 - alpha/2) * sqrt(harrell_C_gbsg5$var),
    
    Uno_C_gbsg5$concordance,
    Uno_C_gbsg5$concordance - 
      qnorm(1 - alpha/2) * sqrt(Uno_C_gbsg5$var),
    Uno_C_gbsg5$concordance + 
      qnorm(1 - alpha/2) * sqrt(Uno_C_gbsg5$var)
  ), 
  nrow = 2,
  ncol = 3, 
  byrow = T,
  dimnames = list(c("Harrell C", "Uno C"),
                  c("Estimate", "2.5 %", "97.5 %"))
)

res_C

# Uno's time dependent AUC

Uno_gbsg5 <-
  timeROC(
    T = gbsg5$ryear, 
    delta = gbsg5$rfs,
    marker = gbsg5$PI,
    cause = 1, 
    weighting = "marginal", 
    times = 4.95,
    iid = TRUE
  )

Uno_AUC_res <- c(
  "Uno AUC" = unname(Uno_gbsg5$AUC[2]),
  "2.5 %" = unname(Uno_gbsg5$AUC["t=4.95"] -
                     qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.95"]),
  "97. 5 %" = unname(Uno_gbsg5$AUC["t=4.95"] +
                       qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.95"])
)

Uno_AUC_res



# Calibration -----------------------------------------

# Observed / Expected ration
t_horizon <- 5

# Observed
obj <- summary(survfit(
  Surv(ryear, rfs) ~ 1, 
  data = gbsg5),
  times = t_horizon)

obs_t <- 1 - obj$surv

# Expected
exp_t <- mean(gbsg5$pred5)

OE_t <- obs_t / exp_t

alpha <- .05
OE_summary <- c(
  "OE" = OE_t,
  "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
  "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
)

OE_summary



# Calibration plot ----------------------------------
gbsg5$pred.cll <- log(-log(1 - gbsg5$pred5))


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
  "pred" = gbsg5$pred5
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


# Overall performance ---------------------------------------
# COMMENT: I will wait for brier function created by Terry
# For now I programmed the Brier score function
brier_score <- function(tfup, status, thorizon, survival) {
  db <- data.frame(
    time = tfup,
    status = status
  )
  
  db$predsurv <- survival
  db$predmort <- 1 - db$predsurv
  
  db <- db %>% mutate(cat = case_when(
    time <= thorizon & status == 1 ~ 1,
    time > thorizon ~ 2,
    time <= thorizon & status == 0 ~ 3
  ))
  
  
  sfit <- db %>% survfit(Surv(time, status == 0) ~ 1, data = .)
  sfit_df <- data.frame(
    time = c(0, sfit$time), surv = c(1, sfit$surv),
    weight = c(1, 1 / sfit$surv)
  )
  
  db2 <- db %>%
    left_join(sfit_df, by = "time") %>%
    select(time, status, predsurv, predmort, cat, weight)
  
  db2$weight[db2$time > thorizon] <- max(db2$weight[db2$weight != max(db2$weight)])
  # summary(data2$weight)
  # tail(data2)
  
  db2$weight[db2$cat == 3] <- 0
  db2$contrib[db2$cat == 1] <- (-db$predsurv[db2$cat == 1])**2
  db2$contrib[db2$cat == 2] <- (1 - db$predsurv[db2$cat == 2])**2
  db2$contrib[db2$cat == 3] <- 0
  db2$bs <- db2$contrib * db2$weight
  
  brier <- (1 / sum(db2$weight)) * sum(db2$bs)
  
  # IPA
  # Null model
  cox_null <- db %>% coxph(Surv(time, status) ~ 1, data = ., x = T, y = T)
  
  db_null <- db2
  db_null$predsurv_null <- db_null %>% predictSurvProb(cox_null, times = thorizon, newdata = .)
  
  sfit_null <- db %>% survfit(Surv(time, status) ~ 1, data = .)
  sfit_null_df <- data.frame(
    time = c(0, sfit_null$time),
    surv = c(1, sfit_null$surv)
  )
  
  db_null2 <- db_null %>%
    left_join(sfit_null_df, by = "time") %>%
    select(time, status, predsurv_null, predmort, weight, cat, surv)
  
  db_null2$weight[db_null2$cat == 3] <- 0
  db_null2$contrib[db_null2$cat == 1] <-
    (-db_null2$predsurv_null[db_null2$cat == 1])**2
  db_null2$contrib[db_null2$cat == 2] <-
    (1 - db_null2$predsurv_null[db_null2$cat == 2])**2
  db_null2$contrib[db_null2$cat == 3] <- 0
  db_null2$bs <- db_null2$contrib * db_null2$weight
  
  brier_null <- (1 / sum(db_null2$weight)) * sum(db_null2$bs)
  IPA <- 1 - (brier / brier_null)
  
  res <- c(brier, brier_null, IPA)
  names(res) <- c("Brier", "Null Brier", "IPA")
  return(res)
}


brier_gbsg5 <-
  brier_score(
    tfup = gbsg5$ryear, status = gbsg5$rfs,
    thorizon = 4.95, survival = 1 - gbsg5$pred5
  )

res_ov <- c("Brier" = brier_gbsg5["Brier"],
            "Scaled Brier" = brier_gbsg5["IPA"])

res_ov



# Clinical utility --------------------------------

# Minimal version (better to use stdca function in the repository):
# source("Functions/stdca.R")

# 1. Set grid of thresholds
thresholds <- seq(0, 1.0, by = 0.01)

# 2. Calculate Aelen johansen for all patients exceeding threshold (i.e. treat-all)
horizon <- 5
survfit_all <- summary(
  survfit(Surv(ryear, rfs) ~ 1, data = gbsg5), 
  times = horizon
)
f_all <- 1 - survfit_all$surv

# 3. Calculate Net Benefit across all thresholds
list_nb <- lapply(thresholds, function(ps) {
  
  # Treat all
  NB_all <- f_all - (1 - f_all) * (ps / (1 - ps))
  
  # Based on threshold
  p_exceed <- mean(gbsg5$pred5 > ps)
  survfit_among_exceed <- try(
    summary(
      survfit(Surv(ryear, rfs) ~ 1, data = gbsg5[gbsg5$pred5 > ps, ]), 
      times = horizon
    ), silent = TRUE
  )
  
  # If a) no more observations above threshold, or b) among subset exceeding..
  # ..no indiv has event time >= horizon, then NB = 0
  if (class(survfit_among_exceed) == "try-error") {
    NB <- 0
  } else {
    f_given_exceed <- 1 - survfit_among_exceed$surv
    TP <- f_given_exceed * p_exceed
    FP <- (1 - f_given_exceed) * p_exceed
    NB <- TP - FP * (ps / (1 - ps))
  }
  
  # Return together
  df_res <- data.frame("threshold" = ps, "NB" = NB, "treat_all" = NB_all)
  return(df_res)
})


# Combine into data frame
df_nb <- do.call(rbind.data.frame, list_nb)
head(df_nb)

# Make basic decision curve plot
par(
  xaxs = "i", 
  yaxs = "i", 
  las = 1, 
  mar = c(6.1, 5.8, 4.1, 2.1), 
  mgp = c(4.25, 1, 0)
)
plot(
  df_nb$threshold, 
  df_nb$NB,
  type = "l", 
  lwd = 2,
  ylim = c(-0.1, 0.6),
  xlim = c(0, 1), 
  xlab = "",
  ylab = "Net Benefit",
  bty = "n", 
)
lines(df_nb$threshold, df_nb$treat_all, type = "l", col = "darkgray", lwd = 2)
abline(h = 0, lty = 2, lwd = 2)
legend(
  "topright", 
  c("Treat all", "Treat none", "Prediction model"),
  lwd = c(2, 2, 2), 
  lty = c(1, 2, 1), 
  col = c("darkgray", "black", "black"), 
  bty = "n"
)
mtext("Threshold probability", 1, line = 2)
title("Validation data")

# NOTES ---------------------------
# 1. To run the apparent validation find "gbsg5" with "rott5"
# from paragraph "Discrimination" on. 
# 2. To run the model with the PGR as additional biomarker
# find "efit1" with "efit1_pgr"
# from paragraph "Discrimination" on. 

# When use apparent validation, be careful to use the correct labels in the plot!



