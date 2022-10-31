# Libraries and options ----------------------------------
# General packages
pkgs <- c("survival", "rms", 
          "timeROC", "riskRegression")
vapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE, quietly = TRUE)
}, FUN.VALUE = logical(length = 1L))

options(show.signif.stars = FALSE)  # display statistical intelligence
palette("Okabe-Ito")  # optionally use color-blind friendly color palette  (needs R 4.0)

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

# Absolute risk calculation -------------
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


### Validation of the original model ------------------

# Discrimination ---------------------------------------
## Validation data
# Harrell's C
harrell_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ PI, 
                               gbsg5, 
                               reverse = TRUE)
# Uno's C
Uno_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ PI, 
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
  timeROC::timeROC(
    T = gbsg5$ryear, 
    delta = gbsg5$rfs,
    marker = gbsg5$PI,
    cause = 1, 
    weighting = "marginal", 
    times = 4.99,
    iid = TRUE
  )

Uno_AUC_res <- c(
  "Uno AUC" = unname(Uno_gbsg5$AUC[2]),
  "2.5 %" = unname(Uno_gbsg5$AUC["t=4.99"] -
                     qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.99"]),
  "97. 5 %" = unname(Uno_gbsg5$AUC["t=4.99"] +
                       qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.99"])
)

Uno_AUC_res


# Calibration -----------------------------------------

# Observed / Expected ratio
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
vcal <- rms::cph(Surv(ryear, rfs) ~ rcs(pred.cll, 3),
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
  "pred" = gbsg5$pred5,
  
  "cll" = gbsg5$pred.cll,
  
  "id" = 1:nrow(gbsg5)
)

dat_cal <- dat_cal[order(dat_cal$pred), ]

dev.new() 
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

# Numerical measures
absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)

numsum_cph <- c(
  "ICI" = mean(absdiff_cph),
  setNames(quantile(absdiff_cph, c(0.5, 0.9)), c("E50", "E90")),
  "Emax" = max(absdiff_cph)
)
numsum_cph

# calibratoin slope (fixed time point)-------------------------------------

gval <- coxph(Surv(ryear, rfs) ~ PI, data=gbsg5)

calslope_summary <- c(
  "calibration slope" = gval$coef,
  "2.5 %"  = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
  "97.5 %" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
)

calslope_summary


# Overall performance ---------------------------------------
score_gbsg5 <-
  riskRegression::Score(list(gbsg5$pred5),
                        formula = Surv(ryear, rfs) ~ 1, 
                        data = gbsg5, 
                        conf.int = TRUE, 
                        times = 4.99,
                        cens.model = "km", 
                        metrics = "brier",
                        summary = "ipa"
)

score_gbsg5$Brier


# Clinical utility --------------------------------

# 1. Set grid of thresholds
thresholds <- seq(0, 1.0, by = 0.01)

# 2. Calculate observed risk  for all patients exceeding threshold (i.e. treat-all)
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

#read off at 23% threshold
df_nb[df_nb$threshold == 0.23,]

# Make basic decision curve plot

# Decision curves plot

# Smoothed decision curve
smooth_nb <- smooth(df_nb$NB, twiceit = TRUE)

dev.new() 
par(xaxs = "i", yaxs = "i", las = 1)
plot(df_nb$threshold,
     smooth_nb,
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
lines(df_nb$threshold, 
      df_nb$treat_all, 
      type = "l", 
      lwd = 3, 
      col = 2)
abline(h = 0, type = "l", lwd = 3, lty = 4, col = 8)
legend("topright",
       c(
         "Treat All",
         "Original model",
         "Treat None"
       ),
       lty = c(1, 2, 4), lwd = 3, 
       col = c(2, 4, 8),
       bty = "n"
)
title("Validation set", cex = 1.5)


### Validation of the extended model including PGR ------------------

# Discrimination ---------------------------------------
## Validation data
# Harrell's C
harrell_C_gbsg5 <- concordance(Surv(ryear, rfs) ~ PI_pgr, 
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
  timeROC::timeROC(
    T = gbsg5$ryear, 
    delta = gbsg5$rfs,
    marker = gbsg5$PI_pgr,
    cause = 1, 
    weighting = "marginal", 
    times = 4.99,
    iid = TRUE
  )

Uno_AUC_res <- c(
  "Uno AUC" = unname(Uno_gbsg5$AUC[2]),
  "2.5 %" = unname(Uno_gbsg5$AUC["t=4.99"] -
                     qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.99"]),
  "97. 5 %" = unname(Uno_gbsg5$AUC["t=4.99"] +
                       qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.99"])
)

Uno_AUC_res


# Calibration -----------------------------------------

# Observed / Expected ratio
t_horizon <- 5

# Observed
obj <- summary(survfit(
  Surv(ryear, rfs) ~ 1, 
  data = gbsg5),
  times = t_horizon)

obs_t <- 1 - obj$surv

# Expected
exp_t <- mean(gbsg5$pred5_pgr)

OE_t <- obs_t / exp_t

alpha <- .05
OE_summary <- c(
  "OE" = OE_t,
  "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
  "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
)

OE_summary


# Calibration plot ----------------------------------
gbsg5$pred.cll <- log(-log(1 - gbsg5$pred5_pgr))


# Estimate actual risk
vcal <- rms::cph(Surv(ryear, rfs) ~ rcs(pred.cll, 3),
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
  "pred" = gbsg5$pred5_pgr,
  
  "cll" = gbsg5$pred.cll,
  
  "id" = 1:nrow(gbsg5)
)

dat_cal <- dat_cal[order(dat_cal$pred), ]

dev.new() 
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

# Numerical measures
absdiff_cph <- abs(dat_cal$pred - dat_cal$obs)

numsum_cph <- c(
  "ICI" = mean(absdiff_cph),
  setNames(quantile(absdiff_cph, c(0.5, 0.9)), c("E50", "E90")),
  "Emax" = max(absdiff_cph)
)
numsum_cph

# calibratoin slope (fixed time point)-------------------------------------

gval <- coxph(Surv(ryear, rfs) ~ PI_pgr, data=gbsg5)

calslope_summary <- c(
  "calibration slope" = gval$coef,
  "2.5 %"  = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
  "97.5 %" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
)

calslope_summary


# Overall performance ---------------------------------------
score_gbsg5 <-
  riskRegression::Score(list(gbsg5$pred5_pgr),
                        formula = Surv(ryear, rfs) ~ 1, 
                        data = gbsg5, 
                        conf.int = TRUE, 
                        times = 4.99,
                        cens.model = "km", 
                        metrics = "brier",
                        summary = "ipa"
  )

score_gbsg5$Brier


# Clinical utility --------------------------------

# 1. Set grid of thresholds
thresholds <- seq(0, 1.0, by = 0.01)

# 2. Calculate observed risk  for all patients exceeding threshold (i.e. treat-all)
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
  p_exceed <- mean(gbsg5$pred5_pgr > ps)
  survfit_among_exceed <- try(
    summary(
      survfit(Surv(ryear, rfs) ~ 1, data = gbsg5[gbsg5$pred5_pgr > ps, ]), 
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

#read off at 23% threshold
df_nb[df_nb$threshold == 0.23,]

# Make basic decision curve plot

# Decision curves plot

# Smoothed decision curve
smooth_nb <- smooth(df_nb$NB, twiceit = TRUE)

dev.new() 
par(xaxs = "i", yaxs = "i", las = 1)
plot(df_nb$threshold,
     smooth_nb,
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
lines(df_nb$threshold, 
      df_nb$treat_all, 
      type = "l", 
      lwd = 3, 
      col = 2)
abline(h = 0, type = "l", lwd = 3, lty = 4, col = 8)
legend("topright",
       c(
         "Treat All",
         "Original model + PGR",
         "Treat None"
       ),
       lty = c(1, 2, 4), lwd = 3, 
       col = c(2, 4, 8),
       bty = "n"
)
title("Validation set", cex = 1.5)

##### ---------



