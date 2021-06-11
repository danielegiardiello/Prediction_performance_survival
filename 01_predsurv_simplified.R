


# Libraries and options ----------------------------------

library(survival)
library(rms)
library(pec)
library(riskRegression)
library(timeROC)

options(show.signif.stars = FALSE)  # display statistical intelligence
palette("Okabe-Ito")  # color-blind friendly  (needs R 4.0)



# Data and recoding ----------------------------------
# Development data

rotterdam$ryear <- rotterdam$rtime/365.25  # time in years
rotterdam$rfs <- with(rotterdam, pmax(recur, death))

# variables used in the analysis
pgr99 <- quantile(rotterdam$pgr, .99) # there is a large outlier of 5000
rotterdam$pgr2 <- pmin(rotterdam$pgr, pgr99) # Winsorized value
rotterdam$csize <- rotterdam$size           # categorized size
rotterdam$cnode <- cut(rotterdam$nodes, 
                       c(-1,0, 3, 50),
                       c("0", "1-3", ">3"))   # categorized node

# Save in the data the restricted cubic spline term using Hmisc::rcspline.eval() package
rcs3_pgr <- rcspline.eval(rotterdam$pgr2, 
                          knots = c(0, 41, 486))
attr(rcs3_pgr, "dim") <- NULL
attr(rcs3_pgr, "knots") <- NULL
rotterdam$pgr3 <- rcs3_pgr

# Validation data
gbsg$ryear <- gbsg$rfstime/365.25
gbsg$rfs   <- gbsg$status           # the GBSG data contains RFS
gbsg$cnode <- cut(gbsg$nodes, 
                  c(-1,0, 3, 51),
                  c("0", "1-3", ">3"))   # categorized node
gbsg$csize <- cut(gbsg$size,  
                  c(-1, 20, 50, 500), #categorized size
                  c("<=20", "20-50", ">50"))
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



# Model development -------------------------------------

efit1 <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade,
               data = rott5, 
               x = T, 
               y = T)

# COMMENT: I will wait for brier function created by Terry
score_gbsg5 <-
  Score(list("Cox development" = efit1),
    formula = Surv(ryear, rfs) ~ 1, 
    data = gbsg5, 
    conf.int = TRUE, 
    times = 4.95,
    cens.model = "km", 
    metrics = "brier",
    summary = "ipa"
  )

score_gbsg5$Brier$score 

# Discrimination ---------------------------------------

# Add linear predictor in the validation set
gbsg5$lp <- predict(efit1, newdata = gbsg5)


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
  Uno_C_vdata$concordance + 
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
    marker = gbsg5$lp,
    cause = 1, 
    weighting = "marginal", 
    times = 4.95,
    iid = TRUE
  )

Uno_AUC_res <- c(
  "Uno AUC" = unname(Uno_gbsg5$AUC[2]),
  "2.5 %" = unname(Uno_gbsg5$AUC["t=4.95"] -
    qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=4.95"]),
  "97. 5 %" = unname(Uno_vdata1$AUC["t=4.95"] +
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

# Predicted risk 
gbsg5$pred <- 1 - predictSurvProb(efit1, 
                                  newdata = gbsg5,
                                  times = t_horizon)
# Expected
exp_t <- mean(gbsg5$pred)

OE_t <- obs_t / exp_t

alpha <- .05
OE_summary <- c(
  "OE" = OE_t,
  "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
  "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
)

OE_summary



# Calibration plot ----------------------------------
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


# Clinical utility --------------------------------

source("Functions/stdca.R")
gbsg5 <- as.data.frame(gbsg5)
dca_gsbg5 <- stdca(
  data = gbsg5, outcome = "rfs", ttoutcome = "ryear",
  timepoint = 5, predictors = "pred", xstop = 1.0,
  ymin = -0.01, graph = FALSE
)

# Decision curves plot
par(xaxs = "i", yaxs = "i", las = 1)
plot(dca_gsbg5$net.benefit$threshold,
     dca_gsbg5$net.benefit$pred,
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
lines(dca_gsbg5$net.benefit$threshold,
      dca_gsbg5$net.benefit$none, 
      type = "l", 
      lwd = 2, 
      lty = 4)
lines(dca_gsbg5$net.benefit$threshold,
      dca_gsbg5$net.benefit$all,
      type = "l", 
      lwd = 2, 
      col = "darkgray")


