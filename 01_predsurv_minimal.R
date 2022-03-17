


# Libraries and options ----------------------------------

# General packages
pkgs <- c("survival", "pec", "rms", "timeROC", "riskRegression")
vapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE, quietly = TRUE)
}, FUN.VALUE = logical(length = 1L))

options(show.signif.stars = FALSE)  # display statistical intelligence
palette("Okabe-Ito")  # color-blind friendly  (needs R 4.0)


# Data and recoding ----------------------------------
# Development data

rotterdam$ryear <- rotterdam$rtime/365.25  # time in years
rotterdam$rfs <- with(rotterdam, pmax(recur, death)) #The variable rfs is a status indicator, 0 = alive without relapse, 1 = death or relapse.

# Fix the outcome for 43 patients who have died but 
# censored at time of recurrence which was less than death time. 
# The actual death time should be used rather than the earlier censored recurrence time.

rotterdam$ryear[rotterdam$rfs == 1 & 
                  rotterdam$recur == 0 & 
                  rotterdam$death == 1 & 
                  (rotterdam$rtime < rotterdam$dtime)] <- 
  
rotterdam$dtime[rotterdam$rfs == 1 &
                rotterdam$recur == 0 & 
                rotterdam$death == 1 & 
                (rotterdam$rtime < rotterdam$dtime)]/365.25  

# variables used in the analysis
pgr99 <- quantile(rotterdam$pgr, .99, type = 1) # there is a large outlier of 5000, used type=1 to get same result as in SAS
rotterdam$pgr2 <- pmin(rotterdam$pgr, pgr99) # Winsorized value
rotterdam$csize <- rotterdam$size           # categorized size
rotterdam$cnode <- cut(rotterdam$nodes, 
                       c(-1,0, 3, 51),
                       c("0", "1-3", ">3"))   # categorized node
rotterdam$grade3 <- as.factor(rotterdam$grade)
levels(rotterdam$grade3) <- c("1-2", "3")

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
gbsg$grade3 <- as.factor(gbsg$grade)
levels(gbsg$grade3) <- c("1-2", "1-2", "3")

# Restricted cubic spline for PGR
rcs3_pgr <- rcspline.eval(gbsg$pgr2, knots = c(0, 41, 486))
attr(rcs3_pgr, "dim") <- NULL
attr(rcs3_pgr, "knots") <- NULL
gbsg$pgr3 <- rcs3_pgr


# Much of the analysis will focus on the first 5 years: create
#  data sets that are censored at 5
temp <- survSplit(Surv(ryear, rfs) ~ ., data = rotterdam, cut = 5,
                  episode="epoch")
rott5 <- subset(temp, epoch == 1)  # only the first 5 years
temp <- survSplit(Surv(ryear, rfs) ~ ., data = gbsg, cut = 5,
                  episode ="epoch")
gbsg5 <- subset(temp, epoch == 1)

# Relevel
rott5$cnode <- relevel(rotterdam$cnode, "0")
gbsg5$cnode <- relevel(gbsg$cnode, "0")


# Model development -------------------------------------

efit1 <- coxph(Surv(ryear, rfs) ~ csize + cnode + grade3,
               data = rott5, 
               x = T, 
               y = T)

# The model with additional PGR marker
efit1_pgr  <- update(efit1, . ~ . + pgr2 + pgr3)


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
    marker = gbsg5$lp,
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

# calibration slope (fixed time point)-------------------------------------
gval <- coxph(Surv(ryear, rfs) ~ lp, data=gbsg5)

calslope_summary <- c(
  "calibration slope" = gval$coef,
  "2.5 %"  = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
  "97.5 %" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
)

calslope_summary



# Overall performance ---------------------------------------
score_gbsg5 <-
  Score(list("cox" = efit1),
        formula = Surv(ryear, rfs) ~ 1, 
        data = gbsg5, 
        conf.int = TRUE, 
        times = 4.99,
        cens.model = "km", 
        metrics = "brier",
        summary = "ipa"
  )

score_gbsg5$Brier$score 


# Clinical utility --------------------------------

# 1. Set grid of thresholds
thresholds <- seq(0, 1.0, by = 0.01)

# 2. Calculate observed risk for all patients exceeding threshold (i.e. treat-all)
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
  p_exceed <- mean(gbsg5$pred > ps)
  survfit_among_exceed <- try(
    summary(
      survfit(Surv(ryear, rfs) ~ 1, data = gbsg5[gbsg5$pred > ps, ]), 
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

#read off at 23% threshold
df_nb[df_nb$threshold == 0.23,]

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
# 1. To run the apparent validation find "gbsg5" and replace with "rott5"
# from paragraph "Discrimination" on. 
# 2. To run the model with the PGR as additional biomarker
# find "efit1" and replace with with "efit1_pgr"
# from paragraph "Discrimination" on. 

# When use apparent validation, be careful to use the correct labels in the plot!



