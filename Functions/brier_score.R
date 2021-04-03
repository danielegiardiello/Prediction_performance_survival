
# Function to calculate Brier Score without development data
# Author: Daniele Giardiello

# Brier Score calculation
# Brier score function
# tfup = follow-up time in the validation data
# status = indicator variable (0=censored, 1=event) in the validation data
# thorizon = fixed time horizon t (in this example 4.95 years)
# survival = predicted survival at time t (this example 5 years)
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
