#' Calculate  R2 and Royston's D
#'
#' @param lp linear predictor (X*beta e.g. survival::predict())
#' @param time follow-up time in the validation data
#' @param status indicator variable (0=censored, 1=event) in the validation data
#' @param ties correction for ties (default is TRUE)
#'
#' @return
#'
#' @author Terry Therneau, Daniele Giardiello
#'
#' @examples


# R2
# Function to calculate R2 and also Royston's D
# although we suggest time-dependent C and AUC
Rsq <- function(lp, time, status, ties = TRUE) {
  pi <- 3.141592653589793 # a user might have overwritten the constant
  phat <- lp
  y2 <- Surv(time, status)
  n <- length(phat)

  if (ties && any(duplicated(phat))) {
    z <- qnorm((1:n - 3 / 8) / (n + .25))
    # per the paper, take the average z over any ties
    index1 <- match(phat, sort(unique(phat)))
    index2 <- rank(phat, ties = "first")
    z2 <- tapply(z[index2], index1, mean)
    qhat <- z2[index1]
  }
  else {
    qhat <- qnorm((rank(phat) - 3 / 8) / (n + .25))
  } # simple case of no ties

  rfit <- coxph(y2 ~ qhat)
  beta <- unname(coef(rfit))
  D <- beta * sqrt(8 / pi)
  R2 <- beta^2 / ((pi^2 / 6) + beta^2)

  c(D = D, R2 = R2) # return vector
}
