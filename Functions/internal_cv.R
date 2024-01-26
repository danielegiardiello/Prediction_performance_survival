#' Calculate optimism-corrected bootstrap internal validation for AUC, Brier and Scaled Brier Score
#' @param db data to calculate the optimism-corrected bootstrap
#' @param B number of bootstrap sample (default 10)
#' @param time follow-up time
#' @param status indicator variable (0=censored, 1=event)
#' @param formula_model formula for the model (Cox model)
#' @param pred.time time horizon as predictor
#' @param formula_ipcw formula to calculate inverse probability censoring weighting
#' @seed set the seed for bootstrap computation (default 2024)
#'
#' @return
#'
#' @author Daniele Giardiello
#'
#' @examples

# Use pacman to check whether packages are installed, if not load
if (!require("pacman")) install.packages("pacman")
library(pacman)

pacman::p_load(
  rio,
  survival,
  rms,
  pec,
  tidyverse,
  timeROC,
  riskRegression
)


bootstrap_cv <- function(db, B = 10,
                         time,
                         status,
                         formula_model,
                         pred.time,
                         formula_ipcw,
                         seed = 2024) {
  frm_model <- as.formula(formula_model)
  frm_ipcw <- as.formula(formula_ipcw)
  db$id <- 1:nrow(db)


  # Duplicate data
  db_ext <- db %>% dplyr::slice(rep(dplyr::row_number(), B))
  db_ext$.rep <- with(db_ext, ave(seq_along(id), id, FUN = seq_along)) # add an index identifying the replications

  db_tbl <- db_ext %>%
    dplyr::group_by(.rep) %>%
    tidyr::nest() %>%
    dplyr::rename(
      orig_data = data,
      id_boot = .rep
    )

  # Create bootstrap samples
  sample_boot <- function(db, B) {
    set.seed(seed)
    db_boot <- matrix(NA, nrow = nrow(db) * B, ncol = ncol(db))
    sample_row <- list()
    for (j in 1:B) {
      sample_row[[j]] <- sample(nrow(db), size = nrow(db), replace = TRUE)
    }
    sample_row <- unlist(sample_row)
    db_boot <- db[sample_row, ]
    db_boot$id_boot <- sort(rep(1:B, nrow(db)))
    db_boot <- db_boot %>%
      dplyr::group_by(id_boot) %>%
      tidyr::nest() %>%
      dplyr::rename(boot_data = data)
    return(db_boot)
  }

  # Join original data and the bootstrap data in a nested tibble
  a <- sample_boot(db, B)
  b <- a %>% left_join(db_tbl)

  # Create optimism-corrected performance measures
  b <- b %>% 
    dplyr::mutate(
      cox_boot = purrr::map(
      boot_data,
      ~ coxph(frm_model, data = ., x = T, y = T)
    ),
    
    cox_apparent = purrr::map(
      orig_data,
      ~ coxph(frm_model, data = ., x = T, y = T)
    ),
    
    # Discrimination assessment ------ 
    # Discrimination - time range (Harrell & Uno c-statistic)
    
    Harrell_C_app =
      purrr::map2_dbl(
      orig_data, cox_apparent,
      ~ survival::concordance(Surv(.x[[time]], .x[[status]])
                              ~ predict(.y, newdata = .x),
                              reverse = TRUE)$concordance
    ),
    
    Harrell_C_orig =
      purrr::map2_dbl(
        orig_data, cox_boot,
        ~ survival::concordance(Surv(.x[[time]], .x[[status]]) 
                                ~ predict(.y, newdata = .x),
                                reverse = TRUE)$concordance
      ),
    
    Harrell_C_boot =
      purrr::map2_dbl(
        boot_data, cox_boot,
        ~ survival::concordance(Surv(.x[[time]],.x[[status]])
                                ~ predict(.y, newdata = .x),
                                reverse = TRUE)$concordance
      ),
    
    Harrell_C_diff = 
      purrr::map2_dbl(
        Harrell_C_boot, Harrell_C_orig,
        function(a, b) {
          a - b
      }
    ),
    
    Uno_C_app =
      purrr::map2_dbl(
        orig_data, cox_apparent,
        ~ survival::concordance(Surv(.x[[time]], .x[[status]])
                                ~ predict(.y, newdata = .x),
                                reverse = TRUE,
                                timewt = "n/G2")$concordance
      ),
    
    Uno_C_orig =
      purrr::map2_dbl(
        orig_data, cox_boot,
        ~ survival::concordance(Surv(.x[[time]], .x[[status]])
                                ~ predict(.y, newdata = .x),
                                reverse = TRUE,
                                timewt = "n/G2")$concordance
      ),
    
    Uno_C_boot =
      purrr::map2_dbl(
        boot_data, cox_boot,
        ~ survival::concordance(Surv(.x[[time]],.x[[status]])
                                ~ predict(.y, newdata = .x),
                                reverse = TRUE,
                                timewt = "n/G2")$concordance
      ),
    
    Uno_C_diff = 
      purrr::map2_dbl(
        Uno_C_boot, Uno_C_orig,
        function(a, b) {
          a - b
      }
    ),
    
    # Discrimination - fixed time horizon (Uno's AUC)
    
    AUC_app = map2_dbl(
      orig_data, cox_apparent,
      ~ timeROC::timeROC(
        T = .x[[time]], delta = .x[[status]],
        marker = predict(.y, newdata = .x),
        cause = 1, weighting = "marginal", times = pred.time,
        iid = FALSE
      )$AUC[[2]]
    ),
    
    AUC_orig = map2_dbl(
      orig_data, cox_boot,
      ~ timeROC::timeROC(
        T = .x[[time]], delta = .x[[status]],
        marker = predict(.y, newdata = .x),
        cause = 1, weighting = "marginal", times = pred.time,
        iid = FALSE
      )$AUC[[2]]
    ),
    
    AUC_boot = purrr::map2_dbl(
      boot_data, cox_boot,
      ~ timeROC::timeROC(
        T = .x[[time]], delta = .x[[status]],
        marker = predict(.y, newdata = .x),
        cause = 1, weighting = "marginal", times = pred.time,
        iid = FALSE
      )$AUC[[2]]
    ),
    
    AUC_diff = purrr::map2_dbl(
      AUC_boot, AUC_orig,
      function(a, b) {
        a - b
      }
    ),
    
    # Overall performance - Brier Score & Scaled Brier 
    Score_app = purrr::map2(
      orig_data, cox_apparent,
      ~ riskRegression::Score(list("Cox" = .y),
                              formula = frm_ipcw,
                              data = .x, times = pred.time,
                              cens.model = "km",
                              metrics = "brier",
                              summary = "ipa"
                              )$Brier[[1]]
    ),
    
    Brier_app = purrr::map_dbl(Score_app, ~ .x$Brier[[2]]),
    IPA_app = purrr::map_dbl(Score_app, ~ .x$IPA[[2]]),
    Score_orig = purrr::map2(
      orig_data, cox_boot,
      ~ riskRegression::Score(list("Cox" = .y),
                              formula = frm_ipcw,
                              data = .x, times = pred.time,
                              cens.model = "km",
                              metrics = "brier",
                              summary = "ipa"
                              )$Brier[[1]]
    ),
    
    Brier_orig = purrr::map_dbl(Score_orig, ~ .x$Brier[[2]]),
    IPA_orig = purrr::map_dbl(Score_orig, ~ .x$IPA[[2]]),
    Score_boot = purrr::map2(
      boot_data, cox_boot,
      ~ riskRegression::Score(list("Cox" = .y),
                              formula = frm_ipcw,
                              data = .x, times = pred.time,
                              cens.model = "km",
                              metrics = "brier",
                              summary = "ipa"
                              )$Brier[[1]]
    ),
    
    Brier_boot = purrr::map_dbl(Score_boot, ~ .x$Brier[[2]]),
    IPA_boot = purrr::map_dbl(Score_boot, ~ .x$IPA[[2]]),
    Brier_diff = purrr::map2_dbl(
      Brier_boot, Brier_orig,
      function(a, b) {
        a - b
      }
    ),
    
    IPA_diff = purrr::map2_dbl(
      IPA_boot, IPA_orig,
      function(a, b) {
        a - b
      }
    ),
    
    
  # Calibration assessment --------
  
  obs_apparent = purrr::map_dbl(orig_data, ~ 
                                  summary(survival::survfit(Surv(.x[[time]], .x[[status]]) ~ 1,
                                                            data = .), time = pred.time)$surv), 
  # obs_apparent and obs_orig should be the same in this case
  obs_boot = purrr::map_dbl(boot_data, ~ 
                              summary(survival::survfit(Surv(.x[[time]], .x[[status]]) ~ 1,
                                                        data = .), time = pred.time)$surv),
  
  pred_apparent = purrr::map2(orig_data, 
                              cox_apparent, 
                              ~ riskRegression::predictRisk(.y,
                                                            newdata = .x,
                                                            time = pred.time)),
  
  pred_orig = purrr::map2(orig_data, 
                          cox_boot, 
                          ~ riskRegression::predictRisk(.y,
                                                        newdata = .x,
                                                        time = pred.time)),
  
  pred_boot = purrr::map2(boot_data, 
                          cox_boot, 
                          ~ riskRegression::predictRisk(.y,
                                                        newdata = .x,
                                                        time = pred.time)),
  
  OE_apparent = purrr::map2_dbl(obs_apparent, 
                                pred_apparent, 
                                ~ (1 - .x) / mean(.y)),
  
  OE_orig = purrr::map2_dbl(obs_apparent, 
                            pred_orig, 
                            ~ (1 - .x) / mean(.y)),
  
  OE_boot = purrr::map2_dbl(obs_boot, 
                            pred_boot, 
                            ~ (1 - .x) / mean(.y)),
  
  OE_diff = purrr::map2_dbl(OE_boot, OE_orig, 
                            function(a, b) {
                              a - b
                            }
  ),
  
  
  # mean calibration - time range assessment
  
  p_apparent = purrr::map2(orig_data, 
                           cox_apparent,
                           ~ predict(.y,
                                     newdata = .x,
                                     type = "expected")),
  
  p_orig = purrr::map2(orig_data, 
                       cox_boot,
                       ~ predict(.y,
                                 newdata = .x,
                                 type = "expected")),
  
  p_boot = purrr::map2(boot_data, 
                       cox_boot,
                       ~ predict(.y,
                                 newdata = .x,
                                 type = "expected")),
  
  intercept_app = purrr::map2_dbl(orig_data, 
                                  p_apparent,
                                  ~ exp(glm(.x[[status]] ~ offset(log(.y)), 
                                            family = poisson,
                                            data = .x,
                                            subset = (.y > 0))$coefficients)),
  
  
  intercept_orig = purrr::map2_dbl(orig_data, 
                                   p_orig,
                                   ~ exp(glm(.x[[status]] ~ offset(log(.y)), 
                                             family = poisson,
                                             data = .x,
                                             subset = (.y > 0))$coefficients)),
  
  intercept_boot = purrr::map2_dbl(boot_data, 
                                   p_boot,
                                   ~ exp(glm(.x[[status]] ~ offset(log(.y)), 
                                             family = poisson,
                                             data = .x,
                                             subset = (.y > 0))$coefficients)),
  
  
  intercept_diff = purrr::map2_dbl(intercept_boot, 
                                   intercept_orig,
                                   function(a, b) {
                                     a - b
                                   }),
  
  
  # weak calibration - fixed t
  
  lp_apparent = purrr::map2(orig_data, 
                            cox_apparent,
                            ~ predict(.y, newdata = .x, type = "lp")),
  
  lp_orig = purrr::map2(orig_data, 
                        cox_boot,
                        ~ predict(.y, newdata = .x, type = "lp")),
  
  lp_boot = purrr::map2(boot_data, 
                        cox_boot,
                        ~ predict(.y, newdata = .x, type = "lp")),
  
  
  slope_apparent = purrr::map2_dbl(orig_data, 
                                   lp_apparent,
                                   ~ coxph(Surv(.x[[time]], .x[[status]]) ~ .y)$coefficients), 
  
  slope_orig =purrr:: map2_dbl(orig_data, 
                               lp_orig,
                               ~ coxph(Surv(.x[[time]], .x[[status]]) ~ .y)$coefficients),
  
  slope_boot = purrr::map2_dbl(boot_data, 
                               lp_boot,
                               ~ coxph(Surv(.x[[time]], .x[[status]]) ~ .y)$coefficients),
  
  slope_diff = purrr::map2_dbl(slope_boot, 
                               slope_orig,
                               function(a, b) {
                                 a - b
                               }),
  
  # weak calibration - time range
  slope_range_apparent = purrr::pmap_dbl(list(lp = lp_apparent,
                                              p = p_apparent,
                                              db = orig_data),
                                         function(lp, p, db){
                                           slope = glm(db[[status]] ~ lp + offset(log(p) - lp), 
                                                       data = as.data.frame(db),
                                                       family = poisson,
                                                       subset = (p > 0))$coefficients[2]
                                         }),
  
  slope_range_orig = purrr::pmap_dbl(list(lp = lp_orig,
                                          p = p_orig,
                                          db = orig_data),
                                     function(lp, p, db){
                                       slope = glm(db[[status]] ~  lp + offset(log(p) - lp), 
                                                   data = as.data.frame(db),
                                                   family = poisson,
                                                   subset = (p > 0))$coefficients[2]
                                     }),
  
  slope_range_boot = purrr::pmap_dbl(list(lp = lp_boot, 
                                          p = p_boot, 
                                          db = boot_data),
                                     function(lp, p, db){
                                       slope = glm(db[[status]] ~  lp + offset(log(p) - lp), 
                                                   data = as.data.frame(db),
                                                   family = poisson,
                                                   subset = (p > 0))$coefficients[2]
                                     }),
  
  slope_range_diff =purrr::map2_dbl(slope_range_boot, slope_range_orig,
                                    function(a, b){
                                      a - b
                                    })
  
  
  )
  
  # Generate output
  # Discrimation & overall performance (Brier) ---
  AUC_corrected <- b$AUC_app[1] - mean(b$AUC_diff)
  Brier_corrected <- b$Brier_app[1] - mean(b$Brier_diff)
  IPA_corrected <- b$IPA_app[1] - mean(b$IPA_diff)
  Harrell_C_corrected <- b$Harrell_C_app[1] - mean(b$Harrell_C_diff)
  
  # Calibration ---
  OE_corrected <- b$OE_apparent[1] - mean(b$OE_diff)
  intercept_corrected <- b$intercept_app[1] - mean(b$intercept_diff)
  slope_corrected <- b$slope_apparent[1] - mean(b$slope_diff)
  slope_range_corrected <- b$slope_range_apparent[1] - mean(b$slope_range_diff)
  
  Uno_C_corrected <- b$Uno_C_app[1] - mean(b$Uno_C_diff)
  res <- c("AUC corrected" = AUC_corrected, 
           "Brier corrected" = Brier_corrected, 
           "IPA corrected" = IPA_corrected, 
           "Harrell C corrected" = Harrell_C_corrected,
           "Uno C corrected" = Uno_C_corrected,
           
           "Mean calibration fixed time point - OE" = OE_corrected,
           "Mean calibration corrected - time range" = intercept_corrected,
           
           "Weak calibration corrected (slope) - fixed time point" = slope_corrected,
           "Weak calibration corrected (slope) - time range" = slope_range_corrected
           )
  return(cbind(res))
}
