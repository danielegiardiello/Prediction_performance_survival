# Install libraries
# using python please use, for example: 
# pip install pandas 

# In R / R studio
# pkgs <- c("reticulate")
# vapply(pkgs, function(pkg) {
#   if (!require(pkg, character.only = TRUE)) install.packages(pkg)
#   require(pkg, character.only = TRUE, quietly = TRUE)
# }, FUN.VALUE = logical(length = 1L))

# py_install("pandas","numpy", "scipy", "statsmodels", 
#            "matplotlib", "sklearn", "lifelines", "warnings")


# Load libraries and data
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)
warnings.filterwarnings("ignore", category = RuntimeWarning) # suppressing warnings
import pandas as pd
import numpy as np
import scipy as sp
import statsmodels.api as smf
import matplotlib.pyplot as plt
import sklearn as sk
import lifelines as lf

# Get work directory
# os.getcwd()
url_rdata = "https://raw.githubusercontent.com/danielegiardiello/Prediction_performance_survival/main/Data/rotterdam.csv"
url_vdata = "https://raw.githubusercontent.com/danielegiardiello/Prediction_performance_survival/main/Data/gbsg.csv"
# NOTE: go to 
# "https://github.com/danielegiardiello/ValLogRegMod/blob/main/Data/vdata.csv"
# then click" Raw" button to the upper right corner of the file preview.
# Copy and paste the url link to have the raw gitHub version of the data
rotterdam = pd.read_csv(url_rdata)
gbsg = pd.read_csv(url_vdata)
# Inspect data:
# print(rotterdam.head(5)) # print the first five rows
# print(gbsg.head(5)) # print the first five rows
# rdata.info() # inspect data as in R str()
# vdata.info() # inspect data as in R str()

## Data manipulation ----
# Development data --
rotterdam["ryear"] = rotterdam.rtime / 365.25
rotterdam["rfs"] = np.maximum(rotterdam.recur, rotterdam.death)

# Prolong the follow-up among those death 
# to the last follow-up info about death
def death_time(df):
  if (df["rfs"] == 1) & (df["recur"] == 0) & (df["death"] == 1) & (df["rtime"] < df["dtime"]):
    return df["dtime"] / 365.25
  else:
    return df["ryear"]
  
rotterdam['ryear'] = rotterdam.apply(death_time, axis = 1)
      
# Converting categorical variables to dummies
# rotterdam.info()
# Rename size and grade
rotterdam = rotterdam.rename(columns = {"size" : "csize",
                                        "grade" : "cgrade"})
rotterdam = pd.get_dummies(data = rotterdam,
                           columns = ["csize", "cgrade"])
rotterdam = rotterdam.rename(columns = {"cgrade_2" : "cgrade_1-2"})
                           
from scipy.stats.mstats import winsorize
rotterdam["pgr2"] = winsorize(rotterdam["pgr"], limits = .01)
rotterdam["nodes2"] = winsorize(rotterdam["nodes"], limits = .01)

# Create function to calculate restricted cubic splines with three knots
def rcs_3(x, data):
    res_x = np.zeros((len(data), 1))
    qknots = [.1, .5, .9]
    knots = np.round(np.quantile(x, q = qknots), 0)
    res_x[:, 0] = (np.power(np.clip((x - knots[0]), a_min = 0, a_max = None), 3) - np.power(np.clip((x - knots[1]), a_min = 0, a_max = None), 3) *((knots[2] - knots[0])/(knots[2] - knots[1])) + np.power(np.clip((x - knots[2]), a_min = 0, a_max = None), 3) * ((knots[1] - knots[0])/(knots[2] - knots[1]))) / ((knots[2] - knots[0])**2)
    return(res_x)
# NOTE: to be extended for 4,5, 6 and 7 knots

# Create variable for knots for PGR
rotterdam["pgr3"] = rcs_3(rotterdam.pgr2, rotterdam)
rotterdam["nodes3"] = rcs_3(rotterdam.nodes2, rotterdam)


# Validation data ---
gbsg["ryear"] = gbsg["rfstime"] / 365.25
gbsg["rfs"]  = gbsg["status"]          # the GBSG data contains RFS

# Create categorical predictors of size and grade
csize = pd.cut(gbsg["size"], 
               bins = [-1, 20, 50, 5000],
               labels =  ["<=20", "20-50", ">50"])
csize = csize.rename("csize")       
               
cgrade = pd.cut(gbsg["grade"],
                bins = [0, 2, 5],
                labels = ["1-2", "3"])
cgrade = cgrade.rename("cgrade")
               
gbsg = pd.concat([gbsg, csize, cgrade], axis = 1)

gbsg = pd.get_dummies(data = gbsg,
                     columns = ["csize", "cgrade"])
                     
# Winsorization PGR in the validation data based on values of development data
gbsg["pgr2"] = np.minimum(gbsg["pgr"], 1360)
gbsg["nodes2"] = np.minimum(gbsg["nodes"], 19)

# Spline PGR in the validation data based on knots of the development data
def rcs_3_eval(x, data, knots):
    res_x = np.zeros((len(data), 1))
    res_x[:, 0] = (np.power(np.clip((x - knots[0]), a_min = 0, a_max = None), 3) - np.power(np.clip((x - knots[1]), a_min = 0, a_max = None), 3) *((knots[2] - knots[0])/(knots[2] - knots[1])) + np.power(np.clip((x - knots[2]), a_min = 0, a_max = None), 3) * ((knots[1] - knots[0])/(knots[2] - knots[1]))) / ((knots[2] - knots[0])**2)
    return(res_x)
 
gbsg["pgr3"] = rcs_3_eval(gbsg["pgr2"], data = gbsg, knots = [0, 41, 486])
gbsg["nodes3"] = rcs_3_eval(gbsg["nodes2"], data = gbsg, knots = [0, 1, 9])


# Sanity check
# rotterdam.info()
# gbsg.info() 

# Much of the analysis will focus on the first 5 years: create
#  data sets that are censored at 5
# ryear, rfs

# Development data
rott5 = rotterdam
rott5["rfs"] = np.where( (rott5["rfs"] == 1) & (rott5["ryear"] > 5), 0, rott5["rfs"])   
rott5["ryear"] = np.where(rott5["ryear"] > 5, 5, rott5["ryear"])

# Validation data 
gbsg5 = gbsg
gbsg5["rfs"] = np.where( (gbsg5["rfs"] == 1) & (gbsg5["ryear"] > 5), 0, gbsg5["rfs"])   
gbsg5["ryear"] = np.where(gbsg5["ryear"] > 5, 5, gbsg5["ryear"])


# Select only variable needed
# Development data 
rott5_sel = rott5[["pid", "rfs", "ryear", 
                   "csize_<=20", "csize_20-50", "csize_>50", 
                   "cgrade_1-2", "cgrade_3", 
                   "nodes2", "nodes3", 
                   "pgr2", "pgr3"]]
                   
# Validation data 
gbsg5_sel = gbsg5[["pid", "rfs", "ryear", 
                   "csize_<=20", "csize_20-50", "csize_>50", 
                   "cgrade_1-2", "cgrade_3", 
                   "nodes2", "nodes3", 
                   "pgr2", "pgr3"]]

# gbsg5_sel.info()
# rott5_sel.info()

## Fitting Cox proportional hazard model ------------------
# Cox model using lifelines library
# Cox model using also sklearn library
# See: https://scikit-survival.readthedocs.io/en/stable/user_guide/00-introduction.html#

from lifelines import CoxPHFitter
rott5_m01 = rott5_sel[["rfs", "ryear",
                       "csize_20-50", "csize_>50",
                       "cgrade_3",
                       "nodes2", "nodes3"]]
                       
gbsg5_m01 = gbsg5_sel[["rfs", "ryear",
                       "csize_20-50", "csize_>50",
                       "cgrade_3",
                       "nodes2", "nodes3"]]                       
                       
cph = CoxPHFitter()
cph_01 = cph.fit(rott5_m01, duration_col = 'ryear', event_col = 'rfs')
cph_01.print_summary()  # access the individual results using cph.summary


# Model with PGR 
# NOTE:  we did not assess the prediction performance of the extended model
#        in the minimal code
rott5_m02 = rott5_sel[["rfs", "ryear",
                       "csize_20-50", "csize_>50",
                       "cgrade_3",
                       "nodes2", "nodes3",
                       "pgr2", "pgr3"]]
                       
gbsg5_m02 = gbsg5_sel[["rfs", "ryear",
                       "csize_20-50", "csize_>50",
                       "cgrade_3",
                       "nodes2", "nodes3",
                       "pgr2", "pgr3"]]
                       
cph = CoxPHFitter()
cph_02 = cph.fit(rott5_m02, duration_col = 'ryear', event_col = 'rfs')
cph_02.print_summary()  # access the individual results using cph.summary

###############
#  Cox model using sklearn library ----------------
from sksurv.linear_model import CoxPHSurvivalAnalysis

# Import performance measures as c-index, time-dependent AUC, Brier Score
from sksurv.metrics import (
    concordance_index_censored, # Harrell c-index
    concordance_index_ipcw, # Uno's c-index
    cumulative_dynamic_auc,
    brier_score,
)

# Transform outcomes to boolean
rott5_rfs = rott5_m01.rfs.astype(bool)
gbsg5_rfs = gbsg5_m01.rfs.astype(bool)

# Creating structured array - development set
rfs = np.array(rott5_rfs)
ryear = np.array(rott5.ryear)
wtype = np.dtype([('rfs', rfs.dtype),('ryear', ryear.dtype)])
rott5_y = np.empty(len(rfs),dtype = wtype)
rott5_y['rfs']= rfs
rott5_y['ryear']= ryear

del rfs # deleting objects
del ryear

# Creating structured array - validation set
rfs = np.array(gbsg5_rfs)
ryear = np.array(gbsg5.ryear)
wtype = np.dtype([('rfs', rfs.dtype),('ryear', ryear.dtype)])
gbsg5_y = np.empty(len(rfs), dtype = wtype)
gbsg5_y['rfs']= rfs
gbsg5_y['ryear']= ryear

del rfs # deleting objects
del ryear

# Run Cox model using sksurv
from sksurv.linear_model import CoxPHSurvivalAnalysis
estimator = CoxPHSurvivalAnalysis()
rott5_m01_x = rott5_m01.drop(["rfs", "ryear"], axis = 1)
gbsg5_m01_x = gbsg5_m01.drop(["rfs", "ryear"], axis = 1)
cph_01_sk = estimator.fit(rott5_m01_x, rott5_y)
# cph_01_sk.coef_ # coefficients using sklearn
# cph_01.params_  # coefficients using lifelines

###  Discrimination ---------------------------------------

## Harrell c-index -----------------
# cindex_dev_01 = cph_01.score(rott5_m01, scoring_method = "concordance_index")
# cindex_dev_02 = cph_02.score(rott5_m02, scoring_method = "concordance_index")
cindex_Harrell_val_01 = cph_01.score(gbsg5_m01, scoring_method = "concordance_index")
# cindex_val_02 = cph_02.score(gbsg5_m02, scoring_method = "concordance_index")


## Uno c-index
# Calculating prognostic index in the validation data using lifelines package
coeff = cph_01.params_
cov = gbsg5_m01[["csize_20-50", "csize_>50",
                 "cgrade_3",
                 "nodes2", "nodes3"]]
lp_val = np.matmul(cov, coeff)

# using sklearn
pi_gbsg5 = cph_01_sk.predict(gbsg5_m01_x) # calculating prognostic index (PI) using scikit

gbsg5_m01 = gbsg5_m01.assign(pi = pi_gbsg5) # Add PI to validatation set

cindex_Uno_val_01 = concordance_index_ipcw(rott5_y, gbsg5_y, estimate = pi_gbsg5, tau = None)

# Time dependent AUC:
# See: sksurv.metrics.cumulative_dynamic_au
auc_val01 = cumulative_dynamic_auc(rott5_y, gbsg5_y, pi_gbsg5, times = 4.99)[1]

# Bootstrap percentile for confidence intervals
B = 2000
# b_dev_m01 = {}
# b_dev_m02 = {}
b_val_m01 = {}
b_val_m01_x = {}
# b_val_m02 = {}
# b_cindex_dev01 = [0] * B
# b_cindex_dev02 = [0] * B
b_cindex_Harrell_val01 = [0] * B
b_cindex_Uno_val01 = [0] * B
# b_cindex_val02 = [0] * B
b_pi_val01 = [0] * B
b_val_rfs = [0] * B
b_rfs = [0] * B
b_ryear = [0] * B
b_val_y = [0] * B
b_cd_auc_val01 = [0] * B

for j in range(B):
  
  # b_dev_m01[j] = sk.utils.resample(rott5_m01,
  #     replace = True, 
  #     n_samples = len(rott5_m01)) # bootstrapping development data
  #     
  # b_dev_m02[j] = sk.utils.resample(rott5_m02,
  #     replace = True, 
  #     n_samples = len(rott5_m02)) # bootstrapping development data with PGR
      
  b_val_m01[j] = sk.utils.resample(gbsg5_m01,
      replace = True, 
      n_samples = len(gbsg5_m01)) # bootstrapping validation data 
      
  b_val_m01_x[j] = b_val_m01[j].drop(["rfs", "ryear"], axis = 1)
      
  # b_val_m02[j] = sk.utils.resample(gbsg5_m02,
  #     replace = True, 
  #     n_samples = len(gbsg5_m02)) # bootstrapping validation data with PGR
      
  # Bootstrap c-index
  # Harrell c-index
  # b_cindex_dev01[j] = cph_01.score(b_dev_m01[j], scoring_method = "concordance_index")
  # b_cindex_dev02[j] = cph_01.score(b_dev_m02[j], scoring_method = "concordance_index")
  b_cindex_Harrell_val01[j] = cph_01.score(b_val_m01[j], scoring_method = "concordance_index")
  # b_cindex_val02[j] = cph_01.score(b_val_m02[j], scoring_method = "concordance_index")

  # Uno c-index
  # Create bootstrapped structured array
  b_val_rfs[j] = b_val_m01[j].rfs.astype(bool)

  # Creating structured array - development set
  b_rfs[j] = np.array(b_val_rfs[j])
  b_ryear[j] = np.array(b_val_m01[j].ryear)
  wtype = np.dtype([('rfs', b_rfs[0].dtype),('ryear', b_ryear[0].dtype)])
  b_val_y[j] = np.empty(len(b_rfs[0]),dtype = wtype)
  b_val_y[j]['rfs']= b_rfs[j]
  b_val_y[j]['ryear']= b_ryear[j]

  b_cindex_Uno_val01[j] = concordance_index_ipcw(rott5_y, b_val_y[j], 
                                                 estimate = np.array(b_val_m01[j].pi), 
                                                 tau = None)[0]
  b_cd_auc_val01[j] = cumulative_dynamic_auc(rott5_y, b_val_y[j], np.array(b_val_m01[j].pi), times = 4.99)[1]
                                                 
                                                 
  
# Save results
k = 2
res_discr = np.reshape(
  (
   round(cindex_Harrell_val_01, k) ,
   round(np.percentile(b_cindex_Harrell_val01, q = 2.5), k),
   round(np.percentile(b_cindex_Harrell_val01, q = 97.5), k),
   
   round(cindex_Uno_val_01[0], k) ,
   round(np.percentile(b_cindex_Uno_val01, q = 2.5), k),
   round(np.percentile(b_cindex_Uno_val01, q = 97.5), k),
   
   round(auc_val01, k),
   round(np.percentile(b_cd_auc_val01, q = 2.5), k),
   round(np.percentile(b_cd_auc_val01, q = 97.5), k),
   ),
   (3, 3)
)

res_discr = pd.DataFrame(res_discr, 
                         columns = ["Estimate", "2.5 %", "97.5 %"],
                         index = ["Harrell C-index", "Uno C-index", "Uno AUC"])
res_discr


# Calibration --------------

# Mean calibration
# Estimated observed survival using Kaplan-Meier estimator

# Using lifelines
from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter()

kmf.fit(gbsg5_m01.ryear, event_observed = gbsg5_m01.rfs,
timeline = np.array((1, 5))) # observed survival using KM 
obs_t = 1 - kmf.survival_function_.iloc[1, 0]

# Using sklearn
# from sksurv.nonparametric import SurvivalFunctionEstimator, kaplan_meier_estimator
# kaplan_meier_estimator(gbsg5_y["rfs"], gbsg5_y["ryear"])

# Expected mortality
gbsg5_pred_t5 =  1 - cph_01.predict_survival_function(gbsg5_m01, times = 5) # using lifelines
exp_t = gbsg5_pred_t5.mean(axis = 1)

OE_t = obs_t / exp_t
alpha = 0.05
OE_t_lower = OE_t * np.exp(-sp.stats.norm.ppf(1 - alpha / 2) * np.sqrt( 1 / sum(gbsg5.rfs)))
OE_t_upper = OE_t * np.exp(sp.stats.norm.ppf(1 - alpha / 2) * np.sqrt( 1 / sum(gbsg5.rfs)))

# Weak calibration (calibration slope)
# using lifelines (using sklearn seems univariable models cannot be fitted)
cph = CoxPHFitter()
calslope_val = gbsg5_m01[["rfs", "ryear", "pi"]]
cph_calslope = cph.fit(calslope_val, duration_col = 'ryear', event_col = 'rfs')
# cph_calslope.print_summary()  # results

# Results of mean and weak calibration (O/E ratio and calibration slope)
k = 2
res_cal = np.reshape(
  (
    round(OE_t, k),
    round(OE_t_lower, k),
    round(OE_t_upper, k),
  
    round(cph_calslope.summary["coef"], k),
    round(cph_calslope.summary["coef lower 95%"], k),
    round(cph_calslope.summary["coef upper 95%"], k)
    ),
  
  (2, 3)
)

res_cal = pd.DataFrame(res_cal, 
                      columns = ["Estimate", "2.5 %", "97.5 %"],
                      index = ["O/E ratio", "Calibration slope"])
res_cal


# Moderate calibration (calibration-in-the-small)

# Calibration plot 
# Method used: The actual probability is estimated using a 'secondary' Cox PH regression model
# using the predicted probabilities as a covariate.

gbsg5["pred_val"] =  1 - cph_01.predict_survival_function(gbsg5_m01, times = 5).iloc[0, :]
gbsg5["pred_val_cll"] = np.log(-np.log(1 - gbsg5["pred_val"])) # cloglog transformation

# Add restricted cubic splines with three knots
def rcs_3(x, data):
    res_x = np.zeros((len(data), 1))
    qknots = [.1, .5, .9]
    knots = np.quantile(x, q = qknots)
    res_x[:, 0] = (np.power(np.clip((x - knots[0]), a_min = 0, a_max = None), 3) - np.power(np.clip((x - knots[1]), a_min = 0, a_max = None), 3) *((knots[2] - knots[0])/(knots[2] - knots[1])) + np.power(np.clip((x - knots[2]), a_min = 0, a_max = None), 3) * ((knots[1] - knots[0])/(knots[2] - knots[1]))) / ((knots[2] - knots[0])**2)
    return(res_x)
  
gbsg5["pred_val_cll3"] = rcs_3(gbsg5.pred_val_cll, gbsg5) # calculate rcs with 3 knots

# Create dataframe to run the secondary Cox model
pred_val_cal = pd.DataFrame({'ryear' : gbsg5.ryear,
                             'rfs' : gbsg5.rfs,
                             'pred_val_cll' : gbsg5.pred_val_cll,
                             'pred_val_cll3' : gbsg5.pred_val_cll3})
cph = CoxPHFitter()
cph_cal = cph.fit(pred_val_cal, duration_col = 'ryear', event_col = 'rfs')
cph_cal.print_summary()  # access the individual results using cph.summary

# Estimate the observed values
obs_val =   1 - cph_cal.predict_survival_function(gbsg5, times = 5).iloc[0, :]

# Bootstrapping observed values 
# NOTE: percentile bootstrap is an alternative to estimate the confidence intervals
# for the estimated observed values, standard error formula would be ideal
B = 2000
b_gbsg5 = {}
b_pred_val_cal = {}
b_cph_cal = {}
b_obs_val = pd.DataFrame(np.empty((len(gbsg5), B)))
for j in range(B):

  b_pred_val_cal[j] = sk.utils.resample(pred_val_cal,
      replace = True, 
      n_samples = len(pred_val_cal)) 
      
  b_cph_cal = cph.fit(b_pred_val_cal[j], duration_col = 'ryear', event_col = 'rfs')
 
  b_obs_val.iloc [:, j] = 1 - b_cph_cal.predict_survival_function(gbsg5, times = 5).iloc[0, :]  
  
# Create df for calibration plot based on secondary log reg
alpha = 0.05
df_cal = pd.DataFrame({
    'obs' : obs_val,
    'pred' : gbsg5.pred_val,
    'lower_95' : np.quantile(b_obs_val, q = alpha / 2, axis = 1),
    'upper_95' : np.quantile(b_obs_val, q = 1 - alpha / 2, axis = 1)
})

# Sorting
df_cal = df_cal.sort_values(by = ['pred'])

# Calibration plots
# Calibration plots
p1 = plt.plot(df_cal.pred, df_cal.obs, "-", 
         label = "Calibration curve based on secondary Cox model", 
         color = "black")
plt.axline(xy1 = (0, 0), xy2 = (1, 1), linestyle = "--", color = "r", 
          label = "Perfect calibration")
p3 = plt.plot(df_cal.pred, df_cal.lower_95, "--", 
         label = "95% confidence intervals", color = "black")
plt.legend(loc = "upper left")
p4 = plt.plot(df_cal.pred, df_cal.upper_95, "--", 
         label = "95% confidence intervals", color = "black")
plt.xlabel("Predicted risk from developed model")
plt.ylabel("Predicted risk from refitted model")
plt.title("Calibration plot")
plt.show()
plt.clf()
plt.cla()
plt.close('all')

# Calibration metrics based on a secondary logistic regression model
cal_metrics = pd.DataFrame(
  {'ICI' : np.mean(abs(df_cal.obs - df_cal.pred)),
   'E50' : np.median(abs(df_cal.obs - df_cal.pred)),
   'E90' : np.quantile(abs(df_cal.obs - df_cal.pred), 
                       0.9, 
                       interpolation = 'midpoint')}, 
  index = [0]
)
cal_metrics

# Overall performance measures --------------
# Brier Score
from sksurv.metrics import brier_score
survs =  cph_01.predict_survival_function(gbsg5_m01, times = 4.99).iloc[0, :]
brier_01 = brier_score(gbsg5_y, gbsg5_y, estimate = survs, times = 4.99)[1]


# Empty Cox model 
kmf.fit(rott5_m01.ryear, event_observed = rott5_m01.rfs,
timeline = np.array((1, 5))) # observed survival using KM in the development data

survs_null = kmf.survival_function_.iloc[1, 0]
survs_null = pd.Series([survs_null] * len(gbsg5))
brier_null = brier_score(gbsg5_y, gbsg5_y, estimate = survs_null, times = 4.99)[1]

# Scaled brier score -
scaled_brier_01 = 1 - (brier_01 / brier_null)

# Create dataframe with all we need for bootstrap Brier Score
df_brier = pd.DataFrame({'ryear' : gbsg5_m01.ryear,
                         'rfs' : gbsg5_m01.rfs,
                         'surv_null' : survs_null,
                         'surv' : survs})

# Bootstrap percentile for confidence intervals
B = 2000
b_valbrier_m01 = {}
b_val_rfs = [0] * B
b_rfs = [0] * B
b_ryear = [0] * B
b_val_y = [0] * B
b_brier_01 = [0] * B
b_brier_null = [0] * B
b_scaled_brier_01 = [0] * B

for j in range(B):
      
  b_valbrier_m01[j] = sk.utils.resample(df_brier,
      replace = True, 
      n_samples = len(df_brier)) # bootstrapping validation data 
    
      
  # Create bootstrapped structured array
  b_val_rfs[j] = b_valbrier_m01[j].rfs.astype(bool)

  # Creating structured array - development set
  b_rfs[j] = np.array(b_val_rfs[j])
  b_ryear[j] = np.array(b_valbrier_m01[j].ryear)
  wtype = np.dtype([('rfs', b_rfs[0].dtype),('ryear', b_ryear[0].dtype)])
  b_val_y[j] = np.empty(len(b_rfs[0]),dtype = wtype)
  b_val_y[j]['rfs']= b_rfs[j]
  b_val_y[j]['ryear']= b_ryear[j]

  # Bootstrap Brier Score
  b_brier_01[j] = brier_score(b_val_y[j], b_val_y[j], estimate = b_valbrier_m01[j].surv, times = 4.99)[1]
  b_brier_null[j] = brier_score(b_val_y[j], b_val_y[j], estimate = b_valbrier_m01[j].surv_null, times = 4.99)[1]
  b_scaled_brier_01[j] = 1 - (b_brier_01[j] / b_brier_null[j])    
  

# Overall performance results
overall_metrics = np.reshape(
  (brier_01[0],
   np.percentile(b_brier_01, q = 2.5),
   np.percentile(b_brier_01, q = 97.5), 
   
   scaled_brier_01[0],
   np.percentile(b_scaled_brier_01, q = 2.5),
   np.percentile(b_scaled_brier_01, q = 97.5)),
   
   (2, 3)
)

overall_metrics = pd.DataFrame(overall_metrics, 
                               columns = ["Estimate", "2.5 %", "97.5 %"], 
                               index = ["Brier Score", "Scaled Brier"])
overall_metrics

# Clinical utility ------

# Using lifelines
from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter()

val_out = pd.DataFrame({
  "pid" : gbsg5.pid,
  "ryear" : gbsg5.ryear,
  "rfs" : gbsg5.rfs,
  "pred_val" : gbsg5.pred_val})

epsilon = .001 # to select at least an observation
thresholds = np.arange(0, float(val_out[["pred_val"]].max() - epsilon), step = 0.01)
leng = np.arange(0, len(thresholds), step = 1)
f_all = 1 - kmf.fit(val_out.ryear, event_observed = val_out.rfs,
                   timeline = np.array((1, 5))).survival_function_.iloc[1, 0]
NB_all = [0] * len(thresholds)
NB = [0] * len(thresholds)

for j in leng:
  NB_all[j] = f_all - (1 - f_all) * (thresholds[j] / (1 - thresholds[j]))
  tdata = val_out[val_out["pred_val"] > thresholds[j]]
  f_given_exceed = 1 - kmf.fit(val_out[val_out["pred_val"] > thresholds[j]].ryear, 
                               event_observed = val_out[val_out["pred_val"] > thresholds[j]].rfs,
                               timeline = np.array((1, 5))).survival_function_.iloc[1, 0]
                               # NOTE: we can also use tdata directly
  TP = f_given_exceed * (len(tdata) / len(val_out)) 
  FP = (1 - f_given_exceed) * (len(tdata) / len(val_out))
  NB[j]= TP - (FP  * (thresholds[j] / (1 - thresholds[j])))
  
  
# Create dataframe
df_dca = pd.DataFrame({
  'threshold' : thresholds,
  'NB_all' : NB_all,
  'NB' : NB
  }
)

# Plot decision curves
plt.plot(df_dca.threshold, 
         df_dca.NB, "--", 
         color = "green", 
         label = "Prediction model", 
         linewidth = 2.0)
plt.plot(df_dca.threshold, 
         df_dca.NB_all, 
         color = "goldenrod", 
         label = "Treat all",
         linewidth = 2.0)
plt.xlim([0, 1])
plt.ylim([-0.1, 0.6])
plt.xlabel("Threshold probability in %")
plt.ylabel("Net Benefit")
plt.title("Decision curve - validation data")
plt.axhline(y = 0, 
            linestyle = 'dashdot', 
            color = 'pink', 
            label = "Treat none",
            linewidth = 2.0)
plt.legend(loc = "upper right")
plt.show()
plt.show()
plt.clf()
plt.cla()
plt.close('all')

# Next steps: 


