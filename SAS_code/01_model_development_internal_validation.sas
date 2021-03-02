*******************************************************************************;
* Program:			STRATOS Cox model development and internal validation.sas *;
* Author:			David McLernon        					 				  *;
* Date: 			26 FEB 2021          					 				  *;
* Purpose: 			Model development and internal validation				  *;
*					Bootstrap to examine optimism  			 				  *;
* Note:				This programme is not currently automated. It is coded	  *;
*					based on the case study in the article. Therefore, 		  *;
*					adapting this code for your own study will require careful*;
*					editing according to your data							  *;
*******************************************************************************;	

* NB Edit folder locations to your own;
OPTIONS MPRINT MLOGIC NONUMBER NODATE SOURCE2 MAXMEMQUERY=MAX /*(fmtsearch = (formats)*/; 
LIBNAME STRATOS 'C:\Users\sme544\Documents\STRATOS';

*macros required;
%INCLUDE 'C:\Users\sme544\Documents\STRATOS\RCSPLINE macro.SAS';
%INCLUDE 'C:\Users\sme544\Documents\STRATOS\STDCA.SAS';


****** DATA CODING, FUNCTIONAL FORM, AND ASSUMPTION CHECKS;

* The model we base our case study on is the Nottingham Prognostic Index;
* I = 0.2 X size + stage + grade;
* size = max tumour diameter in cms where 0 is <=20, 1=21-50, 2 is >50;
* stage = number of pos lymph nodes where 0=zero, 1=1 to 3, 2=4+ pos nodes;
* grade = tumour grade 0=Grade 1 or 2, 1=Grade 3;


*input data to fit new Nottingham PI and code predictors to be consistent with both development and validation datasets;
*(ideally you do not want to categorise predictors but we have no choice given the limitations of the data we are using);
DATA RO(RENAME=(_D=STATUS));
	SET STRATOS.Rotterdam_br_ca;
	NODESCAT=0;
	if 1<=NODES<=3 then NODESCAT=1;
	if NODES>3 then NODESCAT=2;
	GRADE = GRADE-2;
	SURVTIME=_T;
	*Winzorise PGR to the 99th percentile to deal with very large influential values;
	IF PGR>1360 THEN PGR=1360;
	KEEP PID SIZE NODES NODESCAT GRADE _T SURVTIME _D MENO AGE PGR ER;
	proc freq; tables nodes SIZE NODESCAT GRADE STATUS;
RUN;

*** descriptive statistics for age and PGR;
PROC UNIVARIATE DATA=RO;
	VAR AGE PGR;
RUN;

*get median follow-up using reverse kaplan-meier method;
PROC LIFETEST DATA=RO METHOD=pl ATRISK;
        TIME SURVTIME*STATUS(1);
RUN;

*The ZPH command option requests diagnostics on weight Schoenfeld residuals to check proportional hazards assumption;
*model shows strong evidence of non-proportional hazards for SIZE;
PROC PHREG DATA=RO ZPH(GLOBAL TRANSFORM=LOG);
	CLASS SIZE (REF=FIRST) NODESCAT (REF='1') GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE / TIES=EFRON RL;
RUN;

*Administrative censor at 5 years since this is our prediction horizon; 
DATA ROTT;
	SET RO;
	*administrative censoring at 5 years;
	IF _T > 5 THEN STATUS=0;
	IF SURVTIME > 5 THEN SURVTIME=5;
RUN;


*CHECK FUNCTIONAL FORM OF PGR;
* Code PGR as a restricted cubic spline with 3 knots;
*First calculate the 10th, 50th and 90th percentiles for knots;
PROC UNIVARIATE DATA = ROTT;
	VAR PGR;
	OUTPUT OUT=KNOTS PCTLPRE=P_PGR PCTLPTS= 10 50 90;
RUN;

PROC PRINT DATA=KNOTS; RUN;

/* here we find the following values:
P_PGR10 P_PGR50 P_PGR90 
0 		41 		486 
*/

*Use Frank Harrell's RCSPLINE macro for calculating the spline terms - found here: http://biostat.mc.vanderbilt.edu/wiki/Main/SasMacros;
DATA FFPGR;
	SET ROTT;
	%RCSPLINE(PGR, 0, 41, 486);
RUN;

*make copy to keep for making predictions in PROC PHREG;
DATA FFPGR1;
	SET FFPGR;
RUN;

*The following code will produce Suppl Fig 1 - restricted cubic spline plot for PGR;

*Fit Cox model with PGR terms and save predictions at 5 years;
PROC PHREG DATA=FFPGR;
	MODEL SURVTIME*STATUS(0)=PGR PGR1/ TIES=EFRON;
	TEST1: TEST PGR, PGR1;
	TEST2: TEST PGR+PGR1=0;
	BASELINE COVARIATES=FFPGR1 OUT=PGRVAL SURVIVAL=PREDPGR LOWER=PREDPGRL UPPER=PREDPGRUP TIMELIST=5;
RUN;

*- Manipulate data so that we can plot the diagonal ref line and probability of death (not survival);
DATA PGRVAL1;
	SET PGRVAL;
	IF PID=3 THEN DIAG1=0;
	IF PID=3 THEN DIAG2=0;
	IF PID=7 THEN DIAG1=100;
	IF PID=7 THEN DIAG2=100;
	PREDPGR_DTH=1-PREDPGR;
	PREDPGR_LOW=1-PREDPGRL;
	PREDPGR_UPP=1-PREDPGRUP;
	KEEP PID PGR PREDPGR_DTH PREDPGR_LOW PREDPGR_UPP DIAG1 DIAG2;
RUN;

PROC SORT DATA=PGRVAL1;
	BY PGR;
RUN;

title;
FOOTNOTE;
*- this allows editing of the .sge file!;
ODS LISTING SGE=ON STYLE=PRINTER IMAGE_DPI=300 GPATH='C:\';
ODS GRAPHICS ON / RESET=ALL NOBORDER OUTPUTFMT=TIFF /*WIDTH=4IN*/ IMAGENAME="Spline of PGR" ANTIALIAS=OFF/*ANTIALIASMAX=*/;

*- Plot the restricted cubic spline of PGR with predicted probability of death within 5 years;
PROC SGPLOT DATA=PGRVAL1 NOAUTOLEGEND;
XAXIS       	LABEL="PGR" VALUES=(0 TO 1400 BY 100) /*edit labels and value range as appropriate here and below*/;
YAXIS        	LABEL="Predicted probability of mortality" VALUES=(0 TO 0.5 BY 0.1);
  BAND X=PGR LOWER=PREDPGR_LOW UPPER=PREDPGR_UPP / NOFILL LINEATTRS=(COLOR=BLACK PATTERN=MEDIUMDASH THICKNESS=3) NOEXTEND OUTLINE;
  SERIES Y=PREDPGR_DTH X=PGR / LINEATTRS=(COLOR=BLACK THICKNESS=3) ;
run;

ods graphics off;


****************************************FIT SIMPLE MODEL FIRST (without PGR)******************************************************;

*store original model;
PROC PHREG DATA=ROTT;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE   / TIES=EFRON RL;
	STORE SimpModel;
	*BASELINE COVARIATES=ROTVAL OUT=ROTTVALRD XBETA=XB /*TIMELIST=5 SURVIVAL=FIVEYRSURV*/;
	OUTPUT OUT=ROTTX XBETA=XB;
RUN;

*** Estimating Baseline Survival Function under PH;
* - I calculated the baseline survival backwards using predicted probs and PI for 2 different patients and got S0=0.82267 (S0=0.75566136 WHEN NODES REF=1) which would;
* - suggest that the lowest values are not used
*use cards;
data inrisk;
	input SIZE NODESCAT GRADE;
	cards;
1 0 0
;

proc phreg data=ROTT;
CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
model SURVTIME*STATUS(0)= SIZE NODESCAT GRADE / TIES=EFRON RL ;
baseline covariates=inrisk out=outph survival=ps timepoint=(1,2,3,4,5) / method=breslow;
run;

proc print data=outph;
run;

/*Obs SIZE NODESCAT GRADE SURVTIME ps 
1 1 0 0 1 0.97067 
2 1 0 0 2 0.92195 
3 1 0 0 3 0.88126 
4 1 0 0 4 0.85000 
5 1 0 0 5 0.82267 
 */



************** FIRST, RESAMPLE 500 VERSIONS OF THE EXTERNAL DATASET WITH REPLACEMENT TO ALLOW CALCULATION OF 95% CI AND INTERNAL VALIDATION ********************;


*some bootstrap code to edit;

SASFILE FFPGR LOAD; /* a way of loading the dataset into RAM - speeds it up */

PROC SURVEYSELECT DATA=FFPGR OUT=OUTBOOT /* Use PROC SURVEYSELECT and provide input and output dataset names */
SEED=4817 /* can enter 0 for it to select a random seed but remember to type it in here from the output otherwise cannot replicate results */
METHOD=URS /* Unrestricted Random Sampling - simple random sampling */
SAMPRATE=1 /* can accept proportions or percentages but we want n to be size of original database so =1 (or 100) */
OUTHITS /* with replacement */
REP=500; /* number of bootstrap samples */
RUN;

SASFILE FFPGR CLOSE; /* closes frees RAM buffers when done */

ODS LISTING CLOSE; /* turns off ODS listing so no printing of all output. Better than using NOPRINT as it doesn't allow storage of data in output dataset at end */

*save for later use in PGR programme;
DATA STRATOS.OUTBOOT;
	SET OUTBOOT;
RUN;


*APPLY BOOTSTRAP MODEL TO BOOTSTRAPPED DATASET AND THE ORIGINAL DATASET;

PROC PHREG DATA=OUTBOOT NOPRINT;
	*WHERE REPLICATE in (1,2) /*useful line to use when testing out*/;
	BY REPLICATE;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE   / TIES=EFRON RL;
	STORE SimpModelBoot;
	*baseline statement applies each of the 500 models to the original dataset - ext validation;
	BASELINE COVARIATES=ROTT OUT=ROTTBOOT XBETA=XB TIMELIST=5 /*SURVIVAL=FIVEYRSURV*/;
	*output statement applies each of the 500 models to the corresponding dataset - apparent;
	OUTPUT OUT=ROTTAPP XBETA=XB;
RUN;

PROC SORT DATA=ROTTBOOT;
	BY PID REPLICATE;
RUN;

PROC SORT DATA=ROTT;
	BY PID;
RUN;

*select first xbeta per person per replicate and merge original data ; 
DATA ROTTBOOT1;
	SET ROTTBOOT;
	BY PID REPLICATE;
	IF FIRST.REPLICATE;
	KEEP REPLICATE PID XB;
RUN;

DATA ROTTBOOTRD;
	MERGE ROTTBOOT1 ROTT;
	BY PID;
	PROC SORT; BY REPLICATE PID;
RUN;




*****GLOBAL ASSESSMENT FOR OVERALL PERFORMANCE;



***** WE FIRST PRESENT (FOR INTEREST) HOW TO CALCULATE COX-SNELL, NAGELKERKE, AND SCHEMPER AND HENDERSON R SQUARED;
***calculate Generalised R squared - Cox and Snell 1989, Magee 1990, Allison;
*Global LR chi-sq stat for model with all 3 preds Gsq =  Likelihood Ratio = 465.8;
* Cox and Snell Rsq = 1 - exp(-Gsq/n) = 1 - exp (-466/2982)= 14.5;
* Schemper and Henderson = 10.7% (outputted from below model when you block the ODS OUTPUT line);
*The EV option in the PROC statement provides Schemper-Henderson measure (Schemper and Henderson 2000) of the proportion of variation that is explained by a Cox regression;
ods select none;
*fit model and extract likelihood ratio chi-square statistic for testing global=0;
PROC PHREG DATA=ROTT ev;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE / TIES=EFRON RL;
	ODS OUTPUT GlobalTests=GlobalChiSq CensoredSummary=Total FitStatistics=Null;
RUN;
ods select all;

DATA GlobalChiSq1;
	SET GlobalChiSq;
	IF TEST NOT IN ('Likelihood Ratio') THEN DELETE;
	FLAG=1;
	KEEP FLAG ChiSq;
RUN;

DATA Total1;
	SET Total;
	FLAG=1;
	KEEP FLAG Total;
RUN;

DATA Null1;
	SET Null;
	IF Criterion NOT IN ('-2 LOG L') THEN DELETE;
	FLAG=1;
	KEEP FLAG WithOutCovariates;
RUN;

*dataset contains Cox Snell and Nagelkerke's R squared;
DATA RSQUARE;
	MERGE GlobalChiSq1 Total1 Null1;
	BY FLAG;
	RSQ_CS = 1 - EXP((-1)*ChiSq/Total);
	RSQ_N = RSQ_CS/(1-(EXP((-1)*WithOutCovariates/Total)));
	KEEP RSQ_CS RSQ_N;
RUN;

** AS IN PAPER, WE CALCULATE Royston and Sauerbrei’s R SQUARED D AND BOOTSTRAP PERFORMANCE;

*Royston and Sauerbrei’s D;
/*1. To compute D, first the Cox PH model is fitted. 
  2. Then the prognostic index of the model, XB, is transformed to give standard normal order rank statistics (rankits - formed using Blom’s approximation). 
  3. The rankits are multiplied by a factor of SQRT(8/pi) to give Zi (i = 1, n subjects). 
  4. Finally a Cox PH model is fitted to these values; D is the coefficient of Z, say a*, from this second model. 
NOTE: Royston and Sauerbrei (2004) showed that D most accurately measures separation of survival curves when the underlying prognostic index values, XB, are normally distributed. 
The regression on the Z in the second model is then linear and a* is an approximately unbiased estimate of a. They explained that when the XB are not normally distributed, 
linearity in the second model breaks down. D still measures separation because a* in the second model still estimates a, but with bias.*/

**Estimate D for apparent validation;
 *2. Then the prognostic index of the model, XB, is transformed to give standard normal order rank statistics (rankits - formed using Blom’s approximation); 

PROC UNIVARIATE DATA=ROTTX;
	HISTOGRAM;
	VAR XB;
RUN;

PROC RANK DATA=ROTTX NORMAL=BLOM OUT= ROTT_d;
	VAR XB;
RUN;

PROC UNIVARIATE DATA=ROTT_d;
	HISTOGRAM;
	VAR XB;
RUN;

proc freq data=rott_d;
run;

*  3. The rankits are multiplied by a factor of SQRT(8/pi) to give Zi (i = 1, n subjects);
DATA X;
	SET ROTT_d;
	PIE = CONSTANT("pi");
	Z=XB/(SQRT(8/PIE));
RUN;

PROC UNIVARIATE DATA=x;
	HISTOGRAM;
	VAR z;
RUN;

*  4. Finally a Cox PH model is fitted to these values, D is the coefficient of Z, say a*, from this second model;
*Royston's D is parameter estimate of Z, D = 1.05755;
ODS SELECT NONE;
PROC PHREG DATA=X;
	MODEL SURVTIME*STATUS(0)=Z / TIES=EFRON;
	ODS OUTPUT ParameterEstimates=PAREST;
RUN;
ODS SELECT ALL;

*Calculate RsqD from D;
DATA ROYSTON(RENAME=(ESTIMATE=D));
	SET PAREST;
	PIE = CONSTANT("pi");
	*R2D=((ESTIMATE**2)/(8/PI))/(((ESTIMATE**2)/(8/PI))+(1));
	R2D=((ESTIMATE**2)/(8/PIE))/(((ESTIMATE**2)/(8/PIE))+((PIE**2)/6));
	VALIDATE='Apparent';
	ind=1;
	KEEP VALIDATE ESTIMATE R2D IND;
	PROC PRINT; TITLE "Royston's D and RsqD";
RUN;



*** - GET 95% CI;

DATA GETXB;
	SET ROTTX;
	KEEP PID XB;
	PROC SORT;
		BY PID;
RUN;

PROC SORT DATA=OUTBOOT;
	BY PID;
RUN;

DATA OUTBOOT1;
	MERGE OUTBOOT GETXB;
	BY PID;
RUN;

PROC SORT DATA=OUTBOOT1;
	BY REPLICATE PID;
RUN;

PROC RANK DATA=OUTBOOT1 NORMAL=BLOM OUT= ROTT_dAX;
	BY REPLICATE;
	VAR XB;
RUN;

DATA XAX;
	SET ROTT_dAX;
	PIE = CONSTANT("pi");
	Z=XB/(SQRT(8/PIE));
RUN;

ODS SELECT NONE;
PROC PHREG DATA=XAX;
	BY REPLICATE;
	MODEL SURVTIME*STATUS(0)=Z / TIES=EFRON;
	ODS OUTPUT ParameterEstimates=PARESTAX;
RUN;
ODS SELECT ALL;

DATA ROYSTONAX(RENAME=(ESTIMATE=AppD));
	SET PARESTAX;
	PIE = CONSTANT("pi");
	AppR2D=((ESTIMATE**2)/(8/PIE))/(((ESTIMATE**2)/(8/PIE))+((PIE**2)/6));
	KEEP REPLICATE ESTIMATE AppR2D;
RUN;

PROC UNIVARIATE DATA=ROYSTONAX NOPRINT;
	VAR AppD AppR2D;
	OUTPUT OUT = CONFINT PCTLPTS=2.5 97.5 PCTLPRE= AppD_ AppR2D_ PCTLNAME=LOWER95 UPPER95;
RUN;

DATA CONFINT1;
	SET CONFINT;
	IND=1;
RUN;

DATA ROYSTONAX2;
	MERGE ROYSTON CONFINT1;
	BY IND;
	PROC PRINT; TITLE "Royston's D and R2D with 95% CI";
RUN;




*** - NOW WE ASSESS INTERNAL VALIDATION WITH BOOTSTRAP;
* APPARENT BOOTSTRAP;

PROC RANK DATA=ROTTAPP NORMAL=BLOM OUT= ROTT_dA;
	BY REPLICATE;
	VAR XB;
RUN;

DATA XA;
	SET ROTT_dA;
	PIE = CONSTANT("pi");
	Z=XB/(SQRT(8/PIE));
RUN;

ODS SELECT NONE;
PROC PHREG DATA=XA;
	BY REPLICATE;
	MODEL SURVTIME*STATUS(0)=Z / TIES=EFRON;
	ODS OUTPUT ParameterEstimates=PARESTA;
RUN;
ODS SELECT ALL;

DATA ROYSTONA(RENAME=(ESTIMATE=AppD));
	SET PARESTA;
	PIE = CONSTANT("pi");
	AppR2D=((ESTIMATE**2)/(8/PIE))/(((ESTIMATE**2)/(8/PIE))+((PIE**2)/6));
	KEEP REPLICATE ESTIMATE AppR2D;
RUN;

**EXTERNAL BOOTSTRAP;
PROC RANK DATA=ROTTBOOTRD NORMAL=BLOM OUT= ROTVALRDD;
	BY REPLICATE;
	VAR XB;
RUN;

DATA XD;
	SET ROTVALRDD;
	PIE = CONSTANT("pi");
	Z=XB/(SQRT(8/PIE));
RUN;

*Royston's D is parameter estimate of Z, D = 0.95979;
ODS SELECT NONE;
PROC PHREG DATA=XD;
	BY REPLICATE;
	MODEL SURVTIME*STATUS(0)=Z / TIES=EFRON;
	ODS OUTPUT ParameterEstimates=PAREST;
RUN;
ODS SELECT ALL;

DATA ROYSTONEXT(RENAME=(ESTIMATE=BootD));
	SET PAREST;
	PIE = CONSTANT("pi");
	BootR2D=((ESTIMATE**2)/(8/PIE))/(((ESTIMATE**2)/(8/PIE))+((PIE**2)/6));
	KEEP REPLICATE ESTIMATE BootR2D;
RUN;

*merge;
DATA ROYSTONBOOT;
	MERGE ROYSTONa ROYSTONEXT;
	BY REPLICATE;
	OptD = AppD - BootD;
	OptR2D = AppR2D - BootR2D;
	PROC PRINT; TITLE "Royston's D and RsqD Bootstrap values";
RUN;




*****FIXED TIME POINT ASSESSMENT FOR OVERALL PERFORMANCE;



*Brier - see Graf et al 1999;
*scaled brier score = 1-brier/briermax, briermax=mean(1-surv)= bmax * (1-bmax)2 + (1-bmax) * bmax2;
/*Weights - Kaplan Meier with censor as event;
3 groups - Those who have the event up to fixed event time of interest, and those who go beyond fixed t (could be event or event free), also those censored up to;
fixed t (only first 2 groups contribute to BS but all to weights;
Then group 1 calc -surv^2, and group 2 calc (1-surv)^2 where surv is probability of surv at t*; 
Check weight for group 2 (G(t*)) should be came for all patients in that group;*/

title ' ';

****First for apparent validation i.e. model development;
*calculate weights for apparent validation;

PROC LIFETEST DATA=ROTTX METHOD=pl ATRISK OUTSURV=OUTKM /*NOPRINT*/;
        TIME SURVTIME*STATUS(1);
RUN;

*CODE THE 3 GROUPINGS AND DUPLICATE SURVIVAL TIME AS SAS WILL REMOVE THE OFFICIAL SURVIVAL TIME VARIABLE IN BASELINE STATEMENT;
DATA ROTT_B;
	SET ROTTX;
	IF SURVTIME<=4.95 AND STATUS=1 THEN CAT=1;
	IF SURVTIME>4.95 THEN CAT=2;
	IF SURVTIME<=4.95 AND STATUS=0 THEN CAT=3;
	TIME=SURVTIME;
RUN;

*NOW ESTIMATE SURVIVAL AT 5 YEARS IN DEVELOPMENT DATASET;
PROC PHREG DATA=ROTT;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE / TIES=EFRON RL;
	*APPARENT;
	BASELINE COVARIATES=ROTT_B OUT=ROTT_BS TIMELIST=5 SURVIVAL=FIVEYRSURV/ method=breslow;
RUN;

*code up the fixed time of 5 years;
DATA ROTT_BS1(RENAME=(TIME=SURVTIME));
	SET ROTT_BS;
	TIME1=TIME;
	IF TIME1>5 THEN TIME1=5;
	DROP SURVTIME;
	PROC SORT; BY TIME1;
RUN;

*MERGE THE KAPLAN-MEIER WEIGHTS TO THE APPROPRIATE TIMES;
DATA OUTKM1(RENAME=(SURVTIME=TIME1));
	SET OUTKM;
	IF SURVTIME>5 THEN DELETE;
	WEIGHT=1/SURVIVAL;
	KEEP SURVTIME WEIGHT;
RUN;


DATA ROTT_BS2;
	MERGE ROTT_BS1 OUTKM1;
	BY TIME1;
	IF pID=. THEN DELETE;
RUN;

DATA ROTT_BS3;
	SET ROTT_BS2;
	RETAIN _WEIGHT;
	IF NOT MISSING(WEIGHT) THEN _WEIGHT=WEIGHT;
	ELSE WEIGHT=_WEIGHT;
	IF CAT=3 THEN WEIGHT=0;
	IF TIME1=0 THEN DELETE;
	IF CAT=1 THEN CONTRIB=(-FIVEYRSURV)**2;
	IF CAT=2 THEN CONTRIB=(1-FIVEYRSURV)**2;
	IF CAT=3 THEN CONTRIB=0;
	BS=CONTRIB*WEIGHT;
	DROP _WEIGHT;
RUN;

*ESTIMATE BRIER SCORE;
PROC UNIVARIATE DATA=ROTT_BS3 NOPRINT;
	VAR BS WEIGHT;
	OUTPUT OUT=SUMS SUM=SBS SWEIGHT;
	PROC PRINT; 
RUN;

DATA SUMS;
	RETAIN SWEIGHT SBS BRIER;
	SET SUMS;
	BRIER = (1/SWEIGHT)*SBS;
	SWEIGHT=LEFT(SWEIGHT);
	IND=1;
	TITLE 'Brier score';
	PROC PRINT;
RUN;

PROC FREQ DATA=ROTT_BS3;
	TABLE CAT;
RUN;


******- BOOTSTRAP THE 95% CI FOR BRIER SCORE;
****First for apparent validation i.e. model development;
*calculate weights for apparent validation;
PROC SORT DATA=OUTBOOT;
	BY REPLICATE PID;
RUN;

PROC LIFETEST DATA=OUTBOOT METHOD=pl OUTSURV=OUTKM_ NOPRINT;
		BY REPLICATE;
        TIME SURVTIME*STATUS(1);
RUN;

DATA OUTBOOTX;
	SET OUTBOOT;
	KEEP REPLICATE PID;
	PROC SORT; BY PID;
RUN;

*MERGE BOOTSTRAP SET TO THE ROTT_BS1 DATASET FROM EARLIER WITH PREDICTIONS AT 5 YEARS; 
PROC SORT DATA=ROTT_BS1;
	BY PID;
RUN;

DATA ROTT_BS1_;
	MERGE OUTBOOTX ROTT_BS1;
	BY PID;
RUN;

PROC SORT DATA=ROTT_BS1_;
	BY REPLICATE TIME1;
RUN;

*MERGE THE KAPLAN-MEIER WEIGHTS TO THE APPROPRIATE TIMES;
DATA OUTKM1_(RENAME=(SURVTIME=TIME1));
	SET OUTKM_;
	IF SURVTIME>5 THEN DELETE;
	WEIGHT=1/SURVIVAL;
	KEEP REPLICATE SURVTIME WEIGHT;
RUN;

DATA ROTT_BS2_;
	MERGE ROTT_BS1_ OUTKM1_;
	BY REPLICATE TIME1;
	IF pID=. THEN DELETE;
RUN;

DATA ROTT_BS3_;
	SET ROTT_BS2_;
	RETAIN _WEIGHT;
	IF NOT MISSING(WEIGHT) THEN _WEIGHT=WEIGHT;
	ELSE WEIGHT=_WEIGHT;
	IF CAT=3 THEN WEIGHT=0;
	IF TIME1=0 THEN DELETE;
	IF CAT=1 THEN CONTRIB=(-FIVEYRSURV)**2;
	IF CAT=2 THEN CONTRIB=(1-FIVEYRSURV)**2;
	IF CAT=3 THEN CONTRIB=0;
	BS=CONTRIB*WEIGHT;
	DROP _WEIGHT;
RUN;

*ESTIMATE BRIER SCORE;
PROC UNIVARIATE DATA=ROTT_BS3_ NOPRINT;
	BY REPLICATE;
	VAR BS WEIGHT;
	OUTPUT OUT=SUMS_ SUM=SBS SWEIGHT;
	PROC PRINT; 
RUN;

DATA SUMS_;
	RETAIN SWEIGHT SBS BRIER;
	SET SUMS_;
	BRIER = (1/SWEIGHT)*SBS;
	SWEIGHT=LEFT(SWEIGHT);
	TITLE 'Brier score';
	PROC PRINT;
RUN;

*95% CIs are presented with the scaled Brier results later;



*** NOW BRIER BOOTSTRAP PERFORMANCE;

PROC LIFETEST DATA=ROTTAPP METHOD=pl ATRISK /*PLOTS=(S, LS, LLS)*/ OUTSURV=OUTKMA noprint;
		BY REPLICATE;
        TIME SURVTIME*STATUS(1);
RUN;

*for BOOTSTRAPPED external validation;
PROC LIFETEST DATA=ROTT METHOD=pl ATRISK /*PLOTS=(S, LS, LLS)*/ OUTSURV=OUTKM_EXT noprint;
        TIME SURVTIME*STATUS(1);
RUN;

*CODE THE 3 GROUPINGS AND DUPLICATE SURVIVAL TIME AS SAS WILL REMOVE THE OFFICIAL SURVIVAL TIME VARIABLE IN BASELINE STATEMENT;
DATA ROTT_BA;
	SET ROTTAPP;
	IF SURVTIME<=4.95 AND STATUS=1 THEN CAT=1;
	IF SURVTIME>4.95 THEN CAT=2;
	IF SURVTIME<=4.95 AND STATUS=0 THEN CAT=3;
	TIME=SURVTIME;
RUN;

DATA ROTT_B_EX;
	SET ROTT;
	IF SURVTIME<=4.95 AND STATUS=1 THEN CAT=1;
	IF SURVTIME>4.95 THEN CAT=2;
	IF SURVTIME<=4.95 AND STATUS=0 THEN CAT=3;
	TIME=SURVTIME;
RUN;

*NOW ESTIMATE SURVIVAL AT 5 YEARS;
PROC PHREG DATA=ROTTAPP;
	BY REPLICATE;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE   / TIES=EFRON RL;
	*APPARENT;
	BASELINE COVARIATES=ROTT_BA OUT=ROTT_BSA TIMELIST=5 SURVIVAL=FIVEYRSURV/ method=breslow;
RUN;

PROC PHREG DATA=ROTTAPP;
	BY REPLICATE;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE   / TIES=EFRON RL;
	*EXTERNAL;
	BASELINE COVARIATES=ROTT_B_EX OUT=ROTTVAL_BS TIMELIST=5 SURVIVAL=FIVEYRSURV/ method=breslow;
RUN;

*apparent brier;
*MERGE THE KAPLAN-MEIER WEIGHTS TO THE APPROPRIATE TIMES;
DATA ROTT_BS1A(RENAME=(TIME=SURVTIME));
	SET ROTT_BSA;
	TIME1=TIME;
	IF TIME1>5 THEN TIME1=5;
	DROP SURVTIME;
	PROC SORT; BY REPLICATE TIME1;
RUN;

DATA OUTKM1A(RENAME=(SURVTIME=TIME1));
	SET OUTKMA;
	*IF _CENSOR_=1 THEN DELETE;
	IF SURVTIME>5 THEN DELETE;
	WEIGHT=1/SURVIVAL;
	KEEP REPLICATE SURVTIME WEIGHT;
RUN;

DATA ROTT_BS2A;
	MERGE ROTT_BS1A OUTKM1A;
	BY REPLICATE TIME1;
	IF PID=. THEN DELETE;
RUN;

DATA ROTT_BS3A;
	SET ROTT_BS2A;
	RETAIN _WEIGHT;
	IF NOT MISSING(WEIGHT) THEN _WEIGHT=WEIGHT;
	ELSE WEIGHT=_WEIGHT;
	IF CAT=3 THEN WEIGHT=0;
	IF TIME1=0 THEN DELETE;
	IF CAT=1 THEN CONTRIB=(-FIVEYRSURV)**2;
	IF CAT=2 THEN CONTRIB=(1-FIVEYRSURV)**2;
	IF CAT=3 THEN CONTRIB=0;
	BS=CONTRIB*WEIGHT;
	DROP _WEIGHT;
RUN;

*ESTIMATE BRIER SCORE;
PROC UNIVARIATE DATA=ROTT_BS3A NOPRINT;
	BY REPLICATE;
	VAR BS WEIGHT;
	OUTPUT OUT=SUMSA SUM=AppSBS SWEIGHT;
	*PROC PRINT; 
RUN;

DATA SUMSA;
	RETAIN SWEIGHT AppSBS AppBRIER;
	SET SUMSA;
	AppBRIER = (1/SWEIGHT)*AppSBS;
	SWEIGHT=LEFT(SWEIGHT);
	*TITLE 'Brier score';
	*PROC PRINT;
RUN;

*external validation brier;
*MERGE THE KAPLAN-MEIER WEIGHTS TO THE APPROPRIATE TIMES;
DATA ROTTVAL_BS1(RENAME=(TIME=SURVTIME));
	SET ROTTVAL_BS;
	TIME1=TIME;
	IF TIME1>5 THEN TIME1=5;
	DROP SURVTIME;
	PROC SORT; BY TIME1;
RUN;

DATA OUTKM_EXT1(RENAME=(SURVTIME=TIME1));
	SET OUTKM_EXT;
	*IF _CENSOR_=1 THEN DELETE;
	IF SURVTIME>5 THEN DELETE;
	WEIGHT=1/SURVIVAL;
	KEEP SURVTIME WEIGHT;
RUN;

DATA ROTTVAL_BS2;
	MERGE ROTTVAL_BS1 OUTKM_EXT1;
	BY TIME1;
	IF PID=. THEN DELETE;
RUN;

DATA ROTTVAL_BS3;
	SET ROTTVAL_BS2;
	RETAIN _WEIGHT;
	IF NOT MISSING(WEIGHT) THEN _WEIGHT=WEIGHT;
	ELSE WEIGHT=_WEIGHT;
	IF CAT=3 THEN WEIGHT=0;
	IF TIME1=0 THEN DELETE;
	IF CAT=1 THEN CONTRIB=(-FIVEYRSURV)**2;
	IF CAT=2 THEN CONTRIB=(1-FIVEYRSURV)**2;
	IF CAT=3 THEN CONTRIB=0;
	BS=CONTRIB*WEIGHT;
	DROP _WEIGHT;
RUN;

*ESTIMATE BRIER SCORE;
proc sort data= ROTTVAL_BS3;
	BY REPLICATE PID;
RUN;

PROC UNIVARIATE DATA=ROTTVAL_BS3 NOPRINT;
	BY REPLICATE;
	VAR BS WEIGHT;
	OUTPUT OUT=SUMSEX SUM=BootSBS SWEIGHT;
	*PROC PRINT; 
RUN;

DATA SUMSEX;
	RETAIN SWEIGHT BootSBS BootBRIER;
	SET SUMSEX;
	BootBRIER = (1/SWEIGHT)*BootSBS;
	SWEIGHT=LEFT(SWEIGHT);
	*TITLE 'External Brier score';
	*PROC PRINT;
RUN;







********************************Scaled Brier;

/*Scaled Brier = 1 - (model Brier score/null model Brier score), where null model Brier score is null cox; 
100% is perfect, <0 is useless, higher better, hamrful models <0;*/

*apparent Scaled Brier first;
*NOW ESTIMATE SURVIVAL AT 5 YEARS for NULL model;
PROC PHREG DATA=ROTTX;
	MODEL SURVTIME*STATUS(0)= / TIES=EFRON RL;
	*APPARENT;
	BASELINE COVARIATES=ROTT_B OUT=ROTT_BSnull TIMELIST=5 SURVIVAL=FIVEYRSURV_null;
RUN;

DATA ROTT_BS1null;
	SET ROTT_BSnull;
	KEEP PID FIVEYRSURV_null;
	PROC SORT; BY PID;
RUN;

PROC SORT DATA=ROTT_BS3;
	BY PID;
RUN;

*apparent brier for null model;
*MERGE THE NULL MODEL SURVIVAL PROBABILITIES TO BRIER SCORE DATASET FROM EARLIER;
DATA ROTT_BS4;
	MERGE ROTT_BS3 ROTT_BS1null;
	BY PID;
RUN;

DATA ROTT_BS5;
	SET ROTT_BS4;
	IF CAT=1 THEN CONTRIB_NULL=(-FIVEYRSURV_null)**2;
	IF CAT=2 THEN CONTRIB_NULL=(1-FIVEYRSURV_null)**2;
	IF CAT=3 THEN CONTRIB_NULL=0;
	BS_NULL=CONTRIB_NULL*WEIGHT;
	DROP FIVEYRSURV_null;
RUN;

*ESTIMATE BRIER SCORE FOR NULL MODEL;
PROC UNIVARIATE DATA=ROTT_BS5 NOPRINT;
	VAR BS_NULL WEIGHT;
	OUTPUT OUT=SUMNULL SUM=SBS_NULL SWEIGHT;
	PROC PRINT; 
RUN;

DATA SUMNULL;
	RETAIN SWEIGHT SBS_NULL;
	SET SUMNULL;
	SWEIGHT=LEFT(SWEIGHT);
RUN;

*CALCULATE Scaled Brier;
DATA SCALEDB;
	MERGE SUMNULL SUMS;
	BY SWEIGHT;
	NULL_BRIER = (1/SWEIGHT)*SBS_NULL;
	SCALED_B = 1-(BRIER/NULL_BRIER);
	TITLE 'Brier score and Scaled Brier';
	PROC PRINT;
RUN;


***** 95% CI;


*NOW ESTIMATE SURVIVAL AT 5 YEARS for NULL model;
DATA L2ROTT_B_EXX;
	SET OUTBOOT;
	IF SURVTIME<=4.95 AND STATUS=1 THEN CAT=1;
	IF SURVTIME>4.95 THEN CAT=2;
	IF SURVTIME<=4.95 AND STATUS=0 THEN CAT=3;
	*TIME=SURVTIME;
	PROC SORT; BY REPLICATE SURVTIME;
RUN;

PROC PHREG DATA=OUTBOOT NOPRINT;
	BY REPLICATE;
	MODEL SURVTIME*STATUS(0)= / TIES=EFRON RL;
	BASELINE COVARIATES=L2ROTT_B_EXX OUT=ROTT_BSnull_ TIMELIST=5 SURVIVAL=FIVEYRSURV_null;
RUN;

DATA ROTT_BS1null_;
	SET ROTT_BSnull_;
	KEEP REPLICATE PID FIVEYRSURV_null;
	PROC SORT; BY REPLICATE PID;
RUN;

PROC SORT DATA=ROTT_BS3_;
	BY REPLICATE PID;
RUN;

*apparent brier for null model;
*MERGE THE NULL MODEL SURVIVAL PROBABILITIES TO BRIER SCORE DATASET FROM EARLIER;
DATA ROTT_BS4_;
	MERGE ROTT_BS3_ ROTT_BS1null_;
	BY REPLICATE PID;
RUN;

DATA ROTT_BS5_;
	SET ROTT_BS4_;
	IF CAT=1 THEN CONTRIB_NULL=(-FIVEYRSURV_null)**2;
	IF CAT=2 THEN CONTRIB_NULL=(1-FIVEYRSURV_null)**2;
	IF CAT=3 THEN CONTRIB_NULL=0;
	BS_NULL=CONTRIB_NULL*WEIGHT;
	DROP FIVEYRSURV_null;
RUN;

*ESTIMATE BRIER SCORE FOR NULL MODEL;
PROC UNIVARIATE DATA=ROTT_BS5_ NOPRINT;
	BY REPLICATE;
	VAR BS_NULL WEIGHT;
	OUTPUT OUT=SUMNULL_ SUM=SBS_NULL SWEIGHT;
	PROC PRINT; 
RUN;

DATA SUMNULL_;
	RETAIN SWEIGHT SBS_NULL;
	SET SUMNULL_;
	SWEIGHT=LEFT(SWEIGHT);
RUN;

*CALCULATE Scaled Brier;
DATA SCALEDB_;
	MERGE SUMNULL_ SUMS_;
	BY REPLICATE;
	NULL_BRIER = (1/SWEIGHT)*SBS_NULL;
	SCALED_B = 1-(BRIER/NULL_BRIER);
	*PROC PRINT;
RUN;

PROC UNIVARIATE DATA=SCALEDB_ NOPRINT;
	VAR BRIER SCALED_B;
	OUTPUT OUT = CONFINTR_ PCTLPTS=2.5 97.5 PCTLPRE= BRIER_ SCALEDB_ PCTLNAME=LOWER95 UPPER95;
RUN;

DATA CONFINTR1_;
	SET CONFINTR_;
	IND=1;
RUN;

DATA BRIERAX2_;
 	RETAIN BRIER BRIER_LOWER95 BRIER_UPPER95 SCALED_B SCALEDB_LOWER95 SCALEDB_UPPER95;
	MERGE SCALEDB CONFINTR1_;
	BY IND;
	DROP IND SWEIGHT SBS;
	TITLE 'External Brier score and Scaled Brier score with 95% CI';
	PROC PRINT;
RUN;





**** Bootstrap Scaled Brier performance;


PROC PHREG DATA=ROTTAPP NOPRINT;
	BY REPLICATE;
	MODEL SURVTIME*STATUS(0)= / TIES=EFRON RL;
	BASELINE COVARIATES=ROTT_BA OUT=OUTKM_NULL TIMELIST=5 SURVIVAL=FIVEYRSURV_null;
RUN;

*set the 5 year survival probability in the apparent bootstrap dataset;
DATA OUTKM_NULL1;
	SET OUTKM_NULL;
	KEEP REPLICATE PID FIVEYRSURV_null;
	PROC SORT; BY REPLICATE PID;
RUN;

proc sort data= ROTT_BS3A;
	BY REPLICATE PID;
RUN;

DATA ROTT_BS4A;
	MERGE ROTT_BS3A OUTKM_NULL1;
	BY REPLICATE PID;
RUN;

DATA ROTT_BS5A;
	SET ROTT_BS4A;
	/*RETAIN _SURVIVAL;
	IF NOT MISSING(SURVIVAL) THEN _SURVIVAL=SURVIVAL;
	ELSE SURVIVAL=_SURVIVAL;
	IF TIME1=0 THEN DELETE;*/
	IF CAT=1 THEN CONTRIB_NULL=(-FIVEYRSURV_null)**2;
	IF CAT=2 THEN CONTRIB_NULL=(1-FIVEYRSURV_null)**2;
	IF CAT=3 THEN CONTRIB_NULL=0;
	BS_NULL=CONTRIB_NULL*WEIGHT;
	DROP FIVEYRSURV_null;
RUN;


*ESTIMATE BRIER SCORE;
PROC UNIVARIATE DATA=ROTT_BS5A NOPRINT;
	BY REPLICATE;
	VAR BS_NULL WEIGHT;
	OUTPUT OUT=SUMNULLA SUM=AppSBS_NULL SWEIGHT;
	*PROC PRINT; 
RUN;

DATA SUMNULLA;
	RETAIN SWEIGHT AppSBS_NULL;
	SET SUMNULLA;
	SWEIGHT=LEFT(SWEIGHT);
RUN;


DATA SCALEDBA;
	MERGE SUMSA SUMNULLA;
	BY REPLICATE;
	NULL_AppBRIER = (1/SWEIGHT)*AppSBS_NULL;
	AppScaledB = 1-(AppBRIER/NULL_AppBRIER);
	*TITLE 'Brier score and IPA';
	*PROC PRINT;
RUN;



*external;
PROC PHREG DATA=ROTTAPP;
	BY REPLICATE;
	MODEL SURVTIME*STATUS(0)= / TIES=EFRON RL;
	BASELINE COVARIATES=ROTT_B_EX OUT=OUTKM_NULLEX TIMELIST=5 SURVIVAL=FIVEYRSURV_null;
RUN;

DATA OUTKM_NULLEX1;
	SET OUTKM_NULLEX;
	KEEP REPLICATE PID FIVEYRSURV_null;
	PROC SORT; BY REPLICATE PID;
RUN;

PROC SORT DATA=ROTTVAL_BS3;
	BY REPLICATE PID;
RUN;

DATA ROTTVAL_BS4;
	MERGE ROTTVAL_BS3 OUTKM_NULLEX1;
	BY REPLICATE PID;
RUN;

DATA ROTTVAL_BS5;
	SET ROTTVAL_BS4;
	/*RETAIN _SURVIVAL;
	IF NOT MISSING(SURVIVAL) THEN _SURVIVAL=SURVIVAL;
	ELSE SURVIVAL=_SURVIVAL;
	IF TIME1=0 THEN DELETE;*/
	IF CAT=1 THEN CONTRIB_NULL=(-FIVEYRSURV_null)**2;
	IF CAT=2 THEN CONTRIB_NULL=(1-FIVEYRSURV_null)**2;
	IF CAT=3 THEN CONTRIB_NULL=0;
	BS_NULL=CONTRIB_NULL*WEIGHT;
	DROP FIVEYRSURV_null;
RUN;


*ESTIMATE BRIER SCORE;
PROC SORT DATA=ROTTVAL_BS5;
	BY REPLICATE PID;
RUN;

PROC UNIVARIATE DATA=ROTTVAL_BS5 NOPRINT;
	BY REPLICATE;
	VAR BS_NULL WEIGHT;
	OUTPUT OUT=SUMNULLEX SUM=BootSBS_NULL SWEIGHT;
	PROC PRINT; 
RUN;

DATA ScaledB2;
	MERGE SUMSEX SUMNULLEX;
	BY REPLICATE;
	NULL_BootBRIER = (1/SWEIGHT)*BootSBS_NULL;
	BootScaledB = 1-(BootBRIER/NULL_BootBRIER);
RUN;

DATA FINALIPA;
	MERGE ScaledBA ScaledB2;
	BY REPLICATE;
	OptBrier = AppBrier - BootBrier;
	OptScaledB = AppScaledB - BootScaledB;
TITLE 'Optimism Brier score and Scaled Brier score';
	PROC PRINT;
RUN;




*Discrimination;
/*time-dep c: Heagerty and others - https://support.sas.com/resources/papers/proceedings17/SAS0462-2017.pdf
Methods of Estimating Time-Dependent ROC Curves in PROC PHREG
Option Method Reference
IPCW Inverse probability of censoring weighting Uno et al. (2007)
KM Conditional Kaplan-Meier Heagerty, Lumley, and Pepe (2000)
NNE Nearest neighbors Heagerty, Lumley, and Pepe (2000)
RECURSIVE Recursive method Chambless and Diao (2006)*/;

***** GLOBAL ASSESSMENT OF DISCRIMINATION - C-STATISTIC;

*apparent discrimination;
*ODS LISTING SGE=ON STYLE=PRINTER IMAGE_DPI=300 GPATH='C:\Users\sme544\Documents\STRATOS\Figures';
*ODS GRAPHICS ON / RESET=ALL NOBORDER OUTPUTFMT=TIFF /*WIDTH=4IN*/ IMAGENAME="Concord_ROTT_App2" ANTIALIAS=OFF/*ANTIALIASMAX=*/;

*THESE CALCULATE C WHICH WE DO NOT PRESENT IN PAPER;
*CALCULATE HARRELL C BUT IF REPLACE WITH UNO YOU GET UNO'S C - NEED TAU TO EQUAL EVENT TIME IF FIXED OTHERWISE USES ALL;
*Plots here provide ROC curves at 5 and 10 years survival;
PROC PHREG DATA=ROTT CONCORDANCE=/*HARRELL(SE)*/ UNO(SE SEED=8754 ITER=50) TAU=5 /*PLOTS(OVERLAY=INDIVIDUAL)=ROC ROCOPTIONS(AT=5)*/;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE / TIES=EFRON RL;
RUN;

*ODS GRAPHICS OFF;

*GET BOOTSTRAP PERFORMANCE FOR UNO C;

*CALCULATE Harrell's C - NEED TAU TO EQUAL EVENT TIME IF FIXED OTERWISE USES ALL;
PROC PHREG DATA=OUTBOOT CONCORDANCE=HARRELL(SE) /*UNO(SE SEED=8754 ITER=50)*/ TAU=5  /*PLOTS(OVERLAY=INDIVIDUAL)=ROC ROCOPTIONS(AT=5 10)*/ ;
	*WHERE REPLICATE IN (1,2);
	BY REPLICATE;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE  / TIES=EFRON RL ;
	ODS OUTPUT CONCORDANCE=AppHar(RENAME=(Estimate=AppHar stderr=AppSEH) DROP=SOURCE);
RUN;

*calculate Uno C;
PROC PHREG DATA=OUTBOOT /*CONCORDANCE=HARRELL(SE)*/ CONCORDANCE=UNO(SE SEED=8754 ITER=50) TAU=5 /*PLOTS(OVERLAY=INDIVIDUAL)=ROC ROCOPTIONS(AT=5 10)*/ ;
	*WHERE REPLICATE IN (1,2);
	BY REPLICATE;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE  / TIES=EFRON RL ;
	ODS OUTPUT CONCORDANCE=AppUno(RENAME=(Estimate=AppUno stderr=AppSEU) DROP=SOURCE);
RUN;

*EXTERNAL BOOT;
SASFILE ROTT LOAD; /* a way of loading the dataset into RAM - speeds it up */

*THIS SHOULD REPLICATE THE DATASET 500 TIMES WITHOUT REPLACEMENT - TRICK TO GET 500 COPIES - NEEDED FOR CONCORDANCE TO WORK BELOW;
PROC SURVEYSELECT DATA=ROTT OUT=OUTROTT /* Use PROC SURVEYSELECT and provide input and output dataset names */
SEED=853794 /* can enter 0 for it to select a random seed but remember to type it in here from the output otherwise cannot replicate results */
METHOD=SRS /* Unrestricted Random Sampling - simple random sampling */
SAMPRATE=1 /* can accept proportions or percentages but we want n to be size of orginal database so =1 (or 100) */
REP=500; /* number of bootstrap samples */
RUN;

SASFILE ROTT CLOSE; /* closes frees RAM buffers when done */

ODS LISTING CLOSE; /* turns off ODS listing so no printing of all output. Better than using NOPRINT as it doesn't allow storage of data in output dataset at end */

PROC PHREG DATA=OUTROTT CONCORDANCE=HARRELL(SE) /*CONCORDANCE=UNO(SE SEED=8754 ITER=50)*/ TAU=5 /*PLOTS(OVERLAY=INDIVIDUAL)=ROC ROCOPTIONS(AT=5 10)*/;
	*WHERE REPLICATE IN (1,2);
	BY REPLICATE;	
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE  / TIES=EFRON RL NOFIT;
	ROC SOURCE=SimpModelBoot;
	*STORE Stratos.SimpModel;
	*assess functional form (siz p=0.035) and PH ok;
	*ASSESS VAR=(SIZECM) PH / RESAMPLE CRPANEL;
	*ASSESS PH / RESAMPLE CRPANEL;
	*OUTPUT OUT=ROTTX XBETA=XB;
	ODS OUTPUT CONCORDANCE=BootHar(RENAME=(Estimate=BootHar stderr=BootSEH) drop=SOURCE);
RUN;

PROC PHREG DATA=OUTROTT /*CONCORDANCE=HARRELL(SE)*/ CONCORDANCE=UNO(SE SEED=8754 ITER=50) TAU=5 /*PLOTS(OVERLAY=INDIVIDUAL)=ROC ROCOPTIONS(AT=5 10)*/;
	*WHERE REPLICATE IN (1,2);
	BY REPLICATE;	
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE  / TIES=EFRON RL NOFIT;
	ROC SOURCE=SimpModelBoot;
	*STORE Stratos.SimpModel;
	*assess functional form (siz p=0.035) and PH ok;
	*ASSESS VAR=(SIZECM) PH / RESAMPLE CRPANEL;
	*ASSESS PH / RESAMPLE CRPANEL;
	*OUTPUT OUT=ROTTX XBETA=XB;
	ODS OUTPUT CONCORDANCE=BootUno(RENAME=(Estimate=BootUno stderr=BootSEU) drop=SOURCE);
RUN;

DATA HARCON;
	MERGE AppHar BootHar;
	BY REPLICATE;
	OptHarC = AppHar - BootHar;
RUN;

DATA UNOCON;
	MERGE AppUno BootUno;
	BY REPLICATE;
	OptUnoC = AppUno - BootUno;
RUN;

***************** FIXED TIME POINT ASSESSMENT OF DISCRIMINATION - TIME-DEPENDENT AUROC AT 5 YEARS;

*time dependent AUC;
*The PLOTS=AUC option in the PROC PHREG statement plots the AUC curve. The ROCOPTIONS option in the
PROC PHREG statement enables you to specify the inverse probability of censoring weighting (IPCW) method (UNO) to
compute the ROC curves, and the CL suboption requests pointwise confidence limits for the AUC curve. The IAUC
option computes and displays the integrated AUC over time;
*IPCW OR UNO ARE SAME - (CL SEED ITER OPTIONS ONLY APPLY FOR THIS AUC);

*Apparent;
PROC PHREG DATA=ROTT TAU=5 ROCOPTIONS(AUC AT=4.95 METHOD=/*RECURSIVE*/ /*NNE*/ /*KM*/ IPCW (CL SEED=134) IAUC /*OUTAUC=IPCWAUC *//*AUCDIFF*/);
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE / TIES=EFRON RL NOFIT;
	ROC 'NPI' SOURCE=SimpModel;
RUN;

*GET BOOTSTRAP PERFORMANCE FOR UNO AUC;

PROC PHREG DATA=outboot /*PLOTS=AUC*/ TAU=5 ROCOPTIONS(AUC AT=4.95 METHOD=/*RECURSIVE*/ /*NNE*/ /*KM*/ IPCW (CL SEED=134) /*IAUC OUTAUC=IPCWAUC*/);
	*WHERE REPLICATE IN (1,2);
	BY REPLICATE;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE  / TIES=EFRON RL;
	*STORE Stratos.SimpModel;
	*assess functional form (siz p=0.035) and PH ok;
	*ASSESS VAR=(SIZECM) PH / RESAMPLE CRPANEL;
	*ASSESS PH / RESAMPLE CRPANEL;
	*OUTPUT OUT=ROTTX XBETA=XB;
	ODS OUTPUT AUC=AppTDUno(RENAME=(Estimate=AppTDUno stderr=AppTDSE) DROP=SOURCEID UPPER LOWER source) IAUC=AppIAUC;
RUN;

*External boot;
PROC PHREG DATA=outrott /*PLOTS=AUC*/ TAU=5 ROCOPTIONS(AUC AT=4.95 METHOD=/*RECURSIVE*/ /*NNE*/ /*KM*/ IPCW (CL SEED=134) /*IAUC OUTAUC=IPCWAUC*/);
	*WHERE REPLICATE IN (1,2);
	BY REPLICATE;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE  / TIES=EFRON RL NOFIT;
	ROC SOURCE=SimpModelBoot;
	*assess functional form (siz p=0.035) and PH ok;
	*ASSESS VAR=(SIZECM) PH / RESAMPLE CRPANEL;
	*ASSESS PH / RESAMPLE CRPANEL;
	*OUTPUT OUT=ROTTX XBETA=XB;
	ODS OUTPUT AUC=BootTDUno(RENAME=(Estimate=BootTDUno stderr=BootTDSE) DROP=SOURCE sourceid UPPER LOWER);
RUN;

DATA TDUNOCON;
	MERGE AppTDUno BootTDUno;
	BY REPLICATE;
	OptUnoAUC = AppTDUno - BootTDUno;
RUN;




************** Conduct DCA analysis;

*- see %stdca macro - get from https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/decision-curve-analysis;

*let's add PGR as a new marker;
DATA BER;
	SET FFPGR;
	SURVTIME1=SURVTIME;
RUN;

PROC PHREG DATA=ROTT;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE/ TIES=EFRON RL;	
	BASELINE OUT=ORIGEXT COVARIATES=BER SURVIVAL=FIVEYR timelist=5;
RUN;

PROC PHREG DATA=FFPGR;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE pgr pgr1/ TIES=EFRON RL;	
	BASELINE OUT=NEWEXT COVARIATES=ORIGEXT SURVIVAL=FIVEYR_PGR timelist=5;
RUN;

DATA NEWEXT1;
	SET NEWEXT;
	RISK_ORIG=1-FIVEYR;
	RISK_NEW=1-FIVEYR_PGR;
RUN;

%PUT _user_;
GOPTIONS RESET = ALL;

symbol1 i=join c=green;
symbol2 i=join c=red;
symbol3 i=join c=blue;
symbol4 i=join c=darkred;
symbol5 i=join c=gray;

%STDCA(data=NEWEXT1, out=survivalmult, outcome=STATUS, ttoutcome=SURVTIME1, timepoint=5, predictors=RISK_ORIG);
%STDCA(data=NEWEXT1, out=survivalmult_new, outcome=STATUS, ttoutcome=SURVTIME1, timepoint=5, predictors=RISK_NEW);

*Sort by threshold variable;
PROC SORT DATA=survivalmult OUT=kmsort;
	BY threshold;
RUN;

*Rename the variables so that we know they are the Kaplan Meier estimates;
DATA kmsort; SET kmsort(RENAME=(RISK_ORIG=kmmodel all=kmall));
	LABEL kmmodel="Original model: Pr(Recurrence) at 5 years";
	LABEL kmall="Treat All";
RUN;

*Sort by threshold variable;
PROC SORT DATA=survivalmult_new OUT=crsort;
	BY threshold;
RUN;

*Rename the variables so that we know they are the Competing Risk estimates;
DATA crsort; SET crsort(RENAME=(RISK_NEW=crmodel all=crall));
	LABEL crmodel="Original + PGR: Pr(Recurrence) at 5 years";
RUN;

*Merge Kaplan-Meier and Competing Risk data using threshold probabilities;
DATA crsort;
	MERGE kmsort crsort;
	BY threshold;
RUN;

data stratos.AppNB;
	set crsort;
run;

title;
FOOTNOTE;
*- this allows editing of the .sge file!;
ODS LISTING SGE=ON STYLE=PRINTER IMAGE_DPI=300 GPATH='T:\People\d.mclernon\STRATOS\STRATOS\Figures';
ODS GRAPHICS ON / RESET=ALL NOBORDER OUTPUTFMT=TIFF /*WIDTH=4IN*/ IMAGENAME="dca with spline PGR" ANTIALIAS=OFF/*ANTIALIASMAX=*/;

*create graph (decision curve) with treat none, treat all, Kaplan-Meier method, and Competing Risks method;
PROC SGPLOT DATA=crsort;
	YAXIS VALUES=(-0.05 to 0.45 by 0.05) LABEL="Net Benefit";
	XAXIS VALUES=(0.0 to 1 by 0.1) LABEL="Threshold Probability";
	KEYLEGEND "all" "orig" "new" "none" / DOWN=4 POSITION=BOTTOM;
	SERIES Y=kmall X=threshold / LINEATTRS=(COLOR=BLACK THICKNESS=2 PATTERN=SOLID) NAME="all" LEGENDLABEL="Treat All";
		 /*crall*threshold*/
	SERIES Y=kmmodel X=threshold / LINEATTRS=(COLOR=GREEN THICKNESS=2 PATTERN=SOLID) NAME="orig" LEGENDLABEL="Original model: Pr(Rec) at 5 years";
	SERIES Y=crmodel X=threshold / LINEATTRS=(COLOR=BLUE THICKNESS=2 PATTERN=SOLID) NAME="new" LEGENDLABEL="Original model + PGR: Pr(Rec) at 5 years";
	SERIES Y=none X=threshold / LINEATTRS=(COLOR=RED THICKNESS=2 PATTERN=SOLID)  NAME="none" LEGENDLABEL="None";
RUN;

ODS GRAPHICS OFF;




********* NOW PUT ALL TOGETHER TO GET THE OPTIMISM IN PERFORMANCE & 95% CI FOR APPARENT PERFORMANCE **************;
DATA ALL1;
	MERGE ROYSTONBOOT FINALIPA HARCON UNOCON TDUNOCON;
	BY REPLICATE;
RUN;

DATA STRATOS.ALL1;
	SET ALL1;
RUN;


*calculate the average optimism in fit;
proc univariate data= ALL1;
	var OptD OptR2D OptBrier OptScaledB OptHarC OptUnoC OptUnoAUC;
	output out=avgopt mean=Mean_OptD Mean_OptR2D Mean_OptBrier Mean_OptScaleB Mean_OptHarC Mean_OptUnoC Mean_OptUnoAUC;
run;

*-get c-stat from final model;
*enter the aucs for each model;
DATA RESULTSDATA;
	INFILE DATALINES DELIMITER=',';
	INPUT OrigD OrigR2D OrigBrier OrigScBr OrigHarC OrigUnoC OrigUnoAUC;
	DATALINES;
	 1.058, 0.211, 0.205, 0.148, 0.678, 0.677, 0.718
;
data resultsdata_new;
	set resultsdata;
	ind=1;
run;

data avgopt1;
	set avgopt;
	ind=1;
run;

*- calculate the difference between the c from final model and the average optimism in fit to get a nearly unbiased;
*- estimate of the expected value of the external predictive discrimination of the process which generated Capp; 
*- in other words InternVal is an honest estimate of internal validity penalising for overfitting;
data correctedperf;
	merge resultsdata_new avgopt1;
	by ind;
	InternValD=OrigD - Mean_OptD;
	InternValR2D=OrigR2D - Mean_OptR2D;
	InternValBrier=OrigBrier - Mean_OptBrier;
	InternValSBrier=OrigScBr - Mean_OptScaleB;
	InternValHarC=OrigHarC - Mean_OptHarC;
	InternValUnoC=OrigUnoC - Mean_OptUnoC;
	InternValUnoAUC=OrigUnoAUC - Mean_OptUnoAUC;
	drop ind;
run;

title;
title 'Estimate of internal validity penalising for overfitting';
proc print data=correctedperf;
run;

PROC EXPORT DATA=correctedperf
  	OUTFILE= "C:\Users\sme544\Documents\STRATOS\Bootstrap_Val" 
   	DBMS=XLSX REPLACE;
RUN;
