***************************************************************************************;
* Programme: 	STRATOS External validation when original model dataset available.sas  ;
* Author: 		David McLernon														   ;
* Date: 		26th Feb 2021														   ;
* Purpose: 		STRATOS case study for survival prediction model paper				   ;
*				External validation when you have access to original model 			   ;
*				development	data													   ;
* Note:			This programme is not currently automated. It is coded	  			   ;
*				based on the case study in the article. Therefore, 		 		 	   ;
*				adapting this code for your own study will require careful			   ;
*				editing according to your data							  			   ;
***************************************************************************************;	

* NB Edit folder locations to your own;

OPTIONS MPRINT MLOGIC NONUMBER NODATE SOURCE2 MAXMEMQUERY=MAX /*(fmtsearch = (formats)*/ spool; 
LIBNAME STRATOS 'C:\Users\sme544\Documents\STRATOS';

*macros required;
%include 'C:\Users\sme544\Documents\STRATOS\RCSPLINE macro.SAS';


****** DATA CODING, FUNCTIONAL FORM, AND ASSUMPTION CHECKS;

*The model we use is the Nottingham Prognostic Index;
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

*Administrative censor at 5 years since this is our prediction horizon; 
DATA FFPGR;
	SET RO;
	*administrative censoring at 5 years;
	IF _T > 5 THEN STATUS=0;
	IF SURVTIME > 5 THEN SURVTIME=5;
	* RCS for PGR;
	%RCSPLINE(PGR, 0, 41, 486);
RUN;

***** ALL DATA CHECKS, ASSUMPTION CHECKS AND PUNCTIONAL FORM DONE IN MODEL DEVELOPMENT PROGRAMME - " ";
DATA FFPGR1;
	SET FFPGR;
RUN;



********************** Import external validation dataset ***************************************;

*code up predictors;
DATA GBSG_VA(RENAME=( _D=STATUS NODES1=NODESCAT MENO=MENO_IMP SIZECAT=SIZE));
	SET STRATOS.gbsg_br_ca;
	%RCSPLINE(PGR, 0, 41, 486);
	IF SIZE LE 20 THEN SIZECAT=1;
	IF 21 <= SIZE <= 50 THEN SIZECAT=2;
	IF SIZE >50 THEN SIZECAT=3;
	nodes_imp=nodes;
	nodes1=0;
	if 1<=nodes<=3 then nodes1=1;
	if nodes>3 then nodes1=2;
	grade = grade-1;
	IF GRADE IN (0,1) THEN GRADE=0;
	IF GRADE = 2 THEN GRADE=1;
	SURVTIME=_T;
	keep ID SIZECAT NODES1 nodes_imp GRADE _T _D meno age PGR PGR1 ER SURVTIME;
RUN;

*** Descriptive statistics;
PROC FREQ DATA=GBSG_VA;; 
	TABLES SIZE NODESCAT GRADE STATUS;
RUN;

PROC UNIVARIATE DATA=GBSG_VA;
	VAR AGE PGR;
RUN;

*use reverse Kaplan-Meier to obtain median follow-up time;
PROC LIFETEST DATA=GBSG_VA METHOD=pl ATRISK;
        TIME SURVTIME*STATUS(1);
RUN;

*Administrative censor at 5 years as development cohort; 
DATA GBSG_VAL;
	SET GBSG_VA;
	*administrative censoring at 5 years;
	IF _T > 5 THEN STATUS=0;
	IF SURVTIME > 5 THEN SURVTIME=5;
	*- keep a copy of survtime for calibration step;
	SURVTIME1=SURVTIME;
RUN;

*** fit new cox model on validation dataset;
PROC PHREG DATA=GBSG_VAL ZPH(GLOBAL TRANSFORM=LOG) EV;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE / TIES=EFRON RL;
RUN;

* Code PGR as a restricted cubic spline with 3 knots;
*First calculate the 10th, 50th and 90th percentiles for knots;
PROC UNIVARIATE DATA = GBSG_VAL;
	VAR PGR;
	OUTPUT OUT=KNOTS PCTLPRE=P_PGR PCTLPTS= 10 50 90;
RUN;

PROC PRINT DATA=KNOTS; RUN;

/* here we find the following values:
P_PGR10 P_PGR50 P_PGR90 
0       32.5      312 
*/

DATA GBSG_VAL_NEW;
	SET GBSG_VAL;
	DROP PGR1;
RUN;

DATA GBSG_VAL_NEW1;
	SET GBSG_VAL_NEW;
	%RCSPLINE(PGR, 0, 32.5, 312);
RUN;	

/*- CALCULATE THE 1ST AND 3RD QUARTILES FOR INTERQUARTILE HAZARD RATIO ESTIMATION BELOW*/;
PROC UNIVARIATE DATA=GBSG_VAL_NEW1;
	VAR PGR PGR1;
	OUTPUT OUT=PERCENTILES Q3=Q3PGR Q3PGR1 Q1=Q1PGR Q1PGR1;
RUN;
PROC PRINT DATA=PERCENTILES;
RUN;
/*Q3PGR Q3PGR1    Q1PGR       Q1PGR1 
 132 	12.3310  	7 			.003523586   */

PROC PHREG DATA=GBSG_VAL_NEW1 ZPH(GLOBAL TRANSFORM=LOG) EV;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE PGR PGR1/ TIES=EFRON RL;
	* IQ HR - put as negative to get the reverse comparison i.e. 7 v 132;
	CONTRAST 'Test PGR' PGR -125 PGR1 -12.327476414 / estimate=both;
RUN;


**************************************EXTERNAL VALIDATION ******************************************;


****************************************SIMPLE MODEL FIRST (without PGR)******************************************************;

* Fit the NPI model to the development data, and calculate linear predictor of this existing model for patients in the external dataset;
ODS LISTING SGE=ON STYLE=PRINTER IMAGE_DPI=300 GPATH='C:\Users\sme544\Documents\STRATOS';
ODS GRAPHICS ON / RESET=ALL NOBORDER OUTPUTFMT=TIFF /*WIDTH=4IN*/ IMAGENAME="FORM AND PH ROTT testint" ANTIALIAS=OFF/*ANTIALIASMAX=*/;

*The ZPH command option requests diagnostics on weight Schoenfeld residuals to check proportional hazards assumption;
*EV request the Schemper-Henderson R squared;
PROC PHREG DATA=FFPGR ZPH(GLOBAL TRANSFORM=LOG) EV OUTEST=COEFF;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE / TIES=EFRON RL;
	*Store allows SAS to save the model estimates for applying later on to external data for AUROC;
	STORE SimpModel;
	*assess functional form;
	ASSESS PH / RESAMPLE CRPANEL;
	*Calculate model (developed on Rotterdam data) linear predictor to external German data;
	BASELINE COVARIATES=GBSG_VAL OUT=GBSG_VALRD XBETA=XB /*TIMELIST=5 SURVIVAL=FIVEYRSURV*/;
	*Calculate model (developed on Rotterdam data) linear predictor for patients in original Rotterdam dataset;
	OUTPUT OUT=ROTTX XBETA=XB;
RUN;

ODS GRAPHICS OFF;

*** Estimating Baseline Survival Function under PH;
/*proc phreg data=ROTT;
CLASS SIZE (REF=FIRST) NODESCAT (REF='1') GRADE (REF=FIRST);
model SURVTIME*STATUS(0)= / TIES=EFRON RL ;
run;*/

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


PROC SORT DATA=GBSG_VALRD;
	BY ID;
RUN;

PROC SORT DATA=GBSG_VAL;
	BY ID;
RUN;

DATA GBSG_VALRD1;
	SET GBSG_VALRD;
	BY ID;
	IF FIRST.ID;
	KEEP ID XB;
RUN;

*merge linear predictor to external dataset;
DATA GBSG_VALRD2;
	MERGE GBSG_VAL GBSG_VALRD1;
	BY ID;
RUN;





**************** EXTERNAL VALIDATION WHEN YOU HAVE ALL THE DEVELOPMENT DATA ****************;


****************GLOBAL ASSESSMENT FOR OVERALL PERFORMANCE;


******* AS IN PAPER, WE CALCULATE Royston and Sauerbrei’s R SQUARED D ;

/*1. To compute D, first the Cox PH model is fitted - (already done above in ROTTX dataset). 
  2. Then the prognostic index of the model, XB, is transformed to give standard normal order rank statistics (rankits - formed using Blom’s approximation). 
  3. The rankits are multiplied by a factor of SQRT(8/pi) to give Zi (i = 1, n subjects). 
  4. Finally a Cox PH model is fitted to these values; D is the coefficient of Z, say a*, from this second model. 
NOTE: Royston and Sauerbrei (2004) showed that D most accurately measures separation of survival curves when the underlying prognostic index values, XB, are normally distributed. 
The regression on the Z in the second model is then linear and a* is an approximately unbiased estimate of a. They explained that when the XB are not normally distributed, 
linearity in the second model breaks down. D still measures separation because a* in the second model still estimates a, but with bias.*/


**EXTERNAL validation - calculate D and RsqD - use external dataset with linear predictor from original model applied to it;

 *2. Then the prognostic index of the model, XB, is transformed to give standard normal order rank statistics (rankits - formed using Blom’s approximation); 
PROC UNIVARIATE DATA=GBSG_VALRD2;
	HISTOGRAM;
	VAR XB;
RUN;

PROC RANK DATA=GBSG_VALRD2 NORMAL=BLOM OUT= GBSG_VALRDD;
	VAR XB;
RUN;

PROC UNIVARIATE DATA=GBSG_VALRDD;
	HISTOGRAM;
	VAR XB;
RUN;

*  3. The rankits are multiplied by a factor of SQRT(8/pi) to give Zi (i = 1, n subjects);
DATA XD;
	SET GBSG_VALRDD;
	PIE = CONSTANT("pi");
	Z=XB/(SQRT(8/PIE));
RUN;

PROC UNIVARIATE DATA=XD;
	HISTOGRAM;
	VAR z;
RUN;

*  4. Finally a Cox PH model is fitted to these values, D is the coefficient of Z, say a*, from this second model;
*Royston's D is parameter estimate of Z, D = 0.833;
ODS SELECT NONE;
PROC PHREG DATA=XD;
	MODEL SURVTIME*STATUS(0)=Z / TIES=EFRON;
	ODS OUTPUT ParameterEstimates=PAREST;
RUN;
ODS SELECT ALL;

DATA ROYSTONEXT(RENAME=(ESTIMATE=D));
	SET PAREST;
	PIE = CONSTANT("pi");
	R2D=((ESTIMATE**2)/(8/PIE))/(((ESTIMATE**2)/(8/PIE))+((PIE**2)/6));
	VALIDATE='External';
	KEEP VALIDATE ESTIMATE R2D;
	PROC PRINT; TITLE "Royston's D and RsqD";
RUN;






*CALCULATE Brier score - see Graf et al 1999;
/*Weights - Kaplan Meier with censor as event;
3 groups - 1. Those who have the event up to fixed event time of interest, and 2. those who go beyond fixed t (could be event or event free), 3. those censored up to fixed t (only first 2 groups contribute to BS but all to weights);
Then group 1 calc -surv^2, and group 2 calc (1-surv)^2 where surv is probability of surv at t*; 
Check weight for group 2 (G(t*)) - should be same for all patients in that group;*/
title ' ';

*External validation - Calculate Brier score;

*calculate weights for external validation;
PROC LIFETEST DATA=GBSG_VAL METHOD=pl ATRISK OUTSURV=OUTKM_EXT /*NOPRINT*/;
        TIME SURVTIME*STATUS(1);
RUN;

*CODE THE 3 GROUPINGS AND DUPLICATE SURVIVAL TIME AS SAS WILL REMOVE THE OFFICIAL SURVIVAL TIME VARIABLE IN BASELINE STATEMENT;
*in external;
DATA ROTT_B_EX;
	SET GBSG_VAL;
	IF SURVTIME<=4.95 AND STATUS=1 THEN CAT=1;
	IF SURVTIME>4.95 THEN CAT=2;
	IF SURVTIME<=4.95 AND STATUS=0 THEN CAT=3;
	TIME=SURVTIME;
RUN;

*NOW ESTIMATE SURVIVAL AT 5 YEARS IN VALIDATION DATASET;
PROC PHREG DATA=ROTTX;
	CLASS SIZE (REF=FIRST) NODESCAT (REF='1') GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE / TIES=EFRON RL;
	*EXTERNAL;
	BASELINE COVARIATES=ROTT_B_EX OUT=ROTTVAL_BS TIMELIST=4.95 SURVIVAL=FIVEYRSURV XBETA=XB;
RUN;

*MERGE THE KAPLAN-MEIER WEIGHTS TO THE APPROPRIATE TIMES;
DATA ROTTVAL_BS1(RENAME=(TIME=SURVTIME));
	SET ROTTVAL_BS;
	TIME1=TIME;
	IF TIME1>4.95 THEN TIME1=4.95;
	DROP SURVTIME;
	PROC SORT; BY TIME1;
RUN;

DATA OUTKM_EXT1(RENAME=(SURVTIME=TIME1));
	SET OUTKM_EXT;
	*IF _CENSOR_=1 THEN DELETE;
	IF SURVTIME>4.95 THEN DELETE;
	WEIGHT=1/SURVIVAL;
	*KEEP SURVTIME WEIGHT;
RUN;

DATA OUTKM_EXT2;
	SET OUTKM_EXT1;
	S=LAG(TIME1);
	IF TIME1=S THEN DELETE;
	DROP S;
RUN;

DATA ROTTVAL_BS2;
	MERGE ROTTVAL_BS1 OUTKM_EXT1;
	BY TIME1;
	IF ID=. THEN DELETE;
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
PROC UNIVARIATE DATA=ROTTVAL_BS3 NOPRINT;
	VAR BS WEIGHT;
	OUTPUT OUT=SUMSEX SUM=SBS SWEIGHT;
	PROC PRINT; 
RUN;

DATA SUMSEX;
	RETAIN SWEIGHT SBS BRIER;
	SET SUMSEX;
	BRIER = (1/SWEIGHT)*SBS;
	SWEIGHT=LEFT(SWEIGHT);
	TITLE 'External Brier score';
	PROC PRINT;
RUN;



****** Scaled Brier FOR External validation;
/*Scaled Brier = 1 - (model Brier score/null model Brier score), where null model Brier score is null cox; 
100% is perfect, <0 is useless, higher better, hamrful models <0;*/

*NOW ESTIMATE SURVIVAL AT 5 YEARS for NULL model;

PROC PHREG DATA=GBSG_VAL;
	MODEL SURVTIME*STATUS(0)= / TIES=EFRON RL;
	*APPARENT;
	BASELINE COVARIATES=ROTT_B_EX OUT=ROTTVAL_BSnull TIMELIST=5 SURVIVAL=FIVEYRSURV_null;
RUN;

*Brier for null model;
DATA ROTTVAL_BS1null;
	SET ROTTVAL_BSnull;
	KEEP ID FIVEYRSURV_null;
	PROC SORT; BY ID;
RUN;

PROC SORT DATA=ROTTVAL_BS3;
	BY ID;
RUN;

*MERGE THE NULL SURVIVAL PROBS TO THE BRIER SCORE DATASET FROM EARLIER;
DATA ROTTVAL_BS4;
	MERGE ROTTVAL_BS3 ROTTVAL_BS1null;
	BY ID;
RUN;

DATA ROTTVAL_BS5;
	SET ROTTVAL_BS4;
	IF CAT=1 THEN CONTRIB_NULL=(-FIVEYRSURV_null)**2;
	IF CAT=2 THEN CONTRIB_NULL=(1-FIVEYRSURV_null)**2;
	IF CAT=3 THEN CONTRIB_NULL=0;
	BS_NULL=CONTRIB_NULL*WEIGHT;
	DROP FIVEYRSURV_null;
RUN;

*ESTIMATE BRIER SCORE FOR NULL MODEL;
PROC UNIVARIATE DATA=ROTTVAL_BS5 NOPRINT;
	VAR BS_NULL WEIGHT;
	OUTPUT OUT=SUMNULLEX SUM=SBS_NULL SWEIGHT;
	PROC PRINT; 
RUN;

DATA SUMNULLEX;
	RETAIN SWEIGHT SBS_NULL;
	SET SUMNULLEX;
	SWEIGHT=LEFT(SWEIGHT);
RUN;

*CALCULATE Scaled Brier;
DATA SCALEDB;
	MERGE SUMSEX SUMNULLEX;
	BY SWEIGHT;
	NULL_BRIER = (1/SWEIGHT)*SBS_NULL;
	SCALED_B = 1-(BRIER/NULL_BRIER);
	TITLE 'External Brier score and Scaled Brier';
	PROC PRINT;
RUN;




********************************************************************* ABOVE REPEATED FOR ADDED MARKER ***************************************************;

/*- CALCULATE THE 1ST AND 3RD QUARTILES FOR INTERQUARTILE HAZARD RATIO ESTIMATION BELOW*/;
PROC UNIVARIATE DATA=FFPGR;
	VAR PGR PGR1;
	OUTPUT OUT=PERCENTILES Q3=Q3PGR Q3PGR1 Q1=Q1PGR Q1PGR1;
RUN;
PROC PRINT DATA=PERCENTILES;
RUN;
/*Q3PGR Q3PGR1    Q1PGR   Q1PGR1 
 198   14.9704   4       0.000270961 */


ODS LISTING SGE=ON STYLE=PRINTER IMAGE_DPI=300 GPATH='T:\People\d.mclernon\STRATOS\STRATOS\Figures';
ODS GRAPHICS ON / RESET=ALL NOBORDER OUTPUTFMT=TIFF /*WIDTH=4IN*/ IMAGENAME="PH ROTT PGR as spline" ANTIALIAS=OFF/*ANTIALIASMAX=*/;

PROC PHREG DATA=FFPGR ZPH(GLOBAL TRANSFORM=LOG) EV OUTEST=COEFF_PGR;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE PGR PGR1/ TIES=EFRON RL;
	* IQ HR - put as negative to get the reverse comparison i.e. 4 v 198;
	CONTRAST 'Test PGR' PGR -194 PGR1 -14.970129039 / estimate=both;
	STORE SimpModelPGR;
	*assess functional form;
	ASSESS PH / RESAMPLE CRPANEL;
	BASELINE COVARIATES=GBSG_VAL OUT=GBSG_VALRDa XBETA=XB /*TIMELIST=5 SURVIVAL=FIVEYRSURV*/;
	OUTPUT OUT=ROTTXa XBETA=XB;
RUN;

ODS GRAPHICS OFF;

*** Estimating Baseline Survival Function under PH;

*use cards;
data inrisks;
	input SIZE NODESCAT GRADE PGR PGR1;
	cards;
1 1 0 0 0
;

proc phreg data=FFPGR;
CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
model SURVTIME*STATUS(0)= SIZE NODESCAT GRADE PGR PGR1 / TIES=EFRON RL ;
baseline covariates=inrisks out=outph1 survival=ps timepoint=(1,2,3,4,5) / method=breslow;
run;

proc print data=outph1;
run;
/*Obs SIZE NODESCAT GRADE PGR PGR1 SURVTIME ps 
1 1 0 0 0 0 1 0.96277 
2 1 0 0 0 0 2 0.90116 
3 1 0 0 0 0 3 0.85028 
4 1 0 0 0 0 4 0.81161 
5 1 0 0 0 0 5 0.77815 

 */


DATA STRATOS.COEFF_PGR;
	SET COEFF_PGR;
	DROP _TYPE_ _TIES_ _STATUS_ _NAME_ _LNLIKE_;
RUN;


PROC SORT DATA=GBSG_VALRDa;
	BY ID;
RUN;

PROC SORT DATA=GBSG_VAL;
	BY ID;
RUN;

DATA GBSG_VALRD1a;
	SET GBSG_VALRDa;
	BY ID;
	IF FIRST.ID;
	KEEP ID XB;
RUN;

DATA GBSG_VALRD2a;
	MERGE GBSG_VAL GBSG_VALRD1a;
	BY ID;
RUN;


***  check model misspecification;
PROC PHREG DATA=GBSG_VALRD2a;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE PGR PGR1/ OFFSET=XB TIES=EFRON RL;
RUN;

* PH test;
PROC PHREG DATA=GBSG_VALRD2a ZPH(GLOBAL TRANSFORM=LOG);
	MODEL SURVTIME*STATUS(0)=XB / TIES=EFRON RL;
	ASSESS PH / RESAMPLE CRPANEL;
RUN;




******************APPENDIX 7 - WHAT TO DO IF DEVELOPMENT DATASET OR BASELINE HAZARD ARE NOT AVAILABLE ;

****** ASSUME ORIGINAL MODEL DEVELOPMENT PAPER HAS KAPLAN-MEIER PLOTS FOR RISK GROUPS;
****** Define Risk Groups for K-M - CREATES PLOT IN APPENDIX 7, FIG S7;



*** Define Risk Groups for K-M;
* Mean of XB in dev data is 0.272779;
PROC UNIVARIATE DATA=ROTTXa;
	VAR XB;
RUN;
*Centre PI on average risk by subtracting mean;
DATA ROTTRISKa;
	SET ROTTXa;
	XBC = XB - 0.272779;
run;

PROC UNIVARIATE DATA=ROTTRISKa;
	VAR XBc;
RUN;

* generate risk groups;
%macro generate_percentiles(ptile1,ptile2,ptile3,mean); 
/* Output desired percentile values */                         
proc univariate data=rottxa;                                               
   var xb;                                                       
   output out=test1a mean=&mean pctlpts= &ptile1 &ptile2 &ptile3 pctlpre=xb_;                               
run;                                                                 
 
data test1a;
	set test1a;
	z=1;
run;

data rottxa_;
	set rottxa;
	z=1;
run;

data GBSG_VALRD2a_;
	set GBSG_VALRD2a;
	z=1;
run;

data rottxa_1;                                                             
   merge rottxa_ test1a;        
	by z; 
	xbc=xb-&mean;
	riskgp=0;
    if (xb_&ptile1 < xb <= xb_&ptile2) then riskgp=1;
    if (xb_&ptile2 < xb <= xb_&ptile3) then riskgp=2;
	if (xb > xb_&ptile3) then riskgp=3;
	drop z;
run;  

data GBSG_VALRD2a_1;                                                             
   merge GBSG_VALRD2a_ test1a;        
	by z; 
	xbc=xb-&mean;
	riskgp=0;
    if (xb_&ptile1 < xb <= xb_&ptile2) then riskgp=1;
    if (xb_&ptile2 < xb <= xb_&ptile3) then riskgp=2;
	if (xb > xb_&ptile3) then riskgp=3;
	drop z;
run; 
  
%mend;
 
options mprint mlogic symbolgen;
%generate_percentiles(25,50,75,xbmean);

*********** - plot kaplan-meier of risk groups in development dataset and validation dataset;
*mean (SD) of xb centered is 0 (0.66083);
PROC UNIVARIATE DATA=rottxa_1;
	VAR XBC;
RUN;

PROC LIFETEST DATA=rottxa_1 METHOD=pl PLOTS=(S, LS, LLS) OUTSURV=OUTKMa;
        TIME SURVTIME*STATUS(0);
        STRATA RISKGP/TEST=(LOGRANK);	
		*ods output Quartiles=q;
RUN;

DATA OUTKM1a; 
	SET OUTKMa; 
	*LTIME=LOG(TIMELB); 
	*LLS=LOG((-1)*(LOG(SURVIVAL)));
	IF _CENSOR_=1 THEN DELETE; 
	NUM+1;
	KEEP RISKGP SURVIVAL SURVTIME NUM;
RUN;

*- plot kaplan-meier in validation dataset;
*mean (SD) of xb centered is 0.23287906 (0.52679941);
PROC UNIVARIATE DATA=GBSG_VALRD2a_1;
	VAR XBC;
RUN;

*what is the percentage in each risk group - highest in intermediate high risk, lowest in high risk;
/*
riskgp Frequency Percent Cumulative
Frequency Cumulative
Percent 
0 60 8.75 60 8.75 
1 147 21.43 207 30.17 
2 253 36.88 460 67.06 
3 226 32.94 686 100.00 


*/
PROC FREQ DATA=GBSG_VALRD2a_1;
	TABLE RISKGP;
RUN;

PROC LIFETEST DATA=GBSG_VALRD2a_1 METHOD=pl PLOTS=(S, LS, LLS) OUTSURV=OUTKMva;
        TIME SURVTIME*STATUS(0);
        STRATA RISKGP/TEST=(LOGRANK);	
		*ods output Quartiles=q;
RUN;

DATA OUTKMv1a(RENAME=(RISKGP=RISKGPV SURVIVAL=SURVIVALV SURVTIME=SURVTIMEV)); 
	SET OUTKMva; 
	*LTIME=LOG(TIMELB); 
	*LLS=LOG((-1)*(LOG(SURVIVAL)));
	IF _CENSOR_=1 THEN DELETE; 
	NUM+1;
	KEEP RISKGP SURVIVAL SURVTIME NUM;
RUN;

PROC FORMAT;
	VALUE DEV 1='Development dataset';
	VALUE VAL 1='Validation dataset';
RUN;

DATA OVERLAYKMa;
	MERGE OUTKM1a OUTKMV1a;
	BY NUM;
	FORMAT GP1 DEV.;
	FORMAT GPV1 VAL.;
	*TRICK TO FIX LEGEND IN PLOT;
	IF RISKGP=0 THEN GP1=1;
	IF RISKGP=1 THEN GP2=1;
	IF RISKGP=2 THEN GP3=1;
	IF RISKGP=3 THEN GP4=1;

	IF RISKGPV=0 THEN GPV1=1;
	IF RISKGPV=1 THEN GPV2=1;
	IF RISKGPV=2 THEN GPV3=1;
	IF RISKGPV=3 THEN GPV4=1;
RUN;

*CREATE KAPLAN MEIER PLOT OF RISK GROUPS FOR DEVELOPMENT AND VALIDATION DATASETS;
title;
FOOTNOTE;
*- this allows editing of the .sge file!;
ODS LISTING SGE=ON STYLE=PRINTER IMAGE_DPI=300 GPATH='T:\People\d.mclernon\STRATOS\STRATOS\Figures';
ODS GRAPHICS ON / RESET=ALL NOBORDER OUTPUTFMT=TIFF /*WIDTH=4IN*/ IMAGENAME="Overlay KM PGR"  ANTIALIAS=OFF/*ANTIALIASMAX=*/;

PROC SGPLOT DATA=OVERLAYKMa;
	YAXIS VALUES=(0 to 1 by 0.2) LABEL="Event-free survival probability";
	XAXIS VALUES=(0 to 5 by 1) LABEL="Years";
	KEYLEGEND "all" "orig" / ACROSS=1 LOCATION=INSIDE POSITION=TOPRIGHT;
	STEP Y=SURVIVAL X=SURVTIME / LINEATTRS=(COLOR=BLACK THICKNESS=2 PATTERN=SOLID) NAME="all" GROUP=GP1 NOMISSINGGROUP;
	STEP Y=SURVIVAL X=SURVTIME / LINEATTRS=(COLOR=BLACK THICKNESS=2 PATTERN=SOLID) GROUP=GP2 NOMISSINGGROUP;
	STEP Y=SURVIVAL X=SURVTIME / LINEATTRS=(COLOR=BLACK THICKNESS=2 PATTERN=SOLID) GROUP=GP3 NOMISSINGGROUP;
	STEP Y=SURVIVAL X=SURVTIME / LINEATTRS=(COLOR=BLACK THICKNESS=2 PATTERN=SOLID) GROUP=GP4 NOMISSINGGROUP;
	STEP Y=SURVIVALv X=SURVTIMEv / LINEATTRS=(COLOR=BLACK THICKNESS=2 PATTERN=DASH) NAME="orig" GROUP=GPV1 NOMISSINGGROUP;
	STEP Y=SURVIVALv X=SURVTIMEv / LINEATTRS=(COLOR=BLACK THICKNESS=2 PATTERN=DASH) GROUP=GPV2 NOMISSINGGROUP;
	STEP Y=SURVIVALv X=SURVTIMEv / LINEATTRS=(COLOR=BLACK THICKNESS=2 PATTERN=DASH) GROUP=GPV3 NOMISSINGGROUP;
	STEP Y=SURVIVALv X=SURVTIMEv / LINEATTRS=(COLOR=BLACK THICKNESS=2 PATTERN=DASH) GROUP=GPV4 NOMISSINGGROUP;
RUN;

ODS GRAPHICS OFF;




**************** EXTERNAL VALIDATION WHEN YOU HAVE ALL THE DEVELOPMENT DATA FOR PGR MARKER ADDED ****************;


****************GLOBAL ASSESSMENT FOR OVERALL PERFORMANCE;


******* AS IN PAPER, WE CALCULATE Royston and Sauerbrei’s R SQUARED D ;

*Royston and Sauerbrei’s D;
/*1. To compute D, first the Cox PH model is fitted. 
  2. Then the prognostic index of the model, XB, is transformed to give standard normal order rank statistics (rankits - formed using Blom’s approximation). 
  3. The rankits are multiplied by a factor of SQRT(8/pi) to give Zi (i = 1, n subjects). 
  4. Finally a Cox PH model is fitted to these values; D is the coefficient of Z, say a*, from this second model. 
NOTE: Royston and Sauerbrei (2004) showed that D most accurately measures separation of survival curves when the underlying prognostic index values, XB, are normally distributed. 
The regression on the Z in the second model is then linear and a* is an approximately unbiased estimate of a. They explained that when the XB are not normally distributed, 
linearity in the second model breaks down. D still measures separation because a* in the second model still estimates a, but with bias.*/


**EXTERNAL;
PROC UNIVARIATE DATA=GBSG_VALRD2a;
	HISTOGRAM;
	VAR XB;
RUN;

PROC RANK DATA=GBSG_VALRD2a NORMAL=BLOM OUT= GBSG_VALRDDa;
	VAR XB;
RUN;

PROC UNIVARIATE DATA=GBSG_VALRDDa;
	HISTOGRAM;
	VAR XB;
RUN;

DATA XDa;
	SET GBSG_VALRDDa;
	PI = CONSTANT("pi");
	Z=XB/(SQRT(8/PI));
RUN;

PROC UNIVARIATE DATA=XDa;
	HISTOGRAM;
	VAR z;
RUN;

*Royston's D is parameter estimate of Z, D = 0.95979;
ODS SELECT NONE;
PROC PHREG DATA=XDa;
	MODEL SURVTIME*STATUS(0)=Z / TIES=EFRON;
	ODS OUTPUT ParameterEstimates=PAREST;
RUN;
ODS SELECT ALL;

DATA ROYSTONEXTaa(RENAME=(ESTIMATE=D));
	SET PAREST;
	PI = CONSTANT("pi");
	R2D=((ESTIMATE**2)/(8/PI))/(((ESTIMATE**2)/(8/PI))+((PI**2)/6));
	VALIDATE='External';
	KEEP VALIDATE ESTIMATE R2D;
	PROC PRINT; TITLE "Royston's D and RsqD with PGR";
RUN;

 *Royston and Sauerbrei’s R2D for Cox PH model;
*R2D = D2/(8/pi) / ((D2/(8/pi)) + (pi2/6));



*Brier - see Graf et al 1999;
/*Weights - Kaplan Meier with censor as event;
3 groups - Those who have the event up to fixed event time of interest, and those who go beyond fixed t (could be event or event free), also those censored up to;
fixed t (only first 2 groups contribute to BS but all to weights;
Caclculate KM altogether;
Then group 1 calc -surv^2, and group 2 calc (1-surv)^2 where surv is probability of surv at t*; 
Check weight for group 2 (G(t*)) should be came for all patients in that group;*/
title ' ';




*external validation brier;

*for external validation;
PROC LIFETEST DATA=GBSG_VAL METHOD=pl ATRISK /*PLOTS=(S, LS, LLS)*/ OUTSURV=OUTKM_EXTa NOPRINT;
        TIME SURVTIME*STATUS(1);
RUN;

DATA ROTT_B_EXa;
	SET GBSG_VAL;
	IF SURVTIME<=4.95 AND STATUS=1 THEN CAT=1;
	IF SURVTIME>4.95 THEN CAT=2;
	IF SURVTIME<=4.95 AND STATUS=0 THEN CAT=3;
	TIME=SURVTIME;
RUN;

PROC PHREG DATA=ROTTXa;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE PGR PGR1/ TIES=EFRON RL;
	*EXTERNAL;
	BASELINE COVARIATES=ROTT_B_EXa OUT=ROTTVAL_BSa TIMELIST=5 SURVIVAL=FIVEYRSURV;
RUN;


*MERGE THE KAPLAN-MEIER WEIGHTS TO THE APPROPRIATE TIMES;

DATA ROTTVAL_BS1a(RENAME=(TIME=SURVTIME));
	SET ROTTVAL_BSa;
	TIME1=TIME;
	IF TIME1>5 THEN TIME1=5;
	DROP SURVTIME;
	PROC SORT; BY TIME1;
RUN;

DATA OUTKM_EXT1a(RENAME=(SURVTIME=TIME1));
	SET OUTKM_EXTa;
	*IF _CENSOR_=1 THEN DELETE;
	IF SURVTIME>5 THEN DELETE;
	WEIGHT=1/SURVIVAL;
	KEEP SURVTIME WEIGHT;
RUN;

DATA OUTKM_EXT2a;
	SET OUTKM_EXT1a;
	S=LAG(TIME1);
	IF TIME1=S THEN DELETE;
	DROP S;
RUN;

DATA ROTTVAL_BS2a;
	MERGE ROTTVAL_BS1a OUTKM_EXT1a;
	BY TIME1;
	IF ID=. THEN DELETE;
RUN;

DATA ROTTVAL_BS3a;
	SET ROTTVAL_BS2a;
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
PROC UNIVARIATE DATA=ROTTVAL_BS3a NOPRINT;
	VAR BS WEIGHT;
	OUTPUT OUT=SUMSEXa SUM=SBS SWEIGHT;
	PROC PRINT; 
RUN;

DATA SUMSEXa;
	SET SUMSEXa;
	BRIER = (1/SWEIGHT)*SBS;
	SWEIGHT=LEFT(SWEIGHT);
	TITLE 'External Brier score for PGR';
	PROC PRINT;
RUN;

PROC FREQ DATA=ROTTVAL_BS3a;
	TABLE CAT;
RUN;




****** Scaled Brier FOR External validation;
/*Scaled Brier = 1 - (model Brier score/null model Brier score), where null model Brier score is null cox; 
100% is perfect, <0 is useless, higher better, hamrful models <0;*/

PROC SORT DATA=ROTTVAL_BS3a;
	BY ID;
RUN;

DATA ROTTVAL_BS4a;
	MERGE ROTTVAL_BS3a ROTTVAL_BS1null;
	BY ID;
RUN;

DATA ROTTVAL_BS5a;
	SET ROTTVAL_BS4a;
	*RETAIN _SURVIVAL;
	*IF NOT MISSING(SURVIVAL) THEN _SURVIVAL=SURVIVAL;
	*ELSE SURVIVAL=_SURVIVAL;
	*IF TIME1=0 THEN DELETE;
	IF CAT=1 THEN CONTRIB_NULL=(-FIVEYRSURV_null)**2;
	IF CAT=2 THEN CONTRIB_NULL=(1-FIVEYRSURV_null)**2;
	IF CAT=3 THEN CONTRIB_NULL=0;
	BS_NULL=CONTRIB_NULL*WEIGHT;
	DROP FIVEYRSURV_null;
RUN;



*ESTIMATE BRIER SCORE;
PROC UNIVARIATE DATA=ROTTVAL_BS5a NOPRINT;
	VAR BS_NULL WEIGHT;
	OUTPUT OUT=SUMNULLEXa SUM=SBS_NULL SWEIGHT;
	PROC PRINT; 
RUN;

DATA SUMNULLEXa;
	SET SUMNULLEXa;
	SWEIGHT=LEFT(SWEIGHT);
RUN;

DATA ScaledBrierA;
	MERGE SUMSEXa SUMNULLEXa;
	BY SWEIGHT;
	NULL_BRIER = (1/SWEIGHT)*SBS_NULL;
	SCALEDB = 1-(BRIER/NULL_BRIER);
	TITLE 'External Brier score and Scaled Brier with PGR';
	PROC PRINT;
RUN;



********* END OF OVERALL PERFORMANCE FOR ADDED MARKER;







**********************DISCRIMINATION**************************;

/*time-dep c: Heagerty and others - https://support.sas.com/resources/papers/proceedings17/SAS0462-2017.pdf
Methods of Estimating Time-Dependent ROC Curves in PROC PHREG
Option Name Reference:
IPCW Inverse probability of censoring weighting Uno et al. (2007)
KM Conditional Kaplan-Meier Heagerty, Lumley, and Pepe (2000)
NNE Nearest neighbors Heagerty, Lumley, and Pepe (2000)
RECURSIVE Recursive method Chambless and Diao (2006)*/;

***** GLOBAL ASSESSMENT OF DISCRIMINATION - C-STATISTIC;


*apparent discrimination;

*External validation - discrimination;
ODS LISTING SGE=ON STYLE=PRINTER IMAGE_DPI=300 GPATH='C:\Users\sme544\Documents\STRATOS\Figures';
ODS GRAPHICS ON / RESET=ALL NOBORDER OUTPUTFMT=TIFF /*WIDTH=4IN*/ IMAGENAME="Concord_ROTT_ext2" ANTIALIAS=OFF/*ANTIALIASMAX=*/;

*CALCULATE HARRELL C BUT IF REPLACE WITH UNO YOU GET UNO'S C - NEED TAU TO EQUAL EVENT TIME IF FIXED OTERWISE USES ALL;
PROC PHREG DATA=GBSG_VAL /*CONCORDANCE=HARRELL(SE)*/ CONCORDANCE=UNO(SE SEED=8754 ITER=50 /*DIFF*/) TAU=5 PLOTS(OVERLAY=INDIVIDUAL)=ROC ROCOPTIONS(AT=5);
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE PGR PGR1/ TIES=EFRON RL;
	*the following statements call the stored models from earlier and calculate the concordance statistics for model WITH AND WITHOUT PGR;
	ROC 'NPI' SOURCE=SimpModel;
	ROC 'NPI + PGR' SOURCE=SimpModelPGR;
RUN;

ODS GRAPHICS OFF;



*TIME DEPENDENT roc CURVES;

***EXTERNAL DISCRIMINATON;

*The PLOTS=AUC option in the PROC PHREG statement plots the AUC curve. The ROCOPTIONS option in the
PROC PHREG statement enables you to specify the inverse probability of censoring weighting (IPCW) method (UNO) to
compute the ROC curves, and the CL suboption requests pointwise confidence limits for the AUC curve. The IAUC
option computes and displays the integrated AUC over time;
*IPCW OR UNO ARE SAME - (CL SEED ITER OPTIONS ONLY APPLY FOR THIS AUC;

* Heagerty TD AUC;
*NOTE: For the AT= you must specify the time at or just prior to the last event time in your dataset. If you use 5 years (or your last follow-up time this will not work);
PROC PHREG DATA=GBSG_VAL /*PLOTS=AUC*/ /*TAU=5*/ ROCOPTIONS(AUC AT=4.95 METHOD=/*RECURSIVE*/ NNE /*KM*/ /*IPCW (CL SEED=134)*/ /*IAUC*/ /*OUTAUC=IPCWAUC*/);
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE / TIES=EFRON RL NOFIT;
	ROC 'NPI' SOURCE=SimpModel;
	ROC 'NPI + PGR' SOURCE=SimpModelPGR;
RUN;

PROC PHREG DATA=GBSG_VAL /*PLOTS=AUC*/ /*TAU=5*/ ROCOPTIONS(AUC AT=4.95 METHOD=/*RECURSIVE*//*NNE*/ /*KM*/ IPCW (CL SEED=134) /*IAUC*/ /*OUTAUC=IPCWAUC*/);
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE / TIES=EFRON RL NOFIT;
	ROC 'NPI' SOURCE=SimpModel;
	ROC 'NPI + PGR' SOURCE=SimpModelPGR;
RUN;


* rerun Uno and save dataset to plot the AUC yourself;
PROC PHREG DATA=GBSG_VAL /*TAU=5*/ ROCOPTIONS(/*AUC AT=5 AUCDIFF*/ METHOD=/*RECURSIVE*/ /*NNE*/ /*KM*/ IPCW (CL SEED=134) /*IAUC*/ OUTAUC=IPCWAUC2);
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE PGR PGR1/ TIES=EFRON RL NOFIT;
	ROC 'NPI + PGR' SOURCE=SimpModelPGR;
	ROC 'NPI' SOURCE=SimpModel;
RUN;

*plot time-dependent AUC using outauc dataset;
DATA IPCWAUC3;
	SET IPCWAUC2;
	IF SURVTIME >5 THEN DELETE;
RUN;

title;
FOOTNOTE;
*- this allows editing of the .sge file!;
ODS LISTING SGE=ON STYLE=PRINTER IMAGE_DPI=300 GPATH='C:\Users\sme544\Documents\STRATOS\Figures';
ODS GRAPHICS ON / RESET=ALL NOBORDER OUTPUTFMT=TIFF /*WIDTH=4IN*/ IMAGENAME="Time-dep ROC ROTT with PGR EXT uno" ANTIALIAS=OFF/*ANTIALIASMAX=*/;

*PLOT THE TIME-DEPENDENT AUC UP TO 5 YEARS;
PROC SGPLOT DATA=IPCWAUC3 NOAUTOLEGEND ;
XAXIS       	LABEL="Survival time (years)" VALUES=(0 TO 5 BY 1);
YAXIS        	LABEL="AUC" VALUES=(0 TO 1 BY 0.20);
  BAND X=SURVTIME LOWER=_LowerAUC_ UPPER=_UpperAUC_ / transparency=0.5 GROUP=_SOURCE_  /*LINEATTRS=(COLOR=BLACK PATTERN=MEDIUMDASH THICKNESS=3)*/ /*NOEXTEND NOOUTLINE*/;
  SERIES X=SURVTIME Y=_AUC_ / LINEATTRS=(THICKNESS=1) GROUP=_SOURCE_ GROUPLC=_ID_ NAME="Models" LEGENDLABEL="Models";
  REFLINE 0.5 / AXIS=Y LINEATTRS=(COLOR=BLACK THICKNESS=1);
  KEYLEGEND "Models" / BORDER LOCATION=INSIDE POSITION=BOTTOMRIGHT DOWN=2 TITLE="Models";
run;

ods graphics off;

************* END OF DISCRIMINATION CODE ***************************;





*************CALIBRATION - GLOBAL using Crowson's Poisson approach******************;


*- CUMHAZ command provides the expected number of events and XBETA provides linear predictor for original model applied to external cohort;
PROC PHREG DATA=FFPGR;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE/ TIES=EFRON RL;
	BASELINE OUT=GBSG_VAL1 COVARIATES=GBSG_VAL CUMHAZ=P XBETA=lp;
RUN;

*get the original survival time for each person;
DATA PRED1A;
	SET GBSG_VAL1;
	PROC SORT; BY ID SURVTIME1;
RUN;

*calculate log of expected number of events;
DATA PRED1;
	SET PRED1A;
	BY ID SURVTIME;
	*first get the follow-up time for each patient;
	IF LAST.ID;
	*P is the Outcome-Martingale residual;
	LOGP=LOG(P);
	LOGBASE=LOGP-LP;
RUN;

*divide obs by exp and should agree with simple model 1 below;
proc univariate data=pred1 noprint;
	var p STATUS;
	output out=sums sum=sump sumo;
	TITLE ' ';
	proc print;
run;

*MODEL 1 - CALIBRATION IN THE LARGE;
DATA DATA2;
	SET PRED1;
RUN;

*take exponential of Intercept estimate (0.1348) to get SIR (exp(0.1348)) - should agree with ratio of obs/exp above;
PROC GENMOD DATA=DATA2;
	MODEL STATUS=/OFFSET=LOGP DIST=POISSON LINK=LOG;
	OUTPUT OUT=INT PRED=PRED;
RUN;


*MODEL 2 - CALIBRATION SLOPE;
*the coefficient is the slope;
PROC GENMOD DATA=DATA2;
	MODEL STATUS=LP /OFFSET=logbase DIST=POISSON LINK=LOG;
RUN;





********************************** FIXED TIME POINT ASSESSMENT OF CALIBRATION ******************;

********* CALIBRATION-IN-THE-LARGE **********;

*Try O/E using ratio of Kaplan-Meier and Avg predicted prob;

PROC LIFETEST DATA=GBSG_VAL METHOD=KM ATRISK OUTSURV=OUTKM_EXT_C /*NOPRINT*/;
        TIME SURVTIME*STATUS(0);
RUN;

DATA OUTKM_EXT_C1;
	SET OUTKM_EXT_C;
	WHERE _CENSOR_ = 0;
RUN;

*output observed survival at 5 years;
DATA OUTKM_EXT_C2;
	SET OUTKM_EXT_C1;
	BY _CENSOR_;
	IF LAST._CENSOR_;
	n=1;
	KEEP n SURVIVAL;
RUN;

*output predicted;
PROC PHREG DATA=FFPGR;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE/ TIES=EFRON RL;
	BASELINE OUT=GBSG_VAL11 COVARIATES=GBSG_VAL SURVIVAL=PSURV XBETA=lp TIMELIST=5;
RUN;

proc univariate data=GBSG_VAL11 noprint;
	var pSURV;
	output out=mean mean=meanpred;
	TITLE ' ';
	proc print;
run;

DATA MEAN1;
	SET MEAN;
	N=1;
RUN;

* L2 we get O/E of 0.90;
DATA MEAN2;
	MERGE OUTKM_EXT_C2 MEAN1;
	BY N;
	PROP = SURVIVAL/MEANPRED;
RUN;

PROC PRINT DATA=MEAN2;
RUN;


*MODEL 2 - CALIBRATION SLOPE;
PROC PHREG DATA=GBSG_VAL11 ZPH(GLOBAL TRANSFORM=LOG) EV;
	MODEL SURVTIME1*STATUS(0)=LP / TIES=EFRON RL;
RUN;
/*Parameter DF Parameter Estimate Standard Error Chi-Square Pr > ChiSq Hazard Ratio 95% Hazard Ratio Confidence Limits 
  PI  		1 	1.04035 			0.12408 		70.2949 <.0001 			2.830 2.219 3.609 */





**************** - CALIBRATION PLOT *****************;


**calibration plots using Harrell 2020 Stat Med paper, restricted cubic splines;

PROC PHREG DATA=FFPGR;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE/ TIES=EFRON RL;
	BASELINE OUT=GBSG_VALH COVARIATES=GBSG_VAL SURVIVAL=PredProb TIMELIST=5;
RUN;

*get the original survival time for each person;
DATA PREDHar;
	SET GBSG_VALH;
	WHERE SURVTIME1 LE SURVTIME; 
	*Prob death;
	Prob=1-PredProb;
	CLL=LOG(-LOG(1-Prob));
	PROC SORT; BY ID SURVTIME;
RUN;

*CHECK FUNCTIONAL FORM OF CLL;
* Code CLL as a restricted cubic spline with 3 knots;
*First calculate the 10th, 50th and 90th percentiles for knots;
PROC UNIVARIATE DATA = PREDHar;
	VAR CLL;
	OUTPUT OUT=KNOTS PCTLPRE=P_CLL PCTLPTS= 10 50 90;
RUN;

PROC PRINT DATA=KNOTS; RUN;

/* here we find the following values:
P_CLL10 P_CLL50 P_CLL90 
-1.27239 -0.54407 0.26471 

*/

*Use Frank Harrell's RCSPLINE macro for calculating the spline terms - found here: http://biostat.mc.vanderbilt.edu/wiki/Main/SasMacros;
%include 'T:\People\d.mclernon\SAS Progs\RCSPLINE macro.SAS';
DATA PREDHar1;
	SET PREDHar;
	%RCSPLINE(CLL, -1.27239, -0.54407, 0.26471);
RUN;

PROC UNIVARIATE DATA = PREDHar1;
	VAR Prob;
	OUTPUT OUT=gridcox PCTLPRE=P_ PCTLPTS= 1 to 99;
RUN;

PROC TRANSPOSE DATA=GRIDCOX OUT=GRIDCOX1;
	VAR P_1 - P_99;
RUN;

DATA GRIDCOX2(RENAME=(COL1=Prob));
	SET GRIDCOX1;
	CLL=LOG(-LOG(1-COL1));
	*WARNING - TAKE THE BELOW FORMULA FROM LOG WINDOW AFTER RUNNING ABOVE RCSPLINE - THIS WILL BE DIFFERENT ACCORDING TO YOUR DATA!;
	 _kd_= (0.26471 - -1.27239)**.666666666666 ;
	CLL1=max((CLL--1.27239)/_kd_,0)**3+((-0.54407--1.27239)*max((CLL-0.26471)/_kd_,0)**3 -(0.26471--1.27239)*max((CLL--0.54407)/_kd_,0)**3)/(0.26471--0.54407);
RUN;

*Calibration for predictions of 5-year survival probabilities;
PROC PHREG DATA=PREDHar1 ZPH(GLOBAL TRANSFORM=LOG);
	MODEL SURVTIME1*STATUS(0)=CLL CLL1/ TIES=EFRON RL;
	BASELINE OUT=GRIDCOX3 COVARIATES=GRIDCOX2 SURVIVAL=PredProb LOWER=Low UPPER=Upp TIMELIST=5;
RUN;

DATA GRIDCOX4;
	SET GRIDCOX3;
	BY Prob;
	IF FIRST.Prob;
	OBSProb=1-PredProb;
	OBSLow=1-Upp;
	OBSUpp=1-Low;
	IF _NAME_='P_1' THEN DIAG1=0;
	IF _NAME_='P_1' THEN DIAG2=0;
	IF _NAME_='P_15' THEN DIAG1=1;
	IF _NAME_='P_15' THEN DIAG2=1;
RUN;

proc sort data=gridcox4;
by Prob;
run;

*CREATE DENSITY PLOT;
proc univariate data=PREDHar1;
   var Prob;
   histogram Prob / kernel midpoints=(0.005 to 0.995 by 0.01) outkernel=OutHist; /* 100 points in [0,1] */
   *ods select histogram;
run;

data density;
   keep _value_ _density_;
   set OutHist;
run;
 
/* Merge the counts with the predicted probabilities. */
data LogiPlot2;
   set gridcox4(keep=Prob ObsProb OBSLow OBSUpp DIAG1 DIAG2 /*<-ref line*/)
       density;
	*Prob1=1-Prob;
	*Obs1=1-ObsProb;
   *label y0 = "Count"  y1 = "Count";
  * y0=(_density_/2.3)*0.15 - 0.05/*divide by max value of y's from butterflyfringe and scale and offset below 0*/;
   *R=-0.05;
run;

title;
*ODS LISTING SGE=ON STYLE=PRINTER IMAGE_DPI=300 GPATH='T:\People\d.mclernon\STRATOS\STRATOS\Figures';
ODS LISTING SGE=ON STYLE=PRINTER IMAGE_DPI=300 GPATH='C:\Users\sme544\Documents\STRATOS';
ODS GRAPHICS ON / RESET=ALL NOBORDER OUTPUTFMT=TIFF /*WIDTH=4IN*/ IMAGENAME="Calibration plot Ext Harrell with RCS" ANTIALIAS=OFF/*ANTIALIASMAX=*/;

PROC SGPLOT DATA=LogiPlot2 NOAUTOLEGEND ;
XAXIS       	LABEL="Predicted probability of recurrence" VALUES=(0 TO 1 BY 0.20);
YAXIS        	LABEL="Observed recurrence" VALUES=(0 TO 1 BY 0.20);
Y2AXIS;
  SERIES X=Prob Y=ObsProb / LINEATTRS=(COLOR=RED THICKNESS=2) ;
  BAND X=Prob LOWER=OBSLow UPPER=OBSUpp / NOFILL LINEATTRS=(COLOR=BLACK PATTERN=MEDIUMDASH THICKNESS=2) NOEXTEND OUTLINE;
	
  REG X=DIAG1 Y=DIAG2 / LINEATTRS=(PATTERN=SHORTDASH) MARKERATTRS=(SIZE=1PX SYMBOL=CIRCLEFILLED);
  SERIES X=_value_ Y=_Density_ / LINEATTRS=(COLOR=BLACK THICKNESS=1 PATTERN=1) Y2AXIS;
run;

ods graphics off;


*calculation of ICI for 5 year probabilities;
*Calibration for predictions of 5-year survival probabilities;
DATA PREDHarl_Rep;
	SET PREDHar1;
RUN;

PROC PHREG DATA=PREDHar1;
	MODEL SURVTIME1*STATUS(0)=CLL CLL1/ TIES=EFRON RL;
	BASELINE OUT=PREDHarl_Rep1 COVARIATES=PREDHarl_Rep SURVIVAL=PredProb2 TIMELIST=5;
RUN;

DATA ICI;
	SET PREDHarl_Rep1;
	*Observed Prob death;
	Prob2=1-PredProb2;
	Diff=ABS(Prob - Prob2);
RUN;

PROC UNIVARIATE DATA=ICI;
	VAR Diff;
	OUTPUT OUT=ICIS MEAN=ICI MEDIAN=E50 P90=E90;
RUN;




*************CALIBRATION for model with pgr - GLOBAL using Crowson's Poisson approach******************;




*Calibration - Crowson;
*- CUMHAZ provides the expected number of events;
PROC PHREG DATA=FFPGR;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE pgr pgr1/ TIES=EFRON RL;	*STORE Stratos.SimpModel;
	BASELINE OUT=GBSG_VAL1 COVARIATES=GBSG_VAL CUMHAZ=P XBETA=lp;
RUN;

*get the orginal survival time for each person;
DATA PRED1A;
	SET GBSG_VAL1;
	WHERE SURVTIME LE SURVTIME1; 
	PROC SORT; BY ID SURVTIME;
RUN;

*log expected number of events;
DATA PRED1;
	SET PRED1A;
	BY ID SURVTIME;
	*first get the follow-up time for each patient;
	IF LAST.ID;
	LOGP=LOG(P);
	LOGBASE=LOGP-LP;
RUN;

*divide obs by exp and should agree with simple model 1;
proc univariate data=pred1 noprint;
	var p STATUS;
	output out=sums sum=sump sumo;
	TITLE ' ';
	proc print;
run;

*MODEL 1 - IN THE LARGE;
DATA DATA2;
	SET PRED1;
RUN;

*take exponential of observed and get SIR and agrees with ratio of obs/exp above;

PROC GENMOD DATA=DATA2;
	MODEL STATUS=/OFFSET=LOGP DIST=POISSON LINK=LOG;
	OUTPUT OUT=INT PRED=PRED;
RUN;


*MODEL 2 - SLOPE;
*the coefficient is the slope;
PROC GENMOD DATA=DATA2;
	MODEL STATUS=LP /OFFSET=LOGBASE DIST=POISSON LINK=LOG;
RUN;


********************************** FIXED TIME POINT ASSESSMENT OF CALIBRATION ******************;

********* CALIBRATION-IN-THE-LARGE **********;




**calibration plots using Harrell 2020 Stat Med paper, restricted cubic splines;

PROC PHREG DATA=FFPGR;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE pgr pgr1/ TIES=EFRON RL;
	BASELINE OUT=GBSG_VALH COVARIATES=GBSG_VAL SURVIVAL=PredProb TIMELIST=5;
RUN;

*get the original survival time for each person;
DATA PREDHar;
	SET GBSG_VALH;
	WHERE SURVTIME1 LE SURVTIME; 
	*Prob death;
	Prob=1-PredProb;
	CLL=LOG(-LOG(1-Prob));
	PROC SORT; BY ID SURVTIME;
RUN;

*CHECK FUNCTIONAL FORM OF CLL;
* Code PGR as a restricted cubic spline with 3 knots;
*First calculate the 10th, 50th and 90th percentiles for knots;
PROC UNIVARIATE DATA = PREDHar;
	VAR CLL;
	OUTPUT OUT=KNOTS PCTLPRE=P_CLL PCTLPTS= 10 50 90;
RUN;

PROC PRINT DATA=KNOTS; RUN;

/* here we find the following values:
P_CLL10 P_CLL50 P_CLL90 
-1.27239 -0.54407 0.26471 

*/

*Use Frank Harrell's RCSPLINE macro for calculating the spline terms - found here: http://biostat.mc.vanderbilt.edu/wiki/Main/SasMacros;
%include 'T:\People\d.mclernon\SAS Progs\RCSPLINE macro.SAS';
DATA PREDHar1;
	SET PREDHar;
	%RCSPLINE(CLL, -1.27239, -0.54407, 0.26471);
RUN;

PROC UNIVARIATE DATA = PREDHar1;
	VAR Prob;
	OUTPUT OUT=gridcox PCTLPRE=P_ PCTLPTS= 1 to 99;
RUN;

PROC TRANSPOSE DATA=GRIDCOX OUT=GRIDCOX1;
	VAR P_1 - P_99;
RUN;

DATA GRIDCOX2(RENAME=(COL1=Prob));
	SET GRIDCOX1;
	CLL=LOG(-LOG(1-COL1));
	*WARNING - TAKE THE BELOW FORMULA FROM LOG WINDOW AFTER RUNNING ABOVE RCSPLINE - THIS WILL BE DIFFERENT ACCORDING TO YOUR DATA!;
	 _kd_= (0.26471 - -1.27239)**.666666666666 ;
	CLL1=max((CLL--1.27239)/_kd_,0)**3+((-0.54407--1.27239)*max((CLL-0.26471)/_kd_,0)**3 -(0.26471--1.27239)*max((CLL--0.54407)/_kd_,0)**3)/(0.26471--0.54407);
RUN;

*Calibration for predictions of 5-year survival probabilities;
PROC PHREG DATA=PREDHar1 ZPH(GLOBAL TRANSFORM=LOG);
	MODEL SURVTIME1*STATUS(0)=CLL CLL1/ TIES=EFRON RL;
	BASELINE OUT=GRIDCOX3 COVARIATES=GRIDCOX2 SURVIVAL=PredProb TIMELIST=5;
RUN;

DATA GRIDCOX4;
	SET GRIDCOX3;
	BY Prob;
	IF FIRST.Prob;
	OBSProb=1-PredProb;
	IF _NAME_='P_1' THEN DIAG1=0;
	IF _NAME_='P_1' THEN DIAG2=0;
	IF _NAME_='P_15' THEN DIAG1=1;
	IF _NAME_='P_15' THEN DIAG2=1;
RUN;

proc sort data=gridcox4;
by Prob;
run;

*CREATE DENSITY PLOT;
proc univariate data=PREDHar1;
   var Prob;
   histogram Prob / kernel midpoints=(0.005 to 0.995 by 0.01) outkernel=OutHist; /* 100 points in [0,1] */
   *ods select histogram;
run;

data density;
   keep _value_ _density_;
   set OutHist;
run;
 
/* Merge the counts with the predicted probabilities. */
data LogiPlot2;
   set gridcox4(keep=Prob ObsProb DIAG1 DIAG2 /*<-ref line*/)
       density;
	Prob1=1-Prob;
	Obs1=1-ObsProb;
run;

title;
ODS LISTING SGE=ON STYLE=PRINTER IMAGE_DPI=300 GPATH='T:\People\d.mclernon\STRATOS\STRATOS\Figures';
ODS GRAPHICS ON / RESET=ALL NOBORDER OUTPUTFMT=TIFF /*WIDTH=4IN*/ IMAGENAME="Calbration plot pgr Ext Harrell with RCS" ANTIALIAS=OFF/*ANTIALIASMAX=*/;

PROC SGPLOT DATA=LogiPlot2 NOAUTOLEGEND ;
XAXIS       	LABEL="Predicted probability of recurrence" VALUES=(0 TO 1 BY 0.20);
YAXIS        	LABEL="Observed recurrence" VALUES=(0 TO 1 BY 0.20);
Y2AXIS;
  SERIES X=Prob1 Y=OBS1 / LINEATTRS=(COLOR=RED THICKNESS=2) ;

  REG X=DIAG1 Y=DIAG2 / LINEATTRS=(PATTERN=SHORTDASH) MARKERATTRS=(SIZE=1PX SYMBOL=CIRCLEFILLED);
  SERIES X=_value_ Y=_Density_ / LINEATTRS=(COLOR=BLACK THICKNESS=1 PATTERN=1) Y2AXIS;
run;

ods graphics off;

*calculation of ICI for 5 year probabilities;
*Calibration for predictions of 5-year survival probabilities;
DATA PREDHarl_Rep;
	SET PREDHar1;
RUN;

PROC PHREG DATA=PREDHar1;
	MODEL SURVTIME1*STATUS(0)=CLL CLL1/ TIES=EFRON RL;
	BASELINE OUT=PREDHarl_Rep1 COVARIATES=PREDHarl_Rep SURVIVAL=PredProb2 TIMELIST=5;
RUN;

DATA ICI;
	SET PREDHarl_Rep1;
	*Observed Prob death;
	Prob2=1-PredProb2;
	Diff=ABS(Prob - Prob2);
RUN;

PROC UNIVARIATE DATA=ICI;
	VAR Diff;
	OUTPUT OUT=ICIS MEAN=ICI MEDIAN=E50 P90=E90;
RUN;







************** Conduct DCA analysis - see %stdca macro - get from https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/decision-curve-analysis;


*external;

PROC PHREG DATA=FFPGR;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE/ TIES=EFRON RL;	
	BASELINE OUT=ORIGEXT COVARIATES=GBSG_VAL SURVIVAL=FIVEYR /*CUMHAZ=P*/ timelist=5;
RUN;

PROC PHREG DATA=FFPGR;
	CLASS SIZE (REF=FIRST) NODESCAT (REF=FIRST) GRADE (REF=FIRST);
	MODEL SURVTIME*STATUS(0)=SIZE NODESCAT GRADE PGR PGR1/ TIES=EFRON RL;	*STORE Stratos.SimpModel;
	BASELINE OUT=NEWEXT COVARIATES=ORIGEXT SURVIVAL=FIVEYR_PGR /*CUMHAZ=P*/ timelist=5;
RUN;

DATA NEWEXT2;
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

%INCLUDE 'T:\People\d.mclernon\STRATOS\STRATOS\NRI DCA\STDCA.SAS';
%STDCA(data=NEWEXT2, out=survivalmult, outcome=STATUS, ttoutcome=SURVTIME1, timepoint=5, predictors=RISK_ORIG);
%STDCA(data=NEWEXT2, out=survivalmult_new, outcome=STATUS, ttoutcome=SURVTIME1, timepoint=5, predictors=RISK_NEW);

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

data stratos.ExtNB;
	set crsort;
run;

title;
FOOTNOTE;
*- this allows editing of the .sge file!;
ODS LISTING SGE=ON STYLE=PRINTER IMAGE_DPI=300 GPATH='T:\People\d.mclernon\STRATOS\STRATOS\Figures';
ODS GRAPHICS ON / RESET=ALL NOBORDER OUTPUTFMT=TIFF /*WIDTH=4IN*/ IMAGENAME="dca with spline PGR ext" ANTIALIAS=OFF/*ANTIALIASMAX=*/;

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
