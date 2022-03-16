*******************************************************************************;	
* Program: 			STRATOS External validation when original model dataset    ;
*					available.sas  											   ;
* Author: 			David McLernon											   ;
* Date: 			24th Jan 2022											   ;
* Purpose: 			This is SAS code for the STRATOS paper on validation of	   ;
*					survival risk prediction models. This programme covers:    ;
*					1. external validation when the original dataset used to   ;
*					develop the model is available							   ;
*					2. what to do if development dataset or baseline hazard are; 
*					not available (starts LINE 1430)						   ;
* Note:				This programme is not currently automated. It is coded	   ;
*					based on the case study in the article. Therefore, 		   ;
*					adapting this code for your own study will require careful ;
*					editing according to your data							   ;
*					Readers can skip the sections that create "nice" graphs,   ;
*					these all have surrounding comments of  "Optional Block"   ;
*					and "End Optional Block". Though an essential skill for    ;
*					published papers, such maniuplation is outside the primary ;
*					scope of this work.										   ;
* Data:				The Rotterdam and GBSG data sets can be found on the web   ;
*					in various locations and various formats.  For this 	   ;
*					exercise we have used the versions found in the R package. ;
*					(Older web sites have a habit of disappearing and we 	   ;
*					expect R to be more stable). In R we used 				   ;
*				    library(survival)									       ;
*			        write.csv(rotterdam, row.names=FALSE, file="rotterdam.csv");
*			        write.csv(gbsg,      row.names=FALSE, file="gbsg.csv")	   ;
*				       to create the csv files.								   ;
*******************************************************************************;	
options ls=80 mprint mlogic nonumber nodate source2 spool;

* NB Edit folder locations to your own;
* i call libname stratos but you can change to whatever name and location you 
* wish - you may not even wish to use a libname;
libname stratos 'c:\users\sme544\documents\stratos';

* macros required - download and save and specify the correct location after 
* %include;
* - Frank Harrells RCSPLINE macro for calculating the spline terms - found here: 
* http://biostat.mc.vanderbilt.edu/wiki/Main/SasMacros;
%include 'c:\users\sme544\documents\stratos\rcspline macro.sas';
* - Andrew Vickers DCA SAS macro see: 
* https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/
* decision-curve-analysis;
%include 'c:\users\sme544\documents\stratos\stdca.sas';


****** Read in the Rotterdam datset;
proc import out= rotterdam 
            datafile= "c:\users\sme544\documents\stratos\rotterdam.csv" 
            dbms=csv replace;
run;

****** Read in the German datset;
proc import out= gbsg 
            datafile= "c:\users\sme544\documents\stratos\gbsg.csv" 
            dbms=csv replace;
run;

/*
** Create the categorical variables that we use in the model fit
** grade of 1-2 vs 3
** node categories of 0, 1-3, > 3
** size of  <=20, 21-50, >50
**  The size variable is already categorized in this way, in both data sets
*/

*** descriptive statistics for age and PGR;
proc univariate data=rotterdam;
	var age pgr;
run;

proc format;  * make them print in a nice order;
    value nodef 0 = "0"
                1 = "1-3"
                2 = ">3";
    value sizef 0 = "<=20"
                1 = "20-50"
                2 = ">50";
    value gradef 0 = "1-2"
                 1 = "3";

data r1(rename=(pr=pgr)); set rotterdam;
    format nodescat nodef. sizecat sizef. gradecat gradef.;

    nodescat = 1* (1 <= nodes <=3) + 2*(nodes > 3);
    sizecat =  1* (size= "20-50")  + 2*(size = ">50");
    gradecat = grade -2;  * rotterdam has only grade 2 and 3 subjects;

    * recurrence free survival (RFS) = earlier of recurrence or death;
    if (recur = 1) then do;
        survtime = rtime/365.25;  * days to years;
        status = recur;
    end;
    else do;
        survtime = dtime/365.25;
        status = death;
    end;
    
    * Winzorise PGR to the 99th percentile to deal with large influential
    *  values;
    if (pgr > 1360) then pr = 1360; else pr = pgr;
	drop pgr;
run;

* Descriptive statistics at baseline for Table 1;
proc freq data=r1;
    table nodescat gradecat sizecat;
run;

*Administrative censor at 5 years since this is our prediction horizon; 
data ffpgr;
	set r1;
	*administrative censoring at 5 years;
	if survtime > 5 then status=0;
	if survtime > 5 then survtime=5;
	* rcs for pgr - see development program for knot calculations;
	%rcspline(pgr, 0, 41, 486);
run;

************** Import external validation dataset ****************************;

*code up predictors;
data gbsg_va; set gbsg;
    format nodescat nodef. sizecat sizef. gradecat gradef.;
	*Winzorise PGR to 99th percentile to deal with very large influential 
		values;
	if pgr>1360 then pgr=1360;
	%rcspline(pgr, 0, 41, 486);
    nodescat = 1* (1 <= nodes <=3) + 2*(nodes > 3);
	if size le 20 then sizecat=0;
	if 21 <= size <= 50 then sizecat=1;
	if size >50 then sizecat=2;
	gradecat = grade-1;
	if grade in (1,2) then gradecat=0;
	if grade = 3 then gradecat=1;
    survtime = rfstime/365.25;  * days to years;
	keep pid sizecat nodescat gradecat status age pgr pgr1	survtime;
run;


* Plot the overall Kaplan-Meier, as motivation for using a 5 year cut-off;
proc lifetest data=gbsg_va method=pl plots=(s, ls, lls) outsurv=outkm;
        time survtime*status(0);
		ods exclude ProductLimitEstimates;   * suppress the long table;
		title "Overall Kaplan-Meier";
run;

******* Optional Block ****** ;
** to create a nicer Kaplan-Meier plot;

data outkm1; 
	set outkm; 
	if _censor_=1 then delete; 
	keep survival survtime; 
run;

*Create kaplan meier plot for external dataset;
title;
footnote;
*- this allows editing of the .sge file;
ods listing sge=on style=printer image_dpi=300 gpath='c:';
ods graphics on / reset=all noborder outputfmt=tiff /*width=4in*/ 
	imagename="val km" antialias=off/*antialiasmax=*/;

*create graph;
proc sgplot data=outkm1;
	yaxis values=(0 to 1 by 0.2) label="Recurrence-free survival probability";
	xaxis values=(0 to 7 by 1) label="Years";
	step y=survival x=survtime / lineattrs=(color=blue thickness=2 
		pattern=solid) name="all";
run;

ods graphics off;

******  End Optional Block ************;

*** Descriptive statistics;
proc freq data=gbsg_va;
	tables sizecat nodescat gradecat;
run;

proc univariate data=gbsg_va;
	var age pgr;
run;

*use reverse Kaplan-Meier to obtain median follow-up time;
proc lifetest data=gbsg_va method=pl atrisk;
    time survtime*status(1);
    ods exclude ProductLimitEstimates;   * suppress the long table;
    title "Median follow-up time";
run;

*Administrative censor at 5 years as development cohort; 
data gbsg_val;
	set gbsg_va;
	*administrative censoring at 5 years;
	if survtime > 5 then status=0;
	if survtime > 5 then survtime=5;
	*- keep a copy of survtime for calibration step;
	survtime1=survtime;
run;


******************EXTERNAL VALIDATION *********************;

*******************Simple model first (without PGR)***************************;

* Fit the model to the development data, and calculate linear predictor of 
* this existing model for patients in the external dataset;

proc phreg data=ffpgr;
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat / ties=efron rl;
	*Store allows SAS to save the model estimates for applying later on to 
	*external data;
	store simpmodel;
	*Calculate model (developed on Rotterdam data) linear predictor to external 
	*German dataset;
	baseline covariates=gbsg_val out=gbsg_valrd xbeta=xb;
	*Calculate model (developed on Rotterdam data) linear predictor for 
	*patients in original Rotterdam dataset;
	output out=rottx xbeta=xb;
run;

*merge linear predictor to external dataset;
proc sort data=gbsg_valrd;
	by pid;
run;

proc sort data=gbsg_val;
	by pid;
run;

data gbsg_valrd1;
	set gbsg_valrd;
	by pid;
	if first.pid;
	keep pid xb;
run;

data gbsg_valrd2;
	merge gbsg_val gbsg_valrd1;
	by pid;
run;

******************** Extended model (with PGR) *************************;

* Fit the model to the development data, and calculate linear predictor of 
* this existing model for patients in the external dataset;
proc phreg data=ffpgr;
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat pgr pgr1/ ties=efron rl;
	*Store allows SAS to save the model estimates for applying later on to 
	*external data;
	store simpmodelpgr;
	*Calculate model (developed on Rotterdam data) linear predictor to external 
	*German dataset;
	baseline covariates=gbsg_val out=gbsg_valrda xbeta=xb;
	*Calculate model (developed on Rotterdam data) linear predictor for 
	*patients in original Rotterdam dataset;
	output out=rottxa xbeta=xb;
run;

*merge linear predictor to external dataset;
proc sort data=gbsg_valrda;
	by pid;
run;

proc sort data=gbsg_val;
	by pid;
run;

data gbsg_valrd1a;
	set gbsg_valrda;
	by pid;
	if first.pid;
	keep pid xb;
run;

data gbsg_valrd2a;
	merge gbsg_val gbsg_valrd1a;
	by pid;
run;

****** External validation when you have all the development data *************;

*****TIME RANGE ASSESSMENT OF DISCRIMINATION;

*Harrell's C - Need tau to equal event time of interest;
title 'External Harrells C';
proc phreg data=gbsg_val concordance=harrell(se) tau=5;
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat pgr pgr1/ ties=efron rl;
	*the following statements call the stored models from earlier and 
	calculate the concordance statistics for model with and without pgr;
	roc 'npi' source=simpmodel;
	roc 'npi + pgr' source=simpmodelpgr;
run;

*Uno's C - Need tau to equal event time of interest;
title 'External Unos C';
PROC PHREG DATA=GBSG_VAL CONCORDANCE=UNO(SE SEED=8754 ITER=50) TAU=5;
	CLASS sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat pgr pgr1/ ties=efron rl;
	*the following statements call the stored models from earlier and 
	calculate the concordance statistics for model with and without PGR;
	roc 'npi' source=simpmodel;
	roc 'npi + pgr' source=simpmodelpgr;
RUN;


****FIXED TIME POINT ASSESSMENT OF DISCRIMINATION;

*- The last observed event time in  the german breast cancer dataset is 4.966 
*- so use 4.96 for time point as estimation requires events after time of 
*- interest;
title 'External Unos AUC';
proc phreg data=gbsg_val rocoptions(auc at=4.96 method=ipcw (cl seed=134));
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat / ties=efron rl nofit;
	roc 'npi' source=simpmodel;
	roc 'npi + pgr' source=simpmodelpgr;
run;

******* Optional Block ****** ;
* Rerun Uno and save dataset to plot the Time-dependent AUC yourself;
proc phreg data=gbsg_val rocoptions(method=ipcw (cl seed=134) outauc=ipcwauc2);
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat pgr pgr1/ ties=efron rl 
		nofit;
	roc 'npi + pgr' source=simpmodelpgr;
	roc 'npi' source=simpmodel;
run;

*plot time-dependent AUC using outauc dataset;
data ipcwauc3;
	set ipcwauc2;
	if survtime >5 then delete;
run;

title;
footnote;
ods listing sge=on style=printer image_dpi=300 gpath='c:';
ods graphics on / reset=all noborder outputfmt=tiff imagename="Time-dep ROC" 
	antialias=off;

proc sgplot data=ipcwauc3 noautolegend ;
xaxis       	label="Survival time (years)" values=(0 to 5 by 1);
yaxis        	label="AUC" values=(0 to 1 by 0.20);
  band x=survtime lower=_LowerAUC_ upper=_UpperAUC_ / transparency=0.5 
		group=_SOURCE_;
  series x=survtime y=_auc_ / lineattrs=(thickness=1) group=_source_ 
		grouplc=_id_ name="models" legendlabel="models";
  refline 0.5 / axis=y lineattrs=(color=black thickness=1);
  keylegend "models" / border location=inside position=bottomright down=2 
		title="models";
run;

ods graphics off;

******  End Optional Block ************;

************* END OF DISCRIMINATION CODE ***************************;




*************CALIBRATION - TIME RANGE ASSESSMENT ******************;

***Simple model first;
**Using Crowson's Poisson approach;

*- CUMHAZ command provides the expected number of events and XBETA provides 
	linear predictor for original model applied to external cohort;
proc phreg data=ffpgr;
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat/ ties=efron rl;
	baseline out=gbsg_val1 covariates=gbsg_val cumhaz=p xbeta=lp;
run;

*get the original survival time for each person;
data pred1a;
	set gbsg_val1;
	where survtime le survtime1; 
	proc sort; by pid survtime1;
run;

data testl;
	set gbsg_val1;
	proc sort; by survtime;
run;

*calculate log of expected number of events;
data pred1;
	set pred1a;
	by id survtime;
	*first get the follow-up time for each patient;
	if last.pid;
	*p is the outcome-martingale residual;
	logp=log(p);
	logbase=logp-lp;
run;

*divide observed by expected and should agree with simple model 1 below;
proc univariate data=pred1 noprint;
	var p status;
	output out=sums sum=sump sumo;
	TITLE 'Sum of observed and expected events within 5 years';
	proc print;
run;

*Model 1 - Mean calibration (O/E);
data data2;
	set pred1;
run;

*take exponential of Intercept estimate (0.0544) to get SIR (exp(0.0544)) 
- should agree with ratio of obs/exp above;
proc genmod data=data2;
	model status=/offset=logp dist=poisson link=log;
	output out=int pred=pred;
run;

******* Optional Block ****** ;

*For plot (Fig 1B in paper) get all observed times up to the event or censoring 
time of each patient;
data plotall;
	set pred1a;
	proc sort; by pid survtime;
run;

* get the last observed time per patient and save status;
data lastob;
	set pred1a;
	by pid survtime;
	*first get the follow-up time for each patient;
	if last.pid;
	stat=status;
	keep pid stat survtime;
run;

*merge back and then put status for all observed times before last equal to 
zero;
data plotall1;
	merge plotall lastob;
	by pid survtime;
run;

data plotall2;
	set plotall1;
	if stat=. then stat=0;
	proc sort; by pid survtime;
run;

*make macro for every cumulative month to 5 years;
%macro plotoe(month);
	%do i=1 %to &month;
		data timepts;
			set plotall2;
			if survtime*12 gt &i then delete;
		run;

		data timepts2;
			set timepts;
			by pid survtime;
			if last.pid;
			logp=log(p);
		run;

		proc genmod data=timepts2;
			model stat=/offset=logp dist=poisson link=log;
			ods output parameterestimates=param;
		run;

		*get number obs and exp at each time;
		data sumsall2;
			set param;
			month=&i;
			if Parameter='Scale' then delete;
			Est=exp(Estimate);
			Low=exp(LowerWaldCL);
			Upp=exp(UpperWaldCL);
			keep month Est Low Upp;
		run;

		proc append base=sumoe1 data=sumsall2 force;
		run;
	%end;
%mend;

%plotoe(60);

data sumoe2;
	set sumoe1;
	R=1;
run;

title;
ods listing sge=on style=printer image_dpi=300 gpath='c:';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="calbration plot o/e basic" antialias=off;

proc sgplot data=sumoe2 noautolegend ;
xaxis       	label="Month of follow-up" values=(0 to 60 by 6);
yaxis        	label="Number of observed events / Number of expected events" 
	Values=(0 to 1.6 by 0.1) /*offsetmin=0.19*/;
y2axis;
  series x=month y=est / lineattrs=(color=red thickness=2) ;
  refline r / axis=y  lineattrs=(color=black thickness=1);
  band x=month lower=low upper=upp / nofill lineattrs=(color=black 
	pattern=mediumdash thickness=2) noextend outline;
run;

ods graphics off;


******  End Optional Block ************;

*Weak assessment - calibration slope;
*the coefficient is the slope;
proc genmod data=data2;
	model status=lp /offset=logbase dist=poisson link=log;
run;

*Moderate assessment;
** Check functional form of lp;
* Code LP as a restricted cubic spline with 3 knots;
* First calculate the 10th, 50th and 90th percentiles for knots;
proc univariate data = data2 noprint;
	var lp;
	output out=knots pctlpre=p_lp pctlpts= 10 50 90;
run;

proc print data=knots; run;

/* here we find the following values:
 0.35998 1.06273 1.82091 
*/

data test;
	set data2;
	%rcspline(lp, 0.35998, 1.06273, 1.82091);
run;

*Fit Poisson model with spline terms and save predictions;
proc genmod data=test;
	model status=lp lp1 /offset=logp dist=poisson link=log;
	output out=test1 pred=predobs lower=low upper=upp;
run;

*Sort and code diagonal ref line;
data test2;
	set test1;
	if pid=1 then diag1=0;
	if pid=1 then diag2=0;
	if pid=6 then diag1=3;
	if pid=6 then diag2=3;
	proc sort; by p;
run;

*Select unique predictions;
data test3;
	set test2;
	by p;
	if first.p;
run;

*Plot flexible calibration curve;
title;
ods listing sge=on style=printer image_dpi=300 gpath='C:';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="Calibration plot Poisson" ANTIALIAS=OFF;

proc sgplot data=test3 noautolegend ;
	xaxis       	label="Cumulative hazard from Cox model" 
		values=(0 to 2.4 by 0.20);
	yaxis        	label="Cumulative hazard from Poisson model" 
		values=(0 to 2.4 by 0.20);
  	loess x=p y=predobs / nomarkers lineattrs=(color=red thickness=2) ;
  	loess x=p y=low / nomarkers lineattrs=(color=black 
		pattern=mediumdash thickness=2);
  	loess x=p y=upp / nomarkers lineattrs=(color=black 
		pattern=mediumdash thickness=2);
	reg x=diag1 y=diag2 / lineattrs=(pattern=shortdash) markerattrs=(size=1px 
		symbol=circlefilled);
run;

ods graphics off;




**************** FIXED TIME POINT ASSESSMENT OF CALIBRATION ******************;

*** Mean / Calibration-in-the-large **;

*O/E using ratio of Kaplan-Meier and Avg predicted risk;
proc lifetest data=gbsg_val method=km atrisk outsurv=outkm_ext_c noprint;
        time survtime*status(0);
run;

data outkm_ext_c1;
	set outkm_ext_c;
	where _censor_ = 0;
run;

*Output observed survival at 5 years;
data outkm_ext_c2;
	set outkm_ext_c1;
	by _censor_;
	if last._censor_;
	n=1;
	keep n survival;
run;

*Output predicted;
proc phreg data=ffpgr;
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat/ ties=efron rl;
	baseline out=gbsg_val11 covariates=gbsg_val survival=psurv xbeta=lp 
		timelist=5;
run;

proc univariate data=gbsg_val11 noprint;
	var psurv;
	output out=mean mean=meanpred;
run;

data mean1;
	set mean;
	n=1;
run;

data mean2;
	merge outkm_ext_c2 mean1;
	by n;
	prop = (1-survival)/(1-meanpred);
run;

title 'Calibration-in-large using ratio of Kaplan-Meier & Avg predicted risk';
proc print data=mean2;
run;


**Weak - Calibration slope;

proc phreg data=ffpgr;
	class sizecat (ref=first) nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat/ ties=efron rl;
	baseline out=gbsg_valh covariates=gbsg_val survival=predprob xbeta=xb timelist=5;
run;

data t1;
	set gbsg_valh;
	xb_=xb;
run;

proc phreg data=t1;
	model survtime1*status(0) = xb  /  ties=efron rl;
run;

*miscalibration;
proc phreg data=t1;
	model survtime1*status(0) = xb  / offset = xb_ ties=efron rl;
run;


**Moderate assessment - Calibration plots using Austin et al 2020 Stat Med paper;

*get the original survival time for each person;
data predhar;
	set gbsg_valh;
	where survtime1 le survtime; 
	*probability of event;
	prob=1-predprob;
	*Complementary log-log of predicted risk;
	cll=log(-log(1-prob));
	proc sort; by pid survtime;
run;

*Check functional form of cll;
* Code CLL as a restricted cubic spline with 3 knots;
*First calculate the 10th, 50th and 90th percentiles for knots;
proc univariate data = predhar;
	var cll;
	output out=knots pctlpre=p_cll pctlpts= 10 50 90;
run;

proc print data=knots; run;

/* here we find the following values:
Obs P_CLL10 P_CLL50 P_CLL90 
1 -1.14858 -0.44583 0.31235 
*/

*Use Frank Harrell's RCSPLINE macro for calculating the spline terms;
data predhar1;
	set predhar;
	* NOTE: you will need to take the spline information from the log 
	window and edit gridcox2 below;
	%rcspline(cll, -1.14858, -0.44583, 0.31235);
run;

proc univariate data = predhar1 noprint;
	var prob;
	output out=gridcox pctlpre=p_ pctlpts= 1 to 99;
run;

proc transpose data=gridcox out=gridcox1;
	var p_1 - p_99;
run;

data gridcox2(rename=(col1=prob));
	set gridcox1;
	cll=log(-log(1-col1));
	*Warning! - take the below formula from the log window after running above 
	rcspline - this will be different according to your data!;
	 _kd_= (0.31235 - -1.14858)**.666666666666 ;
	cll1=max((cll--1.14858)/_kd_,0)**3
		+((-0.44583--1.14858)*max((cll-0.31235)/_kd_,0)**3
		-(0.31235--1.14858)*max((cll--0.44583)/_kd_,0)**3)/(0.31235--0.44583);
run;

*Calibration for predictions of 5-year survival probabilities;
proc phreg data=predhar1 zph(global transform=log);
	model survtime1*status(0)=cll cll1/ ties=efron rl;
	baseline out=gridcox3 covariates=gridcox2 survival=predprob lower=low 
		upper=upp timelist=5;
run;

data gridcox4;
	set gridcox3;
	by prob;
	if first.prob;
	obsprob=1-predprob;
	obslow=1-upp;
	obsupp=1-low;
	if _name_='p_1' then diag1=0;
	if _name_='p_1' then diag2=0;
	if _name_='p_15' then diag1=1;
	if _name_='p_15' then diag2=1;
run;

proc sort data=gridcox4;
	by Prob;
run;

*Create density plot;
proc univariate data=PREDHar1 noprint;
   var Prob;
   histogram Prob / kernel midpoints=(0.005 to 0.995 by 0.01) outkernel=OutHist; 
run;

data density;
   keep _value_ _density_;
   set OutHist;
run;
 
/* Merge the counts with the predicted probabilities */
data LogiPlot2;
   set gridcox4(keep=Prob ObsProb OBSLow OBSUpp DIAG1 DIAG2 /*<-ref line*/)
       density;
run;

title;
ods listing sge=on style=printer image_dpi=300 gpath='C:';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="Calibration plot Ext Harrell with RCS" ANTIALIAS=OFF;

proc sgplot data=logiplot2 noautolegend ;
	xaxis       	label="Predicted risk from developed model" 
		values=(0 to 1 by 0.20);
	yaxis        	label="Predicted risk from refitted model" 
		values=(0 to 1 by 0.20);
	y2axis			label=" " values=(0 to 20 by 20) display=none;
  	series x=prob y=obsprob / lineattrs=(color=red thickness=2) ;
  	band x=prob lower=obslow upper=obsupp / nofill lineattrs=(color=black 
		pattern=mediumdash thickness=2) noextend outline;
	reg x=diag1 y=diag2 / lineattrs=(pattern=shortdash) markerattrs=(size=1px 
		symbol=circlefilled);
  	series x=_value_ y=_density_ / lineattrs=(color=black thickness=1 pattern=1) 
		y2axis;
run;

ods graphics off;

*Calculation of ICI for 5 year probabilities;
data predharl_rep;
	set predhar1;
run;

proc phreg data=predhar1;
	model survtime1*status(0)=cll cll1/ ties=efron rl;
	baseline out=predharl_rep1 covariates=predharl_rep survival=predprob2 
		timelist=5;
run;

data ici;
	set predharl_rep1;
	*observed prob death;
	prob2=1-predprob2;
	diff=abs(prob - prob2);
run;

proc univariate data=ici noprint;
	var diff;
	output out=icis mean=ici median=e50 p90=e90;
run;

title 'Calibration metrics: ICI, E50, E90';
proc print data=icis;
run;

title;





***EXTENDED MODEL WITH PGR INCLUDED;
**Using Crowson's Poisson approach;

*- CUMHAZ command provides the expected number of events and XBETA provides 
	linear predictor for original model applied to external cohort;
proc phreg data=ffpgr;
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat pgr pgr1/ ties=efron rl;	*store stratos.simpmodel;
	baseline out=gbsg_val1 covariates=gbsg_val cumhaz=p xbeta=lp;
run;

*get the orginal survival time for each person;
data pred1a;
	set gbsg_val1;
	where survtime le survtime1; 
	proc sort; by pid survtime;
run;

*calculate log of expected number of events;
data pred1;
	set pred1a;
	by pid survtime;
	*first get the follow-up time for each patient;
	if last.pid;
	*p is the outcome-martingale residual;
	logp=log(p);
	logbase=logp-lp;
run;

*divide observed by expected and should agree with simple model 1 below;
proc univariate data=pred1 noprint;
	var p status;
	output out=sums sum=sump sumo;
	title 'Sum of observed and expected events within 5 years';
	proc print;
run;

*Model 1 - Mean calibration (O/E);
data data2;
	set pred1;
run;

*take exponential of Intercept estimate (0.0214) to get SIR (exp(0.0544)) 
- should agree with ratio of obs/exp above;
proc genmod data=data2;
	model status=/offset=logp dist=poisson link=log;
	output out=int pred=pred;
run;


******* Optional Block ****** ;

*For plot (Fig 1B in paper) get all observed times up to the event or censoring 
time of each patient;
data plotall;
	set pred1a;
	proc sort; by pid survtime;
run;

* get the last observed time per patient and save status;
DATA lastob;
	set pred1a;
	by pid survtime;
	*first get the follow-up time for each patient;
	if last.pid;
	stat=status;
	keep pid stat survtime;
run;

*merge back and then put status for all observed times before last equal to 
zero;
data plotall1;
	merge plotall lastob;
	by pid survtime;
run;

data plotall2;
	set plotall1;
	if stat=. then stat=0;
	proc sort; by pid survtime;
run;

*make macro for every cumulative month to 5 years;
%macro plotoe(month);
	%do i=1 %to &month;
		data timepts;
			set plotall2;
			if survtime*12 gt &i then delete;
		run;

		data timepts2;
			set timepts;
			by pid survtime;
			if last.pid;
			logp=log(p);
		run;

		proc genmod data=timepts2;
			where p>0;
			model stat=/offset=logp dist=poisson link=log;
			ods output parameterestimates=paramx;
		run;

		*get number obs and exp at each time;
		data sumsall2;
			set paramx;
			month=&i;
			if Parameter='Scale' then delete;
			Est=exp(Estimate);
			Low=exp(LowerWaldCL);
			Upp=exp(UpperWaldCL);
			keep month Est Low Upp;
		run;

		proc append base=sumoex data=sumsall2 force;
		run;
	%end;
%mend;

%plotoe(60);

data sumoe4;
	set sumoex;
	R=1;
run;

title;
ods listing sge=on style=printer image_dpi=300 gpath='c:';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="calbration plot o/e pgr" antialias=off;

proc sgplot data=sumoe4 noautolegend ;
xaxis       	label="Month of follow-up" values=(0 to 60 by 6);
yaxis        	label="Number of observed events / Number of expected events" 
	values=(0 to 1.6 by 0.1) /*offsetmin=0.19*/;
y2axis;
  series x=month y=est / lineattrs=(color=red thickness=2) ;
  refline r / axis=y  lineattrs=(color=black thickness=1);
  band x=month lower=low upper=upp / nofill lineattrs=(color=black 
	pattern=mediumdash thickness=2) noextend outline;
run;

ods graphics off;

******  End Optional Block ************;

*Weak calibration - calibration slope;
*the coefficient is the slope;
proc genmod data=data2;
	model status=lp /offset=logbase dist=poisson link=log;
run;

*Moderate calibration;
** Check functional form of lp;
* Code LP as a restricted cubic spline with 3 knots;
* First calculate the 10th, 50th and 90th percentiles for knots;
proc univariate data = data2 noprint;
	var lp;
	output out=knots pctlpre=p_lp pctlpts= 10 50 90;
run;

proc print data=knots; run;

/* here we find the following values:
 0.27212 0.90521 1.65965 
*/

data test;
	set data2;
	%rcspline(lp, 0.27212, 0.90521, 1.65965);
run;

*Fit Poisson model with spline terms for PI;
proc genmod data=test;
	model status=lp lp1 /offset=logp dist=poisson link=log;
	output out=test1 pred=predobs lower=low upper=upp;
run;

data test2;
	set test1;
	if pid=1 then diag1=0;
	if pid=1 then diag2=0;
	if pid=6 then diag1=3;
	if pid=6 then diag2=3;
	proc sort; by p;
run;

data test3;
	set test2;
	by p;
	if first.p;
run;

*Flexible calibration plot;
title;
ods listing sge=on style=printer image_dpi=300 gpath='C:';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="Calibration plot Poisson PGR" ANTIALIAS=OFF;

proc sgplot data=test3 noautolegend ;
	xaxis       	label="Cumulative hazard from Cox model" 
		values=(0 to 2.4 by 0.20);
	yaxis        	label="Cumulative hazard from Poisson model" 
		values=(0 to 2.4 by 0.20);
  	loess x=p y=predobs / nomarkers lineattrs=(color=red thickness=2) ;
  	loess x=p y=low / nomarkers lineattrs=(color=black 
		pattern=mediumdash thickness=2);
  	loess x=p y=upp / nomarkers lineattrs=(color=black 
		pattern=mediumdash thickness=2);
	reg x=diag1 y=diag2 / lineattrs=(pattern=shortdash) markerattrs=(size=1px 
		symbol=circlefilled);
run;

ods graphics off;





**************** FIXED TIME POINT ASSESSMENT OF CALIBRATION ******************;

*** Calibration-in-the-large **;

*O/E using ratio of Kaplan-Meier and Avg predicted risk;
proc lifetest data=gbsg_val method=km atrisk outsurv=outkm_ext_c noprint;
        time survtime*status(0);
run;

data outkm_ext_c1;
	set outkm_ext_c;
	where _censor_ = 0;
run;

*Output observed survival at 5 years;
data outkm_ext_c2;
	set outkm_ext_c1;
	by _censor_;
	if last._censor_;
	n=1;
	keep n survival;
run;

*Output predicted;
proc phreg data=ffpgr;
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat pgr pgr1/ ties=efron rl;
	baseline out=gbsg_val11 covariates=gbsg_val survival=psurv cumhaz=tch 
		xbeta=lp timelist=5;
run;

proc univariate data=gbsg_val11 noprint;
	var psurv;
	output out=mean mean=meanpred;
run;

data mean1;
	set mean;
	n=1;
run;

data mean2;
	merge outkm_ext_c2 mean1;
	by n;
	prop = (1-survival)/(1-meanpred);
run;

title 'O/E using ratio of Kaplan-Meier and Avg predicted risk';
proc print data=mean2;
run;


*Weak assessment - Calibration slope;

proc phreg data=ffpgr;
	class sizecat (ref=first) nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat pgr pgr1/ ties=efron rl;
	baseline out=gbsg_valh covariates=gbsg_val survival=predprob xbeta=xb timelist=5;
run;

data t1;
	set gbsg_valh;
	xb_=xb;
run;

data t2;
	set t1;
run;

proc phreg data=t1;
	model survtime1*status(0) = xb  /  ties=efron rl;
run;

*miscalibration;
proc phreg data=t1;
	model survtime1*status(0) = xb  / offset = xb_ ties=efron rl;
run;


**Moderate assessment - Calibration plot;

**Calibration plots using Austin et al 2020 Stat Med paper, restricted cubic 
splines;

*get the original survival time for each person;
data predhar;
	set gbsg_valh;
	where survtime1 le survtime; 
	*probability of event;
	prob=1-predprob;
	*Complementary log-log of predicted risk;
	cll=log(-log(1-prob));
	proc sort; by pid survtime;
run;

*Check functional form of cll;
* Code CLL as a restricted cubic spline with 3 knots;
*First calculate the 10th, 50th and 90th percentiles for knots;
proc univariate data = predhar;
	var cll;
	output out=knots pctlpre=p_cll pctlpts= 10 50 90;
run;

proc print data=knots; run;

/* here we find the following values:
P_CLL10 P_CLL50 P_CLL90 
-1.01384 -0.38075 0.37368  
*/

*Use Frank Harrell's RCSPLINE macro for calculating the spline terms;
data predhar1;
	set predhar;
	* NOTE: you will need to take the spline information from the log 
	window and edit gridcox2 below;
	%rcspline(cll, -1.01384, -0.38075, 0.37368);
run;

proc univariate data = predhar1;
	var prob;
	output out=gridcox pctlpre=p_ pctlpts= 1 to 99;
run;

proc transpose data=gridcox out=gridcox1;
	var p_1 - p_99;
run;

data gridcox2(rename=(col1=prob));
	set gridcox1;
	cll=log(-log(1-col1));
	*Warning! - take the below formula from the log window after running above 
	rcspline - this will be different according to your data!;
	 _kd_= (0.37368 - -1.01384)**.666666666666 ;
	 cll1=max((cll--1.01384)/_kd_,0)**3
		+((-0.38075--1.01384)*max((cll-0.37368)/_kd_,0)**3
		-(0.37368--1.01384)*max((cll--0.38075)/_kd_,0)**3)/(0.37368--0.38075);
run;

*Calibration for predictions of 5-year survival probabilities;
proc phreg data=predhar1 zph(global transform=log);
	model survtime1*status(0)=cll cll1/ ties=efron rl;
	baseline out=gridcox3 covariates=gridcox2 survival=predprob lower=low 
		upper=upp timelist=5;
run;

data gridcox4;
	set gridcox3;
	by prob;
	if first.prob;
	obsprob=1-predprob;
	obslow=1-upp;
	obsupp=1-low;
	if _name_='p_1' then diag1=0;
	if _name_='p_1' then diag2=0;
	if _name_='p_15' then diag1=1;
	if _name_='p_15' then diag2=1;
run;

proc sort data=gridcox4;
	by Prob;
run;

*Create density plot;
proc univariate data=PREDHar1 noprint;
   var Prob;
   histogram Prob / kernel midpoints=(0.005 to 0.995 by 0.01) outkernel=OutHist; 
run;

data density;
   keep _value_ _density_;
   set OutHist;
run;
 
/* Merge the counts with the predicted probabilities. */
data LogiPlot2;
   set gridcox4(keep=Prob ObsProb OBSLow OBSUpp DIAG1 DIAG2 /*<-ref line*/)
       density;
run;

title;
ods listing sge=on style=printer image_dpi=300 gpath='c:';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="Calbration plot pgr Ext Harrell with RCS" antialias=off;

proc sgplot data=logiplot2 noautolegend ;
	xaxis       	label="Predicted risk from developed model" 
		values=(0 to 1 by 0.20);
	yaxis        	label="Predicted risk from refitted model" 
		values=(0 to 1 by 0.20);
	y2axis			label=" " values=(0 to 20 by 20) display=none;
  	series x=prob y=obsprob / lineattrs=(color=red thickness=2) ;
  	band x=prob lower=obslow upper=obsupp / nofill lineattrs=(color=black 
		pattern=mediumdash thickness=2) noextend outline;
	reg x=diag1 y=diag2 / lineattrs=(pattern=shortdash) markerattrs=(size=1px 
		symbol=circlefilled);
  	series x=_value_ y=_density_ / lineattrs=(color=black thickness=1 pattern=1) 
		y2axis;
run;

ods graphics off;


*Calculation of ICI for 5 year probabilities;
data predharl_rep;
	set predhar1;
run;

proc phreg data=predhar1;
	model survtime1*status(0)=cll cll1/ ties=efron rl;
	baseline out=predharl_rep1 covariates=predharl_rep survival=predprob2 
		timelist=5;
run;

data ici;
	set predharl_rep1;
	*observed prob death;
	prob2=1-predprob2;
	diff=abs(prob - prob2);
run;

proc univariate data=ici noprint;
	var diff;
	output out=icis mean=ici median=e50 p90=e90;
run;

title 'Calibration metrics: ICI, E50, E90';
proc print data=icis;
run;

title;




*********************** Overall performance;

******* This block calculates the Brier score - see Graf et al 1999;
title ' ';

*Simple model first;
*External validation - Calculate Brier score;

*Calculate weights using Kaplan-Meier;
proc lifetest data=gbsg_val method=pl atrisk outsurv=outkm_ext noprint;
        time survtime*status(1);
run;

*Create 3 groups - Group 1-Those who have the event up to fixed event time of 
*interest, Group 2 - those who go beyond fixed time (could be event or event 
*free), and Group 3- those censored up to fixed time 
*Only first 2 groups contribute to score but all to weights;
data rott_b_ex;
	set gbsg_val;
	*Must make the time just under 5 years since we need some time remaining
	*after timepoint of interest (5 years) and before administrative censoring 
	*(5 years) for it to work;
	if survtime<=4.99 and status=1 then cat=1;
	if survtime>4.99 then cat=2;
	if survtime<=4.99 and status=0 then cat=3;
	*Duplicate survival time as sas will remove the official survival time 
	*variable in baseline statement;
	time=survtime;
run;

*Now estimate survival at 5 years in validation dataset;
proc phreg data=rottx;
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat  / ties=efron rl;
	*external;
	baseline covariates=rott_b_ex out=rottval_bs timelist=5 
		survival=fiveyrsurv xbeta=xb;
run;

*code up the fixed time of 5 years;
data rottval_bs1(rename=(time=survtime));
	set rottval_bs;
	time1=time;
	if time1>5 then time1=5;
	drop survtime;
	proc sort; by time1;
run;

*Merge the kaplan-meier weights to the appropriate times;
data outkm_ext1(rename=(survtime=time1));
	set outkm_ext;
	if survtime>5 then delete;
	weight=1/survival;
run;

data outkm_ext2;
	set outkm_ext1;
	s=lag(time1);
	if time1=s then delete;
	drop s;
run;

data rottval_bs2;
	merge rottval_bs1 outkm_ext1;
	by time1;
	if pid=. then delete;
run;

data rottval_bs3;
	set rottval_bs2;
	retain _weight;
	if not missing(weight) then _weight=weight;
	else weight=_weight;
	if cat=3 then weight=0;
	if time1=0 then delete;
	if cat=1 then contrib=(-fiveyrsurv)**2;
	if cat=2 then contrib=(1-fiveyrsurv)**2;
	if cat=3 then contrib=0;
	bs=contrib*weight;
	drop _weight;
run;

*Estimate brier score;
proc univariate data=rottval_bs3 noprint;
	var bs weight;
	output out=sumsex sum=sbs sweight;
	proc print; 
run;

data sumsex;
	retain sweight sbs brier;
	set sumsex;
	brier = (1/sweight)*sbs;
	sweight=left(sweight);
	title 'External brier score';
	proc print;
run;

******* End of block;

********************************Scaled Brier;
*Scaled Brier = 1 - (model Brier score/null model Brier score), where null 
*model Brier score is null cox; 
*100% is perfect, <0 is useless, higher better, harmful models <0;

*Estimate survival at 5 years for null model;
title 'Null model';

proc phreg data=gbsg_val;
	model survtime*status(0)= / ties=efron rl;
	baseline covariates=rott_b_ex out=rottval_bsnull timelist=5 
		survival=fiveyrsurv_null;
run;

*Brier for null model;
data rottval_bs1null;
	set rottval_bsnull;
	keep pid fiveyrsurv_null;
	proc sort; by pid;
run;

proc sort data=rottval_bs3;
	by pid;
run;

*Merge the null model survival probabilities to brier score dataset from 
*earlier;
data rottval_bs4;
	merge rottval_bs3 rottval_bs1null;
	by pid;
run;

data rottval_bs5;
	set rottval_bs4;
	if cat=1 then contrib_null=(-fiveyrsurv_null)**2;
	if cat=2 then contrib_null=(1-fiveyrsurv_null)**2;
	if cat=3 then contrib_null=0;
	bs_null=contrib_null*weight;
	drop fiveyrsurv_null;
run;

*Estimate brier score for null model;
proc univariate data=rottval_bs5 noprint;
	var bs_null weight;
	output out=sumnullex sum=sbs_null sweight;
	proc print; 
run;

data sumnullex;
	retain sweight sbs_null;
	set sumnullex;
	sweight=left(sweight);
run;

*Calculate Scaled Brier;
data scaledb;
	merge sumsex sumnullex;
	by sweight;
	null_brier = (1/sweight)*sbs_null;
	scaled_b = 1-(brier/null_brier);
	title 'External Brier score and Scaled Brier';
	proc print;
run;


title ' ';

******Start of block;
******Extended model with PGR included;
*External validation - Calculate Brier score;

*Calculate weights using Kaplan-Meier;
proc lifetest data=gbsg_val method=pl atrisk outsurv=outkm_exta noprint;
        time survtime*status(1);
run;

data rott_b_exa;
	set gbsg_val;
	if survtime<=4.99 and status=1 then cat=1;
	if survtime>4.99 then cat=2;
	if survtime<=4.99 and status=0 then cat=3;
	time=survtime;
run;

proc phreg data=rottxa;
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat pgr pgr1/ ties=efron rl;
	*external;
	baseline covariates=rott_b_exa out=rottval_bsa timelist=5 
		survival=fiveyrsurv;
run;

*code up the fixed time of 5 years;
data rottval_bs1a(rename=(time=survtime));
	set rottval_bsa;
	time1=time;
	if time1>5 then time1=5;
	drop survtime;
	proc sort; by time1;
run;

data outkm_ext1a(rename=(survtime=time1));
	set outkm_exta;
	if survtime>5 then delete;
	weight=1/survival;
	keep survtime weight;
run;

data outkm_ext2a;
	set outkm_ext1a;
	s=lag(time1);
	if time1=s then delete;
	drop s;
run;

data rottval_bs2a;
	merge rottval_bs1a outkm_ext1a;
	by time1;
	if pid=. then delete;
run;

data rottval_bs3a;
	set rottval_bs2a;
	retain _weight;
	if not missing(weight) then _weight=weight;
	else weight=_weight;
	if cat=3 then weight=0;
	if time1=0 then delete;
	if cat=1 then contrib=(-fiveyrsurv)**2;
	if cat=2 then contrib=(1-fiveyrsurv)**2;
	if cat=3 then contrib=0;
	bs=contrib*weight;
	drop _weight;
run;

*Estimate brier score;
proc univariate data=rottval_bs3a noprint;
	var bs weight;
	output out=sumsexa sum=sbs sweight;
run;

data sumsexa;
	set sumsexa;
	brier = (1/sweight)*sbs;
	sweight=left(sweight);
	title 'External Brier score for PGR';
	proc print;
run;

******* End of block;

****** Scaled Brier for External validation of model with PGR;
*Scaled Brier = 1 - (model Brier score/null model Brier score), where null 
*model Brier score is null cox; 
*100% is perfect, <0 is useless, higher better, harmful models <0;

proc sort data=rottval_bs3a;
	by pid;
run;

data rottval_bs4a;
	merge rottval_bs3a rottval_bs1null;
	by pid;
run;

data rottval_bs5a;
	set rottval_bs4a;
	if cat=1 then contrib_null=(-fiveyrsurv_null)**2;
	if cat=2 then contrib_null=(1-fiveyrsurv_null)**2;
	if cat=3 then contrib_null=0;
	bs_null=contrib_null*weight;
	drop fiveyrsurv_null;
run;

proc univariate data=rottval_bs5a noprint;
	var bs_null weight;
	output out=sumnullexa sum=sbs_null sweight;
	proc print; 
run;

data sumnullexa;
	set sumnullexa;
	sweight=left(sweight);
run;

data scaledbriera;
	merge sumsexa sumnullexa;
	by sweight;
	null_brier = (1/sweight)*sbs_null;
	scaledb = 1-(brier/null_brier);
	title 'External Brier score and Scaled Brier with PGR';
	proc print;
run;


********* End of overall performance for added marker;


*****Appendix 5 - what to do if development dataset or baseline hazard are not 
available;

*Assume original model development paper has kaplan-meier plots for risk groups;
*Define Risk Groups for K-M - Creates plot in appendix 5, fig s3;

*** Define Risk Groups for K-M;
* Mean of XB in dev data is 0.65232758;
proc univariate data=rottxa;
	var xb;
run;
*Centre PI on average risk by subtracting mean;
data rottriska;
	set rottxa;
	xbc = xb - 0.65232758;
run;

proc univariate data=rottriska;
	var xbc;
run;

* generate risk groups;
%macro generate_percentiles(ptile1,ptile2,ptile3,mean); 
/* Output desired percentile values */                         
proc univariate data=rottxa noprint;                                               
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

data gbsg_valrd2a_1;                                                             
   merge gbsg_valrd2a_ test1a;        
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

*Plot kaplan-meier of risk groups in development dataset and validation dataset;
*mean (SD) of xb centered is 0 (0.64);
proc univariate data=rottxa_1;
	var xbc;
run;

proc lifetest data=rottxa_1 method=pl plots=(s, ls, lls) outsurv=outkma noprint;
        time survtime*status(0);
        strata riskgp/test=(logrank);	
run;

data outkm1a; 
	set outkma; 
	if _censor_=1 then delete; 
	num+1;
	keep riskgp survival survtime num;
run;

*- plot kaplan-meier in validation dataset;
*mean (SD) of xb centered is 0.23620 (0.5046);
proc univariate data=gbsg_valrd2a_1;
	var xbc;
run;

*what is the percentage in each risk group - pasted output below;
proc freq data=gbsg_valrd2a_1;
	table riskgp;
run;
/*
riskgp Frequency Percent Cumulative
Frequency Cumulative
Percent 
0 60 8.75 60 8.75 
1 144 20.99 204 29.74 
2 249 36.30 453 66.03 
3 233 33.97 686 100.00 
*/

*Validation dataset survival by risk group;
proc lifetest data=gbsg_valrd2a_1 method=pl plots=(s, ls, lls) outsurv=outkmva;
        time survtime*status(0);
        strata riskgp/test=(logrank);	
run;

data outkmv1a(rename=(riskgp=riskgpv survival=survivalv survtime=survtimev)); 
	set outkmva; 
	if _censor_=1 then delete; 
	num+1;
	keep riskgp survival survtime num;
run;

proc format;
	value dev 1='Development dataset';
	value val 1='Validation dataset';
run;

data overlaykma;
	merge outkm1a outkmv1a;
	by num;
	format gp1 dev.;
	format gpv1 val.;
	*trick to fix legend in plot;
	if riskgp=0 then gp1=1;
	if riskgp=1 then gp2=1;
	if riskgp=2 then gp3=1;
	if riskgp=3 then gp4=1;

	if riskgpv=0 then gpv1=1;
	if riskgpv=1 then gpv2=1;
	if riskgpv=2 then gpv3=1;
	if riskgpv=3 then gpv4=1;
run;

*CREATE KAPLAN MEIER PLOT OF RISK GROUPS FOR DEVELOPMENT AND VALIDATION DATASETS;
title;
footnote;
*- this allows editing of the .sge file!;
ods listing sge=on style=printer image_dpi=300 gpath='C:';
ods graphics on / reset=all noborder outputfmt=tiff imagename="Overlay KM PGR"  
	antialias=off;

proc sgplot data=overlaykma;
	yaxis values=(0 to 1 by 0.2) label="Recurrence-free survival probability";
	xaxis values=(0 to 5 by 1) label="Years";
	keylegend "all" "orig" / across=1 location=inside position=topright;
	step y=survival x=survtime / lineattrs=(color=black thickness=2 
		pattern=solid) name="all" group=gp1 nomissinggroup;
	step y=survival x=survtime / lineattrs=(color=black thickness=2 
		pattern=solid) group=gp2 nomissinggroup;
	step y=survival x=survtime / lineattrs=(color=black thickness=2 
		pattern=solid) group=gp3 nomissinggroup;
	step y=survival x=survtime / lineattrs=(color=black thickness=2 
		pattern=solid) group=gp4 nomissinggroup;
	step y=survivalv x=survtimev / lineattrs=(color=black thickness=2 
		pattern=dash) name="orig" group=gpv1 nomissinggroup;
	step y=survivalv x=survtimev / lineattrs=(color=black thickness=2 
		pattern=dash) group=gpv2 nomissinggroup;
	step y=survivalv x=survtimev / lineattrs=(color=black thickness=2 
		pattern=dash) group=gpv3 nomissinggroup;
	step y=survivalv x=survtimev / lineattrs=(color=black thickness=2 
		pattern=dash) group=gpv4 nomissinggroup;
run;

ods graphics off;


************** Conduct DCA analysis;
*- see %stdca macro - get from https://www.mskcc.org/departments/epidemiology-
*biostatistics/biostatistics/decision-curve-analysis;

*let's add PGR as a new marker and calculate for basic and extended models in 
external dataset;
proc phreg data=ffpgr;
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat / ties=efron rl;
	baseline out=origext covariates=gbsg_val survival=fiveyr timelist=5;
run;

proc phreg data=ffpgr;
	class sizecat (ref='<=20') nodescat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodescat gradecat pgr pgr1 / ties=efron rl;
	baseline out=newext covariates=origext survival=fiveyr_pgr timelist=5;
run;


data newext2;
	set newext;
	risk_orig=1-fiveyr;
	risk_new=1-fiveyr_pgr;
run;


%put _user_;
goptions reset = all;

symbol1 i=join c=green;
symbol2 i=join c=red;
symbol3 i=join c=blue;
symbol4 i=join c=darkred;
symbol5 i=join c=gray;

%STDCA(data=NEWEXT2, out=survivalmult, outcome=STATUS, ttoutcome=SURVTIME1, 
timepoint=5, predictors=RISK_ORIG);
%STDCA(data=NEWEXT2, out=survivalmult_new, outcome=STATUS, ttoutcome=SURVTIME1, 
timepoint=5, predictors=RISK_NEW);

*Sort by threshold variable;
proc sort data=survivalmult out=kmsort;
	by threshold;
run;

*rename the variables so that we know they are the kaplan meier estimates;
data kmsort; set kmsort(rename=(risk_orig=kmmodel all=kmall));
	label kmmodel="Original model: Prob at 5 years";
	label kmall="Treat All";
RUN;

*Sort by threshold variable;
proc sort data=survivalmult_new out=crsort;
	by threshold;
run;

*Rename the variables so that we know they are the competing model estimates;
data crsort; set crsort(rename=(risk_new=crmodel all=crall));
	label crmodel="Original + PGR: Prob at 5 years";
run;

*merge nb data for original model plus model with pgr using threshold 
*probabilities;
data crsort;
	merge kmsort crsort;
	by threshold;
run;
*******************************************************************************;	

data stratos.extnb;
	set crsort;
run;

******  Optional Block to create nicer DCA curves ************;

title;
footnote;
ods listing sge=on style=printer image_dpi=300 gpath='C:';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="dca with spline PGR ext" antialias=off;

*create decision curve with treat none, treat all, basic model, extended model;
proc sgplot data=crsort;
	yaxis values=(-0.05 to 0.50 by 0.05) label="Net Benefit";
	xaxis values=(0.0 to 1 by 0.1) label="Threshold Probability";
	keylegend "all" "orig" "new" "none" / down=4 position=bottom;
	series y=kmall x=threshold / lineattrs=(color=black thickness=2 
		pattern=solid) name="all" legendlabel="Treat All";
	series y=kmmodel x=threshold / lineattrs=(color=green thickness=2 
		pattern=solid) name="orig" 
		legendlabel="Original model: Pr(Recurrence/death) at 5 years";
	series y=crmodel x=threshold / lineattrs=(color=blue thickness=2 
		pattern=solid) name="new" 
		legendlabel="Original model + PGR: Pr(Recurrence/death) at 5 years";
	series y=none x=threshold / lineattrs=(color=red thickness=2 
		pattern=solid)  name="none" legendlabel="Treat None";
RUN;

ods graphics off;

******  End Optional Block ************;
