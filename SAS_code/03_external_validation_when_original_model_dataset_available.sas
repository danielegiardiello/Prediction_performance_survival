*******************************************************************************;	
* Program: 			STRATOS External validation when original model dataset    ;
*					available.sas  											   ;
* Author: 			David McLernon											   ;
* Date: 			4th Aug 2022											   ;
* Purpose: 			This is SAS code for the STRATOS paper on validation of	   ;
*					survival risk prediction models. This programme covers 	   ;
*					external validation when the original dataset used to      ;
*					develop the model is available							   ;
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
** size of  <=20, 21-50, >50
**  The size variable is already categorized in this way, in both data sets
*/

*** descriptive statistics for age and PGR;
proc univariate data=rotterdam;
	var age pgr;
run;

proc format;  * make them print in a nice order;
    value sizef 0 = "<=20"
                1 = "20-50"
                2 = ">50";
    value gradef 0 = "1-2"
                 1 = "3";

data r1(rename=(pr=pgr nod=nodes)); set rotterdam;
    format sizecat sizef. gradecat gradef.;

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
	if (nodes > 19) then nod = 19; else nod = nodes;
	drop nodes;
run;

* Descriptive statistics at baseline for Table 1;
proc freq data=r1;
    table gradecat sizecat;
run;

*Administrative censor at 5 years since this is our prediction horizon; 
data ffpgr;
	set r1;
	*administrative censoring at 5 years;
	if survtime > 5 then status=0;
	if survtime > 5 then survtime=5;
	* rcs for pgr and nodes - see development program for knot calculations;
	%rcspline(pgr, 0, 41, 486);
	%rcspline(nodes, 0, 1, 9);
run;

************** Import external validation dataset ****************************;
*code up predictors;
data gbsg_va; set gbsg;
    format sizecat sizef. gradecat gradef.;
	*Winzorise PGR to 99th percentile to deal with very large influential 
		values;
	if pgr>1360 then pgr=1360;
	%rcspline(pgr, 0, 41, 486);
	if (nodes > 19) then nodes=19;
	%rcspline(nodes, 0, 1, 9);
	if size le 20 then sizecat=0;
	if 21 <= size <= 50 then sizecat=1;
	if size >50 then sizecat=2;
	gradecat = grade-1;
	if grade in (1,2) then gradecat=0;
	if grade = 3 then gradecat=1;
    survtime = rfstime/365.25;  * days to years;
	label survtime="Survival time (years)";
	keep pid sizecat nodes nodes1 gradecat status age pgr pgr1 survtime;
run;


* Plot the overall Kaplan-Meier, as motivation for using a 5 year cut-off;
title "Kaplan-Meier plot in external dataset";
footnote;
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=png /*width=4in*/ 
imagename="ext km" antialias=off/*antialiasmax=*/ ;

proc lifetest data=gbsg_va method=pl plots=(s(atrisk(atrisktickonly outside)
			=0,1,2,3,4,5,6,7 cb=hw nocensor name=Survival)) outsurv=outkm;
        time survtime*status(0);
		ods exclude ProductLimitEstimates;   * suppress the long table;
		title "Overall Kaplan-Meier";
run;

ods graphics off;

title;


*** Descriptive statistics;
proc freq data=gbsg_va;
	tables sizecat gradecat;
run;

proc univariate data=gbsg_va;
	var nodes pgr;
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


**First, resample 500 versions of the external dataset with replacement to 
allow calculation of 95% CI;

sasfile gbsg_val load;/*a way of loading the dataset into ram - speeds it up*/

proc surveyselect data=gbsg_val out=outboot
seed=4817 /* can enter 0 for it to select a random seed but remember to type 
			it in here from the output otherwise cannot replicate results */
method=urs /* unrestricted random sampling - simple random sampling */
samprate=1 /* can accept proportions or percentages but we want n to be size of 
			original database so =1 (or 100) */
outhits /* with replacement */
rep=500; /* number of bootstrap samples */
run;

sasfile gbsg_val close; /* closes frees ram buffers when done */

ods listing close; /* turns off ods listing so no printing of all output. 
					Better than using NOPRINT as it doesn't allow storage of 
					data in output dataset at end */



******************EXTERNAL VALIDATION *********************;

*******************Simple model first (without PGR)***************************;

* Fit the model to the development data, and calculate linear predictor of 
* this existing model for patients in the external dataset;

proc phreg data=ffpgr;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat gradecat nodes nodes1 / ties=efron rl;
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
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1/ ties=efron rl;
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
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1/ ties=efron rl;
	*the following statements call the stored models from earlier and 
	calculate the concordance statistics for model with and without pgr;
	roc 'npi' source=simpmodel;
	roc 'npi + pgr' source=simpmodelpgr;
run;

*Uno's C - Need tau to equal event time of interest;
title 'External Unos C';
proc phreg data=gbsg_val concordance=uno(se seed=8754 iter=50) tau=5;
	CLASS sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1/ ties=efron rl;
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
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1/ ties=efron rl ;
	roc 'npi' source=simpmodel;
	roc 'npi + pgr' source=simpmodelpgr;
run;

******* Optional Block ****** ;
* Rerun Uno and save dataset to plot the Time-dependent AUC yourself;
proc phreg data=gbsg_val rocoptions(method=ipcw (cl seed=134) outauc=ipcwauc2);
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1/ ties=efron rl 
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
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat/ ties=efron rl;
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

*For plot of O/E at all observed times up to the event or censoring 
time of each patient (Not in the paper or appendix);
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
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="calbration plot o/e basic" antialias=off;

proc sgplot data=sumoe2 noautolegend ;
xaxis       	label="Month of follow-up" values=(0 to 60 by 6);
yaxis        	label="Number of observed events / Number of expected events" 
	Values=(0 to 1.6 by 0.1) /*offsetmin=0.19*/;
y2axis;
  keylegend "ref" "cal" "conf" / location=inside position=bottomright down=3;
  series x=month y=est / lineattrs=(color=red thickness=2) name="cal" 
	legendlabel="O/E based on Poisson approach";
  refline r / axis=y  lineattrs=(color=black thickness=1)
	name="ref" legendlabel="Ideal calibration (O/E=1)";
  band x=month lower=low upper=upp / nofill lineattrs=(color=black 
	pattern=mediumdash thickness=2) noextend outline name="conf" 
		legendlabel="95% confidence interval";
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
* Code LP as a restricted cubic spline with 4 knots (3df);
* First calculate the 5th, 35th, 65th, 95th percentiles for knots;
proc univariate data = data2 noprint;
	var lp;
	output out=knots pctlpre=p_lp pctlpts= /*10 50 90*/ 5 35 65 95;
run;

proc print data=knots; run;

/* here we find the following values:
0.29394 0.88091 1.30518 2.00782 

*/

data test;
	set data2;
	%rcspline(lp, 0.29394, 0.88091, 1.30518, 2.00782);
run;

* generate dataset to allow us to get expected rate over all PI;
data genpi;
	lp=-0.025;
	do i=1 to 101;
		lp + 0.025;
		output;
		end;
run;

*Get the spline formulas from the log window after running test above and paste 
below;
data genpi1;
	set genpi;
 	_kd_= (2.00782 - 0.29394)**.666666666666 ;
	lp1=max((lp-0.29394)/_kd_,0)**3
		+((1.30518-0.29394)*max((lp-2.00782)/_kd_,0)**3
		-(2.00782-0.29394)*max((lp-1.30518)/_kd_,0)**3)/(2.00782-1.30518);
	lp2=max((lp-0.88091)/_kd_,0)**3
		+((1.30518-0.88091)*max((lp-2.00782)/_kd_,0)**3
		-(2.00782-0.88091)*max((lp-1.30518)/_kd_,0)**3)/(2.00782-1.30518);
	*Using a dummy p of 1.0 gives the relative event rate, as compared to the 
	Rotterdam data set as a whole.  (The expecteds add to the total number of 
	deaths, in the reference fit).   So in this case use time=1 as the dummy;
	logp=0;
	drop i _kd_;
run;

data test1;
	set test genpi1;
run;

*Fit Poisson model with spline terms and save predictions;
proc genmod data=test1;
	model status=lp lp1 lp2 /offset=logp dist=poisson link=log;
	output out=test2 pred=predobs lower=low upper=upp;
run;

data test3;
	set test2;
	where pid=.;
	r=1;
run;

*Produce plot in Appendix Fig S2;
title;
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="calibration plot o over e v PI basic" antialias=off;

proc sgplot data=test3 noautolegend ;
xaxis       	label="Prognostic Index" values=(0 to 2.5 by 0.25);
yaxis        	label="Number of observed events / Number of expected events" 
	Values=(0.5 to 3 by 0.5) type=log logbase=10 logstyle=linear ;
  keylegend "ref" "cal" "conf" / location=inside position=top down=3;
  series x=lp y=predobs / lineattrs=(color=red thickness=2) name="cal" 
	legendlabel="O/E based on Poisson approach";
  refline r / axis=y  lineattrs=(color=black thickness=1)
	name="ref" legendlabel="Ideal calibration (O/E=1)";
  band x=lp lower=low upper=upp / nofill lineattrs=(color=black 
	pattern=mediumdash thickness=2) noextend outline name="conf" 
		legendlabel="95% confidence interval";
run;

ods graphics off;



//*Optional block to plot cumulative hazard from Cox model v cumulative hazard
from Poisson model - Not in paper or appendix*//;

*Sort and code diagonal ref line;
data test4;
	set test2;
	if pid=1 then diag1=0;
	if pid=1 then diag2=0;
	if pid=6 then diag1=3;
	if pid=6 then diag2=3;
	if pid=. then delete;
	proc sort; by p;
run;

*Select unique predictions;
data test5;
	set test4;
	by p;
	if first.p;
run;

*Plot flexible calibration curve;
title;
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="Calibration plot Poisson" ANTIALIAS=OFF;

proc sgplot data=test5 noautolegend aspect=1;
	xaxis       	label="Cumulative hazard from Cox model" 
		values=(0 to 2.4 by 0.20);
	yaxis        	label="Cumulative hazard from Poisson model" 
		values=(0 to 2.4 by 0.20);
	keylegend "ref" "cal" "conf" / location=inside position=bottomright;
  	loess x=p y=predobs / nomarkers lineattrs=(color=red thickness=2) name="cal" 
		legendlabel="Calibration curve based on Poisson approach";
  	loess x=p y=low / nomarkers lineattrs=(color=black 
		pattern=mediumdash thickness=2) name="conf" 
		legendlabel="95% confidence interval";
  	loess x=p y=upp / nomarkers lineattrs=(color=black 
		pattern=mediumdash thickness=2);
	reg x=diag1 y=diag2 / lineattrs=(pattern=shortdash) markerattrs=(size=1px 
		symbol=circlefilled) name="ref" legendlabel="Ideal calibration";
run;

ods graphics off;

*End of optional block;

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
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat/ ties=efron rl;
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
	lclprop = prop * exp(-1.96*sqrt(1/285));
	uclprop = prop * exp(+1.96*sqrt(1/285));
run;

title 'Calibration-in-large using ratio of Kaplan-Meier & Avg predicted risk';
proc print data=mean2;
run;


**Weak - Calibration slope;

proc phreg data=ffpgr;
	class sizecat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat/ ties=efron rl;
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
p_cll10 p_cll50 p_cll90 
 -0.98264 -0.44824 0.32253 
*/

*Use Frank Harrell's RCSPLINE macro for calculating the spline terms;
data predhar1;
	set predhar;
	* NOTE: you will need to take the spline information from the log 
	window and edit gridcox2 below;
	%rcspline(cll, -0.98264, -0.44824, 0.32253);
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
	 _kd_= (0.32253 - -0.98264)**.666666666666 ;
	cll1=max((cll--0.98264)/_kd_,0)**3
		+((-0.44824--0.98264)*max((cll-0.32253)/_kd_,0)**3
		-(0.32253--0.98264)*max((cll--0.44824)/_kd_,0)**3)/(0.32253--0.44824);
run;

*Calibration for predictions of 5-year survival probabilities;
ods graphics on;
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
	if _name_='p_9' then diag1=1;
	if _name_='p_9' then diag2=1;
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
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="Calibration plot Ext Harrell with RCS" ANTIALIAS=OFF;

proc sgplot data=logiplot2 noautolegend aspect=1;
	xaxis       	label="Predicted risk from developed model" 
		values=(0 to 1 by 0.20);
	yaxis        	label="Predicted risk from refitted model" 
		values=(0 to 1 by 0.20);
	keylegend "ref" "cal" "conf" "dens" / location=inside position=topleft
		valueattrs=(size=8) down=4 noborder;
	y2axis			label=" " values=(0 to 20 by 20) display=none;
  	series x=prob y=obsprob / lineattrs=(color=red thickness=2) name="cal" 
		legendlabel="Calibration curve based on secondary Cox model";
  	band x=prob lower=obslow upper=obsupp / nofill lineattrs=(color=black 
		pattern=mediumdash thickness=2) noextend outline name="conf" 
		legendlabel="95% confidence interval";
	reg x=diag1 y=diag2 / lineattrs=(pattern=shortdash) markerattrs=(size=1px 
		symbol=circlefilled) name="ref" legendlabel="Ideal calibration";
  	series x=_value_ y=_density_ / lineattrs=(color=black thickness=1 pattern=1) 
		y2axis name="dens" legendlabel="Density function of predicted risk";
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

data icis;
	set icis;
	n=1;
run;

title;



******** BOOTSTRAP 95% CI FOR THE CALIBRATION METRICS ********;

proc phreg data=ffpgr;
	class sizecat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat/ ties=efron rl;
	baseline out=gbsg_valhbt covariates=outboot survival=predprob xbeta=xb timelist=5;
run;

data predharbt;
	set gbsg_valhbt;
	where survtime1 le survtime; 
	*probability of event;
	prob=1-predprob;
	*Complementary log-log of predicted risk;
	cll=log(-log(1-prob));
	proc sort; by replicate pid survtime;
run;

proc univariate data = predharbt noprint;
	by replicate;
	var cll;
	output out=knots pctlpre=p_cll pctlpts= 10 50 90;
run;

data predharbt_;
	merge predharbt knots;
	by replicate;
run;


*This macro runs the RCS on each replicate and appends all 500 runs to one dataset;
%MACRO SPL;
%do i=1 %to 500;
data predharbt_1;
	set predharbt_;
	where replicate=&i;
	%rcspline(cll, p_cll10, p_cll50, p_cll90);
run;

*duplicate dataset so baseline statement works;
data predharbt_rep;
	set predharbt_1;
run;

proc phreg data=predharbt_1 noprint;
	model survtime1*status(0)=cll cll1/ ties=efron rl;
	baseline out=predharbt_rep1 covariates=predharbt_rep survival=predprob2 
		timelist=5;
run;

data ici_;
	set predharbt_rep1;
	*observed prob death;
	prob2=1-predprob2;
	diff=abs(prob - prob2);
run;

proc univariate data=ici_ noprint;
	var diff;
	output out=icis_ mean=ici median=e50 p90=e90;
run;

data icis_1;
	set icis_;
	replicate=&i;
run;

proc append base=icis_2 data=icis_1 force;
run;


%end;
%mend;
%spl;

proc univariate data=icis_2 noprint;
	var ici e50 e90;
	output out = conf_met pctlpts=2.5 97.5 pctlpre= ici_ e50_ e90_ 
		pctlname=lower95 upper95;
run;

data conf_met1;
	set conf_met;
	n=1;
run;

data calmetrics;
 	retain ici ici_lower95 ici_upper95 e50 e50_lower95 e50_upper95 e90 
		e90_lower95 e90_upper95;
	merge icis conf_met1;
	by n;
	drop n;
	title 'External calibration metrics with 95% CI';
	proc print;
run;





***EXTENDED MODEL WITH PGR INCLUDED;
**Using Crowson's Poisson approach;

*- CUMHAZ command provides the expected number of events and XBETA provides 
	linear predictor for original model applied to external cohort;
proc phreg data=ffpgr;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1/ ties=efron rl;
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

*For plot of O/E for all times up to the event or censoring 
time of each patient (Not in paper or appendix);
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
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="calbration plot o/e pgr" antialias=off;

proc sgplot data=sumoe4 noautolegend ;
xaxis       	label="Month of follow-up" values=(0 to 60 by 6);
yaxis        	label="Number of observed events / Number of expected events" 
	values=(0 to 1.6 by 0.1) /*offsetmin=0.19*/;
y2axis;
  keylegend "ref" "cal" "conf" / location=inside position=bottomright down=3;
  series x=month y=est / lineattrs=(color=red thickness=2) name="cal" 
	legendlabel="O/E based on Poisson approach";
  refline r / axis=y  lineattrs=(color=black thickness=1)
	name="ref" legendlabel="Ideal calibration (O/E=1)";
  band x=month lower=low upper=upp / nofill lineattrs=(color=black 
	pattern=mediumdash thickness=2) noextend outline name="conf" 
		legendlabel="95% confidence interval";
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
* Code LP as a restricted cubic spline with 4 knots (3df);
* First calculate the 5th, 35th, 65th and 95th percentiles for knots;
proc univariate data = data2 noprint;
	var lp;
	output out=knots pctlpre=p_lp pctlpts= /*10 50 90*/ 5 35 65 95;
run;

proc print data=knots; run;

/* here we find the following values which we enter into the rcspline 
call below:
p_lp5 		p_lp35 		p_lp65 		p_lp95 
0.090982 	0.67938 	1.11916 	1.86255 
*/

data test;
	set data2;
	%rcspline(lp, 0.090982, 0.67938, 1.11916, 1.86255);
run;

* generate dataset to allow us to get expected rate over all PI;
data genpi;
	lp=-0.5;
	do i=1 to 121;
		lp + 0.025;
		output;
		end;
run;

*Get the spline formulas from the log window after running test above and paste 
below;
data genpi1;
	set genpi;
 	_kd_= (1.86255 - 0.090982)**.666666666666 ;
	lp1=max((lp-0.090982)/_kd_,0)**3
		+((1.11916-0.090982)*max((lp-1.86255)/_kd_,0)**3
 		-(1.86255-0.090982)*max((lp-1.11916)/_kd_,0)**3)/(1.86255-1.11916);
	lp2=max((lp-0.67938)/_kd_,0)**3
		+((1.11916-0.67938)*max((lp-1.86255)/_kd_,0)**3
		-(1.86255-0.67938)*max((lp-1.11916)/_kd_,0)**3)/(1.86255-1.11916);
	*Using a dummy p of 1.0 gives the relative event rate, as compared to the 
	Rotterdam data set as a whole.  (The expecteds add to the total number of 
	deaths, in the reference fit).   So in this case use time=1 as the dummy;
	logp=0;
	drop i _kd_;
run;

data test1;
	set test genpi1;
run;

*Fit Poisson model with spline terms and save predictions;
proc genmod data=test1;
	model status=lp lp1 lp2 /offset=logp dist=poisson link=log;
	output out=test2 pred=predobs lower=low upper=upp;
run;

data test3;
	set test2;
	where pid=.;
	r=1;
run;

*Produce plot in Appendix Fig S4;
title;
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="calibration plot o over e v PI pgr" antialias=off;

proc sgplot data=test3 noautolegend ;
xaxis       	label="Prognostic Index" values=(-0.5 to 2.5 by 0.25);
yaxis        	label="Number of observed events / Number of expected events" 
	Values=(0 to 3.5 by 0.5) type=log logstyle=linear logbase=10;
  keylegend "ref" "cal" "conf" / location=inside position=top down=3;
  series x=lp y=predobs / lineattrs=(color=red thickness=2) name="cal" 
	legendlabel="O/E based on Poisson approach";
  refline r / axis=y  lineattrs=(color=black thickness=1)
	name="ref" legendlabel="Ideal calibration (O/E=1)";
  band x=lp lower=low upper=upp / nofill lineattrs=(color=black 
	pattern=mediumdash thickness=2) noextend outline name="conf" 
		legendlabel="95% confidence interval";
run;

ods graphics off;


//*Optional block to plot cumulative hazard from Cox model v cumulative hazard
from Poisson model - Not in paper or appendix*//;

data test4;
	set test2;
	if pid=1 then diag1=0;
	if pid=1 then diag2=0;
	if pid=6 then diag1=3;
	if pid=6 then diag2=3;
	proc sort; by p;
run;

data test5;
	set test4;
	by p;
	if first.p;
run;

*Flexible calibration plot;
title;
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="Calibration plot Poisson PGR" ANTIALIAS=OFF;

proc sgplot data=test5 noautolegend aspect=1;
	xaxis       	label="Cumulative hazard from Cox model" 
		values=(0 to 2.4 by 0.20);
	yaxis        	label="Cumulative hazard from Poisson model" 
		values=(0 to 2.4 by 0.20);
	keylegend "ref" "cal" "conf" / location=inside position=bottomright;
  	loess x=p y=predobs / nomarkers lineattrs=(color=red thickness=2) name="cal" 
		legendlabel="Calibration curve based on Poisson approach";
  	loess x=p y=low / nomarkers lineattrs=(color=black 
		pattern=mediumdash thickness=2) name="conf" 
		legendlabel="95% confidence interval";
  	loess x=p y=upp / nomarkers lineattrs=(color=black 
		pattern=mediumdash thickness=2);
	reg x=diag1 y=diag2 / lineattrs=(pattern=shortdash) markerattrs=(size=1px 
		symbol=circlefilled) name="ref" legendlabel="Ideal calibration";
run;

ods graphics off;

*** End of optional block;


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
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1/ ties=efron rl;
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
	lclprop = prop * exp(-1.96*sqrt(1/285));
	uclprop = prop * exp(+1.96*sqrt(1/285));
run;

title 'O/E using ratio of Kaplan-Meier and Avg predicted risk';
proc print data=mean2;
run;


*Weak assessment - Calibration slope;

proc phreg data=ffpgr;
	class sizecat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1/ ties=efron rl;
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

/* here we find the following values which you need to enter into rcspline call below:
P_CLL10 P_CLL50 P_CLL90 
-1.05612 -0.39421 0.40613   
*/

*Use Frank Harrell's RCSPLINE macro for calculating the spline terms;
data predhar1;
	set predhar;
	* NOTE: you will need to take the spline information from the log 
	window and edit gridcox2 below;
	%rcspline(cll, -1.05612, -0.39421, 0.40613);
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
	 _kd_= (0.40613 - -1.05612)**.666666666666 ;
	 cll1=max((cll--1.05612)/_kd_,0)**3
		+((-0.39421--1.05612)*max((cll-0.40613)/_kd_,0)**3 
		-(0.40613--1.05612)*max((cll--0.39421)/_kd_,0)**3)/(0.40613--0.39421);
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
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="Calbration plot pgr Ext Harrell with RCS" antialias=off;

proc sgplot data=logiplot2 noautolegend aspect=1;
	xaxis       	label="Predicted risk from developed model" 
		values=(0 to 1 by 0.20);
	yaxis        	label="Predicted risk from refitted model" 
		values=(0 to 1 by 0.20);
	y2axis			label=" " values=(0 to 20 by 20) display=none;
	keylegend "ref" "cal" "conf" "dens" / location=inside position=topleft
		valueattrs=(size=8) down=4 noborder;
  	series x=prob y=obsprob / lineattrs=(color=red thickness=2) name="cal" 
		legendlabel="Calibration curve based on secondary Cox model";
  	band x=prob lower=obslow upper=obsupp / nofill lineattrs=(color=black 
		pattern=mediumdash thickness=2) noextend outline name="conf" 
		legendlabel="95% confidence interval";
	reg x=diag1 y=diag2 / lineattrs=(pattern=shortdash) markerattrs=(size=1px 
		symbol=circlefilled) name="ref" legendlabel="Ideal calibration";
  	series x=_value_ y=_density_ / lineattrs=(color=black thickness=1 pattern=1) 
		y2axis name="dens" legendlabel="Density function of predicted risk";
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

data icis;
	set icis;
	n=1;
run;

title 'Calibration metrics: ICI, E50, E90';
proc print data=icis;
run;

title;




******** BOOTSTRAP 95% CI FOR THE CALIBRATION METRICS ********;

proc phreg data=ffpgr;
	class sizecat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron rl;
	baseline out=gbsg_valhbt covariates=outboot survival=predprob xbeta=xb timelist=5;
run;

data predharbt;
	set gbsg_valhbt;
	where survtime1 le survtime; 
	*probability of event;
	prob=1-predprob;
	*Complementary log-log of predicted risk;
	cll=log(-log(1-prob));
	proc sort; by replicate pid survtime;
run;

proc univariate data = predharbt noprint;
	by replicate;
	var cll;
	output out=knots pctlpre=p_cll pctlpts= 10 50 90;
run;

data predharbt_;
	merge predharbt knots;
	by replicate;
run;


*This macro runs the RCS on each replicate and appends all 500 runs to one dataset;
%MACRO SPL;
%do i=1 %to 500;
data predharbt_1;
	set predharbt_;
	where replicate=&i;
	%rcspline(cll, p_cll10, p_cll50, p_cll90);
run;

*duplicate dataset so baseline statement works;
data predharbt_rep;
	set predharbt_1;
run;

proc phreg data=predharbt_1 noprint;
	model survtime1*status(0)=cll cll1/ ties=efron rl;
	baseline out=predharbt_rep1 covariates=predharbt_rep survival=predprob2 
		timelist=5;
run;

data ici_;
	set predharbt_rep1;
	*observed prob death;
	prob2=1-predprob2;
	diff=abs(prob - prob2);
run;

proc univariate data=ici_ noprint;
	var diff;
	output out=icis_ mean=ici median=e50 p90=e90;
run;

data icis_1;
	set icis_;
	replicate=&i;
run;

proc append base=icis_4 data=icis_1 force;
run;


%end;
%mend;
%spl;

proc univariate data=icis_4 noprint;
	var ici e50 e90;
	output out = conf_met pctlpts=2.5 97.5 pctlpre= ici_ e50_ e90_ 
		pctlname=lower95 upper95;
run;

data conf_met1;
	set conf_met;
	n=1;
run;

data calmetrics;
 	retain ici ici_lower95 ici_upper95 e50 e50_lower95 e50_upper95 e90 
		e90_lower95 e90_upper95;
	merge icis conf_met1;
	by n;
	drop n;
	title 'External calibration metrics with 95% CI';
	proc print;
run;


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
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat  / ties=efron rl;
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
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1/ ties=efron rl;
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



************** Conduct DCA analysis;
*- see %stdca macro - get from https://www.mskcc.org/departments/epidemiology-
*biostatistics/biostatistics/decision-curve-analysis;

*let's add PGR as a new marker and calculate for basic and extended models in 
external dataset;
proc phreg data=ffpgr;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
	baseline out=origext covariates=gbsg_val survival=fiveyr timelist=5;
run;

proc phreg data=ffpgr;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron rl;
	baseline out=newext covariates=origext survival=fiveyr_pgr timelist=5;
run;


data newext2;
	set newext;
	risk_orig=1-fiveyr;
	risk_new=1-fiveyr_pgr;
run;

*count number tp and fp for cut pt of 0.23;
*testpos=1 is number used in Suppl 7 for worked example;
data tp;
	set newext2;
	if risk_orig gt 0.23 then testpos=1;
	else testpos=0;
	if risk_new gt 0.23 then testpos_pgr=1;
	else testpos_pgr=0;
run;

proc freq data=tp;
	table testpos testpos_pgr;
run;

*Kaplan-Meier probability at 5 years - used in Suppl 7 worked example;
proc lifetest data=tp method=pl timelist=5;
	where testpos=1;
    time survtime1*status(0);
run;

proc lifetest data=tp method=pl  timelist=5;
	where testpos_pgr=1;
    time survtime1*status(0);
run;




%put _user_;
goptions reset = all;

symbol1 i=join c=green;
symbol2 i=join c=red;
symbol3 i=join c=blue;
symbol4 i=join c=darkred;
symbol5 i=join c=gray;

*Use the %STDCA macro to calculate net benefit and plot decision curves;
*Note that you can add smooth=yes to get the smoothed curves in the paper;
*Smoothing reduces the visual impact of random noise;
*However, using this option only provides net benefit estimates to 2 decimals;
*rather than 4, so you may wish to run with and without smoothing;
%STDCA(data=NEWEXT2, out=survivalmult, outcome=STATUS, ttoutcome=SURVTIME1, 
	timepoint=5, predictors=RISK_ORIG, xby=0.02, smooth=yes);
%STDCA(data=NEWEXT2, out=survivalmult_new, outcome=STATUS, ttoutcome=SURVTIME1, 
	timepoint=5, predictors=RISK_NEW, xby=0.02, smooth=yes);

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
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="dca with spline PGR ext" antialias=off;

*create decision curve with treat none, treat all, basic model, extended model;
proc sgplot data=crsort;
	yaxis values=(-0.05 to 0.50 by 0.05) label="Net Benefit";
	xaxis values=(0.0 to 1 by 0.1) label="Threshold Probability";
	keylegend "all" "orig" "new" "none" / location=inside down=4 position=topright;
	series y=kmall x=threshold / lineattrs=(color=black thickness=2 
		pattern=mediumdash) name="all" legendlabel="Treat All";
	series y=kmmodel x=threshold / lineattrs=(color=green thickness=2 
		pattern=solid) name="orig" 
		legendlabel="Original model" smoothconnect;
	series y=crmodel x=threshold / lineattrs=(color=blue thickness=2 
		pattern=longdash) name="new" 
		legendlabel="Original model + PGR" smoothconnect;
	series y=none x=threshold / lineattrs=(color=red thickness=2 
		pattern=shortdash)  name="none" legendlabel="Treat None";
RUN;

ods graphics off;

******  End Optional Block ************;
