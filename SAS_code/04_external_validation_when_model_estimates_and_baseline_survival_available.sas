*******************************************************************************;	
* Program: 			STRATOS External validation when model estimates and       ;
*					baseline survival available.sas							   ;
*					available.sas  											   ;
* Author: 			David McLernon											   ;
* Date: 			24th Jan 2022											   ;
* Purpose: 			This is SAS code for the STRATOS paper on validation of	   ;
*					survival risk prediction models. This programme covers:    ;
*					1. external validation when the original dataset used to   ;
*					develop the model is not available but the PI and baseline ;
*					survival are available									   ;
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

****** Read in the German datset (external);
proc import out= gbsg 
            datafile= "c:\users\sme544\documents\stratos\gbsg.csv" 
            dbms=csv replace;
run;

*** Estimating Baseline Survival Function under PH;
** Here S0 at 5 years =  0.80153;
** For model with PGR added, S0 = 0.75852;
** S(t)=S0(t)**exp(ß1X1+ß2X2+ß3X3+...+ßpXp)=S0(t)**exp(PI);

*code up predictors and apply the original model;
data l2gbsg_va;
	set gbsg;
	if size le 20 then sizecat=0;
	if 21 <= size <= 50 then sizecat=1;
	if size >50 then sizecat=2;
    nodescat = 1* (1 <= nodes <=3) + 2*(nodes > 3);
	gradecat = grade-1;
	if grade in (1,2) then gradecat=0;
	if grade = 3 then gradecat=1;
    survtime = rfstime/365.25;  * days to years;
	*Winzorise PGR to the 99th percentile to deal with very large influential 
	values;
	if pgr>1360 then pgr=1360;
	*use formula to code up spline term for PGR;
	PGR1=max(PGR/61.81,0)**3
		+(41*max((PGR-486)/61.81,0)**3 -486*max((PGR-41)/61.81,0)**3)/445;
	*now calculate the prob of survival to 5 years;
	sizepi = 0;
	if sizecat=1 then sizepi=0.383;
	if sizecat=2 then sizepi=0.664;
	nodepi = 0;
	if nodescat=1 then nodepi=0.360;
	if nodescat=2 then nodepi=1.063;
	gradepi = 0;
	if gradecat=1 then gradepi=0.375;
	pi = sizepi + nodepi + gradepi;
	psurv = 0.80153**exp(pi);
	prob = 1- psurv;
	*now calculate the prob of survival to 5 years for model with pgr included;
	sizepi2 = 0;
	if sizecat=1 then sizepi2=0.362;
	if sizecat=2 then sizepi2=0.641;
	nodepi2 = 0;
	if nodescat=1 then nodepi2=0.381;
	if nodescat=2 then nodepi2=1.059;
	gradepi2 = 0;
	if gradecat=1 then gradepi2=0.317;
	pi2 = sizepi2 + nodepi2 + gradepi2 - 0.003*pgr + 0.013*pgr1;
	psurv2 = 0.75852**exp(pi2);
	prob2 = 1-psurv2;
	keep pid sizecat nodescat gradecat age pgr pgr1 survtime status
	sizepi nodepi gradepi pi psurv prob sizepi2 nodepi2 gradepi2 pi2 psurv2 
		prob2;
run;

*** Descriptive statistics;
proc freq data=l2gbsg_va;
	tables sizecat nodescat gradecat status;
run;

proc univariate data=l2gbsg_va;
	var age pgr;
run;

*use reverse Kaplan-Meier to obtain median follow-up time;
proc lifetest data=l2gbsg_va method=pl atrisk;
        time survtime*status(1);
run;

*Administrative censor at 5 years as development cohort; 
data l2gbsg_val;
	set l2gbsg_va;
	*administrative censoring at 5 years;
	if survtime > 5 then status=0;
	if survtime > 5 then survtime=5;
	*- keep a copy of survtime for calibration step;
	survtime1=survtime;

	*CLL is the Outcome-Martingale residual, or expected number of events;
	CLL = LOG(-LOG(PSURV));
	CLL_ = LOG(-LOG(PSURV2));
RUN;


**First, resample 500 versions of the external dataset with replacement to 
allow calculation of 95% CI;

sasfile l2gbsg_val load;/*a way of loading the dataset into ram - speeds it up*/

proc surveyselect data=l2gbsg_val out=outboot
seed=4817 /* can enter 0 for it to select a random seed but remember to type 
			it in here from the output otherwise cannot replicate results */
method=urs /* unrestricted random sampling - simple random sampling */
samprate=1 /* can accept proportions or percentages but we want n to be size of 
			original database so =1 (or 100) */
outhits /* with replacement */
rep=500; /* number of bootstrap samples */
run;

sasfile l2gbsg_val close; /* closes frees ram buffers when done */

ods listing close; /* turns off ods listing so no printing of all output. 
					Better than using NOPRINT as it doesn't allow storage of 
					data in output dataset at end */


*****Assessment of discrimination in external dataset;

*Calculate weights first with inverse Kaplan-Meier;
proc lifetest data=l2gbsg_val method=km atrisk outsurv=l2outkm_ext noprint;
        time survtime*status(1);
run;

data l2rott_b_ex;
	set l2gbsg_val;
	proc sort; by survtime;
run;

data l2outkm_ext1;
	set l2outkm_ext;
	weight=1/survival;
	keep survtime survival weight;
run;

data l2rottval_bs2;
	merge l2rott_b_ex l2outkm_ext1;
	by survtime;
	if pid=. then delete;
run;

data outkm_ext_uno;
	set l2rottval_bs2;
	retain _survival;
	if not missing(survival) then _survival=survival;
	else survival=_survival;
	w=survival;
	wsq=w*w;
	drop _survival;
run;

* Get KM censor probability at time T 4.96 years (chose 4.96 instead of 5 
because no one left after 5 years so won't work - took time of last event);
data ctime1;
	set l2outkm_ext;
	if survtime gt 4.96 then delete;
run;

*get last prob at time t;
data ctime2(rename=(survival=probc));
	if 0 then
		set ctime1 nobs=nobs end=eof;
	set ctime1 point=nobs;
	keep survival;
	output;
	stop;
run;

title 'This is the KM censoring probability at time T';
proc print data=ctime2;
run;


* Get KM survival probability at time T 4.96 years (chose 4.96 instead of 5 
because no one left after 5 years so won't work - took time of last event);
proc lifetest data=outkm_ext_uno method=km atrisk outsurv=stime noprint;
        time survtime*status(0);
run;

data stime1;
	set stime;
	if survtime gt 4.96 then delete;
run;

*get last survival prob at time t;
data stime2(rename=(survival=probt));
	if 0 then
		set stime1 nobs=nobs end=eof;
	set stime1 point=nobs;
	keep survival;
	output;
	stop;
run;

title 'This is the KM survival probability at time T';
proc print data=stime2;
run;


* Matrix code to calculate Uno's C-statistic;
* Should get C=0.639;
* NOTE: You must enter manual information below;
proc iml;
use outkm_ext_uno var {status survtime pi wsq};
read all;
show names;
ch = 0;
dh = 0;
/*MANUALLY ENTER TOTAL NUMBER IN DATASET*/
n=686;
/*MANUALLY ENTER TIME*/
t=4.96;

do i = 1 to n;
   		do j = 1 to n;
	  		if ((survtime[i,1] < survtime[j,1]) & (pi[i,1] > pi[j,1]) 
				& (survtime[i,1] < t)) then ch=ch + (status[i,1]/(wsq[i,1]));
	  		if ((survtime[i,1] < survtime[j,1]) & (pi[i,1] = pi[j,1]) 
			  & (survtime[i,1] < t)) then ch=ch + 0.5*(status[i,1]/(wsq[i,1]));

			if ((survtime[i,1] < survtime[j,1]) & (survtime[i,1] < t)) 
				then dh=dh + (status[i,1]/(wsq[i,1]));
   	 	end;
end;

c_hat = ch /dh;
create cstat var {ch dh c_hat};
append;
close cstat;
quit;

title 'Simple model Uno C';
proc print data = cstat; 
run;


* Matrix code to calculate Uno's Time-dependent AUROC at 5 years;
* Should get C=0.693;
* NOTE: You must enter manual information below;
* The below formula comes from Blanche et al 2013 Biom J paper;

proc iml;
use outkm_ext_uno var {status survtime pi w};
read all;
show names;
ch = 0;
/*MANUALLY ENTER TOTAL NUMBER IN DATASET*/
n=686;
/*MANUALLY ENTER TIME*/
t=4.96;
/*MANUALLY ENTER KM PROB OF CENSOR AT TIME T FROM 'CTIME2' ABOVE*/
probc=0.37674;
/*MANUALLY ENTER KM PROB OF SURVIVAL AT TIME T FROM 'STIME2' ABOVE*/
probt=0.49939;

do i = 1 to n;
   	do j = 1 to n;
	  	if ((survtime[i,1] <= t) & (survtime[j,1] > t) & (pi[i,1] > pi[j,1])) 
			then ch=ch + (status[i,1]/(n*n*probc*(w[i,1])));
	  	if ((survtime[i,1] <= t) & (survtime[j,1] > t) & (pi[i,1] = pi[j,1])) 
			then ch=ch + 0.5*(status[i,1]/(n*n*probc*(w[i,1])));
   	 end;
end;

c_hat = ch /(probt*(1-probt));
create cstat var {ch c_hat};
append;
close cstat;
quit;

title 'Simple model Uno time-dependent AUROC';
proc print data = cstat; 
run;

***** Now for extended model with PGR;
* Matrix code to calculate Uno's C-statistic;
* Should get C=0.667, SAS gives 0.665 when data available;
* NOTE: You must enter manual information below;
proc iml;
use outkm_ext_uno var {status survtime pi2 wsq};
read all;
show names;
ch = 0;
dh = 0;
/*MANUALLY ENTER TOTAL NUMBER IN DATASET*/
n=686;
/*MANUALLY ENTER TIME*/
t=4.96;

do i = 1 to n;
   		do j = 1 to n;
	  		if ((survtime[i,1] < survtime[j,1]) & (pi2[i,1] > pi2[j,1]) 
				& (survtime[i,1] < t)) then ch=ch + (status[i,1]/(wsq[i,1]));
	  		if ((survtime[i,1] < survtime[j,1]) & (pi2[i,1] = pi2[j,1]) 
			  & (survtime[i,1] < t)) then ch=ch + 0.5*(status[i,1]/(wsq[i,1]));

			if ((survtime[i,1] < survtime[j,1]) & (survtime[i,1] < t)) 
				then dh=dh + (status[i,1]/(wsq[i,1]));
   	 	end;
end;

c_hat = ch /dh;
create cstat var {ch dh c_hat};
append;
close cstat;
quit;

title 'Extended model Uno C';
proc print data = cstat; 
run;

* Matrix code to calculate Uno's Time-dependent AUROC at 5 years;
* Should get C=0.722;
* NOTE: You must enter manual information below;
* The below formula comes from Blanche et al 2013 Biom J paper;

proc iml;
use outkm_ext_uno var {status survtime pi2 w};
read all;
show names;
ch = 0;
/*MANUALLY ENTER TOTAL NUMBER IN DATASET*/
n=686;
/*MANUALLY ENTER TIME*/
t=4.96;
/*MANUALLY ENTER KM PROB OF CENSOR AT TIME T FROM 'CTIME2' ABOVE*/
probc=0.37674;
/*MANUALLY ENTER KM PROB OF SURVIVAL AT TIME T FROM 'STIME2' ABOVE*/
probt=0.49939;

do i = 1 to n;
   	do j = 1 to n;
	  	if ((survtime[i,1] <= t) & (survtime[j,1] > t) & (pi2[i,1] > pi2[j,1])) 
			then ch=ch + (status[i,1]/(n*n*probc*(w[i,1])));
	  	if ((survtime[i,1] <= t) & (survtime[j,1] > t) & (pi2[i,1] = pi2[j,1])) 
			then ch=ch + 0.5*(status[i,1]/(n*n*probc*(w[i,1])));
   	 end;
end;

c_hat = ch /(probt*(1-probt));
create cstat var {ch c_hat};
append;
close cstat;
quit;

title 'Extended model Uno time-dependent AUROC';
proc print data = cstat; 
run;

************* END OF DISCRIMINATION CODE ***************************;



**************** FIXED TIME POINT ASSESSMENT OF CALIBRATION ******************;


***Simple model first;

* This is calibration-in-the-large using Kaplan-Meier and average prediction ;
proc lifetest data=l2gbsg_val method=km atrisk outsurv=l2outkm_ext_c;
        time survtime*status(0);
run;

data l2outkm_ext_c1;
	set l2outkm_ext_c;
	where _censor_ = 0;
run;

*output observed survival at 5 years;
data l2outkm_ext_c2;
	set l2outkm_ext_c1;
	by _censor_;
	if last._censor_;
	n=1;
	keep n survival;
run;

*output predicted;
proc univariate data=l2gbsg_val noprint;
	var prob;
	output out=mean mean=meanpred;
	title ' ';
	proc print;
run;

data mean1;
	set mean;
	n=1;
run;

data mean2;
	merge l2outkm_ext_c2 mean1;
	by n;
	prop = (1-survival)/(meanpred);
run;

title 'Calibration-in-large using ratio of Kaplan-Meier & Avg predicted risk';
proc print data=mean2;
run;

**********- Bootstrap the 95% CIs for calibration-in-the-large ***********;

proc lifetest data=outboot method=km atrisk outsurv=l2outkm_ext_c_ noprint;
		by replicate;
        time survtime*status(0);
run;

data l2outkm_ext_c1_;
	set l2outkm_ext_c_;
	where _censor_ = 0;
run;

*output observed survival at 5 years;
data l2outkm_ext_c2_;
	set l2outkm_ext_c1_;
	by replicate _censor_;
	if last._censor_;
	keep replicate survival;
run;

*output predicted;
proc univariate data=outboot noprint;
	by replicate;
	var prob;
	output out=mean_ mean=meanpred;
run;

data mean2_;
	merge l2outkm_ext_c2_ mean_;
	by replicate;
	prop = (1-survival)/(meanpred);
	proc print;
run;

proc univariate data=mean2_ noprint;
	var prop;
	output out=confincal pctlpts=2.5 97.5 pctlpre=prop_ 
		pctlname=lower95 upper95;
run;

data confincal1;
	set confincal;
	n=1;
run;

data callarge;
 	retain prop prop_lower95 prop_upper95;
	merge mean2 confincal1;
	by n;
	drop n survival meanpred;
	title 'External Calibration-in-the-large with 95% CI';
	proc print;
run;



** Calibration slope;
data slope;
	set l2gbsg_val;
	pi_=pi;
run;

proc phreg data=slope;
	model survtime1*status(0) = pi  /  ties=efron rl;
run;

*miscalibration;
proc phreg data=slope;
	model survtime1*status(0) = pi  / offset = pi_ ties=efron rl;
run;




**Calibration plots using Austin et al 2020 Stat Med paper;

*Check functional form of CLL;
* Code CLL as a restricted cubic spline with 3 knots;
*First calculate the 10th, 50th and 90th percentiles for knots;
proc univariate data = l2gbsg_val;
	var cll;
	output out=knots pctlpre=p_cll pctlpts= 10 50 90;
run;

proc print data=knots; run;

/* here we find the following values:
P_CLL10 P_CLL50 P_CLL90 
-1.14854 -0.44554 0.31246 
*/

*Use Frank Harrell's RCSPLINE macro for calculating the spline terms;
data predhar1;
	set l2gbsg_val;
	* NOTE: you will need to take the spline information from the log 
	window and edit gridcox2 below;
	%RCSPLINE(CLL, -1.14854, -0.44554, 0.31246);
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
	 _kd_= (0.31246 - -1.14854)**.666666666666 ;
	cll1=max((cll--1.14854)/_kd_,0)**3
		+((-0.44554--1.14854)*max((cll-0.31246)/_kd_,0)**3
		-(0.31246--1.14854)*max((cll--0.44554)/_kd_,0)**3)/(0.31246--0.44554);
run;

*Calibration for predictions of 5-year survival probabilities;
proc phreg data=predhar1 zph(global transform=log);
	model survtime1*status(0)=cll cll1/ ties=efron rl;
	baseline out=gridcox3 covariates=gridcox2 survival=predprob lower=l 
		upper=u timelist=5;
run;

data gridcox4;
	set gridcox3;
	by prob;
	if first.prob;
	obsprob=1-predprob;
	l1=1-u;
	u1=1-l;
	if _name_='p_1' then diag1=0;
	if _name_='p_1' then diag2=0;
	if _name_='p_15' then diag1=1;
	if _name_='p_15' then diag2=1;
run;

proc sort data=gridcox4;
	by Prob;
run;

*Create density plot;
proc univariate data=predhar1 noprint;
   var prob;
   histogram prob / kernel midpoints=(0.005 to 0.995 by 0.01) outkernel=outhist; 
run;

data density;
   keep _value_ _density_;
   set outhist;
run;
 
/* Merge the counts with the predicted probabilities */
data LogiPlot2;
   set gridcox4(keep=Prob ObsProb DIAG1 DIAG2 L1 U1 /*<-ref line*/)
       density;
run;

title;
ods listing sge=on style=printer image_dpi=300 gpath='c:';
ods graphics on / reset=all noborder outputfmt=tiff  
	imagename="Calbration plot Ext Harrell with RCS" antialias=off;

proc sgplot data=logiplot2 noautolegend ;
	xaxis       	label="Predicted risk from developed model" 
		values=(0 to 1 by 0.20);
	yaxis        	label="Predicted risk from refitted model" 
		values=(0 to 1 by 0.20);
	y2axis			values=(0 to 20 by 1) display=none;
  	series x=prob y=obsprob / lineattrs=(color=red thickness=2) ;
  	band x=prob lower=l1 upper=u1 / nofill lineattrs=(color=black 
		pattern=mediumdash thickness=2) noextend outline;
  	reg x=diag1 y=diag2 / lineattrs=(pattern=shortdash) markerattrs=(size=1px 
		symbol=circlefilled);
  	series x=_value_ y=_density_/lineattrs=(color=black thickness=1 pattern=1)
		y2axis;
run;

ods graphics off;


*Calculation of ICI for 5 year probabilities;
data predharl_rep;
	set predhar1;
run;

proc phreg data=predhar1 noprint;
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
proc print;run;



******** BOOTSTRAP 95% CI FOR THE CALIBRATION METRICS ********;

proc univariate data = outboot noprint;
	by replicate;
	var cll;
	output out=knots pctlpre=p_cll pctlpts= 10 50 90;
run;

data predhar1_;
	merge outboot knots;
	by replicate;
run;


*This macro runs the RCS on each replicate and appends all 500 runs to one dataset;
%MACRO SPL;
%do i=1 %to 500;
data predhar1_1;
	set predhar1_;
	where replicate=&i;
	%rcspline(cll, p_cll10, p_cll50, p_cll90);
run;

*duplicate dataset so baseline statement works;
data predharl_rep_;
	set predhar1_1;
run;

proc phreg data=predhar1_1 noprint;
	model survtime1*status(0)=cll cll1/ ties=efron rl;
	baseline out=predharl_rep1_ covariates=predharl_rep_ survival=predprob2 
		timelist=5;
run;

data ici_;
	set predharl_rep1_;
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






********Extended model calibration;

* This is calibration-in-the-large using Kaplan-Meier and average prediction ;
proc lifetest data=l2gbsg_val method=km atrisk outsurv=l2outkm_ext_c;
        time survtime*status(0);
run;

data l2outkm_ext_c1;
	set l2outkm_ext_c;
	where _censor_ = 0;
run;

*output observed survival at 5 years;
data l2outkm_ext_c2;
	set l2outkm_ext_c1;
	by _censor_;
	if last._censor_;
	n=1;
	keep n survival;
run;

*output predicted;
proc univariate data=l2gbsg_val noprint;
	var prob2;
	output out=mean mean=meanpred;
run;

data mean1;
	set mean;
	n=1;
run;

data mean2_pgr;
	merge l2outkm_ext_c2 mean1;
	by n;
	prop = (1-survival)/(meanpred);
run;

title 'Calibration-in-large using ratio of Kaplan-Meier & Avg predicted risk';
proc print data=mean2_pgr;
run;

**********- Bootstrap the 95% CIs for calibration-in-the-large ***********;

/*
*no need to rerun below if you have already run for model without pgr;
proc lifetest data=outboot method=km atrisk outsurv=l2outkm_ext_c_ noprint;
		by replicate;
        time survtime*status(0);
run;

data l2outkm_ext_c1_;
	set l2outkm_ext_c_;
	where _censor_ = 0;
run;

*output observed survival at 5 years;
data l2outkm_ext_c2_;
	set l2outkm_ext_c1_;
	by replicate _censor_;
	if last._censor_;
	keep replicate survival;
run;
*/

*output predicted;
proc univariate data=OUTBOOT noprint;
	by replicate;
	var prob2;
	output out=mean__ mean=meanpred;
run;

data mean2_pgr_;
	merge l2outkm_ext_c2_ mean__;
	by replicate;
	prop = (1-survival)/(meanpred);
run;

proc univariate data=mean2_pgr_ noprint;
	var prop;
	output out = confincal_ pctlpts=2.5 97.5 pctlpre= prop_ 
		pctlname=lower95 upper95;
run;

data confincal1_;
	set confincal_;
	n=1;
run;

data callarge_;
 	retain prop prop_lower95 prop_upper95;
	merge mean2_pgr confincal1_;
	by n;
	drop n ;
	title 'External Calibration-in-the-large with 95% CI (PGR included)';
	proc print;
run;



** Calibration slope;
data slope1;
	set l2gbsg_val;
	pi2_=pi2;
run;

proc phreg data=slope1;
	model survtime1*status(0) = pi2  /  ties=efron rl;
run;

*miscalibration;
proc phreg data=slope1;
	model survtime1*status(0) = pi2  / offset = pi2_ ties=efron rl;
run;




**Calibration plots using Austin et al 2020 Stat Med paper;

*Check functional form of CLL;
* Code CLL as a restricted cubic spline with 3 knots;
*First calculate the 10th, 50th and 90th percentiles for knots;
proc univariate data = l2gbsg_val;
	var cll_;
	output out=knots pctlpre=p_cll pctlpts= 10 50 90;
run;

proc print data=knots; run;

/* here we find the following values:

P_CLL10 P_CLL50 P_CLL90 
-1.01594 -0.38373 0.37213  
*/
*Use Frank Harrell's RCSPLINE macro for calculating the spline terms;
data predharpgr1;
	set l2gbsg_val;
	* NOTE: you will need to take the spline information from the log 
	window and edit gridcox2 below;
	%rcspline(cll_, -1.01594, -0.38373, 0.37213);
run;

proc univariate data = predharpgr1 noprint;
	var prob2;
	output out=gridcox pctlpre=p_ pctlpts= 1 to 99;
run;

proc transpose data=gridcox out=gridcox1;
	var p_1 - p_99;
run;

data gridcox2(rename=(col1=prob));
	set gridcox1;
	cll_=log(-log(1-col1));
	*Warning! - take the below formula from the log window after running above 
	rcspline - this will be different according to your data!;
 	_kd_= (0.37213 - -1.01594)**.666666666666;
	cll_1=max((cll_--1.01594)/_kd_,0)**3
		+((-0.38373--1.01594)*max((cll_-0.37213)/_kd_,0)**3
		-(0.37213--1.01594)*max((cll_--0.38373)/_kd_,0)**3)/(0.37213--0.38373);
run;

*calibration for predictions of 5-year survival probabilities;
proc phreg data=predharpgr1 zph(global transform=log);
	model survtime1*status(0)=cll_ cll_1/ ties=efron rl;
	baseline out=gridcox3 covariates=gridcox2 survival=predprob lower=l 
		upper=u timelist=5;
run;

data gridcox4;
	set gridcox3;
	by prob;
	if first.prob;
	obsprob=1-predprob;
	l1=1-u;
	u1=1-l;
	if _name_='p_1' then diag1=0;
	if _name_='p_1' then diag2=0;
	if _name_='p_15' then diag1=1;
	if _name_='p_15' then diag2=1;
run;

proc sort data=gridcox4;
	by prob;
run;

*Create density plot;
proc univariate data=predharpgr1;
   var prob2;
   histogram prob2 /kernel midpoints=(0.005 to 0.995 by 0.01) outkernel=OutHist;
run;

data density;
   keep _value_ _density_;
   set OutHist;
run;
 
/* Merge the counts with the predicted probabilities. */
data LogiPlot2;
   set gridcox4(keep=Prob ObsProb DIAG1 DIAG2 L1 U1 /*<-ref line*/)
       density;
run;

title;
ods listing sge=on style=printer image_dpi=300 gpath='C:';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="Calbration plot PGR Ext Harrell with RCS" antialias=off;

proc sgplot data=logiplot2 noautolegend ;
xaxis       	label="Predicted risk from developed model" 
		values=(0 to 1 by 0.20);
yaxis        	label="Predicted risk from refitted model" 
		values=(0 to 1 by 0.20);
y2axis	values=(0 to 20 by 1) display=none;
  series x=prob y=obsprob / lineattrs=(color=red thickness=2) ;
  band x=prob lower=l1 upper=u1 / nofill lineattrs=(color=black 
		pattern=mediumdash thickness=2) noextend outline;
  reg x=diag1 y=diag2 / lineattrs=(pattern=shortdash) 
		markerattrs=(size=1px symbol=circlefilled);
  series x=_value_ y=_density_ / lineattrs=(color=black thickness=1 
		pattern=1) y2axis;
run;

ods graphics off;

*Calculation of ICI for 5 year probabilities;
data predharpgrl_rep;
	set predharpgr1;
run;

proc phreg data=predharpgr1;
	model survtime1*status(0)=cll_ cll_1/ ties=efron rl;
	baseline out=predharpgrl_rep1 covariates=predharpgrl_rep survival=predprob2 
		timelist=5;
run;

data ici_pgr;
	set predharpgrl_rep1;
	*observed prob death;
	prob3=1-predprob2;
	diff=abs(prob2 - prob3);
run;

proc univariate data=ici_pgr noprint;
	var diff;
	output out=ici_pgrs mean=ici median=e50 p90=e90;
run;

data ici_pgrs;
	set ici_pgrs;
	n=1;
run;

proc print; run;


******** BOOTSTRAP 95% CI FOR THE CALIBRATION METRICS ********;
proc univariate data = outboot noprint;
	by replicate;
	var cll_;
	output out=knots pctlpre=p_cll pctlpts= 10 50 90;
run;

data predharpgr1_;
	merge outboot knots;
	by replicate;
run;

*This macro runs the RCS on each replicate and appends to PREDHarl_2;
%MACRO SPL1;
%do i=1 %to 500;
data predhar1_1;
	set predharpgr1_;
	where replicate=&i;
	%rcspline(cll_, p_cll10, p_cll50, p_cll90);
run;

*duplicate dataset so baseline statement works;
data predharl_rep_;
	set predhar1_1;
run;

proc phreg data=predhar1_1 noprint;
	model survtime1*status(0)=cll_ cll_1/ ties=efron rl;
	baseline out=predharl_rep1_ covariates=predharl_rep_ survival=predprob2 
		timelist=5;
run;

data ici_;
	set predharl_rep1_;
	*observed prob death;
	prob3=1-predprob2;
	diff=abs(prob2 - prob3);
run;

proc univariate data=ici_ noprint;
	var diff;
	output out=icis_ mean=ici median=e50 p90=e90;
run;

data icis_1;
	set icis_;
	replicate=&i;
run;

proc append base=icis_3 data=icis_1 force;
run;

%end;
%mend;
%SPL1;

proc univariate data=icis_3 noprint;
	var ici e50 e90;
	output out = conf_met_ pctlpts=2.5 97.5 pctlpre= ici_ e50_ e90_ 
		pctlname=lower95 upper95;
run;

data conf_met_1;
	set conf_met_;
	n=1;
run;

data calmetrics_;
 	retain ici ici_lower95 ici_upper95 e50 e50_lower95 e50_upper95 e90 
		e90_lower95 e90_upper95;
	merge ici_pgrs conf_met_1;
	by n;
	drop n;
	title 'External Calibration Metrics with 95% CI for model with PGR added';
	proc print;
run;




*********************** Overall performance;

******* This block calculates the Brier score - see Graf et al 1999;
title ' ';

*Simple model first;
*External validation - Calculate Brier score;

*Calculate weights using Kaplan-Meier;
proc lifetest data=l2gbsg_val method=km atrisk outsurv=l2outkm_ext noprint;
        time survtime*status(1);
run;

*Create 3 groups - Group 1-Those who have the event up to fixed event time of 
*interest, Group 2 - those who go beyond fixed time (could be event or event 
*free), and Group 3- those censored up to fixed time 
*Only first 2 groups contribute to score but all to weights;
data l2rott_b_ex;
	set l2gbsg_val;
	*Must make the time just under 5 years since we need some time remaining
	*after timepoint of interest (5 years) and before administrative censoring 
	*(5 years) for it to work;
	if survtime<=4.99 and status=1 then cat=1;
	if survtime>4.99 then cat=2;
	if survtime<=4.99 and status=0 then cat=3;
	proc sort; by survtime;
run;

data l2outkm_ext1;
	set l2outkm_ext;
	weight=1/survival;
	keep survtime survival weight;
run;

data l2outkm_ext2;
	set l2outkm_ext1;
	s=lag(survtime);
	if survtime=s then delete;
	drop s;
run;

data l2rottval_bs2;
	merge l2rott_b_ex l2outkm_ext1;
	by survtime;
	if pid=. then delete;
run;

data l2rottval_bs3;
	set l2rottval_bs2;
	retain _weight;
	if not missing(weight) then _weight=weight;
	else weight=_weight;
	if cat=3 then weight=0;
	if survtime=0 then delete;
	if cat=1 then contrib=(-psurv)**2;
	if cat=2 then contrib=(1-psurv)**2;
	if cat=3 then contrib=0;
	bs=contrib*weight;
	drop _weight;
run;

*Estimate brier score;
proc univariate data=l2rottval_bs3 noprint;
	var bs weight;
	output out=l2sumsex sum=sbs sweight;
	proc print; 
run;

data l2sumsex;
	retain sweight sbs brier;
	set l2sumsex;
	brier = (1/sweight)*sbs;
	sweight=left(sweight);
	title 'External Brier score';
	proc print;
run;


**** - Calculate bootstrapped 95% ci for brier score;

*calculate weights for external validation;
proc lifetest data=outboot method=km atrisk outsurv=l2outkm_extx noprint;
		by replicate;
        time survtime*status(1);
run;

data l2rott_b_exx;
	set outboot;
	if survtime<=4.99 and status=1 then cat=1;
	if survtime>4.99 then cat=2;
	if survtime<=4.99 and status=0 then cat=3;
	proc sort; by replicate survtime;
run;

data l2outkm_ext1x;
	set l2outkm_extx;
	weight=1/survival;
	keep replicate survtime survival weight;
run;

data l2rottval_bs2x;
	merge l2rott_b_exx l2outkm_ext1x;
	by replicate survtime;
	if pid=. then delete;
run;

data l2rottval_bs3x;
	set l2rottval_bs2x;
	retain _weight;
	if not missing(weight) then _weight=weight;
	else weight=_weight;
	if cat=3 then weight=0;
	if survtime=0 then delete;
	if cat=1 then contrib=(-psurv)**2;
	if cat=2 then contrib=(1-psurv)**2;
	if cat=3 then contrib=0;
	bs=contrib*weight;
	drop _weight;
run;

*Estimate brier score;
proc univariate data=l2rottval_bs3x noprint;
	by replicate;
	var bs weight;
	output out=l2sumsexx sum=sbs sweight;
	proc print; 
run;

data l2sumsexx;
	retain sweight sbs brier;
	set l2sumsexx;
	brier = (1/sweight)*sbs;
	sweight=left(sweight);
run;



***Scaled brier for external validation;
proc phreg data=l2gbsg_val;
	model survtime*status(0)= / ties=efron rl;
	baseline covariates=l2rott_b_ex out=l2rottval_bsnull timelist=5 
		survival=fiveyrsurv_null;
run;

*Brier for null model;
data l2rottval_bs1null;
	set l2rottval_bsnull;
	keep pid fiveyrsurv_null;
	proc sort; by pid;
run;

proc sort data=l2rottval_bs3;
	by pid;
run;

*Merge the null survival probs to the brier score dataset from earlier;
data l2rottval_bs4;
	merge l2rottval_bs3 l2rottval_bs1null;
	by pid;
run;

data l2rottval_bs5;
	set l2rottval_bs4;
	if cat=1 then contrib_null=(-fiveyrsurv_null)**2;
	if cat=2 then contrib_null=(1-fiveyrsurv_null)**2;
	if cat=3 then contrib_null=0;
	bs_null=contrib_null*weight;
	drop fiveyrsurv_null;
run;

*Estimate brier score for null model;
proc univariate data=l2rottval_bs5 noprint;
	var bs_null weight;
	output out=l2sumnullex sum=sbs_null sweight;
	proc print; 
run;

data l2sumnullex;
	retain sweight sbs_null;
	set l2sumnullex;
	sweight=left(sweight);
run;

*Calculate scaled brier for external validation;
data l2ipa;
	merge l2sumsex l2sumnullex;
	by sweight;
	null_brier = (1/sweight)*sbs_null;
	scaledb = 1-(brier/null_brier);
	ind=1;
		title 'External Brier score and Scaled Brier score';
	proc print;
run;


**** - Calculate bootstrapped 95% ci for scaled brier;

proc phreg data=outboot noprint;
	by replicate;
	model survtime*status(0)= / ties=efron rl;
	baseline covariates=l2rott_b_exx out=l2rottval_bsnullx timelist=5 
		survival=fiveyrsurv_null;
run;

*Brier for null model;
data l2rottval_bs1nullx;
	set l2rottval_bsnullx;
	keep replicate pid fiveyrsurv_null;
	proc sort; by replicate pid;
run;

proc sort data=l2rottval_bs3x;
	by replicate pid;
run;

*Merge the bull survival probs to the brier score dataset from earlier;
data l2rottval_bs4x;
	merge l2rottval_bs3x l2rottval_bs1nullx;
	by replicate pid;
run;

data l2rottval_bs5x;
	set l2rottval_bs4x;
	if cat=1 then contrib_null=(-fiveyrsurv_null)**2;
	if cat=2 then contrib_null=(1-fiveyrsurv_null)**2;
	if cat=3 then contrib_null=0;
	bs_null=contrib_null*weight;
	drop fiveyrsurv_null;
run;

*estimate brier score for null model;
proc univariate data=l2rottval_bs5x noprint;
	by replicate;
	var bs_null weight;
	output out=l2sumnullexx sum=sbs_null sweight;
	proc print; 
run;

data l2sumnullexx;
	retain sweight sbs_null;
	set l2sumnullexx;
	sweight=left(sweight);
run;

*calculate scaled brier for external validation;
data l2ipax;
	merge l2sumsexx l2sumnullexx;
	by replicate sweight;
	null_brier = (1/sweight)*sbs_null;
	scaledb = 1-(brier/null_brier);
	*proc print;
run;

proc univariate data=l2ipax noprint;
	var brier scaledb;
	output out = confintr pctlpts=2.5 97.5 pctlpre= brier_ scaledb_ 
		pctlname=lower95 upper95;
run;

data confintr1;
	set confintr;
	ind=1;
run;

data brierax2;
 	retain brier brier_lower95 brier_upper95 scaledb scaledb_lower95 
		scaledb_upper95;
	merge l2ipa confintr1;
	by ind;
	drop ind sweight sbs sbs_null null_brier;
	title 'External Brier score and Scaled Brier score with 95% CI';
	proc print;
run;

************* Above repeated for added marker ********************;

*Simple model first;
*External validation - Calculate Brier score;

*Calculate weights using Kaplan-Meier;
proc lifetest data=l2gbsg_val method=pl atrisk outsurv=l2outkm_exta noprint;
        time survtime*status(1);
run;

*Create 3 groups - Group 1-Those who have the event up to fixed event time of 
*interest, Group 2 - those who go beyond fixed time (could be event or event 
*free), and Group 3- those censored up to fixed time 
*Only first 2 groups contribute to score but all to weights;
data l2rott_b_exa;
	set l2gbsg_val;
	if survtime<=4.99 and status=1 then cat=1;
	if survtime>4.99 then cat=2;
	if survtime<=4.99 and status=0 then cat=3;
	proc sort; by survtime;
run;

data l2outkm_ext1a;
	set l2outkm_exta;
	weight=1/survival;
	keep survtime weight;
run;

data l2rottval_bs2a;
	merge l2rott_b_exa l2outkm_ext1a;
	by survtime;
	if pid=. then delete;
run;

data l2rottval_bs3a;
	set l2rottval_bs2a;
	retain _weight;
	if not missing(weight) then _weight=weight;
	else weight=_weight;
	if cat=3 then weight=0;
	if survtime=0 then delete;
	if cat=1 then contrib=(-psurv2)**2;
	if cat=2 then contrib=(1-psurv2)**2;
	if cat=3 then contrib=0;
	bs=contrib*weight;
	drop _weight;
run;

*estimate brier score;
proc univariate data=l2rottval_bs3a noprint;
	var bs weight;
	output out=l2sumsexa sum=sbs sweight;
	proc print; 
run;

*brier exactly same as full dataset;
data l2sumsexa;
	set l2sumsexa;
	brier = (1/sweight)*sbs;
	sweight=left(sweight);
	title 'External brier score for pgr';
	proc print;
run;


**** - Calculate bootstrapped 95% ci for brier score;

*calculate weights for external validation;
proc lifetest data=outboot method=km atrisk outsurv=l2outkm_extx_ noprint;
		by replicate;
        time survtime*status(1);
run;

data l2rott_b_exx_;
	set outboot;
	if survtime<=4.99 and status=1 then cat=1;
	if survtime>4.99 then cat=2;
	if survtime<=4.99 and status=0 then cat=3;
	proc sort; by replicate survtime;
run;

data l2outkm_ext1x_;
	set l2outkm_extx_;
	weight=1/survival;
	keep replicate survtime survival weight;
run;

data l2rottval_bs2x_;
	merge l2rott_b_exx_ l2outkm_ext1x_;
	by replicate survtime;
	if pid=. then delete;
run;

data l2rottval_bs3x_;
	set l2rottval_bs2x_;
	retain _weight;
	if not missing(weight) then _weight=weight;
	else weight=_weight;
	if cat=3 then weight=0;
	if survtime=0 then delete;
	if cat=1 then contrib=(-psurv2)**2;
	if cat=2 then contrib=(1-psurv2)**2;
	if cat=3 then contrib=0;
	bs=contrib*weight;
	drop _weight;
run;

*estimate brier score;
proc univariate data=l2rottval_bs3x_ noprint;
	by replicate;
	var bs weight;
	output out=l2sumsexx_ sum=sbs sweight;
	*proc print; 
run;

data l2sumsexx_;
	retain sweight sbs brier;
	set l2sumsexx_;
	brier = (1/sweight)*sbs;
	sweight=left(sweight);
run;


*** Scaled Brier for model with PGR included;
proc sort data=l2rottval_bs3a;
	by pid;
run;

data l2rottval_bs4a;
	merge l2rottval_bs3a l2rottval_bs1null;
	by pid;
run;

data l2rottval_bs5a;
	set l2rottval_bs4a;
	if cat=1 then contrib_null=(-fiveyrsurv_null)**2;
	if cat=2 then contrib_null=(1-fiveyrsurv_null)**2;
	if cat=3 then contrib_null=0;
	bs_null=contrib_null*weight;
	drop fiveyrsurv_null;
run;

*Estimate brier score for null model;
proc univariate data=l2rottval_bs5a noprint;
	var bs_null weight;
	output out=l2sumnullexa sum=sbs_null sweight;
	proc print; 
run;

data l2sumnullexa;
	set l2sumnullexa;
	sweight=left(sweight);
run;

data l2scaledba;
	merge l2sumsexa l2sumnullexa;
	by sweight;
	null_brier = (1/sweight)*sbs_null;
	scaledb = 1-(brier/null_brier);
	ind=1;
	title 'External Brier score and Scaled Brier with PGR';
	proc print;
run;



**** - Calculate bootstrapped 95% ci for scaled brier;

proc phreg data=outboot noprint;
	by replicate;
	model survtime*status(0)= / ties=efron rl;
	baseline covariates=l2rott_b_exx_ out=l2rottval_bsnullx_ timelist=5 
		survival=fiveyrsurv_null;
run;

*Brier for null model;
data l2rottval_bs1nullx_;
	set l2rottval_bsnullx_;
	keep replicate pid fiveyrsurv_null;
	proc sort; by replicate pid;
run;

proc sort data=l2rottval_bs3x_;
	by replicate pid;
run;

*Merge the null survival probs to the brier score dataset from earlier;
data l2rottval_bs4x_;
	merge l2rottval_bs3x_ l2rottval_bs1nullx_;
	by replicate pid;
run;

data l2rottval_bs5x_;
	set l2rottval_bs4x_;
	if cat=1 then contrib_null=(-fiveyrsurv_null)**2;
	if cat=2 then contrib_null=(1-fiveyrsurv_null)**2;
	if cat=3 then contrib_null=0;
	bs_null=contrib_null*weight;
	drop fiveyrsurv_null;
run;

*estimate brier score for null model;
proc univariate data=l2rottval_bs5x_ noprint;
	by replicate;
	var bs_null weight;
	output out=l2sumnullexx_ sum=sbs_null sweight;
	*proc print; 
run;

data l2sumnullexx_;
	retain sweight sbs_null;
	set l2sumnullexx_;
	sweight=left(sweight);
run;

*Calculate scaled brier for external validation;
data l2ipax_;
	merge l2sumsexx_ l2sumnullexx_;
	by replicate sweight;
	null_brier = (1/sweight)*sbs_null;
	scaledb = 1-(brier/null_brier);
	*proc print;
run;

proc univariate data=l2ipax_ noprint;
	var brier scaledb;
	output out = confintr_ pctlpts=2.5 97.5 pctlpre= brier_ scaledb_ 
		pctlname=lower95 upper95;
run;

data confintr1_;
	set confintr_;
	ind=1;
run;

data brierax2;
 	retain brier brier_lower95 brier_upper95 scaledb scaledb_lower95 
		scaledb_upper95;
	merge l2scaledba confintr1_;
	by ind;
	drop ind sweight sbs sbs_null null_brier;
	title 'External Brier score and scaled Brier score with 95% CI';
	proc print;
run;




*** DCA;

*external;

%PUT _user_;
GOPTIONS RESET = ALL;

symbol1 i=join c=green;
symbol2 i=join c=red;
symbol3 i=join c=blue;
symbol4 i=join c=darkred;
symbol5 i=join c=gray;

%STDCA(data=L2GBSG_VAL, out=survivalmult, outcome=STATUS, ttoutcome=SURVTIME, timepoint=5, predictors=PROB);
%STDCA(data=L2GBSG_VAL, out=survivalmult_new, outcome=STATUS, ttoutcome=SURVTIME, timepoint=5, predictors=PROB2);

*Sort by threshold variable;
proc sort data=survivalmult out=kmsort;
	by threshold;
run;

*rename the variables so that we know they are the kaplan meier estimates;
data kmsort; set kmsort(rename=(prob=kmmodel all=kmall));
	label kmmodel="Original model: Pr(Recurrence) at 5 years";
	label kmall="Treat All";
RUN;

*Sort by threshold variable;
proc sort data=survivalmult_new out=crsort;
	by threshold;
run;

*rename the variables so that we know they are the competing risk estimates;
data crsort; set crsort(rename=(prob2=crmodel all=crall));
	label crmodel="Original + PGR: Pr(Recurrence) at 5 years";
RUN;

*Merge Kaplan-Meier and Competing Risk data using threshold probabilities;
data crsort;
	merge kmsort crsort;
	by threshold;
run;

title;
footnote;
*- this allows editing of the .sge file!;
ods listing sge=on style=printer image_dpi=300 gpath='c:';
ods graphics on / reset=all noborder outputfmt=tiff 
		imagename="dca with spline pgr ext" antialias=off/*antialiasmax=*/;

*create graph (decision curve) with treat none, treat all, simple model
and extended model;
proc sgplot data=crsort;
	yaxis values=(-0.05 to 0.55 by 0.05) label="Net Benefit";
	xaxis values=(0.0 to 1 by 0.1) label="Threshold Probability";
	keylegend "all" "orig" "new" "none" / down=4 position=bottom;
	series y=kmall x=threshold / lineattrs=(color=black thickness=2 
		pattern=solid) name="all" legendlabel="Treat All";
	series y=kmmodel x=threshold / lineattrs=(color=green thickness=2 
		pattern=solid) name="orig" 
		legendlabel="Original model: Pr(Rec) at 5 years";
	series y=crmodel x=threshold / lineattrs=(color=blue thickness=2 
		pattern=solid) name="new" 
		legendlabel="Original model + PGR: Pr(Rec) at 5 years";
	series y=none x=threshold / lineattrs=(color=red thickness=2 
		pattern=solid)  name="none" 
		legendlabel="None";
RUN;

ods graphics off;
