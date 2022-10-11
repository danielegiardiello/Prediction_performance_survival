*******************************************************************************;
* Program:			STRATOS Model development and internal validation 		   ;
*					with PGR final.sas 	   									   ;
* Author:			David McLernon        					 				   ;
* Date: 			27th July 2022          					 			   ;
* Purpose: 			This is SAS code for the STRATOS paper on validation of	   ;
*					survival risk prediction models. This programme covers     ;
*					model development internal validation 					   ;
* Note:				This programme is not currently automated. It is coded	   ;
*					based on the case study in the article. Therefore, 		   ;
*					adapting this code for your own study will require careful ;
*					editing according to your data							   ;
*					Readers can skip the sections that create "nice" graphs,   ;
*					these all have surrounding comments of  "Optional Block"   ;
*					and "End Optional Block". Though an essential skill for    ;
*					published papers, such manipulation is outside the primary ;
*					scope of this work.										   ;
* Data:				The Rotterdam and GBSG data sets can be found on the web   ;
*					in various locations and various formats.  For this 	   ;
*					exercise we have used the versions found in the R package. ;
*					(Older web sites have a habit of disappearing and we 	   ;
*					expect R to be more stable). In R we used 				   ;
*				   	library(survival)									       ;
*			       	write.csv(rotterdam, row.names=FALSE, file="rotterdam.csv");
*			       	write.csv(gbsg,      row.names=FALSE, file="gbsg.csv")	   ;
*				       to create the csv files.								   ;
* Updates:			Updated on 4th Oct 2022 with bootstrap corrected mean and  ;
*					weak calibration at fixed time point assessment and over   ;
*					the time range											   ;
*******************************************************************************;	
** Lines up to 260 are the same as in the STRATOS Model development and internal 
** validation final.sas programme;

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

/*
** Create the categorical variables that we use in the model fit
** grade of 1-2 vs 3
** size of  <=20, 21-50, >50
** The size variable is already categorized in this way, in both data sets
*/

*** descriptive statistics for age and PGR;
proc univariate data=rotterdam;
	var age pgr nodes;
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
    
    * Winzorise PGR and nodes to the 99th percentile to deal with large influential
    *  values, found from output of univariate procedure above;
    if (pgr > 1360) then pr = 1360; else pr = pgr;
	drop pgr;
	if (nodes > 19) then nod = 19; else nod = nodes;
	drop nodes;
run;

** Check functional form of nodes;
* Code nodes as a restricted cubic spline with 3 knots;
* First calculate the 10th, 50th and 90th percentiles for knots;
proc univariate data = r1 noprint;
	var nodes;
	output out=knots pctlpre=p_node pctlpts= 10 50 90;
run;

proc print data=knots; run;

/* in output window we find the following values:
p_node10 p_node50 p_node90 
 0 			1 			9 
*/

*Use Frank Harrell's RCSPLINE macro for calculating the spline terms - found 
*here: http://biostat.mc.vanderbilt.edu/wiki/Main/SasMacros;
* instead of %rcspline you could also use following statement within PHREG:
EFFECT spl = spline(nodes / details naturalcubic basis=tpf(noint) 
	knotmethod=percentilelist(10 50 90));

data r11;
	set r1;
	%rcspline(nodes, 0, 1, 9);
	label survtime= 'Survival time (years)';
run;


* Descriptive statistics at baseline for Table 1;
proc freq data=r11;
    table gradecat sizecat;
run;


* Plot the overall Kaplan-Meier;
* - n=1275 events by 5 years;
title "Kaplan-Meier plot";
footnote;
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=png /*width=4in*/ 
imagename="dev km" antialias=off/*antialiasmax=*/ ;

proc lifetest data=r11 method=pl plots=(s(atrisk(atrisktickonly outside)
			=0,5,10,15,20 cb=hw nocensor name=Survival)) outsurv=outkm;
        time survtime*status(0);
		ods exclude ProductLimitEstimates;   * suppress the long table;
		title "Overall Kaplan-Meier";
run;

ods graphics off;

title;



*get median follow-up using reverse kaplan-meier method;
proc lifetest data=r11 method=pl atrisk;
    time survtime*status(1);
    ods exclude ProductLimitEstimates;   * suppress the long table;
    title "Median follow-up time";
run;


*Administrative censor at 5 years since this is our prediction horizon; 
data rott;
	set r11;
	*administrative censoring at 5 years;
	if survtime > 5 then status=0;
	if survtime > 5 then survtime=5;
run;

** Check functional form of pgr;
* Code PGR as a restricted cubic spline with 3 knots;
* First calculate the 10th, 50th and 90th percentiles for knots;
proc univariate data = rott noprint;
	var pgr;
	output out=knots pctlpre=p_pgr pctlpts= 10 50 90;
run;

proc print data=knots; run;

/* here we find the following values:
P_PGR10 P_PGR50 P_PGR90 
0 		41 		486 
*/

*Use Frank Harrell's RCSPLINE macro for calculating the spline terms - found 
*here: http://biostat.mc.vanderbilt.edu/wiki/Main/SasMacros;
data ffpgr;
	set rott;
	*make copy of survival time for calibration bootstrap later;
	survtime1=survtime;
	%rcspline(pgr, 0, 41, 486);
run;

******* Optional Block ****** ;
*The following code will produce Suppl Fig 1 - restricted cubic spline plot for 
*PGR;
data ffpgr1;
	set ffpgr;
run;

title 'Fit Cox model with PGR terms and save predictions at 5 years';
proc phreg data=ffpgr;
	model survtime*status(0)=pgr pgr1/ ties=efron;
	test1: test pgr, pgr1;
	test2: test pgr+pgr1=0;
	baseline covariates=ffpgr1 out=pgrval survival=predpgr lower=predpgrl 
	upper=predpgrup timelist=5;
run;

* Manipulate data so that we can plot the diagonal ref line and probability 
* of death (not survival);
data pgrval1;
	set pgrval;
	if pid=3 then diag1=0;
	if pid=3 then diag2=0;
	if pid=7 then diag1=100;
	if pid=7 then diag2=100;
	predpgr_dth=1-predpgr;
	predpgr_low=1-predpgrl;
	predpgr_upp=1-predpgrup;
	keep pid pgr predpgr_dth predpgr_low predpgr_upp diag1 diag2;
run;

proc sort data=pgrval1;
	by pgr;
run;

title;
footnote;
*- this allows editing of the .sge file!;
ods listing sge=on style=printer image_dpi=300 
gpath='C:';
ods graphics on / reset=all noborder outputfmt=tiff /*width=4in*/ 
imagename="Spline of PGR" antialias=off/*antialiasmax=*/;

* Plot the restricted cubic spline of PGR with predicted probability of death 
* within 5 years;
proc sgplot data=pgrval1 noautolegend;
	xaxis       	label="PGR" values=(0 to 1400 by 100) /*edit labels and value 
		range as appropriate here and below*/;
	yaxis        	label="Predicted probability of recurrence or death" 
				values=(0 to 0.6 by 0.1);
  	band x=pgr lower=predpgr_low upper=predpgr_upp / nofill lineattrs=
		(color=black pattern=mediumdash thickness=3) noextend outline;
  	series y=predpgr_dth x=pgr / lineattrs=(color=black thickness=3) ;
run;

ods graphics off;

******  End Optional Block ************;


******************** FIT MODEL WITH PGR MARKER INCLUDED ***********************;

* Calculate the 1st and 3rd quartiles for interquartile hazard ratio 
* estimation below;
proc univariate data=ffpgr noprint;
	var pgr pgr1;
	output out=percentiles q3=q3pgr q3pgr1 q1=q1pgr q1pgr1;
run;
proc print data=percentiles;
run;
/*Q3PGR Q3PGR1    Q1PGR   Q1PGR1 
 198   14.9704   4       0.000270961 */

* Calculate the 1st and 3rd quartiles for interquartile hazard ratio 
* estimation below;
proc univariate data=ffpgr noprint;
	var nodes nodes1;
	output out=percentiles q3=q3nod q3nod1 q1=q1nod q1nod1;
run;
proc print data=percentiles;
run;
/*q3nod q3nod1 	q1nod 	q1nod1 
	4 	0.41512 0 		0 
 */


*** Estimating Baseline Survival Function under PH;
* - Take lowest values to calculate baseline survival;
data inrisks;
	set ffpgr;
	where sizecat=0 & nodes=0 & gradecat=0 & PGR=0 & PGR1=0;
	n+1;
	if n>1 then delete;
	keep sizecat nodes nodes1 gradecat PGR PGR1;
run;

Title 'Extended model with proportional hazards assessment';
*Note: Can block off saving the plots by blocking ods lines below;
ods listing sge=on style=printer image_dpi=300 gpath='c:';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="Extended model with proportional hazards assessment" 
	antialias=off;

proc phreg data=ffpgr /*zph(global transform=log) ev*/;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat gradecat nodes nodes1 pgr pgr1/ ties=efron rl;
	* Interquartile Hazard Ratio - calculate below by subtracting Q3PGR from 
	* Q1PGR above, and Q3PGR1 from Q1PGR1;
	contrast 'Test PGR' PGR -194 PGR1 -14.970129039 / estimate=both;
	contrast 'Test Nodes' nodes -4 nodes1 -0.41512 / estimate=both;
	* store model for easy validation later;
	store SimpModelPGR;
	*assess functional form;
	assess ph / resample crpanel;
	* This statement estimates baseline survival at yearly times;
	baseline covariates=inrisks out=outph1 survival=ps timepoint=(1,2,3,4,5) / 
		method=breslow;
	* Save the PI to the original dataset;
	output out=rottxa xbeta=xb;
run;

ods graphics off;

proc print data=outph1;
run;

/* 
Baseline survival copied from output window..
sizecat gradecat pgr nodes nodes1 pgr1 survtime ps 
<=20 	1-2 	0 	0 		0 		0 	1 		0.96123 
<=20 	1-2 	0 	0 		0 		0 	2 		0.89709 
<=20 	1-2 	0 	0 		0 		0 	3 		0.84235 
<=20 	1-2 	0 	0 		0 		0 	4 		0.79835 
<=20 	1-2 	0 	0 		0 		0 	5 		0.76055 
*/

*Call bootstrap resampled dataset from previous programme (without pgr);
data outboot;
	set stratos.outboot;
run;

proc phreg data=outboot noprint;
	*where replicate in (1,2) /*useful line to use when testing out*/;
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron rl;
	store simpmodelboot;
	*baseline statement applies each of the 500 models to the original dataset 
	external bootstrap validation;
	baseline covariates=ffpgr out=fpgrboot xbeta=xb timelist=5;
	*output statement applies each of the 500 models to the corresponding 
	dataset - apparent bootstrap validation;
	output out=fpgrapp xbeta=xb;
run;


*****TIME RANGE ASSESSMENT OF DISCRIMINATION;

*Apparent discrimination;

*Harrell's C - Need tau to equal event time of interest;
title 'Apparent Harrells C';
proc phreg data=ffpgr concordance=harrell(se) tau=5;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron rl;
run;

*Uno's C - Need tau to equal event time of interest;
title 'Apparent Unos C';
proc phreg data=ffpgr concordance=uno(se seed=8754 iter=50) tau=5;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron rl;
run;

*Internal validation: Bootstrap to assess optimism in performance;
*Only run the below for the C statistic you wish to use;

*get apparent bootstrap performance for Harrell's C in 500 datasets;
proc phreg data=outboot concordance=harrell(se) tau=5;
	*This line is useful to test bootstrap on first 2 resampled datasets only;
	*where replicate in (1,2);
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron rl;
	ods output concordance=apphar(rename=(estimate=apphar stderr=appseh) 
		drop=source);
run;

*Get apparent bootstrap performance for Uno's C in 500 datasets;
proc phreg data=outboot concordance=uno(se seed=8754 iter=50) tau=5;
	*This line is useful to test bootstrap on first 2 resampled datasets only;
	*where replicate in (1,2);
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron rl;
	ods output concordance=appuno(rename=(estimate=appuno stderr=appseu) 
		drop=source);
run;

*Next run the external bootstrap performance on original dataset;
sasfile ffpgr load; /* a way of loading the dataset into ram - speeds it up */

*This should replicate the dataset 500 times without replacement - trick to get 
*500 copies - needed for concordance to work below;
proc surveyselect data=ffpgr out=outrott
seed=853794
method=srs /* simple random sampling without replacement */
samprate=1 /* can accept proportions or percentages but we want n to be size of 
			orginal database so =1 (or 100) */
rep=500; /* number of bootstrap samples */
run;

sasfile ffpgr close; /* closes frees ram buffers when done */

ods listing close; 

*Get external bootstrap performance for Harrell's C in 500 copies of Rotterdam 
*dataset;
proc phreg data=outrott concordance=harrell(se) tau=5;
	*where replicate in (1,2);
	by replicate;	
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron rl 
			nofit;
	roc source=simpmodelboot;
	ods output concordance=boothar(rename=(estimate=boothar stderr=bootseh) 
			drop=source);
run;

*get external bootstrap performance for Uno's c in 500 copies of rotterdam 
	dataset;
proc phreg data=outrott concordance=uno(se seed=8754 iter=50) tau=5;
	*where replicate in (1,2);
	by replicate;	
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron rl 
			nofit;
	roc source=simpmodelboot;
	ods output concordance=bootuno(rename=(estimate=bootuno stderr=bootseu) 
		drop=source);
run;

*calculate optimism for Harrell's C per replicate;
data harcon;
	merge apphar boothar;
	by replicate;
	optharc = apphar - boothar;
run;

*calculate optimism for Uno's C per replicate;
data unocon;
	merge appuno bootuno;
	by replicate;
	optunoc = appuno - bootuno;
run;

****FIXED TIME POINT ASSESSMENT OF DISCRIMINATION;
title 'Apparent Unos AUC';
proc phreg data=ffpgr tau=5 rocoptions(auc at=4.98 method=ipcw (cl seed=134));
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron rl ;
	*roc 'npi + pgr' source=simpmodelpgr;
run;

*Internal validation: Bootstrap performance;
*Get apparent bootstrap performance for Uno's AUC in 500 datasets;
proc phreg data=outboot tau=5 rocoptions(auc at=4.98 method=ipcw (cl seed=134));
	*where replicate in (1,2);
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 /ties=efron rl;
	ods output auc=apptduno(rename=(estimate=apptduno stderr=apptdse) 
		drop=sourceid upper lower source);
run;

*Get external bootstrap performance for Uno's AUC in 500 copies of rotterdam 
	dataset;
proc phreg data=outrott tau=5 rocoptions(auc at=4.98 method=ipcw (cl seed=134));
	*where replicate in (1,2);
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron rl 
		nofit;
	roc source=simpmodelboot;
	ods output auc=boottduno(rename=(estimate=boottduno stderr=boottdse) 
		drop=source sourceid upper lower);
run;

*calculate optimism for Uno's AUC;
data tdunocon;
	merge apptduno boottduno;
	by replicate;
	optunoauc = apptduno - boottduno;
run;








*************** This block carries out internal validation assessment of;
*************** calibration using bootstrap;
*************** Optimism of time range calibration ******************;

*** Simple model first;
*** Using Crowson's Poisson approach;

*- CUMHAZ command provides the expected number of events and XBETA provides 
	linear predictor for original model applied to external cohort;

*- For calibration we need to determine performance of each model developed;
*- on the bootstrap dataset applied to the original data. The average of the;
*- 500 performance estimates is the optimism;
*- First mean calibration;

* This needs to be run in chunks of 100 as we get memory problems otherwise;
data outbootsub;
	set outboot;
	if 1<=replicate<=50 then gp=1;
	else if 51<=replicate<=100 then gp=2;
	else if 101<=replicate<=150 then gp=3;
	else if 151<=replicate<=200 then gp=4;
	else if 201<=replicate<=250 then gp=5;
	else if 251<=replicate<=300 then gp=6;
	else if 301<=replicate<=350 then gp=7;
	else if 351<=replicate<=400 then gp=8;
	else if 401<=replicate<=450 then gp=9;
	else gp=10;
run;

%macro calib(subset);
	%do i=1 %to &subset;

		proc phreg data=outbootsub;
			where gp=&i;
			by replicate;
			class sizecat (ref='<=20') gradecat (ref=first);
			model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1/ ties=efron rl;
			baseline out=ffpgrboot covariates=ffpgr cumhaz=p xbeta=lp;
		run;

		*get the original survival time for each person;
		data pred1abt;
			set ffpgrboot;
			where survtime le survtime1; 
			proc sort; by replicate pid survtime;
		run;

		*calculate log of expected number of events;
		data pred1bt;
			set pred1abt;
			by replicate pid survtime;
			*first get the follow-up time for each patient;
			if last.pid;
			*p is the outcome-martingale residual;
			logp=log(p);
			logbase=logp-lp;
		run;

		*Model 1 - Mean calibration (O/E);
		*take exponential of Intercept estimate;
		proc genmod data=pred1bt;
			by replicate;
			model status=/offset=logp dist=poisson link=log;
			ods output parameterestimates=par;
		run;

		data meancalbt(rename=(estimate=tr_meancal));
			set par;
			if parameter="Scale" then delete;
			keep replicate estimate;
		run;

		*Weak assessment - calibration slope;
		*the coefficient is the slope;
		proc genmod data=pred1bt;
			by replicate;
			model status=lp /offset=logbase dist=poisson link=log;
			ods output parameterestimates=par1;
		run;

		data weakcalbt(rename=(estimate=tr_weakcal));
			set par1;
			if parameter ne "lp" then delete;
			keep replicate estimate;
		run;

		*Combine mean and weak time range calibration;
		data timeranc;
			merge meancalbt weakcalbt;
			by replicate;
		run;

		*delete large datasets to save memory;
		proc datasets library=work nolist; 
			delete ffpgrboot pred1abt pred1bt; 
		quit;

		proc append base=timeranc2 data=timeranc force;
		run;

	%end;

%mend;

%calib(10);



*************** Optimism of Fixed time point calibration *****************;

*** Mean / Calibration-in-the-large **;

**- O/E using ratio of (1-Kaplan-Meier) and Avg predicted risk at 5 years;
**- Lets first get the Kalan-Meier for 5 years in original dataset;

proc lifetest data=ffpgr method=km atrisk outsurv=outkm_ext_c noprint;
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
	keep survival;
run;

*replicate it 500 times for easy merge later with average predictions;
data outkm_ext_c3(rename=(i=replicate));
	set outkm_ext_c2;
	do i=1 to 500;
	output;
	end;
run;

*Output predicted for each of the 500 bootstrap datasets applied to the original;
*development data;
proc phreg data=outboot;
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1/ ties=efron rl;
	baseline out=ffpgrbt covariates=ffpgr survival=psurv xbeta=lp 
		timelist=5;
run;

*Extract the predictions and get mean;
proc univariate data=ffpgrbt noprint;
	by replicate;
	var psurv;
	output out=mean mean=meanpred;
run;

*Merge observed and predicted means, calculate proportion and 95% CI;
data mean2;
	merge outkm_ext_c3 mean;
	by replicate;
	f_meancal = (1-survival)/(1-meanpred);
run;

title 'Calibration-in-large using ratio of Kaplan-Meier & Avg predicted risk';
proc print data=mean2;
run;


**Weak - Calibration slope;

proc phreg data=ffpgrbt;
	by replicate;
	model survtime1*status(0) = lp /  ties=efron rl;
	ods output parameterestimates=par2;
run;

data f_weakcalbt(rename=(estimate=f_weakcal));
	set par2;
	keep replicate estimate;
run;

*Combine mean and weak fixed time calibration;
data timefixc;
	merge mean2 f_weakcalbt;
	by replicate;
run;








************************************ Overall performance;

******* This block calculates the Brier score - see Graf et al 1999;
title ' ';

****First for apparent validation i.e. model development;
*calculate weights using Kaplan-Meier for apparent validation;
proc lifetest data=rottxa method=pl atrisk outsurv=outkm noprint;
        time survtime*status(1);
run;

*Create 3 groups - Group 1-Those who have the event up to fixed event time of 
*interest, Group 2 - those who go beyond fixed time (could be event or event 
*free), and Group 3- those censored up to fixed time 
*Only first 2 groups contribute to score but all to weights;
data rott_b;
	set rottxa;
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

*Now estimate survival at 5 years in development dataset;
title 'Basic model output';
proc phreg data=ffpgr;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron 
		rl;
	*apparent;
	baseline covariates=rott_b out=rott_bs timelist=5 survival=fiveyrsurv/ 
		method=breslow;
run;

*code up the fixed time of 5 years;
data rott_bs1(rename=(time=survtime));
	set rott_bs;
	time1=time;
	if time1>5 then time1=5;
	drop survtime;
	proc sort; by time1;
run;

*Merge the kaplan-meier weights to the appropriate times;
data outkm1(rename=(survtime=time1));
	set outkm;
	if survtime>5 then delete;
	weight=1/survival;
	keep survtime weight;
run;

data rott_bs2;
	merge rott_bs1 outkm1;
	by time1;
	if pid=. then delete;
run;

*Then for group 1 calculate -surv^2, and group 2 calculate (1-surv)^2 where surv 
*is probability of surv at t; 
data rott_bs3;
	set rott_bs2;
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
proc univariate data=rott_bs3 noprint;
	var bs weight;
	output out=sums sum=sbs sweight;
	proc print; 
run;

data sums;
	retain sweight sbs brier;
	set sums;
	brier = (1/sweight)*sbs;
	sweight=left(sweight);
	ind=1;
	title 'Apparent Brier score with PGR';
	proc print;
run;

******* End of block;


********* This block bootstraps the 95% CI for the Brier score;
****First calculate weights in boostrapped datasets;
proc sort data=outboot;
	by replicate pid;
run;

proc lifetest data=outboot method=pl outsurv=outkm_ noprint;
		by replicate;
        time survtime*status(1);
run;

data outbootx;
	set outboot;
	keep replicate pid;
	proc sort; by pid;
run;

*Merge bootstrap set to the rott_bs1 dataset from earlier which contains the  
*predictions at 5 years; 
proc sort data=rott_bs1;
	by pid;
run;

data rott_bs1_;
	merge outbootx rott_bs1;
	by pid;
run;

proc sort data=rott_bs1_;
	by replicate time1;
run;

*Merge the kaplan-meier weights to the appropriate times;
data outkm1_(rename=(survtime=time1));
	set outkm_;
	if survtime>5 then delete;
	weight=1/survival;
	keep replicate survtime weight;
run;

data rott_bs2_;
	merge rott_bs1_ outkm1_;
	by replicate time1;
	if pid=. then delete;
run;

*Create groups as before;
data rott_bs3_;
	set rott_bs2_;
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
proc univariate data=rott_bs3_ noprint;
	by replicate;
	var bs weight;
	output out=sums_ sum=sbs sweight;
	proc print; 
run;

data sums_;
	retain sweight sbs brier;
	set sums_;
	brier = (1/sweight)*sbs;
	sweight=left(sweight);
	title 'Brier score';
	proc print;
run;

*Note that the 95% CIs will be presented with the scaled Brier results later;
********* End of block;


************ This block performs the bootstrapped internal validation for the Brier score;


*** Now do the Brier bootstrap apparent performance;
proc lifetest data=fpgrapp method=pl outsurv=outkma noprint;
		by replicate;
        time survtime*status(1);
run;

*for bootstrap external validation;
proc lifetest data=ffpgr method=pl outsurv=outkm_ext noprint;
        time survtime*status(1);
run;

*Create 3 groups - Group 1-Those who have the event up to fixed event time of 
*interest, Group 2 - those who go beyond fixed time (could be event or event 
*free), and Group 3- those censored up to fixed time 
*Only first 2 groups contribute to score but all to weights;
data rott_ba;
	set fpgrapp;
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

*And do same for original data;
data rott_b_ex;
	set ffpgr;
	if survtime<=4.99 and status=1 then cat=1;
	if survtime>4.99 then cat=2;
	if survtime<=4.99 and status=0 then cat=3;
	time=survtime;
run;

*Estimate survival at 5 years for apparent bootstrap;
proc phreg data=fpgrapp noprint;
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron 
		rl;
	baseline covariates=rott_ba out=rott_bsa timelist=5 survival=fiveyrsurv/ 
		method=breslow;
run;

*Estimate survival at 5 years for external bootstrap;
proc phreg data=fpgrapp noprint;
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1 / ties=efron 
		rl;
	baseline covariates=rott_b_ex out=rottval_bs timelist=5 survival=fiveyrsurv/ 
		method=breslow;
run;

*Apparent bootstrap for Brier;
*Merge the kaplan-meier weights to the appropriate times;
data rott_bs1a(rename=(time=survtime));
	set rott_bsa;
	time1=time;
	if time1>5 then time1=5;
	drop survtime;
	proc sort; by replicate time1;
run;

data outkm1a(rename=(survtime=time1));
	set outkma;
	if survtime>5 then delete;
	weight=1/survival;
	keep replicate survtime weight;
run;

data rott_bs2a;
	merge rott_bs1a outkm1a;
	by replicate time1;
	if pid=. then delete;
run;

data rott_bs3a;
	set rott_bs2a;
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
proc univariate data=rott_bs3a noprint;
	by replicate;
	var bs weight;
	output out=sumsa sum=appsbs sweight;
run;

data sumsa;
	retain sweight appsbs appbrier;
	set sumsa;
	appbrier = (1/sweight)*appsbs;
	sweight=left(sweight);
run;

*External bootstrap for brier;
*Merge the kaplan-meier weights to the appropriate times;
data rottval_bs1(rename=(time=survtime));
	set rottval_bs;
	time1=time;
	if time1>5 then time1=5;
	drop survtime;
	proc sort; by time1;
run;

data outkm_ext1(rename=(survtime=time1));
	set outkm_ext;
	if survtime>5 then delete;
	weight=1/survival;
	keep survtime weight;
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
proc sort data= rottval_bs3;
	by replicate pid;
run;

proc univariate data=rottval_bs3 noprint;
	by replicate;
	var bs weight;
	output out=sumsex sum=bootsbs sweight;
run;

data sumsex;
	retain sweight bootsbs bootbrier;
	set sumsex;
	bootbrier = (1/sweight)*bootsbs;
	sweight=left(sweight);
run;


************ End of block;



********************************Scaled Brier;
*Scaled Brier = 1 - (model Brier score/null model Brier score), where null 
*model Brier score is null cox; 
*100% is perfect, <0 is useless, higher better, harmful models <0;

*Apparent Scaled Brier first;
*Estimate survival at 5 years for null model;
title 'Null model';
proc phreg data=rottxa;
	model survtime*status(0)= / ties=efron rl;
	baseline covariates=rott_b out=rott_bsnull timelist=5 
		survival=fiveyrsurv_null;
run;

data rott_bs1null;
	set rott_bsnull;
	keep pid fiveyrsurv_null;
	proc sort; by pid;
run;

proc sort data=rott_bs3;
	by pid;
run;

*apparent brier for null model;
*Merge the null model survival probabilities to brier score dataset from 
*earlier;
data rott_bs4;
	merge rott_bs3 rott_bs1null;
	by pid;
run;

data rott_bs5;
	set rott_bs4;
	if cat=1 then contrib_null=(-fiveyrsurv_null)**2;
	if cat=2 then contrib_null=(1-fiveyrsurv_null)**2;
	if cat=3 then contrib_null=0;
	bs_null=contrib_null*weight;
	drop fiveyrsurv_null;
run;

*Estimate brier score for null model;
proc univariate data=rott_bs5 noprint;
	var bs_null weight;
	output out=sumnull sum=sbs_null sweight;
	proc print; 
run;

data sumnull;
	retain sweight sbs_null;
	set sumnull;
	sweight=left(sweight);
run;

*Calculate Scaled Brier;
data scaledb;
	merge sumnull sums;
	by sweight;
	null_brier = (1/sweight)*sbs_null;
	scaled_b = 1-(brier/null_brier);
	title 'Brier score and Scaled Brier with PGR';
	proc print;
run;



******** This block calculates the Bootstrapped 95% CI for Scaled Brier Score;

*Now estimate survival at 5 years for null model;
data l2rott_b_exx;
	set outboot;
	if survtime<=4.99 and status=1 then cat=1;
	if survtime>4.99 then cat=2;
	if survtime<=4.99 and status=0 then cat=3;
	proc sort; by replicate survtime;
run;

proc phreg data=outboot;
	by replicate;
	model survtime*status(0)= / ties=efron rl;
	baseline covariates=l2rott_b_exx out=rott_bsnull_ timelist=5 
		survival=fiveyrsurv_null;
run;

data rott_bs1null_;
	set rott_bsnull_;
	keep replicate pid fiveyrsurv_null;
	proc sort; by replicate pid;
run;

proc sort data=rott_bs3_;
	by replicate pid;
run;

*apparent brier for null model;
*Merge  null model survival probabilities to brier score dataset from earlier;
data rott_bs4_;
	merge rott_bs3_ rott_bs1null_;
	by replicate pid;
run;

data rott_bs5_;
	set rott_bs4_;
	if cat=1 then contrib_null=(-fiveyrsurv_null)**2;
	if cat=2 then contrib_null=(1-fiveyrsurv_null)**2;
	if cat=3 then contrib_null=0;
	bs_null=contrib_null*weight;
	drop fiveyrsurv_null;
run;

*Estimate brier score for null model;
proc univariate data=rott_bs5_ noprint;
	by replicate;
	var bs_null weight;
	output out=sumnull_ sum=sbs_null sweight;
	proc print; 
run;

data sumnull_;
	retain sweight sbs_null;
	set sumnull_;
	sweight=left(sweight);
run;

*Calculate Scaled Brier;
data scaledb_;
	merge sumnull_ sums_;
	by replicate;
	null_brier = (1/sweight)*sbs_null;
	scaled_b = 1-(brier/null_brier);
run;

proc univariate data=scaledb_ noprint;
	var brier scaled_b;
	output out = confintr_ pctlpts=2.5 97.5 pctlpre= brier_ scaledb_ 
		pctlname=lower95 upper95;
run;

data confintr1_;
	set confintr_;
	ind=1;
run;

***** End of block;


***** This contains the Brier score and Scaled Brier score with their  
***** bootstrapped 95% CI;

data brierax2_;
 	retain brier brier_lower95 brier_upper95 scaled_b scaledb_lower95 scaledb_upper95;
	merge scaledb confintr1_;
	by ind;
	drop ind sweight sbs;
	title 'Brier score and Scaled Brier score with 95% CI';
	proc print;
run;




******** Block to do Bootstrap to assess optimism in performance of the 
******** Scaled Brier;

* Begin with apparent bootstrap performance;
proc phreg data=fpgrapp noprint;
	by replicate;
	model survtime*status(0)= / ties=efron rl;
	baseline covariates=rott_ba out=outkm_null timelist=5 
		survival=fiveyrsurv_null;
run;

*set the 5 year survival probability in the apparent bootstrap dataset;
data outkm_null1;
	set outkm_null;
	keep replicate pid fiveyrsurv_null;
	proc sort; by replicate pid;
run;

proc sort data= rott_bs3a;
	by replicate pid;
run;

data rott_bs4a;
	merge rott_bs3a outkm_null1;
	by replicate pid;
run;

data rott_bs5a;
	set rott_bs4a;
	if cat=1 then contrib_null=(-fiveyrsurv_null)**2;
	if cat=2 then contrib_null=(1-fiveyrsurv_null)**2;
	if cat=3 then contrib_null=0;
	bs_null=contrib_null*weight;
	drop fiveyrsurv_null;
run;

*Estimate brier score;
proc univariate data=rott_bs5a noprint;
	by replicate;
	var bs_null weight;
	output out=sumnulla sum=appsbs_null sweight;
run;

data sumnulla;
	retain sweight appsbs_null;
	set sumnulla;
	sweight=left(sweight);
run;

data scaledba;
	merge sumsa sumnulla;
	by replicate;
	null_appbrier = (1/sweight)*appsbs_null;
	appscaledb = 1-(appbrier/null_appbrier);
run;

*External bootstrap performance;
proc phreg data=fpgrapp noprint;
	by replicate;
	model survtime*status(0)= / ties=efron rl;
	baseline covariates=rott_b_ex out=outkm_nullex timelist=5 
		survival=fiveyrsurv_null;
run;

data outkm_nullex1;
	set outkm_nullex;
	keep replicate pid fiveyrsurv_null;
	proc sort; by replicate pid;
run;

proc sort data=rottval_bs3;
	by replicate pid;
run;

data rottval_bs4;
	merge rottval_bs3 outkm_nullex1;
	by replicate pid;
run;

data rottval_bs5;
	set rottval_bs4;
	if cat=1 then contrib_null=(-fiveyrsurv_null)**2;
	if cat=2 then contrib_null=(1-fiveyrsurv_null)**2;
	if cat=3 then contrib_null=0;
	bs_null=contrib_null*weight;
	drop fiveyrsurv_null;
run;

*Estimate brier score;
proc sort data=rottval_bs5;
	by replicate pid;
run;

proc univariate data=rottval_bs5 noprint;
	by replicate;
	var bs_null weight;
	output out=sumnullex sum=bootsbs_null sweight;
	proc print; 
run;

data scaledb2;
	merge sumsex sumnullex;
	by replicate;
	null_bootbrier = (1/sweight)*bootsbs_null;
	bootscaledb = 1-(bootbrier/null_bootbrier);
run;

******** End block;

* Contains optimism corrected Brier and Scaled Brier scores per replicate;
data finalipa;
	merge scaledba scaledb2;
	by replicate;
	optbrier = appbrier - bootbrier;
	optscaledb = appscaledb - bootscaledb;
title 'Optimism Brier score and Scaled Brier score with PGR';
	proc print;
run;



** NOW PUT ALL TOGETHER TO GET THE OPTIMISM IN PERFORMANCE & 95% CI FOR 
** APPARENT PERFORMANCE;
data all1;
	*If you did not calculate all of the performance measures you will need to 
		delete the related dataset from the merge statement;
	merge finalipa harcon unocon tdunocon timefixc timeranc2;
	by replicate;
run;

data stratos.all_pgr1;
	set all1;
run;

*calculate the average optimism in fit;
proc univariate data= all1 noprint;
	var OptBrier OptScaledB OptHarC OptUnoC OptUnoAUC f_meancal f_weakcal 
		tr_meancal tr_weakcal;
	output out=avgopt mean=Mean_OptBrier Mean_OptScaleB Mean_OptHarC 
		Mean_OptUnoC Mean_OptUnoAUC Mean_OptFMean Mean_OptFWeak Mean_OptTMean 
		Mean_OptTWeak;
run;

**** Manually enter the performance measures for each model;
data resultsdata;
	infile datalines delimiter=',';
	input OrigBrier OrigScBr OrigHarC OrigUnoC OrigUnoAUC OrigFMean OrigFWeak 
		OrigTMean OrigTWeak;
	datalines;
	0.206, 0.161, 0.689, 0.688, 0.727, 1, 1, 1, 1
;
data resultsdata_new;
	set resultsdata;
	ind=1;
run;

data avgopt1;
	set avgopt;
	ind=1;
run;

*Calculate the difference between the measure from final model and the average 
*optimism in fit to get a nearly unbiased estimate of the expected value of the 
*external performance measure, in other words, internal validation is an honest 
*estimate of internal validity penalising for overfitting;
data correctedperf;
	merge resultsdata_new avgopt1;
	by ind;
	InternValBrier=OrigBrier - Mean_OptBrier;
	InternValSBrier=OrigScBr - Mean_OptScaleB;
	InternValHarC=OrigHarC - Mean_OptHarC;
	InternValUnoC=OrigUnoC - Mean_OptUnoC;
	InternValUnoAUC=OrigUnoAUC - Mean_OptUnoAUC;
	InternValFMean=OrigFMean - Mean_OptFMean;
	InternValFWeak=OrigFWeak - Mean_OptFWeak;
	InternValTMean=OrigTMean - EXP(Mean_OptTMean);
	InternValTWeak=OrigTWeak - Mean_OptTWeak;
	drop ind;
run;

title;
title 'Estimate of internal validity penalising for overfitting for model with PGR';
proc print data=correctedperf;
run;

proc export data=correctedperf
  	outfile= "C:\Users\sme544\Bootstrap_Val_PGR_Calibration" 
   	dbms=xlsx replace;
run;
