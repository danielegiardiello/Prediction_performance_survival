*******************************************************************************;
* Program:			STRATOS Model development and internal validation.sas 	   ;
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
*				   library(survival)									       ;
*			       write.csv(rotterdam, row.names=FALSE, file="rotterdam.csv") ;
*			       write.csv(gbsg,      row.names=FALSE, file="gbsg.csv")	   ;
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

/*
** Create the categorical variables that we use in the model fit
** grade of 1-2 vs 3
** size of  <=20, 21-50, >50
** The size variable is already categorized in this way, in both data sets
*/

*** descriptive statistics for age, nodes and PGR;
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
run;


******* Optional Block ****** ;
*The following code will produce Suppl Fig 1 - restricted cubic spline plot for 
*nodes;
data r111;
	set r11;
run;

title 'Fit Cox model with nodes terms and save predictions at 5 years';
proc phreg data=r11;
	model survtime*status(0)=nodes nodes1/ ties=efron;
	test1: test nodes, nodes1;
	test2: test nodes+nodes1=0;
	baseline covariates=r111 out=nodval survival=prednod lower=prednodl 
	upper=prednodup timelist=5;
run;

* Manipulate data so that we can plot the diagonal ref line and probability 
* of death (not survival);
data nodval1;
	set nodval;
	if pid=3 then diag1=0;
	if pid=3 then diag2=0;
	if pid=7 then diag1=100;
	if pid=7 then diag2=100;
	prednod_dth=1-prednod;
	prednod_low=1-prednodl;
	prednod_upp=1-prednodup;
	keep pid nodes prednod_dth prednod_low prednod_upp diag1 diag2;
run;

proc sort data=nodval1;
	by nodes;
run;

title;
footnote;
*- this allows editing of the .sge file!;
ods listing sge=on style=printer image_dpi=300 
gpath='C:';
ods graphics on / reset=all noborder outputfmt=tiff /*width=4in*/ 
imagename="Spline of Nodes" antialias=off/*antialiasmax=*/;

* Plot the restricted cubic spline of nodes with predicted probability of death 
* within 5 years;
proc sgplot data=nodval1 noautolegend;
	xaxis       	label="Number of nodes" values=(0 to 20 by 2) /*edit labels and value 
		range as appropriate here and below*/;
	yaxis        	label="Predicted probability of recurrence or death" 
				values=(0 to 1 by 0.1);
  	band x=nodes lower=prednod_low upper=prednod_upp / nofill lineattrs=
		(color=black pattern=mediumdash thickness=3) noextend outline;
  	series y=prednod_dth x=nodes / lineattrs=(color=black thickness=3) ;
run;

ods graphics off;

******  End Optional Block ************;




* Descriptive statistics at baseline for Table 1;
proc freq data=r11;
    table gradecat sizecat;
run;

* Plot the overall Kaplan-Meier, as motivation for using a 5 year cut-off;
* - n=1275 events by 5 years;
proc lifetest data=r11 method=pl plots=(s, ls, lls) outsurv=outkm;
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

*Create kaplan meier plot for development dataset;
title;
footnote;
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=png /*width=4in*/ 
imagename="dev km" antialias=off/*antialiasmax=*/;

* create graph; 
proc sgplot data=outkm1;
	yaxis values=(0 to 1 by 0.2) label="Recurrence-free survival probability";
	xaxis values=(0 to 16 by 1) label="Years";
	step y=survival x=survtime / lineattrs=(color=blue thickness=2 pattern=solid) 
	name="all";
run;

ods graphics off;

******  End Optional Block ************;


*get median follow-up using reverse kaplan-meier method;
proc lifetest data=r11 method=pl atrisk;
    time survtime*status(1);
    ods exclude ProductLimitEstimates;   * suppress the long table;
    title "Median follow-up time";
run;

/*
** The ZPH command option requests diagnostics using scaled Schoenfeld
**  residuals, to check proportional hazards assumption;
**  The model shows strong evidence of non-proportional hazards for all.
*/
title;
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=tiff 
	imagename="PH" ANTIALIAS=OFF;
proc phreg data=r11 zph(global transform=log fit=spline);
	class sizecat (ref=first) gradecat (ref=first);
	model survtime*status(0)=sizecat gradecat nodes nodes1 / ties=efron rl;
run;
ods graphics off;
 
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
	%rcspline(pgr, 0, 41, 486);
run;


******* Optional Block ****** ;
*The following code will produce Suppl Fig - restricted cubic spline plot for 
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

ods graphics on;
proc phreg data=/*r11*/ ffpgr /*zph(global transform=identity fit=loess)*/;
	ltimepgr = pgr*log(survtime);
	ltimepgr1 = pgr1*log(survtime);
	class sizecat (ref=first) /*nodescat (ref=first)*/ gradecat (ref=first);
	model survtime*status(0)=sizecat gradecat nodes nodes1 pgr pgr1 pgr*ltimepgr pgr1*ltimepgr1/ ties=efron rl;
run;



********************FIT BASIC MODEL (without PGR)*******************;
*** Estimating Baseline Survival Function under PH;
* - Take lowest values to calculate baseline survival;
data inrisk;
	set rott;
	where sizecat=0 & nodes=0 & gradecat=0;
	n+1;
	if n>1 then delete;
	keep sizecat nodes nodes1 gradecat;
run;

*Fit model;
title 'Basic model';

* Calculate the 1st and 3rd quartiles for interquartile hazard ratio 
* estimation below;
proc univariate data=rott noprint;
	var nodes nodes1;
	output out=percentiles q3=q3nod q3nod1 q1=q1nod q1nod1;
run;
proc print data=percentiles;
run;
/*q3nod q3nod1 	q1nod 	q1nod1 
	4 	0.41512 0 		0 
 */

Title 'Basic model';

proc phreg data=ffpgr;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat gradecat nodes nodes1 / ties=efron rl;
	* Interquartile Hazard Ratio - calculate below by subtracting Q3nod from 
	* Q1nod above, and Q3nod1 from Q1nod1;
	contrast 'Test Nodes' nodes -4 nodes1 -0.41512 / estimate=both;
	* store model for easy validation later;
	store simpmodel;
	* This statement estimates baseline survival at yearly times;
	baseline covariates=inrisk out=outph survival=ps timepoint=(1,2,3,4,5) /
		method=breslow;
	* Save the PI to the original dataset;
	output out=rottx xbeta=xb;
run;

proc print data=outph;
run;

*we find following results for baseline survival at yearly timepoints;
/*sizecat gradecat nodes nodes1 survtime 	ps */
/*<=20 		1-2 	0 		0 		1 		0.96872 
/* <=20 	1-2 	0 		0 		2 		0.91676 
/* <=20 	1-2 	0 		0 		3 		0.87193 
/* <=20 	1-2 	0 		0 		4 		0.83548 
/* <=20 	1-2 	0 		0 		5 		0.80385 
*/
******Resample 500 versions of the internal dataset with replacement to allow 
******calculation of 95% CI and internal validation;


*some bootstrap code to edit;

sasfile ffpgr load; /* a way of loading the dataset into ram - speeds it up */

proc surveyselect data=ffpgr out=outboot 
seed=4817 /* can enter 0 for it to select a random seed but remember to type 
				it in here from the output otherwise cannot replicate results */
method=urs /* unrestricted random sampling - simple random sampling */
samprate=1 /* can accept proportions or percentages but we want n to be size 
				of original database so =1 (or 100) */
outhits /* with replacement */
rep=500; /* number of bootstrap samples */
run;

sasfile ffpgr close; /* closes frees ram buffers when done */

ods listing close; /* turns off ods listing so no printing of all output. 
					Better than using noprint as it doesn't allow storage of 
					data in output dataset at end */

*save for later use in PGR programme;
data stratos.outboot;
	set outboot;
run;


*Apply bootstrap model to bootstrapped dataset and the original dataset;

proc phreg data=outboot noprint;
	*where replicate in (1,2) /*useful line to use when testing out*/;
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat   / ties=efron rl;
	store simpmodelboot;
	*baseline statement applies each of the 500 models to the original dataset 
	external bootstrap validation;
	baseline covariates=rott out=rottboot xbeta=xb timelist=5;
	*output statement applies each of the 500 models to the corresponding 
	dataset - apparent bootstrap validation;
	output out=rottapp xbeta=xb;
run;



*****TIME RANGE ASSESSMENT OF DISCRIMINATION;

*Apparent discrimination;

*Harrell's C - Need tau to equal event time of interest;
title 'Apparent Harrells C';
proc phreg data=rott concordance=harrell(se) tau=5;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
run;

*Uno's C - Need tau to equal event time of interest;
title 'Apparent Unos C';
proc phreg data=rott concordance=uno(se seed=8754 iter=50) tau=5;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
run;

*Internal validation: Bootstrap to assess optimism in performance;
*Only run the below for the C statistic you wish to use;

*Get apparent bootstrap performance for Harrell's C in 500 datasets;
proc phreg data=outboot concordance=harrell(se) tau=5;
	*This line is useful to test bootstrap on first 2 resampled datasets only;
	*where replicate in (1,2);
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
	ods output concordance=apphar(rename=(estimate=apphar stderr=appseh) 
		drop=source);
run;

*Get apparent bootstrap performance for Uno's C in 500 datasets;
proc phreg data=outboot concordance=uno(se seed=8754 iter=50) tau=5;
	*This line is useful to test bootstrap on first 2 resampled datasets only;
	*where replicate in (1,2);
	by replicate;
	class sizecat (ref='<=20')gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
	ods output concordance=appuno(rename=(estimate=appuno stderr=appseu) 
		drop=source);
run;

*Next run the external bootstrap performance on original dataset;
sasfile rott load; /* a way of loading the dataset into ram - speeds it up */

*This should replicate the dataset 500 times without replacement - trick to get 
*500 copies - needed for concordance to work below;
proc surveyselect data=rott out=outrott
seed=853794 
method=srs /* simple random sampling without replacement */
samprate=1 /* can accept proportions or percentages but we want n to be size of 
			orginal database so =1 (or 100) */
rep=500; /* number of bootstrap samples */
run;

sasfile rott close; /* closes frees ram buffers when done */

ods listing close; 

*Get external bootstrap performance for Harrell's C in 500 copies of Rotterdam 
*dataset;
proc phreg data=outrott concordance=harrell(se) tau=5;
	*where replicate in (1,2);
	by replicate;	
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
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
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
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
proc phreg data=rott tau=5 rocoptions(auc at=4.98 method=ipcw (cl seed=134));
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
	*roc 'npi' source=simpmodel;
run;

*Internal validation: Bootstrap performance;
*Get apparent bootstrap performance for Uno's AUC in 500 datasets;
proc phreg data=outboot tau=5 rocoptions(auc at=4.98 method=ipcw (cl seed=134));
	*where replicate in (1,2);
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
	ods output auc=apptduno(rename=(estimate=apptduno stderr=apptdse) 
		drop=sourceid upper lower source);
run;

*Get external bootstrap performance for Uno's AUC in 500 copies of rotterdam 
	dataset;
proc phreg data=outrott tau=5 rocoptions(auc at=4.98 method=ipcw (cl seed=134));
	*where replicate in (1,2);
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
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



************************************ Overall performance;

******* This block calculates the Brier score - see Graf et al 1999;
title ' ';

****First for apparent validation i.e. model development;
*calculate weights using Kaplan-Meier for apparent validation;
proc lifetest data=rottx method=pl atrisk outsurv=outkm noprint;
        time survtime*status(1);
run;

*Create 3 groups - Group 1-Those who have the event up to fixed event time of 
*interest, Group 2 - those who go beyond fixed time (could be event or event 
*free), and Group 3- those censored up to fixed time 
*Only first 2 groups contribute to score but all to weights;
data rott_b;
	set rottx;
	*Must make the time just under 5 years since we need some time remaining
	*after timepoint of interest (5 years) and		 before administrative censoring 
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
proc phreg data=rott;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
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
	title 'Brier score';
	proc print;
run;

******* End of block;


********* This block bootstraps the 95% CI for the Brier score;
****First calculate weights in bootstrapped datasets;
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
proc lifetest data=rottapp method=pl outsurv=outkma noprint;
		by replicate;
        time survtime*status(1);
run;

*for bootstrap external validation;
proc lifetest data=rott method=pl outsurv=outkm_ext noprint;
        time survtime*status(1);
run;

*Create 3 groups - Group 1-Those who have the event up to fixed event time of 
*interest, Group 2 - those who go beyond fixed time (could be event or event 
*free), and Group 3- those censored up to fixed time 
*Only first 2 groups contribute to score but all to weights;
data rott_ba;
	set rottapp;
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
	set rott;
	if survtime<=4.99 and status=1 then cat=1;
	if survtime>4.99 then cat=2;
	if survtime<=4.99 and status=0 then cat=3;
	time=survtime;
run;

*Estimate survival at 5 years for apparent bootstrap;
proc phreg data=rottapp noprint;
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
	baseline covariates=rott_ba out=rott_bsa timelist=5 survival=fiveyrsurv/ 
		method=breslow;
run;

*Estimate survival at 5 years for external bootstrap;
proc phreg data=rottapp noprint;
	by replicate;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
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

proc phreg data=rottx;
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

*Apparent brier for null model;
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
	title 'Brier score and Scaled Brier';
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

proc phreg data=outboot noprint;
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
 	retain brier brier_lower95 brier_upper95 scaled_b scaledb_lower95 
		scaledb_upper95;
	merge scaledb confintr1_;
	by ind;
	drop ind sweight sbs;
	title 'Brier score and Scaled Brier score with 95% CI';
	proc print;
run;



******** Block to do Bootstrap to assess optimism in performance of the 
******** Scaled Brier;

* Begin with apparent bootstrap performance;
proc phreg data=rottapp noprint;
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
proc phreg data=rottapp noprint;
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
title 'Optimism Brier score and Scaled Brier score';
	proc print;
run;



** NOW PUT ALL TOGETHER TO GET THE OPTIMISM IN PERFORMANCE & 95% CI FOR 
** APPARENT PERFORMANCE;
data all1;
	*If you did not calculate all of the performance measures you will need to 
		delete the related dataset from the merge statement;
	merge finalipa harcon unocon tdunocon;
	by replicate;
run;

data stratos.all1;
	set all1;
run;

*calculate the average optimism in fit;
proc univariate data= all1 noprint;
	var OptBrier OptScaledB OptHarC OptUnoC OptUnoAUC;
	output out=avgopt mean=Mean_OptBrier Mean_OptScaleB Mean_OptHarC 
		Mean_OptUnoC Mean_OptUnoAUC;
run;

**** Manually enter the performance measures for each model;
data resultsdata;
	infile datalines delimiter=',';
	input OrigBrier OrigScBr OrigHarC OrigUnoC OrigUnoAUC;
	datalines;
	0.207, 0.155, 0.682, 0.682, 0.721
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
	drop ind;
run;

title;
title 'Estimate of internal validity penalising for overfitting';
proc print data=correctedperf;
run;

proc export data=correctedperf
  	outfile= "C:\Users\sme544\Bootstrap_Val" 
   	dbms=xlsx replace;
run;


title ' ';

************** Conduct DCA analysis;

*- see %stdca macro - get from https://www.mskcc.org/departments/epidemiology-
*biostatistics/biostatistics/decision-curve-analysis;

*let's add PGR as a new marker and calculate for basic and extended models;
data ber;
	set ffpgr;
	survtime1=survtime;
run;

proc phreg data=rott;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat / ties=efron rl;
	baseline out=origext covariates=ber survival=fiveyr timelist=5;
run;

proc phreg data=ffpgr;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes nodes1 gradecat pgr pgr1/ ties=efron rl;
	baseline out=newext covariates=origext survival=fiveyr_pgr timelist=5;
run;

data newext1;
	set newext;
	risk_orig=1-fiveyr;
	risk_new=1-fiveyr_pgr;
run;

*count number tp and fp for cut pt of 0.23;
data tp;
	set newext1;
	if risk_orig gt 0.23 then testpos=1;
	else testpos=0;
	if risk_new gt 0.23 then testpos_pgr=1;
	else testpos_pgr=0;
run;

proc freq data=tp;
	table testpos*status testpos_pgr*status;
run;



%put _user_;
goptions reset = all;

symbol1 i=join c=green;
symbol2 i=join c=red;
symbol3 i=join c=blue;
symbol4 i=join c=darkred;
symbol5 i=join c=gray;

%STDCA(data=NEWEXT1, out=survivalmult, outcome=STATUS, ttoutcome=SURVTIME1, 
	timepoint=5, predictors=RISK_ORIG);
%STDCA(data=NEWEXT1, out=survivalmult_new, outcome=STATUS, ttoutcome=SURVTIME1, 
	timepoint=5, predictors=RISK_NEW);

*sort by threshold variable;
proc sort data=survivalmult out=kmsort;
	by threshold;
run;

*rename the variables so that we know they are the kaplan meier estimates;
data kmsort; set kmsort(rename=(risk_orig=kmmodel all=kmall));
	LABEL kmmodel="Original model: Pr(Recurrence or death) at 5 years";
	LABEL kmall="Treat All";
RUN;

*sort by threshold variable;
proc sort data=survivalmult_new out=crsort;
	by threshold;
run;

*Rename the variables so that we know they are the competing model estimates;
data crsort; set crsort(rename=(risk_new=crmodel all=crall));
	label crmodel="Original + PGR: Pr(Recurrence or death) at 5 years";
run;

*merge nb data for original model plus model with pgr using threshold 
*probabilities;
data crsort;
	merge kmsort crsort;
	by threshold;
run;

data stratos.appnb;
	set crsort;
run;

******  Optional Block to create nicer DCA curves ************;

title;
footnote;
*- this allows editing of the .sge file!;
ods listing sge=on style=printer image_dpi=300 gpath='C:\Users\sme544';
ods graphics on / reset=all noborder outputfmt=tiff /*width=4in*/ 
	imagename="dca with spline pgr" antialias=off/*antialiasmax=*/;

*create decision curve with treat none, treat all, basic model, extended model;
proc sgplot data=crsort;
	yaxis values=(-0.05 to 0.45 by 0.05) label="Net Benefit";
	xaxis values=(0.0 to 1 by 0.1) label="Threshold Probability";
	keylegend "all" "orig" "new" "none" / location=inside down=4 position=topright;
	series y=kmall x=threshold / lineattrs=(color=black thickness=2 
			pattern=mediumdash) name="all" legendlabel="Treat All";
	series Y=kmmodel X=threshold / lineattrs=(color=green thickness=2 
			pattern=solid) name="orig" 
			legendlabel="Original model";
	series Y=crmodel X=threshold / lineattrs=(color=blue thickness=2 
			pattern=longdash) name="new" 
			legendlabel="Original model + PGR";
	series Y=none X=threshold / lineattrs=(color=red thickness=2 
			pattern=shortdash)  name="none" legendlabel="Treat None";
RUN;

ods graphics off;

******  End Optional Block ************;


