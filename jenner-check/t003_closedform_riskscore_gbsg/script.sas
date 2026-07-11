* Adapted from SAS_code/04_external_validation_when_model_estimates_and_baseline_survival_available.sas
* (lines 47-120: apply the published Cox model equation - baseline
* survival S0 plus the linear predictor coefficients from the
* development paper - directly to the external GBSG cohort, without
* needing the original development data or refitting anything. This
* is the external-validation-when-only-the-model-equation-is-available
* scenario the script is named for).
* Runs against a 30-row de-identified sample of the repository own
* Data/gbsg.csv (15 events / 15 censored) inlined as DATALINES,
* instead of PROC IMPORT from the original author local Windows folder.;

data gbsg;
    infile datalines dsd dlm=' ' truncover;
    input pid age meno size grade nodes pgr er hormon rfstime status;
datalines;
383 65 1 20 1 2 912 606 0 991 1
1009 49 0 45 3 6 0 22 0 375 1
777 49 0 55 2 7 0 0 0 679 1
131 59 1 30 2 8 0 0 1 238 1
757 25 0 29 2 3 0 0 0 343 1
782 45 0 29 2 1 0 0 0 578 1
96 63 1 30 2 7 0 0 1 698 1
433 54 1 23 3 10 13 6 0 307 1
500 35 0 35 2 4 0 0 1 1183 1
1762 42 0 25 2 7 0 2 0 249 1
618 52 0 27 2 5 0 8 0 336 1
512 67 1 65 2 8 0 0 0 745 1
1677 48 0 16 3 10 0 0 0 525 1
179 36 0 30 3 2 0 0 0 956 1
820 32 0 57 3 24 0 13 0 448 1
1007 56 1 20 1 1 2 334 1 1756 0
1525 57 1 35 3 1 1490 209 1 1342 0
694 57 1 55 3 6 22 186 0 859 0
1566 47 0 40 2 1 0 90 0 276 0
1391 40 0 30 3 3 2 16 0 1459 0
468 66 1 20 2 9 1 11 1 1977 0
1081 57 1 7 2 1 1 1 1 2030 0
1549 62 1 35 2 1 2 70 0 733 0
130 65 1 30 2 5 0 36 1 1855 0
1523 45 0 45 3 2 0 0 0 1469 0
1065 65 1 17 1 1 935 200 0 967 0
723 50 0 18 3 1 1 4 0 17 0
496 33 0 19 2 2 0 0 1 2128 0
1570 58 1 52 3 5 0 0 0 1177 0
1249 54 1 30 3 3 1 0 0 1089 0
;
run;

*** Estimating Baseline Survival Function under PH;
** Here baseline is where all predictors=0;
** Here S0 at 5 years =  0.80385;
** For model with PGR added, S0 = 0.76055;
** S(t)=S0(t)**exp(pi1X1+pi2X2+pi3X3+...+piPXP)=S0(t)**exp(PI);

* code up predictors and apply the original model;
data l2gbsg_va;
	set gbsg;
	if size le 20 then sizecat=0;
	if 21 <= size <= 50 then sizecat=1;
	if size >50 then sizecat=2;
	gradecat = grade-1;
	if grade in (1,2) then gradecat=0;
	if grade = 3 then gradecat=1;
    survtime = rfstime/365.25;  * days to years;
	* Winsorise nodes and PGR to the 99th percentile to deal with very large
	* influential values;
	if pgr>1360 then pgr=1360;
	if (nodes > 19) then nodes=19;
	* use formula to code up spline term for PGR;
	PGR1=max(PGR/61.81,0)**3
		+(41*max((PGR-486)/61.81,0)**3 -486*max((PGR-41)/61.81,0)**3)/445;
	* use formula to code up spline term for number of nodes;
	nodes1=max(nodes/4.33,0)**3
		+(max((nodes-9)/4.33,0)**3 -9*max((nodes-1)/4.33,0)**3)/8;
	* now calculate the prob of survival to 5 years;
	sizepi = 0;
	if sizecat=1 then sizepi=0.342;
	if sizecat=2 then sizepi=0.574;
	gradepi = 0;
	if gradecat=1 then gradepi=0.362;
	pi = sizepi + 0.304*nodes - 0.811*nodes1 + gradepi;
	psurv = 0.80385**exp(pi);
	prob = 1- psurv;
	* now calculate the prob of survival to 5 years for model with pgr included;
	sizepi2 = 0;
	if sizecat=1 then sizepi2=0.320;
	if sizecat=2 then sizepi2=0.554;
	gradepi2 = 0;
	if gradecat=1 then gradepi2=0.305;
	pi2 = sizepi2 + 0.305*nodes - 0.820*nodes1 + gradepi2 - 0.003*pgr + 0.013*pgr1;
	psurv2 = 0.76055**exp(pi2);
	prob2 = 1-psurv2;
	keep pid sizecat nodes nodes1 gradecat age pgr pgr1 survtime status
	sizepi gradepi pi psurv prob sizepi2 gradepi2 pi2 psurv2
		prob2;
run;

*** Descriptive statistics;
proc freq data=l2gbsg_va;
	tables sizecat gradecat status;
run;

proc univariate data=l2gbsg_va;
	var age pgr nodes;
run;

*use reverse Kaplan-Meier to obtain median follow-up time;
proc lifetest data=l2gbsg_va method=pl atrisk;
        time survtime*status(1);
run;

*Administrative censor at 5 years as development cohort;
data l2gbsg_val;
	set l2gbsg_va;
	if survtime > 5 then status=0;
	if survtime > 5 then survtime=5;
	survtime1=survtime;

	*CLL is the Outcome-Martingale residual, or expected number of events;
	CLL = LOG(-LOG(PSURV));
	CLL_ = LOG(-LOG(PSURV2));
RUN;

proc print data=l2gbsg_val (obs=10);
	var pid prob prob2 survtime status;
run;
