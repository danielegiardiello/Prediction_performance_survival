* Adapted from SAS_code/03_external_validation_when_original_model_dataset_available.sas
* (lines 47-100, 213-244: fit the Cox model on the development
* (Rotterdam) data, then use PROC PHREG BASELINE COVARIATES= to score
* the fitted linear predictor onto the external (GBSG) validation
* cohort - the external-validation-when-the-original-development-
* -dataset-is-available scenario this script is named for).
* Runs against 30-row de-identified samples of the repository own
* Data/rotterdam.csv and Data/gbsg.csv inlined as DATALINES, instead
* of PROC IMPORT from the original author local Windows folder.;

****** Read in the Rotterdam dataset (development, sampled, inlined);
data rotterdam;
    infile datalines dsd dlm=' ' truncover;
    input pid year age meno size $ grade nodes pgr er hormon chemo
          rtime rtime_status dtime death;
    recur = rtime_status;
datalines;
1336 1990 44 0 <=20 3 0 0 6 0 0 992.00 1 2837.00 1
1715 1987 47 0 20-50 3 0 122 33 0 0 3556.00 1 3858.00 0
2101 1989 66 1 20-50 3 0 1 162 0 0 527.00 1 922.00 1
1991 1987 60 1 20-50 3 3 5 0 0 0 116.00 1 273.00 1
1368 1987 61 1 <=20 3 1 24 192 0 0 742.00 1 3663.00 0
1372 1989 43 0 20-50 3 13 61 114 0 0 294.00 1 991.00 1
2372 1985 46 0 >50 3 1 1230 98 0 0 2727.00 1 2839.00 1
1371 1989 79 1 >50 3 1 431 356 0 0 504.00 1 741.00 1
2106 1984 64 1 20-50 3 2 0 7 0 0 324.00 1 451.00 1
1344 1990 40 0 20-50 3 0 6 0 0 0 420.00 1 1101.00 1
2617 1990 55 1 20-50 3 12 1 1 0 0 151.00 1 332.00 1
1607 1987 65 1 >50 2 15 55 55 0 0 479.00 1 2351.00 1
1346 1986 45 0 20-50 3 0 2 5 0 0 334.00 1 2417.00 1
1454 1984 50 1 20-50 3 0 0 5 0 0 1124.00 1 1849.00 1
1331 1990 67 1 <=20 3 0 0 303 0 0 1472.00 1 3046.00 0
51 1989 55 1 20-50 2 0 101 273 0 0 3709.00 0 3709.00 0
390 1989 75 1 <=20 3 2 363 758 0 0 2712.00 0 2712.00 0
20 1993 62 1 <=20 3 0 22 84 0 0 2498.00 0 2498.00 0
922 1991 33 0 >50 2 1 230 67 0 1 1747.00 0 1747.00 0
52 1992 29 0 20-50 3 4 4 10 0 1 3008.00 0 3008.00 0
48 1991 59 1 20-50 3 2 69 73 1 0 2209.00 0 2209.00 0
38 1982 45 0 20-50 3 0 55 0 0 0 5775.00 0 5775.00 0
456 1990 57 1 <=20 2 0 0 13 0 0 2908.00 0 2908.00 0
783 1983 45 0 20-50 3 0 0 26 0 0 3446.00 0 3446.00 0
915 1991 68 1 20-50 2 0 304 276 0 0 2787.00 0 2787.00 0
1291 1988 57 1 <=20 3 9 2 5 0 0 3892.00 0 3892.00 0
17 1993 61 1 20-50 3 4 0 604 1 0 2734.00 0 2734.00 0
2964 1984 74 1 <=20 3 0 12 152 0 0 2033.00 0 2033.00 1
1292 1991 74 1 20-50 2 0 180 625 0 0 2817.00 0 2817.00 0
1046 1993 56 1 <=20 3 0 0 244 0 0 2524.00 0 2524.00 0
;
run;

****** Read in the German (GBSG) dataset (external validation, sampled, inlined);
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

*** Derive the modeling covariates on the development (Rotterdam) data;
proc format;  * make them print in a nice order;
    value sizef 0 = "<=20"
                1 = "20-50"
                2 = ">50";
    value gradef 0 = "1-2"
                 1 = "3";
run;

data r0(rename=(pr=pgr)); set rotterdam;
    format sizecat sizef. gradecat gradef.;
    sizecat =  1* (size= "20-50")  + 2*(size = ">50");
    gradecat = grade -2;  * rotterdam has only grade 2 and 3 subjects;
    if (recur = 1) then do;
        survtime = rtime/365.25;  * days to years;
        status = recur;
    end;
    else do;
        survtime = dtime/365.25;
        status = death;
    end;
    if (pgr > 1360) then pr = 1360; else pr = pgr;
	drop pgr;
run;

data ffpgr(rename=(nod=nodes)); set r0;
	* Winsorise nodes to the 99th percentile;
	if (nodes > 19) then nod = 19; else nod = nodes;
	drop nodes;
run;

*** Derive the same covariates on the external (GBSG) data;
data gbsg_val;
	set gbsg;
	if size le 20 then sizecat=0;
	if 21 <= size <= 50 then sizecat=1;
	if size >50 then sizecat=2;
	gradecat = grade-1;
	if grade in (1,2) then gradecat=0;
	if grade = 3 then gradecat=1;
    survtime = rfstime/365.25;  * days to years;
	if (nodes > 19) then nodes=19;
run;

* Fit the Cox model on the development data, then use BASELINE
* COVARIATES= to score the linear predictor onto the external GBSG data
* - this is the core of external validation when the original
* development dataset is available. TIMELIST=5 (the papers 5-year
* prediction horizon) is added to the original BASELINE call so it
* returns exactly one scored row per covariate pattern in input order,
* which the following merge relies on;
title 'Fit model on Rotterdam, score onto GBSG';
proc phreg data=ffpgr;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat gradecat nodes / ties=efron rl;
	baseline covariates=gbsg_val out=gbsg_valrd xbeta=xb timelist=5;
run;
title;

* Merge the scored linear predictor back onto the external cohort.
* gbsg_valrd carries one row per gbsg_val row, in the same order, so a
* positional merge recovers pid alongside the scored xb;
data gbsg_valrd2;
	merge gbsg_val gbsg_valrd;
run;

proc print data=gbsg_valrd2 (obs=10);
	var pid sizecat gradecat nodes survtime status xb;
run;
