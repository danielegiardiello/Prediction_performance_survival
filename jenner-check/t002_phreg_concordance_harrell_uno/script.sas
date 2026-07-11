* Adapted from SAS_code/01_model_development_internal_validation.sas
* (lines 48-96, 235-241, 455-469: read data, derive sizecat/gradecat,
* Winsorize nodes, administrative censoring at 5 years, then Harrells C
* and Unos C - the two concordance statistics this repo compares for
* internal validation of the Cox model).
* Runs against a 30-row de-identified sample of the repository own
* Data/rotterdam.csv (15 events / 15 censored) inlined as DATALINES,
* instead of PROC IMPORT from the original author local Windows folder.;

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
	label survtime = 'Survival time (years)';
run;

data r1(rename=(nod=nodes)); set r0;
	* Winsorise nodes to the 99th percentile;
	if (nodes > 19) then nod = 19; else nod = nodes;
	drop nodes;
run;

* Administrative censor at 5 years since this is the prediction horizon;
data rott;
	set r1;
	if survtime > 5 then status=0;
	if survtime > 5 then survtime=5;
run;

* Harrells C - need tau to equal event time of interest;
title 'Apparent Harrells C';
proc phreg data=rott concordance=harrell(se) tau=5;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes gradecat / ties=efron rl;
run;
title;

* Unos C - need tau to equal event time of interest;
title 'Apparent Unos C';
proc phreg data=rott concordance=uno(se seed=8754 iter=50) tau=5;
	class sizecat (ref='<=20') gradecat (ref=first);
	model survtime*status(0)=sizecat nodes gradecat / ties=efron rl;
run;
title;
