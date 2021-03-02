
/************************************************************************************************************
PROGRAM:	STDCA.sas
PROGRAMMER:	Daniel Sjoberg
DATE:		6/10/2013
UPDATED:	2/2/2015 by Emily Vertosick
NOTE:		stdca is an extension to dca for a survival-time endpoint. 
			The program calculates the points on a decision curve and optionally
			plots the decision curve, where <xvars> is the probability of failure 
			at time # associated with a single variable or multivariable
			model. 
*************************************************************************************************************/

%MACRO STDCA(	
	data=,					/*Name of input dataset*/
	out=,					/*Name of output dataset containing calculated net benefit*/
	outcome=, 				/*outcome variable, 1=event, 0=nonevent/censor*/
	ttoutcome=, 			/*time to event/censor*/
	timepoint=,				/*Timepoint of interest (i.e. risk at xx years)*/
	predictors=,			/*List variables with the predicted probabilities separated by a space*/
	xstart=0.01,			/*Low Threshold Value, this is the lowest probability the Net Benefit will be calculated for*/
	xstop=0.99,				/*High Threshold, this is the largest probability the Net Benefit will be calculated for*/
	xby=0.01,				/*By Value for Threshold Values, this is legnth of the interval at which the Net Benefit is calculated.*/
	harm=,					/*list of harms to apply to each predictor*/
	intervention=no,		/*calculate number of interventions avoided (yes/no)*/
	interventionper=100,	/*Number of intervetion per xx patients*/
	probability=,			/*list indicating whether each predictor is a probability*/
	graph=yes,				/*indicates if graph is requested or suppressed (yes/no)*/
	ymin=-0.05,				/*minimum net benefit that will be plotted*/
	interventionmin=0,		/*minimum reduction in interventions that will be plotted*/
	competerisk=no,			/*No for Kaplan-Meier, Yes for Competing Risk*/
	smooth=no,				/*use loess smoothing on decision curve graph (yes/no)*/
	/*GPLOT OPTIONS*/
	vaxis=,
	haxis=,
	legend=,
	plot_options=,
	plot_statements=
	);

DATA _NULL_;
	*removing multiple spaces in a row;
	CALL SYMPUTX("predictors",COMPBL("&predictors."));
	CALL SYMPUTX("harm",COMPBL("&harm."));
	CALL SYMPUTX("probability",UPCASE(COMPBL("&probability.")));

	*saving out number of predictors, harms, and probs specified;
	CALL SYMPUTX("varn",COUNTW("&predictors."," "));
	CALL SYMPUTX("harmn",COUNTW("&harm."," "));
	CALL SYMPUTX("probn",COUNTW("&probability."," "));
RUN;

/*Assigns a macro with a variable name for each predictor.*/
DATA _NULL_;
	%DO predvars=1 %TO &varn.;
		CALL SYMPUTX("var"||strip(put(&predvars.,9.0)),SCAN(COMPBL("&predictors."),&predvars.," "));
	%END;
RUN;

/*These error messages deal with necessary information being missing from the macro call.
This stops the macro if: data= is missing, outcome= is missing, predictors= is missing,
graph= is not "yes" or "no", intervention= is not "yes" or "no", model variable names are "all" or "none",
or if harm or probability is specified but there is not a harm or probability assigned for each predictor
referenced.*/

/*Checking that all required variables are defined*/
%IF %LENGTH(&data.)=0 OR %LENGTH(&outcome.)=0 OR %LENGTH(&ttoutcome.)=0 OR %LENGTH(&predictors.)=0 OR %LENGTH(&timepoint.)=0 %THEN %DO;
	%PUT ERR%STR()OR:  data, outcome, ttoutcome, timepoint and predictors must be specified in the macro call.;
	%GOTO QUIT;
%END;

/*Checking that graph and intervention options are correctly specified*/
%IF %UPCASE(&graph)^=NO & %UPCASE(&graph)^=YES %THEN %DO;
	%PUT ERR%STR()OR:  graph option must be YES or NO;
	%GOTO QUIT;
%END;

%IF %UPCASE(&intervention)^=NO & %UPCASE(&intervention)^=YES %THEN %DO;
	%PUT ERR%STR()OR:  intervention option must be YES or NO;
	%GOTO QUIT;
%END;

/*This code checks that competerisk=YES or NO if it is specified in the macro call.*/
%IF %LENGTH(&competerisk)>0 & %UPCASE(&competerisk)^=NO & %UPCASE(&competerisk)^=YES %THEN %DO;
	%PUT ERR%STR()OR:  if the competerisk option is specified, it must be YES or NO;
	%GOTO QUIT;
%END;

%END;

/*Check that the smooth option is correctly specified*/
%IF %UPCASE(&smooth)^=NO & %UPCASE(&smooth)^=YES & %LENGTH(&smooth)^=0 %THEN %DO;
	%PUT ERR%STR()OR:  smooth option must be YES or NO;
	%GOTO QUIT;
%END;

*If harm or probabilities specified, then the dimension must match predictors;
%IF (&harmn.^=&varn. AND &harmn.^=0) OR (&probn.^=&varn. AND &probn.^=0) %THEN %DO;
	%PUT ERR%STR()OR:  The specified number of harms and indicators of probability must be equal to the number of predictors specified.;
	%GOTO QUIT;
%END;

*Model variable names being checked cannot be equal to "all" or "none";
%DO name=1 %TO &varn.;
 	%IF %SYSFUNC(UPCASE(&&var&name.))=NONE OR %SYSFUNC(UPCASE(&&var&name.))=ALL %THEN %DO;
		%PUT ERR%STR()OR:  Variable names cannot be equal to "all" or "none";
		%GOTO QUIT;
	%END;
%END;

/*This code will generate an error message for incorrect values input for the "xstart", "xstop",
and "xby" options in the macro call. None of these values should be below 0 or above 1.*/
%IF %SYSEVALF(&xstart. < 0) OR %SYSEVALF(&xstart. > 1) OR %SYSEVALF(&xstop. < 0) OR %SYSEVALF(&xstop. > 1)
	OR %SYSEVALF(&xby. < 0) OR %SYSEVALF(&xby. > 1) %THEN %DO;
	%PUT ERR%STR()OR:  Values specified in xstart, xstop and xby options must be greater than 0 and less than 1.;
	%GOTO QUIT;
%END;

/*This code will generate an error message for negative values input for the "timepoint" option.*/
%IF %SYSEVALF(&timepoint. < 0) %THEN %DO;
	%PUT ERR%STR()OR:  Value specified in timepoint option must be greater than 0.;
	%GOTO QUIT;
%END;

/*This code is to check that dataset specified in "out=" option is a valid SAS dataset name, since
	otherwise the macro will run until it needs to save out this dataset and then stop.*/
%IF %LENGTH(&out)>0 & %SYSFUNC(MVALID(work,&out,data))^=1 %THEN %DO;
	%PUT ERR%STR()OR:  The name specified for the outcome dataset (&out) is not a valid SAS dataset name.;
	%GOTO QUIT;
%END;

/*These error messages deal with situations where all necessary information is specified in the
macro call but dataset, outcome and/or predictor variables do not exist as specified. This stops
the macro if: dataset does not exist, outcome variable does not exist in dataset, or any predictor
does not exist in dataset.*/

/*First, this checks that the dataset specified exists.*/
%IF %SYSFUNC(EXIST(&data)) %THEN %DO;

	/*If the dataset does exist, this checks that the outcome and ttoutcome variables exist in this data.*/
	%LET dsid=%SYSFUNC(OPEN(&data,i));
	%LET outcomecheck=%SYSFUNC(VARNUM(&dsid,&outcome));
	%LET ttoutcomecheck=%SYSFUNC(VARNUM(&dsid,&outcome));
	%LET close=%SYSFUNC(CLOSE(&dsid));

	/*If dataset exists but outcome variable is not in the data, print error and exit macro.*/
	%IF &outcomecheck.=0 %THEN %DO;
		%PUT ERR%STR()OR:  The outcome variable &outcome is missing from dataset &data.;
		%GOTO QUIT;
	%END;

	/*If dataset and outcome variable exist, but outcome variable is not in the data, print error and exit macro.*/
	%IF &ttoutcomecheck.=0 %THEN %DO;
		%PUT ERR%STR()OR:  The time to outcome variable &ttoutcome is missing from dataset &data.;
		%GOTO QUIT;
	%END;

	/*If dataset and outcome and ttoutcome variables exist, this checks that all predictor variables exist in this data.*/
	%ELSE %DO check = 1 %TO &varn.;
		%LET dsid=%SYSFUNC(OPEN(&data,i));
		%LET predcheck=%SYSFUNC(VARNUM(&dsid,&&var&check.));
		%LET rc=%SYSFUNC(CLOSE(&dsid));

		/*If dataset and outcome variable exist but any predictor variable is not in the data,
		print error and exit macro.*/
		%IF &predcheck.=0 %THEN %DO;
			%PUT ERR%STR()OR:  The predictor variable &&var&check. is missing from dataset &data.;
			%GOTO QUIT;
		%END;
	%END;
%END;

/*If the dataset does not exist, print error and exit macro.*/
%ELSE %DO;
	%PUT ERR%STR()OR: dataset &data does not exist.;
	%GOTO QUIT;
%END;

/*After checking that all required variables have been specified in the macro call and that the
dataset and outcome and predictor variables referenced all exist, continue with the rest of the
decision curve analysis.*/

*assigning each predictor, harm, and probability an ID and default value if not specified;
DATA _NULL_;
	%DO abc=1 %TO &varn.;
		CALL SYMPUTX("var"||strip(put(&abc.,9.0)),SCAN(COMPBL("&predictors."),&abc.," ") );
		CALL SYMPUTX("harm"||strip(put(&abc.,9.0)),COALESCE(SCAN(COMPBL("&harm."),&abc.," "),0));
		CALL SYMPUTX("prob"||strip(put(&abc.,9.0)),UPCASE(COALESCEC(SCAN(COMPBL("&probability."),&abc.," "),"YES")));
	%END;
RUN;


*deleting missing observations;
DATA stdcamacro_data;
	SET &data;
	IF NOT MISSING(&outcome.);
	IF NOT MISSING(&ttoutcome.);
	%DO abc=1 %TO &varn.;
		IF NOT MISSING(&&var&abc.);
	%END;

	*this variable is needed if the outcome is competing risk, and 
	a variable is included that needs to be converted to a probability 
	with Cox regression;
	outcome_binary=(&outcome.=1);
	KEEP outcome_binary &outcome. &ttoutcome. &predictors.;
RUN;

*creating dataset and macro variables with variable labels;
	PROC CONTENTS DATA=stdcamacro_data  OUT=stdcamacro_contents;
	RUN;
	DATA _NULL_ test;
		SET stdcamacro_contents;
		%DO abc=1 %TO &varn.;
			if STRIP(UPCASE(name))=STRIP(UPCASE("&&var&abc..")) then id=&abc.;
		%END;
		IF NOT MISSING(id);
		CALL SYMPUTX("varlab"||strip(put(id,9.0)),COALESCEC(label,name));
	RUN;

PROC SQL NOPRINT;
	*Getting number of observations;
	SELECT COUNT(*) INTO :n FROM stdcamacro_data;
	*Getting number of events;
	SELECT COUNT(*) INTO :eventn FROM stdcamacro_data WHERE &outcome.=1;
	*Getting number of non-events;
	SELECT COUNT(*) INTO :noneventn FROM stdcamacro_data WHERE &outcome.=0;
QUIT;

*Asserting outcome is coded as 0 or 1;
%IF %SYSFUNC(SUM(&noneventn., &eventn.)) ^= &n. & %UPCASE(&competerisk.)^=YES %THEN %DO;
	%PUT ERR%STR()OR:  &outcome. must be coded as 0 and 1;
	%GOTO QUIT;
%END;

*asserting all inputs are between 0 and 1 OR specified as non-probabilities.  
If not a probability, then converting it to prob with cox regression;
%DO abc=1 %TO &varn.;
	*checking range for probabilites;
	%IF &&prob&abc..=YES %THEN %DO;
	 	PROC SQL NOPRINT;
			SELECT MAX(&&var&abc..) INTO :varmax&abc. FROM stdcamacro_data;
			SELECT MIN(&&var&abc..) INTO :varmin&abc. FROM stdcamacro_data;
		QUIT;

		*any probabilities not between 0 and 1, then printing error;
	 	%IF %SYSEVALF(&&varmax&abc..>1) OR %SYSEVALF(&&varmin&abc..<0) %THEN %DO;
			%PUT ERR%STR()OR:  &&var&abc.. must be between 0 and 1 OR specified as a non-probability in the probability option;
			%GOTO QUIT;
		%END;
	%END;

	*if not probability, converting to prob with cox regression, and replacing original value with prob.;
	%IF &&prob&abc..=NO %THEN %DO;
		*estimating risk at timepoint with Cox model;
		PROC PHREG DATA=stdcamacro_data NOPRINT;
			MODEL &ttoutcome.*outcome_binary(0) = &&var&abc..; 
			*saving out survival estimates;
			BASELINE OUT=stdca_surv_&&var&abc.. COVARIATES=stdcamacro_data (keep=&&var&abc..) SURVIVAL=surv / NOMEAN METHOD=pl;
		RUN;

		PROC SQL NOPRINT UNDO_POLICY=NONE;
			*within each covariate value, estimating risk at timepoint (minimum value on survival curve before timepoint);
			CREATE TABLE stdca_surv_&&var&abc..2 AS
			SELECT DISTINCT 
					&&var&abc.., 1-MIN(surv) as &&var&abc.._cuminc 
			FROM stdca_surv_&&var&abc.. (where=(&ttoutcome.<=&timepoint.))
			GROUP BY &&var&abc..
			;

			*merging in cumulative incidence estimate and replacing previous covariate value with risk value.;
			CREATE TABLE stdcamacro_data (drop=&&var&abc.. rename=(&&var&abc.._cuminc=&&var&abc..)) AS
			SELECT s.*, r.&&var&abc.._cuminc
			FROM stdcamacro_data s
				LEFT JOIN stdca_surv_&&var&abc..2 r
				ON s.&&var&abc..=r.&&var&abc..
			;
		QUIT;
		%PUT WARN%STR()ING:  &&var&abc.. converted to a probability using Cox regression.  Due to linearity and proportional hazards assumptions, miscalibration may occur.;
	%END;
%END;

%IF %UPCASE(&competerisk.)=YES %THEN %DO;
	%stdca_crciest(	data=stdcamacro_data,		/*INPUT DATASET NAME*/
					outcome=&outcome.,			/*OUTCOME VARIABLE NAME (0=Event not observed, 1=Event Observed)*/
					ttoutcome=&ttoutcome.,		/*TIME TO OUTCOME OR CENSOR*/
					timepoint=&timepoint., 		/*TIME POINT CUMULATIVE INCIDENCE ESTIMATE*/
					macroname=pd				/*NAME OF MACRO NAME WITH CUMULATIVE INCIDENCE ESTIMATE*/
							);
%END;
%IF %UPCASE(&competerisk.)^=YES %THEN %DO;
	%stdca_kmciest(	data=stdcamacro_data,		/*INPUT DATASET NAME*/
					outcome=&outcome.,			/*OUTCOME VARIABLE NAME (0=Event not observed, 1=Event Observed)*/
					ttoutcome=&ttoutcome.,		/*TIME TO OUTCOME OR CENSOR*/
					timepoint=&timepoint., 		/*TIME POINT CUMULATIVE INCIDENCE ESTIMATE*/
					macroname=pd				/*NAME OF MACRO NAME WITH CUMULATIVE INCIDENCE ESTIMATE*/
							);
%END;

*This creates a new "xstop" so that when using "xstop" and "xby" options, the lines on the graph extend to the value of
"xstop" even if the last value of "xstart + xby" is greater than "xstop".;

DATA _NULL_;
	CALL SYMPUTX("xstop2",&xstop.);
	IF MOD((&xstop.-&xstart.),&xby.)~=0 THEN DO;
		CALL SYMPUTX("xstop2",&xstop.+&xstart.);
	END;	
RUN;

*creating dataset that is one line per threshold for the treat all and treat none strategies;
DATA stdcamacro_nblong (DROP=t);
	LENGTH model $100.;

	DO t=&xstart. TO &xstop2. BY &xby.;
		threshold=ROUND(t,0.00001);

		*creating the TREAT ALL row;
		model="all";
		nb=&pd. - (1-&pd.)*threshold/(1-threshold);
		output;

		*creating the TREAT NONE row;
		model="none";
		nb=0;
		output;
	END;
RUN;

/*Smoothing using LOESS*/

*ensure stdcamacro_models: datasets are empty;
PROC DATASETS LIB=WORK NOPRINT;
	DELETE stdcamacro_models:;
RUN;
QUIT;

*Looping over predictors and calculating net benefit for each of them.;
%DO abc=1 %TO &varn.;

	*Create stdcamacro_models&abc. dataset to start. These datasets were deleted above so that
	old datasets are not combined but an error is given if the dataset does not exist when
	trying to append data.;
	PROC SQL NOPRINT;
		CREATE TABLE WORK.stdcamacro_models&abc.
		 (model CHARACTER(80), threshold NUMERIC, nb NUMERIC);
	QUIT;

	*Looping over predictors and calculating net benefit for each of them.;
	%DO thresholdid=1 %TO %EVAL(%SYSFUNC(CEIL(%SYSEVALF((&xstop.-&xstart.)/&xby.)))+1);
		%LET threshold=%SYSEVALF((&xstart.-&xby.)+(&xby.*&thresholdid. ));

		*initializing macro vars to no longer exist;
		%SYMDEL px /NOWARN;
		%SYMDEL px_n /NOWARN;
		%SYMDEL px_totn /NOWARN;
		%SYMDEL pdgivenx /NOWARN;
		%SYMDEL nb /NOWARN;

		*calculating number of true and false positives;
	 	PROC SQL NOPRINT;
			SELECT COUNT(*) INTO :px_n FROM stdcamacro_data WHERE &&var&abc..>&threshold.;
			SELECT COUNT(*) INTO :px_totn FROM stdcamacro_data;
			SELECT &px_n./&px_totn. INTO :px FROM stdcamacro_data;
		QUIT;

		%IF &px.=0 %THEN %LET nb=.r;
		%ELSE %DO;
			%IF %UPCASE(&competerisk.)=YES %THEN %DO;
				%stdca_crciest(	data=stdcamacro_data (where=(&&var&abc..>&threshold.)),		/*INPUT DATASET NAME*/
								outcome=&outcome.,			/*OUTCOME VARIABLE NAME (0=Event not observed, 1=Event Observed)*/
								ttoutcome=&ttoutcome.,		/*TIME TO OUTCOME OR CENSOR*/
								timepoint=&timepoint., 		/*TIME POINT CUMULATIVE INCIDENCE ESTIMATE*/
								macroname=pdgivenx				/*NAME OF MACRO NAME WITH CUMULATIVE INCIDENCE ESTIMATE*/
										);
			%END;
			%IF %UPCASE(&competerisk.)^=YES %THEN %DO;
				%stdca_kmciest(	data=stdcamacro_data (where=(&&var&abc..>&threshold.)),		/*INPUT DATASET NAME*/
								outcome=&outcome.,			/*OUTCOME VARIABLE NAME (0=Event not observed, 1=Event Observed)*/
								ttoutcome=&ttoutcome.,		/*TIME TO OUTCOME OR CENSOR*/
								timepoint=&timepoint., 		/*TIME POINT CUMULATIVE INCIDENCE ESTIMATE*/
								macroname=pdgivenx				/*NAME OF MACRO NAME WITH CUMULATIVE INCIDENCE ESTIMATE*/
										);
			%END;

			%IF &pdgivenx.=. %THEN %LET nb=.f;
			%IF &pdgivenx.^=. %THEN %LET nb=&pdgivenx.*&px. - (1-&pdgivenx.)*&px.*&threshold./(1-&threshold.);
		%END;

		*creating one line dataset with nb.;
		DATA stdcamacro_temp;
			length model $100.;
			model="&&var&abc..";
			threshold=ROUND(&threshold.,0.00001);
			nb=&nb.;
		RUN;

		*creating dataset with nb for models only.;
		DATA stdcamacro_models&abc.;
			SET stdcamacro_models&abc. stdcamacro_temp;
		RUN;

		*deleting results dataset;
		PROC DATASETS LIB=WORK NOPRINT;
			DELETE stdcamacro_temp;
		RUN;
		QUIT;

	%END;

	/*After running for all thresholds for each predictor and saving each predictor dataset separately, then smooth*/

	%IF %UPCASE(&smooth.)=YES %THEN %DO;

		PROC LOESS DATA=stdcamacro_models&abc.;
			MODEL nb=threshold / ALL;
			ODS OUTPUT OutputStatistics=smooth_&abc.;
		RUN;

		PROC DATASETS LIB=WORK NOPRINT;
			DELETE stdcamacro_models&abc.;
		RUN;

		DATA stdcamacro_models&abc.(keep=threshold nb model);
			SET smooth_&abc.(rename=(Pred=nb));
			model="&&var&abc..";
			FORMAT nb 5.2;
		RUN;

		PROC SORT DATA=stdcamacro_models&abc. NODUPKEY;
			BY threshold;
		RUN;

	%END;

%END;

/*Merge data from predictors together.*/
DATA stdcamacro_nblong_final;
	SET stdcamacro_nblong stdcamacro_models:;
RUN;

/*STDCA warnings*/

PROC SQL NOPRINT;
	CREATE TABLE stdcamacro_warnings AS
	SELECT DISTINCT model, nb, min(threshold) as threshold
	FROM stdcamacro_nblong_final
	WHERE nb in (.f .r)
	GROUP BY model
	ORDER BY model, threshold
	;

	SELECT COUNT(*) into :warn_tot 
	FROM stdcamacro_nblong_final
	WHERE nb in (.f .r)
	;

QUIT;

DATA _NULL_;
	SET stdcamacro_warnings;
	BY model threshold;

	IF FIRST.threshold; 
	count+1;
	CALL SYMPUTX("warning_model"||strip(put(count,BEST12.)),model);
	CALL SYMPUTX("warning_t"||strip(put(count,BEST12.)),put(threshold,BEST12.));
	CALL SYMPUTX("warning_type"||strip(put(count,BEST12.)),put(nb,BEST12.));

	CALL SYMPUTX("warning_n",put(count,BEST12.));
RUN;

*only print warnings if there is at least one.;
%IF &warn_tot.>0 %THEN %DO;
%DO abc=1 %TO &warning_n.;
	%IF &&warning_type&abc..=R %THEN %PUT 
		WARN%STR()ING:  &&warning_model&abc..: No observations with risk greater than &&warning_t&abc.., and therefore net benefit not calculable in this range.;
	%IF &&warning_type&abc..=F %THEN %PUT 
		WARN%STR()ING:  &&warning_model&abc..: No observations with risk greater than &&warning_t&abc.. that have followup through the timepoint selected, and therefore net benefit not calculable in this range.;
%END;
%END;

*making NB dataset one line per threshold probability;
PROC SORT DATA=stdcamacro_nblong_final;
	BY threshold;
RUN;

PROC TRANSPOSE DATA=stdcamacro_nblong_final OUT=stdcamacro_nbT (DROP=_name_);
	BY threshold;
	ID model;
	VAR nb;
RUN;

DATA stdcamacro_nb &out.;
	SET stdcamacro_nbT;
	WHERE threshold<=&xstop2.;
	*Changed this from "&xstop." to "&xstop2." to account for situations where &xstart. + &xby.*_N_ is greater than &xstop.;

	*applying variable labels;
	label threshold="Threshold Probability";
	label all="Net Benefit: Treat All";
	label none="Net Benefit: Treat None";

	%DO abc=1 %TO &varn.;
		*correcting NB if harms are specified and labelling Net Benefit;
		label &&var&abc..="Net Benefit: &&varlab&abc..";
		%IF %LENGTH(&harm.)>0 %THEN %DO;
			&&var&abc..=&&var&abc.. - &&harm&abc..;
			label &&var&abc..="Net Benefit: &&varlab&abc.. (&&harm&abc.. harm applied)";
		%END;

		
		*transforming to interventions avoided;
		&&var&abc.._i=(&&var&abc.-all)*&interventionper./(threshold/(1-threshold));
		label &&var&abc.._i="Intervention: &&varlab&abc..";
		%IF %LENGTH(&harm.)>0 %THEN %DO;
			label &&var&abc.._i="Intervention: &&varlab&abc.. (&&harm&abc.. harm applied)";
		%END;
		
		*label smoothed net benefit;
		%IF %UPCASE(&smooth)=YES & %LENGTH(&harm.)>0 %THEN %DO;
			label &&var&abc..="Smoothed Net Benefit: &&varlab&abc.. (&&harm&abc.. harm applied)";
		%END;
		%ELSE %IF %UPCASE(&smooth)=YES & %LENGTH(&harm.)=0 %THEN %DO;
			label &&var&abc..="Smoothed Net Benefit: &&varlab&abc..";
		%END;

	%END;

RUN;


*quitting macro if no graph was requested;
%IF %UPCASE(&graph.)=NO %THEN %GOTO QUIT;

***************************************;
********  PLOTTING DCA    *************;
***************************************;
*CREATING VARIABLE LIST FOR GPLOT;
%IF %UPCASE(&intervention.)=NO %THEN %DO;
	%LET plotlist=all none &predictors.;
	%LET ylabel=Net Benefit;
	%LET plotrange=&ymin. <= col1;
%END;
%ELSE %DO;
	%LET ylabel=Net reduction in interventions per &interventionper. patients;
	%LET plotrange=col1 >= &interventionmin.;
	%DO abc=1 %TO &varn.;
		%IF &abc.=1 %THEN %LET plotlist=&&var&abc.._i;
		%ELSE %LET plotlist=&plotlist. &&var&abc.._i;
	%END;
%END;

*transposing data to one line per threshold for model type;
PROC TRANSPOSE DATA=stdcamacro_nb OUT=stdcamacro_plot;
	BY threshold;
	VAR &plotlist.;
RUN;

**** ORDERING CATEGORIES IN GPLOT ****;
*This code keeps the lines / categories ordered in the legend in the order they were entered into the STDCA macro statement.;

*labeling transpose variables;
DATA stdcamacro_plot2;
	SET stdcamacro_plot;
	label	col1="&ylabel."
			_label_="Model Label"
			_name_="Model";

	*setting variables outside plotrange to missing;
	IF NOT (&plotrange.) THEN col1=.;

	*This creates a numeric variable that corresponds to the order that the variables were
	entered into the STDCA macro command.;
	%DO order=1 %TO &varn.+2;
		piece&order.=SCAN("all none &predictors.",&order.," ");
		IF _name_=piece&order. THEN ordernum=&order.;
	%END;

RUN;

*create dataset to hold format for "ordernum" variable for graph;
DATA cntlin(
	KEEP=fmtname start label);
	SET stdcamacro_plot2(RENAME=(_LABEL_=label ordernum=start));
	fmtname="order";
RUN;

*sort format dataset and keep unique observations only;
PROC SORT DATA=cntlin OUT=cntlin NODUPKEYS;
	BY start;
RUN;

*load format for order number variable;
PROC FORMAT CNTLIN=cntlin;
RUN;

*drop unnecessary "piece*" variables and format "ordernum" variable for graph legend;
DATA stdcamacro_plot2(DROP=piece:); SET stdcamacro_plot2;
	FORMAT ordernum $order.;
RUN;

PROC SORT DATA=stdcamacro_plot2;
	BY threshold ordernum;
RUN;

/*Plotting DCA*/

PROC GPLOT DATA=stdcamacro_plot2;
	AXIS1 &vaxis. LABEL=(ANGLE=90) MINOR=NONE;
	AXIS2 &haxis. MINOR=NONE;
	LEGEND1 &legend. LABEL=NONE FRAME;

	PLOT col1*threshold=ordernum / SKIPMISS LEGEND=LEGEND1 VAXIS=AXIS1 HAXIS=AXIS2 &plot_options.;
	SYMBOL INTERPOL=JOIN;

	&plot_statements.;
RUN;
QUIT;

%PUT _USER_;

*location label for quitting macro early;
%QUIT:

/*deleting all macro datasets*/
/*
PROC DATASETS LIB=WORK;
	DELETE stdcamacro_:;
RUN;
QUIT;
*/
%MEND;


/*cumulative incidence estimates via Kaplan Meier*/
%MACRO stdca_kmciest(	data=,		/*INPUT DATASET NAME*/
						outcome=,	/*OUTCOME VARIABLE NAME (0=Event not observed, 1=Event Observed)*/
						ttoutcome=,	/*TIME TO OUTCOME OR CENSOR*/
						timepoint=, /*TIME POINT CUMULATIVE INCIDENCE ESTIMATE*/
						macroname=	/*NAME OF MACRO NAME WITH CUMULATIVE INCIDENCE ESTIMATE*/
						);
	ODS LISTING CLOSE;
	*get KM estimates;
	PROC LIFETEST data=&data. TIMELIST=(&timepoint.);
		TIME &ttoutcome.*&outcome.(0);
		ODS OUTPUT ProductLimitEstimates=stdcamacro_kmestall;
	RUN;
	ODS LISTING;


	*making the macro global to be passed outside of the KM estimate macro;
	%GLOBAL &macroname.;
	PROC SQL NOPRINT;
		SELECT failure INTO :&macroname. FROM stdcamacro_kmestall;
	QUIT;
%MEND;

/*cumulative incidence estimates in presence of competing risks*/
%MACRO stdca_crciest(	data=,		/*INPUT DATASET NAME*/
						outcome=,	/*OUTCOME VARIABLE NAME (0=Event not observed, 1=Event Observed)*/
						ttoutcome=,	/*TIME TO OUTCOME OR CENSOR*/
						timepoint=, /*TIME POINT CUMULATIVE INCIDENCE ESTIMATE*/
						macroname=	/*NAME OF MACRO NAME WITH CUMULATIVE INCIDENCE ESTIMATE*/
						);
	*the CIF macro cannot input a dataset name called 'alldata (where=(time>200))',
	so we must create another dataset first that is the appropriate subset;
	DATA stdca_crciest0;
		SET &data.;
	RUN;

	%CIF(data=stdca_crciest0, out=stdca_crciest, time=&ttoutcome., status=&outcome., event=1, options=NOPLOT);

	PROC SQL NOPRINT;
		SELECT count(*) INTO :len 
		FROM stdca_crciest
		WHERE &ttoutcome.<&timepoint.
		;
		SELECT count(*) INTO :gtn 
		FROM stdca_crciest
		WHERE &ttoutcome.>&timepoint.
		;
		SELECT max(CIF) INTO :crcuminc_temp 
		FROM stdca_crciest
		WHERE &ttoutcome.<=&timepoint.
		;
	QUIT;

	%GLOBAL &macroname.;
	DATA _NULL_;
		IF &len.>0 & &gtn.>0 THEN CALL SYMPUTX("&macroname.",PUT(&crcuminc_temp.,BEST12.));
		ELSE IF &len.=0 THEN CALL SYMPUTX("&macroname.","0");
		ELSE IF &gtn.=0 THEN CALL SYMPUTX("&macroname.",".");
	RUN;
%MEND;

