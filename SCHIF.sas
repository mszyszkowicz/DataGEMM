/*************************************************************************************************
*
* Nonlinear concentraton-response function: model fitting and plotting procedure
*
* Date: Sept 19, 2016
* Version 2.17
*
* Purpose: This macro will fit a series of nonlinear concentration-response functions using 
*		   cox proportional hazards models to identify optimal nonlinear relationship
*
* It will produce three outputs in the file directory specified by users: 
*
*	(1) a graph (.png) showing nonlinear relationship based on an optimal model, and if pre-specified, it will overlay an ensemble model,  
*	    based on either all models examined or the best 3 models, according to their -2LogLik
* 	(2) a summary table (.pdf) listing coefficient, standard error, loglik, and function form of each cox model that was fitted
*   (3) a summary table (.pdf) showing descriptive stats of air pollution exposure variable from the original and the trimmed datasets 
*
* For more details about this modeling approach, please refer to the accompanying paper Masoud, Szyszkowicz,...,Burnett et al 2016 Air Quality, Atmosphere and Health
*
* Should there be any question with this macro, please contact:
* 	Drs Hong Chen (Hong.Chen@oahpp.ca) and Rick Burnett (Rick.Burnett@hc-sc.gc.ca)
*
*
* Please read the following notes before running the macro
*
* 1. All categorical variables should to be prepared as a series of dichotomous variables prior to calling macro
*
* 2. Exposure variable need to be define by using "fitvar" parameter
*
* 3. The following 3 parameters are optional: "strata", "label_exposure", and "label_unit"
*
* 4. To run time-fixed cox model, users need to specify "case" and "time" parameters, 
*		and omit "start" and "stop" parameters
*
* 5. To run time-varying cox model, users need to specify "case", "start", and "stop" parameters,
*		and omit "time" parameter
* 
* 6. If users do not wish to translate data, please specify translate=N or NO
*
* 7. By default, the option "best_model" is set to be 1. This means that the program will display concentration-response curve for best fitting model only.  
*
*    7.1. Set "best_model"=2, the program will display concentration-response curve from ensemble models based on all models examined.
*
*    7.2. Set "best_model"=3, the program will display two concentration-response curves by overlaying best fitting model and ensemble models based on best 3 models examined.
*
*    7.3. Set "best_model"=4, the program will display two concentration-response curves by overlaying best fitting model and ensemble models based on ALL models examined.
*
* 8. By default, the option "lowest" is set to be "minimum". This means that the concentration-response curve covers the full range of the exposure variable. 
*
*    8.1. Set "lowest"=1, the program will display concentration-response curve from 1 to max concentration of the exposure variable.
*
* 9. By default, the option "x_min" is set to be 1. This means that the concentration-response plot has x axis with min=1. 
*
*    9.1. Set "x_min" as any positive number, x axis will start from this user-defined value.
*
*    9.2. Set "x_min" as any NON-numeric character including blank (eg, x_min=, x_min=NO, x_min=Null), x axis will start from the minimum of exposure variable.
*
* 10. By default, the option "export" is set to be NULL. This means that this program will NOT produce any matrix for HRs derived over the range of AP exposure. 
*
*    10.1. Set "export" as YES or Y, this program will produce two csv files: 
*			(1) a table "export_hr.csv" presenting HRs derived from ~1000 curves at a series of locations: 0% to 100% by 1% over the range of AP exposure, and
*			(2) a table "ap_exp.csv" presenting concentrations at 0% to 100% by 1% over the range of AP exposure.
*
*
* To illustrate how to use this macro, 7 examples are given below:
* 
* Example 1: (translate data, fit time-fixed cox models, with strata variable, only show optimal model, c-r curve from min to max)
*
* %fitap(datain=cohort_ABC,perc_trim=1,translate=yes,dataout=allmodels,case=status,time=time,fitvar=no2,covvars=age sex, strata=inst, 
*	label_exposure=NO2, label_unit=ppb, output_path=D:\Working directory\Cohort);
* 
* Example 2: (translate data, fit time-fixed cox models, without strata variable, only show optimal model, c-r curve from min to max)
*
* %fitap(datain=cohort_ABC,perc_trim=1,translate=y,dataout=allmodels,case=status,time=time,fitvar=no2,covvars=age sex, strata=, 
*	label_exposure=NO2, output_path=D:\Working directory\Cohort);
*
* Example 3: (translate data, fit time-fixed cox models, with two strata variables, only show optimal model, c-r curve from min to max)
*
* %fitap(datain=cohort_ABC,perc_trim=5,translate=YES,dataout=allmodels,case=death,time=time,fitvar=pm25,covvars=ses,strata=age sex, 
*	label_exposure=pm25, label_unit=ug/m3, output_path=D:\Working directory\Cohort);
*
* Example 4: (translate data, fit time-varying cox models, with strata variable, only show optimal model, c-r curve from min to max)
*
* %fitap(datain=cohort_ABC,perc_trim=1,dataout=allmodels,case=status,start=T1,stop=T2,fitvar=pm25,covvars=sex ses,strata=age, 
*	label_exposure=pm25, label_unit=ug/m3, output_path=D:\Working directory\Cohort);
*
* Example 5: (no translate data, fit time-varying cox models, with strata variable, only show optimal model, c-r curve from min to max)
*
* %fitap(datain=cohort_ABC,perc_trim=1,translate=N,dataout=allmodels,case=Status,start=TStart,stop=TStop,fitvar=pm25,covvars=sex ses income bmi,strata=age, 
*	label_exposure=pm25, label_unit=ug/m3, output_path=D:\Working directory\Cohort);
*
* Example 6: (no translate data, fit time-varying cox models, with strata variable, c-r curve from min to max, only show ensemble model based on all models examined)
* 
* Note that to display only ensemble model based on all models examined, please define best_model=2
*
* %fitap(datain=cohort_ABC,perc_trim=1,translate=no,dataout=allmodels,case=Status,start=TStart,stop=TStop,fitvar=pm25,covvars=sex ses income bmi,strata=age, 
*	label_exposure=pm25, label_unit=ug/m3, best_model=2, output_path=D:\Working directory\Cohort);
* 
* Example 7: (no translate data, fit time-varying cox models, with strata variable, c-r curve from min to max, overlay optimal model with ensemble model based on all models examined)
* 
* Note that to overlay best fitting model with ensemble model based on all models examined, please define best_model=4
*
* %fitap(datain=cohort_ABC,perc_trim=1,translate=no,dataout=allmodels,case=Status,start=TStart,stop=TStop,fitvar=pm25,covvars=sex ses income bmi,strata=age, 
*	label_exposure=pm25, label_unit=ug/m3, output_path=D:\Working directory\Cohort, best_model=4);
* 
* Example 8: (no translate data, fit time-varying cox models, with strata variable, c-r curve from min to max, overlay optimal model with ensemble model based on best 3 models examined)
* 
* Note that to overlay best fitting model with ensemble model based on best 3 models examined, please define best_model=3
*
* %fitap(datain=cohort_ABC,perc_trim=1,translate=No,dataout=allmodels,case=Status,start=TStart,stop=TStop,fitvar=pm25,covvars=sex ses income bmi,strata=age, 
*	label_exposure=pm25, label_unit=ug/m3, output_path=D:\Working directory\Cohort, best_model=3);
*
* Example 9: (no translate data, fit time-varying cox models, with strata variable, c-r curve from 1 to max, overlay optimal model with ensemble model based on all models examined)
* 
* %fitap(datain=cohort_ABC,perc_trim=1,lowest=1,translate=No,dataout=allmodels,case=Status,start=TStart,stop=TStop,fitvar=pm25,covvars=sex ses income bmi,strata=age, 
*	label_exposure=pm25, label_unit=ug/m3, output_path=D:\Working directory\Cohort, best_model=4);
*
* Example 10: (no translate data, fit time-varying cox models, with strata variable, x axis starts at 10 unit, c-r curve from 1 to max, overlay optimal model with ensemble model based on all models examined)
* 
* %fitap(datain=cohort_ABC,perc_trim=1,lowest=1,translate=No,dataout=allmodels,case=Status,start=TStart,stop=TStop,fitvar=pm25,covvars=sex ses income bmi,strata=age, 
*	label_exposure=pm25, label_unit=ug/m3, x_min=10, output_path=D:\Working directory\Cohort, best_model=4);
*
* Example 11: (translated data, fit time-varying cox models, with strata variable, x axis starts at 10 unit, c-r curve from 1 to max, export HR matrix, overlay optimal model with ensemble model based on all models examined)
* 
* %fitap(datain=cohort_ABC,perc_trim=1,lowest=1,translate=Yes,dataout=allmodels,case=Status,start=TStart,stop=TStop,fitvar=pm25,covvars=sex ses income bmi,strata=age, 
*	label_exposure=pm25, label_unit=ug/m3, x_min=10, output_path=D:\Working directory\Cohort, best_model=4, export=Yes);
*
*************************************************************************************************/



/**************************
* main function to select optimal model and ensemble model
***************************/

option LINESIZE=MAX;  

%macro fitap(datain=,perc_trim=0,lowest=minimum,translate=yes,dataout=,case=,time=,start=,stop=,fitvar=,covvars=,strata=,label_exposure=,label_unit=,x_min=1,best_model=1,export=,output_path=);

/*Clean up any existing global macro variables*/
%symdel low_pct low_wc low_pct_1 low_wc_1 linear_model time start stop LL_min_all sum_LL_all LL_min_3 sum_LL_3 ensemble_model overlay;  

/*lowest must be either min or 1, but not anything else*/
%if %upcase(&lowest.) ne MINIMUM and &lowest. ne 1 %then %do;
	%put "WARNING: lowest must be either MINIMUM or 1"; 
	%abort;
%end;
 
/*Identify user's choice on model overlay*/

/*Show best fitting model only*/
%if &best_model.=1 %then %do;
	%let ensemble_model=;
	%let overlay=;
	/* %put &ensemble_model.; */
	/* %put %length(&ensemble_model.); */
%end;
/*Show only ensemble model based on all models examined*/
%else %if &best_model.=2 %then %do;
	%let ensemble_model=ensemble;
	%let overlay=ALL;
	/* %put &ensemble_model.; */
	/* %put %length(&ensemble_model.); */
%end;
/*Overlay best fitting model with ensemble model based on best 3 models examined*/
%else %if &best_model.=3 %then %do;
	%let ensemble_model=;
	%let overlay=best3;
	/* %put &overlay.; */
	/* %put %length(&overlay.); */
%end;
/*Overlay best fitting model with ensemble model based on ALL models examined*/
%else %if &best_model.>=4 %then %do;
	%let ensemble_model=;
	%let overlay=ALL;
	/* %put &overlay.; */
	/* %put %length(&overlay.); */
%end;

/*Consider non-linear models first*/
%let linear_model = 0;

/*convert any conc < 1 to 1*/
/*
data &datain.; set &datain.;
	if &fitvar.<1 then &fitvar.=1;
run;
*/

/*Trim data*/
%if 0<=&perc_trim.<= 10 %then %do;
	%let trim_l=%SYSEVALF(&perc_trim.);			
	%let trim_r=%SYSEVALF(100-&perc_trim.);
%end;
%else %do; %put "WARNING: perc_trim must be an integer between 0 and 10"; %abort;%end;

%if 0=&perc_trim. %then %do;
	data aftertrim; set &datain.; run;
	proc sql noprint;
	select min(&fitvar.) into: pctll from aftertrim;
	select max(&fitvar.) into: pctlr from aftertrim;
	quit;
%end;
%else %do; 
	proc univariate data=&datain. noprint;
	   var &fitvar.;
	   output out=percentiles pctlpts=&trim_l. &trim_r. pctlpre=ppp;
	run;
	data _null_;set percentiles;call symput('pctll',ppp%left(&trim_l.));call symput('pctlr',ppp%left(&trim_r.));run;
	data aftertrim; set &datain.; where &pctll.<=&fitvar.<=&pctlr.; run;
%end;

/*Descriptive statistics of exp variable in the original and trimmed datasets*/
proc means data=&datain. n nmiss min q1 mean median q3 max;
	var &fitvar.;
	output out=stats_original n=total_obs nmiss=miss_obs min=min_ap q1=q1_ap mean=mean_ap median=median_ap q3=q3_ap max=max_ap; 
run;
proc means data=aftertrim n nmiss min q1 mean median q3 max;
	var &fitvar.;
	output out=stats_aftertrim n=total_obs nmiss=miss_obs min=min_ap q1=q1_ap mean=mean_ap median=median_ap q3=q3_ap max=max_ap; 
run;
data stats_original (drop=_type_ _freq_); set stats_original;
	data_description="original dataset";
	exp_var="&fitvar.";
run;
data stats_aftertrim (drop=_type_ _freq_); set stats_aftertrim;
	data_description="trimmed dataset";
	exp_var="&fitvar.";
run;
data overall_stat; set stats_original stats_aftertrim; run;

/*Translate data*/
%put "check translate";
%put &translate.;
%put &pctll.;
%put &pctlr.;

%let pctll_bk=&pctll.; 	/* save a backup macro variable for &pctll. for translate=N */

%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
	/*%let Tran=1;*/ /* translate z to have the min of 1 */
	%let Tran=0;	/* translate z to have the min of 1 */
%end;
%else %do;
	%let Tran=0;	/* no translate */
	%let pctll=0;	/* no translate */
%end;
%put &Tran.;
%put &pctll.;

/*retain an original copy to be reused at each time when modelpct() is called*/
data aftertrim_backup; set aftertrim; run;

/*set tau=0.1*/
%let set_tau=0.1;
%let model_tau=1;	

/* %LET tau=%SYSEVALF(&set_tau.*(&p100.-&p0.)); */
%LET count=0;

/*run 8 models to decide the model type and pct with smallest coefficient*/
%do mdltp=1 %TO 2;
%do pctt=0 %TO 75 %BY 25;
	%modelpct(modeltype=&mdltp.,pct=&pctt.);
%end;%end;

data fourmodel_1;set mtll_:; run;
proc sort data=fourmodel_1;by WithCovariates;run;
data fourmodel_1;set fourmodel_1; if _N_=1;run;
data _null_;set fourmodel_1; 
call symput('modeltp_a',modeltype); call symput('low_pct_a',pct); call symput('low_wc_a',WithCovariates);
run;

/*set tau=0.2*/
%let set_tau=0.2;
%let model_tau=2;

/* %LET tau=%SYSEVALF(&set_tau.*(&p100.-&p0.)); */
%LET count=0;

/*run 8 models to decide the model type and pct with smallest coefficient*/
%do mdltp=1 %TO 2;
%do pctt=0 %TO 75 %BY 25;
	%modelpct(modeltype=&mdltp.,pct=&pctt.);
%end;%end;

data fourmodel_2;set mtll_:; run;
proc sort data=fourmodel_2;by WithCovariates;run;
data fourmodel_2;set fourmodel_2; if _N_=1;run;
data _null_;set fourmodel_2; 
call symput('modeltp_b',modeltype); call symput('low_pct_b',pct); call symput('low_wc_b',WithCovariates);
run;

/*compare and find the optimal tau*/
%if &low_wc_a.<=low_wc_b. %then %do;
	%let modeltp=&modeltp_a.;
	%let low_pct=&low_pct_a.;
	%let low_wc=&low_wc_a.;
	%let set_tau=0.1;
	/* %LET tau=%SYSEVALF(&set_tau.*(&p100.-&p0.)); */
	%let model_tau=9;
	%put &set_tau.;
	%put &model_tau.;

	/*retain rejected tau and related models*/
	data dataout_reject; set dataout_2; call symput('set_tau_reject',0.2); run;

%end;
%else %do;
	%let modeltp=&modeltp_b.;
	%let low_pct=&low_pct_b.;
	%let low_wc=&low_wc_b.;
	%let set_tau=0.2;
	/* %LET tau=%SYSEVALF(&set_tau.*(&p100.-&p0.)); */
	%let model_tau=9;
	%put &set_tau.;
	%put &model_tau.;

	/*retain rejected tau and related models*/
	data dataout_reject; set dataout_1; call symput('set_tau_reject',0.1); run;

%end;

/*compare -5/+5 percentile around low_pct from above, which is either 0 or 25 or 50 or 75 */
%let low_wc_1=&low_wc.; 
%let low_pct_1=&low_pct.;

%do pct_=(&low_pct.+5) %to (&low_pct.-5) %by -10; 
	%modelpct(modeltype=&modeltp., pct=&pct_.); 
	data _null_; set fits; call symput('new_wc',WithCovariates); call symput('new_pct',pct);	run;
	%if &new_wc.<&low_wc. %then %do; %let low_wc_1=&new_wc.; %let low_pct_1=&new_pct.; %end;
%end;

/* STOP if reaching mu=-15th, 100th, or LL is no longer smaller */
%do %while ( &low_pct_1. >= -10 and &low_pct_1. <= 95 and &low_pct_1. NE &low_pct.);
	%let low_pct_temp = %SYSEVALF(&low_pct_1. + (&low_pct_1. - &low_pct.));	
	%let low_pct=&low_pct_1.;
	%modelpct(modeltype=&modeltp., pct=&low_pct_temp.); 
	data _null_; set fits; call symput('new_wc',WithCovariates); call symput('new_pct',pct);	run;
	%if &new_wc.<&low_wc_1. %then %do; %let low_wc_1=&new_wc.; %let low_pct_1=&new_pct.; %end;
%end;

proc datasets;delete mtll: percentiles fourmodel: fits ;run;quit;

/*drop last run if mu was -20 or 105*/
data &dataout.; set &dataout.; where (pct NE -20); run; 
data &dataout.; set &dataout.; where (pct NE 105); run; 

/*add tau and append the rejected model outputs*/
proc sort data=&dataout.; by iteration; run;
data &dataout.; set &dataout.; tau=&set_tau.; run;

data dataout1_8; set &dataout.; where iteration<=8; run;

data dataout9_n; set &dataout.; where iteration>8; run;
data dataout9_n; set dataout9_n; iteration=iteration+8; run;

proc sort data=dataout_reject;by iteration;run;
data dataout_reject; set dataout_reject; tau=&set_tau_reject.; iteration=iteration+8; run;

data &dataout.; set dataout1_8 dataout_reject; run;
data &dataout.; set &dataout. dataout9_n; run;

/* print out chosen percentage and corresponding coefficient */
data &dataout.;set &dataout.; 
	rename pct=mu;
    rename pctl=z_at_mu;
	drop WithoutCovariates;
run;
proc sort data=&dataout.;by iteration;run;


/* Calculate ensemble weights using 3 models around the best fit, ie., based on best mu with +/- 5th%
   The three model include last 2 models from the search + a 3rd model corresponding mu+5 of last model */

%symdel best_LL final_mu;  

data &dataout. (drop=WithCovariates); set &dataout.; format LL d18.5; LL=WithCovariates; run;
data &dataout.; set &dataout.; rename LL=WithCovariates; run;

/* find optimal model corresonding to minimum LL */
proc sql noprint;
select min(WithCovariates) into: best_LL from &dataout.;
quit;

data &dataout.; set &dataout.; id=1; run;
data qaqc1; set &dataout.; run;
proc sort data=qaqc1; by WithCovariates; run;
data qaqc2; set qaqc1; id=_n_; run;
data qaqc2; set qaqc2; rename WithCovariates=best_LL2; where id=1; run;

data &dataout.; merge &dataout. (in=fro) qaqc2 (keep=id best_LL2); by id; if fro; run;

data &dataout.; set &dataout.; 
	if (WithCovariates = best_LL2) then do;
		best3=1; call symput('final_mu',mu); 
	end;
run;
data &dataout.; set &dataout.; drop id best_LL2; run;

data &dataout.; set &dataout.; final_mu=&final_mu.; final_form=&modeltp.; run;

/* find 2 other alternative models */
proc sql noprint;
select sum(best3) into: count_best_LL from &dataout.;		/* num of models with same smallest LL */
quit;

proc sql noprint;
select max(iteration) into: last_iteration from &dataout.;	/* iteration id for last run */
quit;

%if &count_best_LL.=2 %then %do;							/* a special case where last 2 runs had same smallest LL */

	data &dataout.; set &dataout.; if mu=final_mu and tau=&set_tau. then call symput('best_final_LL_iteration',iteration); run;

	%if &last_iteration. > &best_final_LL_iteration. %then %do;		/* followed by additional run with a larger LL*/
		data &dataout.; set &dataout.; 
			if iteration eq &last_iteration. then best3=3; 
		run;	
	%end;
	%if &last_iteration. eq &best_final_LL_iteration. %then %do;	/* followed by no more additional run*/
		data &dataout.; set &dataout.; 
			if &new_pct. > &low_pct_1. and mu=(final_mu-10) and final_form=modeltype and tau=&set_tau. then best3=3; 		/* ascending*/
			else if &new_pct. < &low_pct_1. and mu=(final_mu+10) and final_form=modeltype and tau=&set_tau. then best3=3;	/* descending*/
		run;
	%end;
	%end;
%if &count_best_LL.=1 %then
	%do;
	data &dataout.; set &dataout.; 
		if final_mu=-15 then do;				/* a special case where last run reached mu=-15 */
			if mu=-10 and tau=&set_tau. then best3=2; 
			if mu=-5 and tau=&set_tau. then best3=3;
			end;
		else if final_mu=100 then do;			/* a special case where last run reached mu=100 */
			if mu=95 and tau=&set_tau. then best3=2; 
			if mu=90 and tau=&set_tau. then best3=3;
			end;
		else do;								/* all other cases */
			if mu=(final_mu-5) and final_form=modeltype and tau=&set_tau. then best3=2; 
			if mu=(final_mu+5) and final_form=modeltype and tau=&set_tau. then best3=3;
			end;
	run;
	%end;

/* derive LL from -2LL */
data &dataout.; set &dataout.; 
	LL = WithCovariates/(-2);
run;

/* calculate ensemble weights for all models */
proc sql noprint;
select max(LL) into: LL_min_all from &dataout.;
quit;

data &dataout.; set &dataout.; 
	LL_diff = exp(LL-&LL_min_all.);
run;
proc sql noprint;
select sum(LL_diff) into: sum_LL_all 
	from &dataout.;
quit;

data &dataout.; set &dataout.; wt=LL_diff/&sum_LL_all.;run;
data &dataout.; set &dataout.; drop LL_diff; run;

/* calculate ensemble weights for 3 final models */
data final3models; set &dataout.; where best3>=1; run;

proc sql noprint;
select max(LL) into: LL_min_3 from final3models;
quit;
data final3models; set final3models; LL_diff = exp(LL-&LL_min_3.); run;
proc sql noprint;
select sum(LL_diff) into: sum_LL_3 from final3models;
quit;

data &dataout.; set &dataout.; 
	if best3>=1 then LL_diff = exp(LL-&LL_min_3.);
	else LL_diff=.;
run;
data &dataout.; set &dataout.; wt_final3=LL_diff/&sum_LL_3.;run;
data &dataout.; set &dataout.; if best3=. then wt_final3=.; drop LL_diff; run;


/* Bootstrap to sample beta based on joint model */

%if %length(&overlay.)>0 %then %do;

	%simulate_z2(indata=&dataout., model="joint", over_lay=&overlay.);
	proc sort data=newap3; by id; run;

	%simulate_beta(indata=&dataout., model="joint", over_lay=&overlay.);
	proc sort data=simdata; by id; run;

	data newap3; merge newap3(in=fro) simdata; by id; if fro; run;

	/* count max num of sim ap data points */
	proc contents data=newap3 out=output_z; run;
	proc sort data=output_z; by varnum; run;
	data _null_; 
	    set output_z;
	    call symputx("maximum",varnum-2);
	run;

	/* beta*transformed(z) for each simulated z data points */
	data newap3 (drop=beta id); set newap3;
		array CC{&maximum.} C1-C&maximum.;
		do i=1 to &maximum.;
			CC{i}=exp(beta*CC{i});
		end;
	run;

	/* derive median and 2.5th% and 97th% */
	%do i=1 %to &maximum.;
		data newap_&i.; set newap3;
			keep C&i.;
		run;
		proc univariate data=newap_&i. noprint;
	   		var C&i.;
	   		output out=beta_distr_&i. pctlpts= 2.5 50 97.5 pctlpre=P;
		run;
		data beta_distr_&i.;set beta_distr_&i.; id=&i.; run;
		%if &i.=1 %then %do;
			data beta_distr; set beta_distr_&i.; run;
		%end;
		%else %do;
			data beta_distr; set beta_distr beta_distr_&i.; run;		/* contain 4 variables: id, p2.5, p50, p97.5 */
		%end;
	%end;

	/* translate=N, rescale &pctll. to its original pctll */
	%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
		%let pctll_sim = &pctll.;		/* no change */
		%let pctlr_sim = &pctlr.;		/* no change */
	%end;
	%else %do;
		%let pctll_sim = &pctll_bk.;	/* scale pctlr_sim to max-min+1 */ 
		%let pctlr_sim = &pctlr.;		/* no change */
	%end;

	%if &lowest.= 1 %then %do;
		%let pctll_sim = 1;				/* force to extend the min to have 1*/
		%let pctlr_sim = &pctlr.;		/* no change */
	%end;

	data newap4(keep = z_sim); 
	   do i = &pctll_sim. to &pctlr_sim. by 0.1;               				
	      z_sim = i;
	      output;
	   end;
	run;
	data newap4; set newap4; id=_n_; run;
	proc sort data=newap4; by id; run;

	proc sort data=beta_distr; by id; run;
	data beta_distr; merge beta_distr(in=fro) newap4; if fro; by id; run;

	%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
		%if &modeltp.=1 %then %do;
			data newap_joint; set beta_distr; 
				rename p50=rr_mean;
				rename p97_5=rr_ucl;
				rename p2_5=rr_lcl;

				ap = z_sim; 				
			run;
		%end;
		%else %if &modeltp.=2 %then %do;
			data newap_joint; set beta_distr; 
				rename p50=rr_mean;
				rename p97_5=rr_ucl;
				rename p2_5=rr_lcl;

				ap = z_sim;  	
			run;
		%end;
	%end;
	%else %do;
		data newap_joint; set beta_distr; 
			rename p50=rr_mean;
			rename p97_5=rr_ucl;
			rename p2_5=rr_lcl;

			ap = z_sim;
		run;
	%end;

	/* export 1000*101 matrix */

	%if %upcase(&export.)= YES or %upcase(&export.)= Y %then %do;

		/* first, find ranking from 0th% to 100th% by 1, and reduce RR matrix to 1000:101 for export */

		data sampleN (drop=i); 
			do i=1 to &maximum.;
				CC=i;
				output;
			end;
		run; 
		proc univariate data=sampleN ;
	   		var CC;
	   		output out=sampleID pctlpts= 0 to 100 by 1 pctlpre=PP;
		run;
		data sampleID2 (drop=i PP0-PP100); set sampleID;
			array P{101} PP0-PP100;	
			do i=1 to 101;	
				CC_P=P{i};
				output;
			end;
		run;
		data sampleID2; set sampleID2;
			CC_P=int(CC_P);
		run;
		data sampleID2; set sampleID2;
			CC=compress('C'||CC_P);
			id=CC_P;
		run;
		proc sql noprint;                              
	 		select CC 
			into :CC1 - :CC101
	 	from sampleID2;
		quit;

		/* %put _user_; */
		%do i=1 %to 101;
			%put &&CC&i.;
			data export_newap_&i.; set newap3;
				keep &&CC&i.;
			run;
			data export_newap_&i.; set export_newap_&i.;
				id=_n_;
			run;
			proc sort data=export_newap_&i.; by id; run;
			%if &i.=1 %then %do;
				data export_newap; set export_newap_&i.; run;
			%end;
			%else %do;
				data export_newap; merge export_newap (in=fro) export_newap_&i.; if fro; by id; run;		
			%end;
		%end;
		data export_newap; set export_newap; drop id; run;
		proc export data=export_newap
	    outfile="&output_path.\export_hr.csv"
	    dbms=csv
	    replace;
		run;

		/* second, create another dataset to contain 101 data points from the range of AP exp */

		proc sort data=sampleID2; by id; run;
		data ap_exp; merge sampleID2 (in=fro) newap4; by id; if fro; run; 
		data ap_exp (drop=id CC_P); set ap_exp;
			CC=compress('CC'||id);
			rename z_sim=AP;
		run;
		proc export data=ap_exp
	    outfile="&output_path.\ap_exp.csv"
	    dbms=csv
	    replace;
		run;

	%end;

%end;


/* Bootstrap to sample beta based on optimal model */

%simulate_z2(indata=&dataout., model="optimal");
proc sort data=newap3; by id; run;

%simulate_beta(indata=&dataout., model="optimal");
proc sort data=simdata; by id; run;

data newap3; merge newap3(in=fro) simdata; by id; if fro; run;

/* count max num of sim ap data points */
proc contents data=newap3 out=output_z; run;
proc sort data=output_z; by varnum; run;
data _null_; 
    set output_z;
    call symputx("maximum",varnum-2);
run;
/* beta*transformed(z) for each simulated z data points */
data newap3 (drop=beta id); set newap3;
	array CC{&maximum.} C1-C&maximum.;
	do i=1 to &maximum.;
		CC{i}=exp(beta*CC{i});
	end;
run;
/* derive median and 2.5th% and 97th% */
%do i=1 %to &maximum.;
	data newap_&i.; set newap3;
		keep C&i.;
	run;
	proc univariate data=newap_&i. noprint;
   		var C&i.;
   		output out=beta_distr_&i. pctlpts= 2.5 50 97.5 pctlpre=P;
	run;
	data beta_distr_&i.;set beta_distr_&i.; id=&i.; run;
	%if &i=1 %then %do;
		data beta_distr; set beta_distr_&i.; run;
	%end;
	%else %do;
		data beta_distr; set beta_distr beta_distr_&i.; run;		/* contain 4 variables: id, p2.5, p50, p97.5 */
	%end;
%end;

/* translate=N, rescale &pctll. to its original pctll */
%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
	%let pctll_sim = &pctll.;		/* no change */
	%let pctlr_sim = &pctlr.;		/* no change */
%end;
%else %do;
	%let pctll_sim = &pctll_bk.;	/* scale pctlr_sim to max-min+1 */ 
	%let pctlr_sim = &pctlr.;		/* no change */
%end;

%if &lowest.= 1 %then %do;
	%let pctll_sim = 1;				/* force to extend the min to have 1*/
	%let pctlr_sim = &pctlr.;		/* no change */
%end;

data newap4(keep = z_sim); 
   do i = &pctll_sim. to &pctlr_sim. by 0.1;               				
      z_sim = i;
      output;
   end;
run;
data newap4; set newap4; id=_n_; run;
proc sort data=newap4; by id; run;

proc sort data=beta_distr; by id; run;
data beta_distr; merge beta_distr(in=fro) newap4; if fro; by id; run;

%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
	%if &modeltp.=1 %then %do;
		data newap_optimal; set beta_distr; 
			rename p50=rr_mean;
			rename p97_5=rr_ucl;
			rename p2_5=rr_lcl;

			ap = z_sim; 				
		run;
	%end;
	%else %if &modeltp.=2 %then %do;
		data newap_optimal; set beta_distr; 
			rename p50=rr_mean;
			rename p97_5=rr_ucl;
			rename p2_5=rr_lcl;

			ap = z_sim;  	
		run;
	%end;
%end;
%else %do;
	data newap_optimal; set beta_distr; 
		rename p50=rr_mean;
		rename p97_5=rr_ucl;
		rename p2_5=rr_lcl;

		ap = z_sim;
	run;
%end;


/* plot joint and optimal models */

%plot_cr(indata_joint=newap_joint, indata_optimal=newap_optimal, expo=&label_exposure., unit=&label_unit., ensemblemodel=&ensemble_model., over_lay=&overlay.); 

/* add a pure linear model with z */

%let linear_model = 1;

%put &linear_model.;
%put &modeltp.;
%put &low_pct_1.;

%modelpct(modeltype=9, pct=&low_pct_1.); 

/* add a pure log(z) model */

%let linear_model = 2;

proc datasets;delete fits ;run;quit;
%modelpct(modeltype=9, pct=&low_pct_1.); 

/* print description of original/trimmed datasets and model fitting, and clean up working directory */

ods pdf file="&output_path.\model_fitting_summary.pdf";

title;
proc print data=overall_stat; 
title "Summary statistics of air pollution exposure variable in the original and trimmed datasets";
run;
title;

title;
data &dataout.;set &dataout. (drop=LL); 
	length model_function $15.;

	if modeltype=1 then model_function='z*logit';
	else if modeltype=2 then model_function='log(z)*logit';
	else if modeltype=3 then model_function='pure linear';
	else if modeltype=4 then model_function='pure log';
	rename WithCovariates=LL;
	drop WithoutCovariates pct pctl best3 final_form final_mu modeltype;
run;
/* proc sort data=&dataout.;by iteration;run; */
data &dataout.; set &dataout.; drop iteration Criterion; run;		/* simplif summary table so as to fit in one page */
data &dataout.; set &dataout.; rename mu=location; run;
data &dataout.; set &dataout.; rename z_at_mu=mu; rename wt_final3=final_wt; run;
/* data &dataout.; set &dataout.; tau=&set_tau.; run;	*/
data prt; set &dataout.; run;
proc print data=prt;
title "Summary description of model fitting";
run; 
title;

ods pdf close;

/* export output table to csv format */
data &dataout.; set &dataout.; rename LL=Minus2LL; run;

proc export data=&dataout.
    outfile="&output_path.\model_fitting_summary.csv"
    dbms=csv
    replace;
run;

/* clean up non-essential datasets */
proc datasets library=work;
	delete sample: sample_all beta_distr final3models 
		   newap newap_joint newap_optimal param prt overall_stat simdata temp_sim temp_sim2 qaqc: percentiles output_z
 		   Z_sim_distr new_datain mtll: fits combined Beta_distr: dataout: newap: sampleID: sampleN export_newap_:;
run;quit;

%Mend fitap;



/**************************
* fit cox PH model with z
**************************/

%macro modelpct(modeltype=,pct=);

/*translate ap data to have the min of 1 or 0, depending on log or linear model*/
/*note that this only applies to translate=Y*/
/*note that if translate=N, then &pctll.=0 and &Tran.=0*/

/*for log model: translate dataset to have the min of 0, if translate=Y*/
%if &modeltype.=2 %then %do;
	data aftertrim; set aftertrim_backup; z=&fitvar.-&pctll.+&Tran.; run;
%end;
/*for linear model: translate dataset to have the min of 0*/
%else %if &modeltype.=1 %then %do;
	data aftertrim; set aftertrim_backup; z=&fitvar.-&pctll.; run;
%end;
/*for pure linear and pure log model: use original dataset*/
%else %if &modeltype.=9 %then %do;
	data aftertrim; set aftertrim_backup; z=&fitvar.; run;
%end;

proc sql noprint;
select max(z) into: p100 from aftertrim;
select min(z) into: p0 from aftertrim;
quit;

proc univariate data=aftertrim noprint;
   var z;
   output out=percentiles pctlpts= 5 to 95 by 5 pctlpre=P;
run;
data percentiles; set percentiles; p_5=&p0.-(p5-&p0.); run;			/* add a variable p_5 denoting 5 pctl below p0*/
data percentiles; set percentiles; p_10=&p0.-(p10-&p0.); run;		/* add a variable p_10 denoting 10 pctl below p0*/
data percentiles; set percentiles; p_15=&p0.-(p15-&p0.); run;		/* add a variable p_15 denoting 15 pctl below p0*/

proc transpose data=percentiles out=percentiles;run;

data _null_;set percentiles;
	call symput(_name_,col1);
run;

/* calc tau based on set_tau */
%LET tau=%SYSEVALF(&set_tau.*(&p100.-&p0.));

/*calc log transformed ap data*/
%if &pct.= -5 %then %do;
	%let mu=&p_5.; 		 					/* if 5 percentile less than p0, then change to _5 */

	data aftertrim; set aftertrim; 
	%if &modeltype.=1 %then %do;
	   APvar=(z)*(1/(1+exp(-(z-&mu.)/&tau.)));
	%end;
	%if &modeltype.=2 %then %do;
		%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
	      APvar= log(1+z)*(1/(1+exp(-(z-&mu.)/&tau.)));
		%end;
		%else %do;
	      APvar= log(1+z)*(1/(1+exp(-(z-&mu.)/&tau.)));
		%end;
	%end;run;
%end;
%else %if &pct. = -10 %then %do;
	%let mu=&p_10.; 		 				/* if 10 percentile less than p0, then change to _10 */

	data aftertrim; set aftertrim; 
	%if &modeltype.=1 %then %do;
	   APvar=(z)*(1/(1+exp(-(z-&mu.)/&tau.)));
	%end;
	%if &modeltype.=2 %then %do;
		%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
	      APvar= log(1+z)*(1/(1+exp(-(z-&mu.)/&tau.)));
		%end;
		%else %do;
	      APvar= log(1+z)*(1/(1+exp(-(z-&mu.)/&tau.)));
		%end;
	%end;run;
%end;
%else %if &pct. < -10 %then %do;
	%let mu=&p_15.; 		 				/* if 15 or more percentile less than p0, then change to _15 */

	data aftertrim; set aftertrim; 
	%if &modeltype.=1 %then %do;
	   APvar=(z)*(1/(1+exp(-(z-&mu.)/&tau.)));
	%end;
	%if &modeltype.=2 %then %do;
		%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
	      APvar= log(1+z)*(1/(1+exp(-(z-&mu.)/&tau.)));
		%end;
		%else %do;
	      APvar= log(1+z)*(1/(1+exp(-(z-&mu.)/&tau.)));
		%end;
	%end;run;
%end;
%else %do;
	%let mu=&&p&pct.;

	data aftertrim; set aftertrim; 
	%if &modeltype.=1 %then %do;
	   APvar=(z)*(1/(1+exp(-(z-&mu.)/&tau.)));
	%end;
	%if &modeltype.=2 %then %do;
		%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
	      APvar= log(1+z)*(1/(1+exp(-(z-&mu.)/&tau.)));
		%end;
		%else %do;
	      APvar= log(1+z)*(1/(1+exp(-(z-&mu.)/&tau.)));
		%end;
	%end;run;
%end;

/* non-linear model*/

%if &linear_model. eq 0 %then %do;
	%if %length(&start.)=0 %then %do;
		%if %length(&time.)=0 %then %do; %put "WARNING: TIME variable is missing"; %abort;%end;

		%if %length(&strata.)=0 %then %do;
			proc phreg data=aftertrim;
			model &time.*&case.(0)=APvar &covvars.;
			ods output FitStatistics=fits(where=(Criterion='-2 LOG L')) 
					ParameterEstimates=param (where=(parameter='APvar'));
			run;
		%end;
		%else %do;
			proc phreg data=aftertrim;
			model &time.*&case.(0)=APvar &covvars.;
			strata &strata.;
			ods output FitStatistics=fits(where=(Criterion='-2 LOG L')) 
					ParameterEstimates=param (where=(parameter='APvar'));
			run;
		%end;
	%end;
	%if %length(&start.)>0 %then %do;
		%if %length(&stop.)=0 %then %do; %put "WARNING: STOP variable is missing"; %abort;%end;

		%if %length(&strata.)=0 %then %do;
			proc phreg data=aftertrim;
			model (&start.,&stop.)*&case.(0)=APvar &covvars.;
			ods output FitStatistics=fits(where=(Criterion='-2 LOG L')) 
					ParameterEstimates=param (where=(parameter='APvar'));
			run;
		%end;
		%else %do;
			proc phreg data=aftertrim;
			model (&start.,&stop.)*&case.(0)=APvar &covvars.;
			strata &strata.;
			ods output FitStatistics=fits(where=(Criterion='-2 LOG L')) 
					ParameterEstimates=param (where=(parameter='APvar'));
			run;
		%end;

	%end;
%end;

/* pure linear z model */

%if &linear_model. eq 1 %then %do;
	%if %length(&start.)=0 %then %do;
		%if %length(&time.)=0 %then %do; %put "WARNING: TIME variable is missing"; %abort;%end;

		%if %length(&strata.)=0 %then %do;
			proc phreg data=aftertrim;
			model &time.*&case.(0)=z &covvars.;
			ods output FitStatistics=fits(where=(Criterion='-2 LOG L')) 
					ParameterEstimates=param (where=(parameter='z'));
			run;
		%end;
		%else %do;
			proc phreg data=aftertrim;
			model &time.*&case.(0)=z &covvars.;
			strata &strata.;
			ods output FitStatistics=fits(where=(Criterion='-2 LOG L')) 
					ParameterEstimates=param (where=(parameter='z'));
			run;
		%end;
	%end;
	%if %length(&start.)>0 %then %do;
		%if %length(&stop.)=0 %then %do; %put "WARNING: STOP variable is missing"; %abort;%end;

		%if %length(&strata.)=0 %then %do;
			proc phreg data=aftertrim;
			model (&start.,&stop.)*&case.(0)=z &covvars.;
			ods output FitStatistics=fits(where=(Criterion='-2 LOG L')) 
					ParameterEstimates=param (where=(parameter='z'));
			run;
		%end;
		%else %do;
			proc phreg data=aftertrim;
			model (&start.,&stop.)*&case.(0)=z &covvars.;
			strata &strata.;
			ods output FitStatistics=fits(where=(Criterion='-2 LOG L')) 
					ParameterEstimates=param (where=(parameter='z'));
			run;
		%end;
	%end;
%end;

/* pure log(z) model */

%if &linear_model. eq 2 %then %do;
	data new_datain; set aftertrim; ap_new=log(z); run;
	%if %length(&start.)=0 %then %do;
		%if %length(&time.)=0 %then %do; %put "WARNING: TIME variable is missing"; %abort;%end;

		%if %length(&strata.)=0 %then %do;
			proc phreg data=new_datain;
			model &time.*&case.(0)=ap_new &covvars.;
			ods output FitStatistics=fits(where=(Criterion='-2 LOG L')) 
					ParameterEstimates=param (where=(parameter='ap_new'));
			run;
		%end;
		%else %do;
			proc phreg data=new_datain;
			model &time.*&case.(0)=ap_new &covvars.;
			strata &strata.;
			ods output FitStatistics=fits(where=(Criterion='-2 LOG L')) 
					ParameterEstimates=param (where=(parameter='ap_new'));
			run;
		%end;
	%end;
	%if %length(&start.)>0 %then %do;
		%if %length(&stop.)=0 %then %do; %put "WARNING: STOP variable is missing"; %abort;%end;

		%if %length(&strata.)=0 %then %do;
			proc phreg data=new_datain;
			model (&start.,&stop.)*&case.(0)=ap_new &covvars.;
			ods output FitStatistics=fits(where=(Criterion='-2 LOG L')) 
					ParameterEstimates=param (where=(parameter='ap_new'));
			run;
		%end;
		%else %do;
			proc phreg data=new_datain;
			model (&start.,&stop.)*&case.(0)=ap_new &covvars.;
			strata &strata.;
			ods output FitStatistics=fits(where=(Criterion='-2 LOG L')) 
					ParameterEstimates=param (where=(parameter='ap_new'));
			run;
		%end;
	%end;
%end;

data param; set param; call symput('new_coef',Estimate); call symput('new_std',StdErr);run;

%if &linear_model. eq 0 %then %do;
data fits;set fits; length modeltype $4.; modeltype=&modeltype.; pct=&pct.; pctl=&mu.; coef=&new_coef.; 
	stderr=&new_std.; 
run;
/*data fits;set fits; if modeltype=1 then pctl=pctl-&Tran.; run;*/
%end;

%if &linear_model. eq 1 %then %do;
data fits;set fits; length modeltype $4.; modeltype=3; pct=.; pctl=.; coef=&new_coef.; stderr=&new_std.; run;
%end;

%if &linear_model. eq 2 %then %do;
data fits;set fits; length modeltype $4.; modeltype=4; pct=.; pctl=.; coef=&new_coef.; stderr=&new_std.; run;
%end;

/* first 8 model runs corresponding to tau=0.1 */
%if &model_tau.=1 %then %do;

	%LET count=%SYSEVALF(&count. + 1);
	data fits;set fits; iteration=&count.; run;

	/* the very first model examined */
	%if &count.=1 %then %do;
	data dataout_1;set fits;run; %put "lastly here here here"; %put &count.; %end;
	%else %do;
	data dataout_1;set dataout_1 fits; run; %end;

	%if &pct.<0 %then %do; %let abs_pct=%sysfunc(abs(&pct.)); data mtll_&modeltype._&abs_pct._;set fits; run;; %end;
	%else %do; data mtll_&modeltype._&pct.;set fits; run; %end;

%end;
/* second 8 model runs corresponding to tau=0.2 */
%else %if &model_tau.=2 %then %do;

	%LET count=%SYSEVALF(&count. + 1);
	data fits;set fits; iteration=&count.; run;

	/* the very first model examined */
	%if &count.=1 %then %do;
	data dataout_2;set fits;run; %put "lastly here here here"; %put &count.; %end;
	%else %do;
	data dataout_2;set dataout_2 fits; run; %end;

	%if &pct.<0 %then %do; %let abs_pct=%sysfunc(abs(&pct.)); data mtll_&modeltype._&abs_pct._;set fits; run;; %end;
	%else %do; data mtll_&modeltype._&pct.;set fits; run; %end;

%end;
/* continue the rest of model runs, based on selected tau, either 0.1 or 0.2 */
%else %if &model_tau.=9 %then %do;

	%LET count=%SYSEVALF(&count. + 1);
	data fits;set fits; iteration=&count.; run;

	/* based on the selected tau, initiate output dataset*/
	%if &count.=9 %then %do;
		%if &set_tau.=0.1 %then %do;
			data &dataout.;set dataout_1;run; 
		%end;
		%else %if &set_tau.=0.2 %then %do;
			data &dataout.;set dataout_2;run; 
		%end;
	%end;
    /* append the 9th model and onwards*/
	%put "Append the 9th model and onwards";
	data &dataout.;set &dataout. fits; run; 

	%if &pct.<0 %then %do; %let abs_pct=%sysfunc(abs(&pct.)); data mtll_&modeltype._&abs_pct._;set fits; run; %end;
	%else %do; data mtll_&modeltype._&pct.;set fits; run; %end;

%end;
run;
%Mend modelpct;



/**************************
* Simulate z over 1000 realizations
**************************/

%macro simulate_z1(modeltype=, nn=);

/* based on the range of trimmed and translated ap_data */
%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
	%if %upcase(&lowest.)= MINIMUM %then %do;
		%let max_ap = %SYSEVALF(&pctlr.-&pctll.+&Tran.); 
		data newap(keep = z_sim); 
		   do i = &Tran. to &max_ap. by 0.1;               				
		      z_sim = i;
		      output;
		   end;
		run;
	%end;
	/* based on the range of 1 and max value of true ap_data */
	%else %if &lowest.= 1 %then %do;
		data newap(keep = z_sim); 
		   do i = 1 to &pctlr. by 0.1;               				
		      z_sim = i;
		      output;
		   end;
		run;
	%end;
	%else %do; %put "WARNING: lowest must be either MINIMUM or 1"; %abort;%end;
%end;
%else %do;	
	%if %upcase(&lowest.)= MINIMUM %then %do;
		/* translate=N, scale &pctlr. by deducing it by original pctll AND &pctll to 1 - only use this with log(z) but NOT log(1+z) */
		/*%let pctlr_simz = %SYSEVALF(&pctlr.-&pctll_bk.+1); */		/* scale pctlr_simz to max-min+1 */ 
		/*%let pctll_simz = 1; */									/* scale pctll_simz to 1 */

		/* translate=N, scale &pctlr. by deducing it by original pctll AND pctlr - not used, otherwise an error on HRs>1 for ap<1 would occur! */
		/* %let pctlr_simz = &pctlr.; 	*/							/* scale pctlr_simz to max */ 
		/* %let pctll_simz = &pctll_bk.;	*/						/* scale pctll_simz to pctll_bk */

		/* translate=N, scale &pctlr. by deducing it by 0 AND pctlr */
		%let pctlr_simz = %SYSEVALF(&pctlr.-&pctll.+&Tran.);		/* scale pctlr_simz to max */ 
		%let pctll_simz = &Tran.;									/* scale pctll_simz to 0 */
	
		data newap(keep = z_sim); 
		   do i = &pctll_simz. to &pctlr_simz. by 0.1;               				
		      z_sim = i;
		      output;
		   end;
		run;
	%end;
	/* based on the range of 1 and max value of true ap_data */
	%else %if &lowest.= 1 %then %do;
		data newap(keep = z_sim); 
		   do i = 1 to &pctlr. by 0.1;               				
		      z_sim = i;
		      output;
		   end;
		run;
	%end;
	%else %do; %put "WARNING: lowest must be either MINIMUM or 1"; %abort;%end;
%end;

/* if linear model, scale the min to 0 */
proc sql noprint;
select min(z_sim) into: Tran2 from newap;
quit;
%if &modeltype.=1 %then %do; 
	data newap; set newap;   
		z_sim=z_sim-&Tran2.;	
	/*	z_sim=z_sim-&Tran.; */	
	run;
%end;

/* prepare nn x _n_ matrix to store z_sim */
proc sort data=newap; by z_sim; run;
data newap; set newap; seqno = _n_ ; run;

proc transpose data=newap out=newap_wide prefix=C;	/* col name is C1, C2, ..., Cxxx */
    id seqno;
    var z_sim;
run;
data newap2 (drop=_name_ i); set newap_wide; 			/* nn*num_z matrix */
	do i=1 to &nn.;
		output;
	end;		
run;

/* extract min and max of z_sim */
proc sql noprint;
select max(z_sim) into: z_p100 from newap;
select min(z_sim) into: z_p0 from newap;
quit;

proc univariate data=newap noprint;
   var z_sim;
   output out=percentiles pctlpts= 0 to 100 by 5 pctlpre=z_p;
run;
data percentiles; set percentiles; z_p_5=&z_p0.-(z_p5-&z_p0.); run;			/* add a variable z_p_5 denoting 5 pctl below z_p0*/
data percentiles; set percentiles; z_p_10=&z_p0.-(z_p10-&z_p0.); run;		/* add a variable z_p_10 denoting 10 pctl below z_p0*/
data percentiles; set percentiles; z_p_15=&z_p0.-(z_p15-&z_p0.); run;		/* add a variable z_p_15 denoting 15 pctl below z_p0*/

proc transpose data=percentiles out=percentiles;run;

%mend simulate_z1; 

%macro simulate_z2(indata=, model=, over_lay=);

/* prepare z for use in simulation, varying depending on perc, set_tau, and funcform */
%if &model.="joint" %then %do;

	/* all models examined */
	%if %length(&over_lay.)>0 %then %do;

		%if %upcase(&over_lay.)=ALL %then %do;
			data temp_sim; set &indata.; run;
			data temp_sim; set temp_sim; seqno = _n_; wt_final3=wt; run;
		%end;
		%else %if %upcase(&over_lay.) ne ALL %then %do;
			data temp_sim; set &indata.; where best3>=1; run;
			data temp_sim; set temp_sim; seqno = _n_; run;
		%end;

		proc sql noprint;
		select max(seqno) into: num_models from temp_sim;		/* total number of models examined */
		quit;

		%let counter=1;

		%do k=1 %to &num_models.;
			data _null_; set temp_sim; 
				call symput('wt_each',wt_final3); call symput('z_tau',tau);
				call symput('z_mu',mu);call symput('mtype',modeltype);
				where seqno=&k.; 
			run;

			/* num of models */
			%let N=%SYSEVALF(1000*&wt_each., floor);

			%if &N.> 0 %then %do;
				/* simulate ap data */
				%simulate_z1(modeltype=&mtype., nn=&N.);
				data _null_;set percentiles;
					call symput(_name_,col1);
					run;
				proc sql noprint;
					select max(seqno) into: num_z from newap;		/* total number of z_sim data points */
				quit;

				%let xxx=&num_z.;

				/* calc tau based on set_tau of each model run */
				%LET z_tau_val=%SYSEVALF(&z_tau.*(&z_p100.-&z_p0.));

				/*calc log transformed ap data*/
				%if &z_mu.= -5 %then %do;
					%let z_mu_val=&z_p_5.; 		 					/* if 5 percentile less than p0, then change to _5 */

					data newap2; set newap2; 
						array CC{&xxx.} C1-C&xxx.;
						%do i=1 %to &xxx.;

							%if &mtype.=1 %then %do;
						   	CC{&i.}=(CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
							%end;
							%if &mtype.=2 %then %do;
								%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
						   			CC{&i.}= log(1+CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
								%end;
								%else %do;
						   			CC{&i.}= log(1+CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
								%end;

							%end;

						%end;					
					run;
				%end;
				%else %if &z_mu. = -10 %then %do;
					%let z_mu_val=&z_p_10.; 		 				/* if 10 percentile less than p0, then change to _10 */

					data newap2; set newap2; 
						array CC{&xxx.} C1-C&xxx.;
						%do i=1 %to &xxx.;

							%if &mtype.=1 %then %do;
							   CC{&i.}=(CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
							%end;
							%if &mtype.=2 %then %do;
								%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
						   			CC{&i.}= log(1+CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
								%end;
								%else %do;
						   			CC{&i.}= log(1+CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
								%end;

							%end;

						%end;					
					run;
				%end;
				%else %if &z_mu. < -10 %then %do;
					%let z_mu_val=&z_p_15.; 		 				/* if 15 or more percentile less than p0, then change to _15 */

					data newap2; set newap2; 
						array CC{&xxx.} C1-C&xxx.;
						%do i=1 %to &xxx.;

							%if &mtype.=1 %then %do;
							   CC{&i.}=(CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
							%end;
							%if &mtype.=2 %then %do;
								%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
						   			CC{&i.}= log(1+CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
								%end;
								%else %do;
						   			CC{&i.}= log(1+CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
								%end;

							%end;
						
						%end;					
					run;
				%end;
				%else %do;
					%let xx=&z_mu.;
					%let z_mu_val=&&z_p&xx.;
					%put &z_mu_val.;

					data newap2; set newap2; 
						array CC{&xxx.} C1-C&xxx.;
						%do i=1 %to &xxx.;
				
							%if &mtype.=1 %then %do;
							   CC{&i.}=(CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
							%end;
							%if &mtype.=2 %then %do;
								%if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
						   			CC{&i.}= log(1+CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
								%end;
								%else %do;
						   			CC{&i.}= log(1+CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
								%end;

							%end;
						
						%end;					
					run;
				%end;

				%if &counter.=1 %then %do;
					data newap3; set newap2; run;
				%end;
				%else %do;
					data newap3; set newap3 newap2; run;
				%end;

				%let counter=%SYSEVALF(&counter.+1);

			%end; /* if &N.> 0 */

		%end; /* do k=1 */

	%end; /*if %length(&over_lay.)>0*/

%end; /*if &model.="joint"*/

/* optimal model */
%if &model.="optimal" %then %do;

	data temp_sim2; set &indata.; where best3=1; run;
	proc sort data=temp_sim2; by best3 iteration; run;			/* in case last 2 models had same LL */
	proc sort data=temp_sim2 NODUPKEY; by best3; run;			/* in case last 2 models had same LL */
	data _null_; set temp_sim2; 
		call symput('wt_each',wt_final3); call symput('z_tau',tau);
		call symput('z_mu',mu);call symput('mtype',modeltype);
	run;

	/* num of models */
	%let N=1000;

	%simulate_z1(modeltype=&mtype., nn=&N.);

	data _null_;set percentiles;
		call symput(_name_,col1);
	run;
	proc sql noprint;
		select max(seqno) into: num_z from newap;		/* total number of z_sim data points */
	quit;

	%let xxx=&num_z.;

	/* calc tau based on set_tau of each model run */
	%LET z_tau_val=%SYSEVALF(&z_tau.*(&z_p100.-&z_p0.));

	/* extract mu */
	%if &z_mu.= -5 %then %do;
		%let z_mu_val=&z_p_5.; 	
	%end;
	%else %if &z_mu.= -10 %then %do;
		%let z_mu_val=&z_p_10.;
	%end;
	%else %if &z_mu.= -15 %then %do;
		%let z_mu_val=&z_p_15.; 
	%end;
	%else %do;
		%let xx=&z_mu.;
		%let z_mu_val=&&z_p&xx.;
	%end;

	data newap3; set newap2; 
		array CC{&xxx.} C1-C&xxx.;
		%do i=1 %to &xxx.;
			%if &mtype.=1 %then %do;
			   CC{&i.}=(CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
			%end;
			%if &mtype.=2 %then %do;
			   %if %upcase(&translate.)= YES or %upcase(&translate.)= Y %then %do;
				  CC{&i.}= log(1+CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
			   %end;
			   %else %do;
			      CC{&i.}= log(1+CC{&i.})*(1/(1+exp(-(CC{&i.}-&z_mu_val.)/&z_tau_val.)));
			   %end;

			%end;
		%end;					
	run;
%end;

data newap3; set newap3; id=_n_; run;
data newap3; set newap3; where 1<=id<=1000; run;	

%mend simulate_z2;


/**************************
* Simulate beta over 1000 realizations
**************************/

%macro simulate_beta(indata=, model=, over_lay=);

%put %upcase(&over_lay.);

%if &model.="joint" %then %do;

	/* all models examined */
	%if %length(&over_lay.)>0 and %upcase(&over_lay.)=ALL %then %do;
		data temp_sim; set &indata.; run;
		data temp_sim; set temp_sim; seqno = _n_; run;

		proc sql noprint;
		select max(seqno) into: num_models from temp_sim;		/* total number of models examined */
		quit;

		%do k=1 %to &num_models.;
			data _null_; set temp_sim; 
				call symput('wt_each',wt); call symput('coef_each',coef);call symput('std_each',stderr);
				where seqno=&k.; 
			run;		

			%let N=%SYSEVALF(1000*&wt_each., floor);
			data sample_&k.(keep=beta);
				call streaminit(4321);
				do i=1 to &N.;
					beta=rand("Normal", &coef_each., &std_each.); 
					output;
				end;
			run;

			%if &k.=1 %then %do; 
				data sample_all; set sample_&k.; run; 
			%end;
			%else %do;
				data sample_all; set sample_all sample_&k.; run; 
			%end;
		%end;

		data simdata; set sample_all; run;
		data simdata; set simdata; id=_n_; run;

	%end;

	/* 3 best models examined */
	%if %length(&over_lay.)>0 and %upcase(&over_lay.) ne ALL %then %do;
		data temp_sim; set &indata.; where best3>=1; run;
		data temp_sim; set temp_sim; seqno = _n_; run;

		proc sql noprint;
		select max(seqno) into: num_models from temp_sim;		/* total number of models examined */
		quit;

		/* first model */
		data _null_; set temp_sim; 
			call symput('wt_each',wt_final3); call symput('coef_each',coef);call symput('std_each',stderr);
			where seqno=1; 
		run;		

		%let N=%SYSEVALF(1000*&wt_each., floor);
		data sample1(keep=beta);
			call streaminit(4321);
			do i=1 to &N.;
				beta=rand("Normal", &coef_each., &std_each.); 
				output;
			end;
		run;
		/* second model */
		data _null_; set temp_sim; 
			call symput('wt_each',wt_final3); call symput('coef_each',coef);call symput('std_each',stderr);
			where seqno=2; 
		run;		

		%let N=%SYSEVALF(1000*&wt_each., floor);
		data sample2(keep=beta);
			call streaminit(4321);
			do i=1 to &N.;
				beta=rand("Normal", &coef_each., &std_each.); 
				output;
			end;
		run;
		/* third model */
	   	data _null_; set temp_sim; 
			call symput('wt_each',wt_final3); call symput('coef_each',coef);call symput('std_each',stderr);
			where seqno=3; 
		run;		

		%let N=%SYSEVALF(1000*&wt_each., floor);
		data sample3(keep=beta);
			call streaminit(4321);
			do i=1 to &N.;
				beta=rand("Normal", &coef_each., &std_each.); 
				output;
			end;
		run;

		data sample_all; set sample1; run; 
		data sample_all; set sample_all sample2; run; 
		data sample_all; set sample_all sample3; run; 

		data simdata; set sample_all; run;
		data simdata; set simdata; id=_n_; run;

	%end;
%end;

/* optimal model */
%if &model.="optimal" %then %do;
	data temp_sim2; set &indata.; where best3=1; run;
	proc sort data=temp_sim2; by best3 iteration; run;			/* in case last 2 models had same LL */
	proc sort data=temp_sim2 NODUPKEY; by best3; run;			/* in case last 2 models had same LL */
	data _null_; set temp_sim2; 
		call symput('wt_each',wt_final3); call symput('coef_each',coef); call symput('std_each',stderr);
	run;
	%let N=1000;
	data sample(keep=beta);
		call streaminit(4321);
		do i=1 to &N.;
			beta=rand("Normal", &coef_each., &std_each.); 
			output;
		end;
	run;
	data sample_all; set sample; run;
	data simdata; set sample_all; id=_n_; run;
%end;

%mend simulate_beta;


/**************************
* plot c-r function
**************************/

%macro plot_cr(indata_joint=, indata_optimal=, expo=, unit=, ensemblemodel=, over_lay=); 

%if %length(&expo.)=0 %then %do; 
	%let label2='Air pollution concentration';
%end;
%else %do;
    %IF %length(&unit.)>0 %then %do; %let label2="&expo. concentration (&unit.)"; %end;
	%IF %length(&unit.)=0 %then %do; %let label2="&expo. concentration"; %end;
%end;
%put &label2.;
%put &unit.;
%put &x_min.;

/* %if %length(&over_lay.)>0 %then %do; */

%if %length(&over_lay.)>0 and %length(&ensemblemodel.)=0 %then %do;

	/* data &indata_joint.; set &indata_joint.; model="Ensemble"; run; */
	proc sort data=&indata_joint.; by ap; run;
	data &indata_joint.; set &indata_joint.; 
		rename rr_mean=joint_rr_mean;
		rename rr_lcl=joint_rr_lcl;
		rename rr_ucl=joint_rr_ucl;
		id=_n_;
	run;

	/* data &indata_optimal.; set &indata_optimal.; model="Optimal"; run; */
	proc sort data=&indata_optimal.; by ap; run;
	data &indata_optimal.; set &indata_optimal.; 
		rename rr_mean=optimal_rr_mean;
		rename rr_lcl=optimal_rr_lcl;
		rename rr_ucl=optimal_rr_ucl;
		id=_n_; 
	run;

	/* data combined; set &indata_joint. &indata_optimal.; run; */
	data combined; merge &indata_joint. (in=fro) &indata_optimal. (keep=id optimal_rr_mean optimal_rr_lcl optimal_rr_ucl); 
		by id; if fro; 
	run;

	ods listing close;
	ods html image_dpi=200 file='fitmodel.html' path="&output_path." style=listing;
	ods graphics / reset noborder width=600px height=400px imagename="Nonlinear_&sysdate.";

	%if &x_min. ge 0 %then %do;
		title1 "Nonlinear C-R Curve";
		proc sgplot data=combined noautolegend;  
		   band x=ap lower=joint_rr_lcl upper=joint_rr_ucl / 
				fillattrs=(color=LTGRAY);
		   series x=ap y=joint_rr_mean / lineattrs=(color=blue);
		   band x=ap lower=optimal_rr_lcl upper=optimal_rr_ucl / 
				NOFILL LINEATTRS=(color=red PATTERN=dash THICKNESS=1);
		   series x=ap y=optimal_rr_mean / lineattrs=(color=red);

		   xaxis label=&label2. min=&x_min.;	
		   yaxis label='Hazard Ratio (95% CI)';
		   refline 1 / LINEATTRS = (color=black);
		run;
	%end;
	%else %do;
		title1 "Nonlinear C-R Curve";
		proc sgplot data=combined noautolegend;  
		   band x=ap lower=joint_rr_lcl upper=joint_rr_ucl / 
				fillattrs=(color=LTGRAY);
		   series x=ap y=joint_rr_mean / lineattrs=(color=blue);
		   band x=ap lower=optimal_rr_lcl upper=optimal_rr_ucl / 
				NOFILL LINEATTRS=(color=red PATTERN=dash THICKNESS=1);
		   series x=ap y=optimal_rr_mean / lineattrs=(color=red);

		   xaxis label=&label2.;	
		   yaxis label='Hazard Ratio (95% CI)';
		   refline 1 / LINEATTRS = (color=black);
		run;
	%end;

	ods html close;
	ods listing;
%end;
%else %do;
	data &indata_optimal.; set &indata_optimal.; model="Optimal"; run;
	data combined; set &indata_optimal.; run;

	%IF %length(&ensemblemodel.)>0 %then %do; 
		data &indata_joint.; set &indata_joint.; model="Ensemble"; run;
		data combined; set &indata_joint.; run;
	%end;

	ods listing close;
	ods html image_dpi=200 file='fitmodel.html' path="&output_path." style=listing;
	ods graphics / reset noborder width=600px height=400px imagename="Nonlinear_&sysdate.";

	%if &x_min. ge 0 %then %do;
		title1 "Nonlinear C-R Curve";
		proc sgplot data=combined noautolegend;  
		   series x=ap y=rr_mean / group=model;
		   band x=ap lower=rr_lcl upper=rr_ucl / group=model transparency=0.75;
		   xaxis label=&label2. min=&x_min.;	
		   yaxis label='Hazard Ratio (95% CI)';
		   refline 1;
		   /* keylegend / location=inside position=topleft across=1; */
		run;
	%end;
	%else %do;
		title1 "Nonlinear C-R Curve";
		proc sgplot data=combined noautolegend;  
		   series x=ap y=rr_mean / group=model;
		   band x=ap lower=rr_lcl upper=rr_ucl / group=model transparency=0.75;
		   xaxis label=&label2.;	
		   yaxis label='Hazard Ratio (95% CI)';
		   refline 1;
		   /* keylegend / location=inside position=topleft across=1; */
		run;
	%end;

	ods html close;
	ods listing;
%end;

%mend plot_cr;












