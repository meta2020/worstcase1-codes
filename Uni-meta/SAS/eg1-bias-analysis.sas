/*--------------------------------------------------------------------*/
/* SIMULATION-BASED WORST-CASE BOUND OF META-ANALYSIS OF INTERVENTION */
/*                                                                    */
/* EXAMPLE 1                                                          */
/*--------------------------------------------------------------------*/

proc datasets library=work kill noprint;quit;
/*CLOSE ALL NOTES*/
ods _all_ close;
options nonotes;

/*IMPORT RAW DATASET*/
/*VARIABLES: LOGOR AND PRECISION (=1/SIGMA)*/
proc import  out=dat_input
	datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example1\eg1-data.csv"
	dbms = csv replace;
	getnames = yes;
run;

/*LOAD MACRO OF CALCULATING BIAS*/
%include "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example1\SAS\macro_bias1.sas";

/*SIMULATION-BASED WORST-CASE BOUND AND STORE RESULTS*/
%let kk = 100; %let rr = 10; %let pl = 0;
ods html file="result1-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example1\SAS\";

%let start_time = %sysfunc(datetime());
%bias1(start=0.1, end=1, by=0.1, repeatn=&rr, dtin = dat_input, k=&kk, seed=2000, plvl = &pl);
%put **** Duration (K=&kk) = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
proc export data=repeat
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example1\SAS\result1-K&kk.-R&rr..csv"
    dbms=csv replace;
run;
ods html close;


