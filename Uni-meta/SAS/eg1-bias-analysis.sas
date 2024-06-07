/*--------------------------------------------------------------------*/
/* SIMULATION-BASED WORST-CASE BOUND OF META-ANALYSIS OF INTERVENTION */
/*                                                                    */
/* EXAMPLE 1                                                          */
/*--------------------------------------------------------------------*/

/*CREATE LOG FILE*/
proc printto log = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Uni-meta\SAS\NEQ\log.txt"; run;

/*CLEAR LIBRARY*/
proc datasets library=work kill noprint;quit;
/*CLOSE ALL NOTES*/
ods _all_ close;
options nonotes;


/*IMPORT RAW DATASET*/
/*VARIABLES: LOGOR AND PRECISION (=1/SIGMA)*/
proc import  out=dat_input
	datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Uni-meta\eg1-data.csv"
	dbms = csv replace;
	getnames = yes;
run;



/*LOAD MACRO OF CALCULATING BIAS*/
%include "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Uni-meta\SAS\macro_bias1.sas";

/*SIMULATION-BASED WORST-CASE BOUND AND STORE RESULTS*/
/*SET HOW MANY TIMES FOR REPEAT*/
%let rr = 1; 
/*SET WHETHER TO PRINT THE OPTMODEL RESULT*/
%let pl = 1;
/*SET SEED*/
%let seedn = 2000;
/*SET HOW MANY RANDOM VARIABLES TO USE*/


%let kk = 1000; 
ods html file="result1-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Uni-meta\SAS\NEQ\";
%let start_time = %sysfunc(datetime());
%bias1(start=0.1, end=1, by=0.1, repeatn=&rr, dtin = dat_input, k=&kk, seed=&seedn, plvl = &pl);
%put **** Duration (K=&kk) = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
proc export data=repeat
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Uni-meta\SAS\NEQ\result1-K&kk.-R&rr..csv"
    dbms=csv replace;
run;
%let kk = 2000; 
ods html file="result1-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Uni-meta\SAS\NEQ\";
%let start_time = %sysfunc(datetime());
%bias1(start=0.1, end=1, by=0.1, repeatn=&rr, dtin = dat_input, k=&kk, seed=&seedn, plvl = &pl);
%put **** Duration (K=&kk) = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
proc export data=repeat
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Uni-meta\SAS\NEQ\result1-K&kk.-R&rr..csv"
    dbms=csv replace;
run;
%let kk = 5000; 
ods html file="result1-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Uni-meta\SAS\NEQ\";
%let start_time = %sysfunc(datetime());
%bias1(start=0.1, end=1, by=0.1, repeatn=&rr, dtin = dat_input, k=&kk, seed=&seedn, plvl = &pl);
%put **** Duration (K=&kk) = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
proc export data=repeat
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Uni-meta\SAS\NEQ\result1-K&kk.-R&rr..csv"
    dbms=csv replace;
run;
%let kk = 20000; 
ods html file="result1-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Uni-meta\SAS\NEQ\";
%let start_time = %sysfunc(datetime());
%bias1(start=0.1, end=1, by=0.1, repeatn=&rr, dtin = dat_input, k=&kk, seed=2000, plvl = &pl);
%put **** Duration (K=&kk) = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
proc export data=repeat
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Uni-meta\SAS\NEQ\result1-K&kk.-R&rr..csv"
    dbms=csv replace;
run;


ods html close;


proc printto; run;
