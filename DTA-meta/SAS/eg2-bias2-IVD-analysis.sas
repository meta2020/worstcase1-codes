/*WORST-CASE BOUND OF THE SROC/SAUC IN META-ANALYSIS OF DIAGNOSTIC TEST ACCURACY */
/*EXAMPLE 2*/

/*CLOSE ALL NOTES*/
proc datasets library=work kill noprint;
quit;
ods _all_ close;
options NONOTES;
/*EXAMPLE 2 IN MAIN: TROPNIN*/
proc import out=data_input 
datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\anadata\eg2-IVD-cc.csv"
dbms = csv replace;
getnames = yes;
run;
proc import out=par_all_wide 
datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\anadata\ml-par-IVD.csv"
dbms = csv replace;
getnames = yes;
run;

/*LOAD MACRO OF CALCULATING BIAS*/
%include "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\macro_bias2.sas";


/*SIMULATION-BASED WORST-CASE BOUND AND STORE RESULTS*/
%let kk = 20000; %let rr = 1; %let pl = 1;
ods html file="IVDresult2-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\";

%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=2000, plvl = &pl, assumption = 1);
%put **** Duration1 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result1;
set result;
as = 1;
run;
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=2000, plvl = &pl, assumption = 2);
%put **** Duration2 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result2;
set result;
as = 2;
run;
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=2000, plvl = &pl, assumption = 3);
%put **** Duration2 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result3;
set result;
as = 3;
run;
data result;
set result1 result2 result3;
run;

proc export data=result
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\IVDresult2-K&kk.-R&rr..csv"
    dbms=csv replace;
run;
ods html close;
