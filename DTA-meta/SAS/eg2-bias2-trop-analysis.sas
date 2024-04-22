/*--------------------------------------------------------------------------------*/
/* WORST-CASE BOUND OF THE SROC/SAUC IN META-ANALYSIS OF DIAGNOSTIC TEST ACCURACY */
/* EXAMPLE 2                                                                      */
/*--------------------------------------------------------------------------------*/

/*CLOSE ALL NOTES*/
proc datasets library=work kill noprint;
quit;
ods _all_ close;
options NONOTES;

/*EXAMPLE 2 IN MAIN: TROPNIN*/
proc import out=data_input 
datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\anadata\eg2-trop-cc.csv"
dbms = csv replace;
getnames = yes;
run;
proc import out=par_all_wide 
datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\anadata\ml-par-trop.csv"
dbms = csv replace;
getnames = yes;
run;

/*LOAD MACRO OF CALCULATING BIAS*/
%include "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\macro_bias2.sas";


/*SIMULATION-BASED WORST-CASE BOUND AND STORE RESULTS*/

%let kk = 1000; %let rr = 1; %let pl = 1;
ods html file="TROPresult2-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\";

%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=2000, plvl = &pl, assumption = 1);
data result1;
set repeat;
as = 1;
run;
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=2000, plvl = &pl, assumption = 2);
data result2;
set repeat;
as = 2;
run;
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=2000, plvl = &pl, assumption = 3);
data result3;
set repeat;
as = 3;
run;
data result;
set result1 result2 result3;
run;
%put **** Duration = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 

proc export data=result
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\TROPresult2-K&kk.-R&rr..csv"
    dbms=csv replace;
run;
