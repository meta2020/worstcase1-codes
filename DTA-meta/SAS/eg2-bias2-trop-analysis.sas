/*-------------------------------------------------------------------------------*/
/*WORST-CASE BOUND OF THE SROC/SAUC IN META-ANALYSIS OF DIAGNOSTIC TEST ACCURACY */
/*                                                                               */
/* EXAMPLE 2: TROPNIN                                                            */
/*-------------------------------------------------------------------------------*/

/*CREATE LOG FILE*/
proc printto log = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\trop\log-trop.txt"; run;

/*CLOSE ALL NOTES*/
proc datasets library=work kill noprint; quit;
ods _all_ close;
options NONOTES;



/*EXAMPLE 2 IN MAIN: TROPNIN*/
proc import out=data_input 
datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\anadata\eg2-trop-cc.csv"
dbms = csv replace;
getnames = yes;
run;
proc import out=par_all_wide 
datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\anadata\ml-par-trop.csv"
dbms = csv replace;
getnames = yes;
run;


/*LOAD MACRO OF CALCULATING BIAS*/
%include "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\macro_bias2.sas";

/*SIMULATION-BASED WORST-CASE BOUND AND STORE RESULTS*/
/*SET HOW MANY TIMES FOR REPEAT*/
%let rr = 1; 
/*SET WHETHER TO PRINT THE OPTMODEL RESULT*/
%let pl = 1;
/*SET SEED*/
%let seedn = 2000;
/*SET HOW MANY RANDOM VARIABLES TO USE*/

%let kk = 1000;
ods html file="TROPresult2-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\trop\";
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=&seedn, plvl = &pl, assumption = 1);
%put **** Duration1 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result1; set repeat; as = 1; run;
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=&seedn, plvl = &pl, assumption = 2);
%put **** Duration2 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result2; set repeat; as = 2; run;
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=&seedn, plvl = &pl, assumption = 3);
%put **** Duration2 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result3; set repeat; as = 3; run;
data result; set result1 result2 result3; run;
proc export data=result
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\trop\TROPresult2-K&kk.-R&rr..csv"
    dbms=csv replace;
run;

%let kk = 2000;
ods html file="TROPresult2-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\trop\";
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=&seedn, plvl = &pl, assumption = 1);
%put **** Duration1 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result1; set repeat; as = 1; run;
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=&seedn, plvl = &pl, assumption = 2);
%put **** Duration2 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result2; set repeat; as = 2; run;
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=&seedn, plvl = &pl, assumption = 3);
%put **** Duration2 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result3; set repeat; as = 3; run;
data result; set result1 result2 result3; run;
proc export data=result
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\trop\TROPresult2-K&kk.-R&rr..csv"
    dbms=csv replace;
run;

%let kk = 5000;
ods html file="TROPresult2-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\trop\";
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=&seedn, plvl = &pl, assumption = 1);
%put **** Duration1 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result1; set repeat; as = 1; run;
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=&seedn, plvl = &pl, assumption = 2);
%put **** Duration2 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result2; set repeat; as = 2; run;
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=&seedn, plvl = &pl, assumption = 3);
%put **** Duration2 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result3; set repeat; as = 3; run;
data result; set result1 result2 result3; run;
proc export data=result
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\trop\TROPresult2-K&kk.-R&rr..csv"
    dbms=csv replace;
run;

%let kk = 20000;
ods html file="TROPresult2-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\trop\";
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=&seedn, plvl = &pl, assumption = 1);
%put **** Duration1 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result1; set repeat; as = 1; run;
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=&seedn, plvl = &pl, assumption = 2);
%put **** Duration2 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result2; set repeat; as = 2; run;
%let start_time = %sysfunc(datetime());
%bias2(start=0.1, end=1, by=0.1, repeatn = &rr, dtin = data_input, dpar = par_all_wide, k=&kk, seed=&seedn, plvl = &pl, assumption = 3);
%put **** Duration2 = %sysfunc(putn(%sysevalf(%sysfunc(datetime()) - &start_time),time8.)) ****; 
data result3; set repeat; as = 3; run;
data result; set result1 result2 result3; run;
proc export data=result
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\DTA-meta\SAS\trop\TROPresult2-K&kk.-R&rr..csv"
    dbms=csv replace;
run;

proc printto; run;
