/*WORST-CASE BOUND OF THE SROC/SAUC IN META-ANALYSIS OF DIAGNOSTIC TEST ACCURACY */
/*EXAMPLE 2*/

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

/*SIMULATION-BASED WORST-CASE BOUND AND STORE RESULTS*/
ods html;
%include "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes-testCI\Example2\macro_bias2.sas";

/*K=2000*/
%bias2(start=0.1, end=1, by=0.1, dtin = data_input, dpar = par_all_wide, k=2000, seed=2000, plvl = 0, assumption = 1);
data result1;
set result;
as = 1;
run;
%bias2(start=0.1, end=1, by=0.1, dtin = data_input, dpar = par_all_wide, k=2000, seed=2000, plvl = 0, assumption = 2);
data result2;
set result;
as = 2;
run;
%bias2(start=0.1, end=1, by=0.1, dtin = data_input, dpar = par_all_wide, k=2000, seed=2000, plvl = 0, assumption = 3);
data result3;
set result;
as = 3;
run;
data result;
set result1 result2 result3;
run;

proc export data=result
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\result2-2000trop.csv"
    dbms=csv replace;
run;


/*K=1000*/
%bias2(start=0.1, end=1, by=0.1, dtin = data_input, dpar = par_all_wide, k=1000, seed=2000, plvl = 0, assumption = 1);
data result1;
set result;
as = 1;
run;
%bias2(start=0.1, end=1, by=0.1, dtin = data_input, dpar = par_all_wide, k=1000, seed=2000, plvl = 0, assumption = 2);
data result2;
set result;
as = 2;
run;
%bias2(start=0.1, end=1, by=0.1, dtin = data_input, dpar = par_all_wide, k=1000, seed=2000, plvl = 0, assumption = 3);
data result3;
set result;
as = 3;
run;
data result;
set result1 result2 result3;
run;

proc export data=result
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\result2-1000trop.csv"
    dbms=csv replace;
run;

/*K=20000*/
%bias2(start=0.1, end=1, by=0.1, dtin = data_input, dpar = par_all_wide, k=20000, seed=2000, plvl = 0, assumption = 1);
data result1;
set result;
as = 1;
run;
%bias2(start=0.1, end=1, by=0.1, dtin = data_input, dpar = par_all_wide, k=20000, seed=2000, plvl = 0, assumption = 2);
data result2;
set result;
as = 2;
run;
%bias2(start=0.1, end=1, by=0.1, dtin = data_input, dpar = par_all_wide, k=20000, seed=2000, plvl = 0, assumption = 3);
data result3;
set result;
as = 3;
run;
data result;
set result1 result2 result3;
run;

proc export data=result
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\result2-20000trop.csv"
    dbms=csv replace;
run;
