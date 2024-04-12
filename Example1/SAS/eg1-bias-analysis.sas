/**/
/*SIMULATION-BASED WORST-CASE BOUND OF META-ANALYSIS OF INTERVENTION*/
/**/
/*EXAMPLE 1*/
/**/

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


/*SIMULATION-BASED WORST-CASE BOUND AND STORE RESULTS*/
ods html;
%include "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example1\SAS\macro_bias1.sas";

/*K=1000*/
%bias1(start=0.1, end=1, by=0.1, dtin = dat_input, k=1000,  seed=2000, plvl = 1);
proc export data=result
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example1\SAS\result1-1000.csv"
    dbms=csv replace;
run;

/*K=2000*/
%bias1(start=0.1, end=1, by=0.1, dtin = dat_input, k=2000,  seed=2000, plvl = 1);
proc export data=result
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example1\SAS\result1-2000.csv"
    dbms=csv replace;
run;

/*K=20000*/
%bias1(start=0.1, end=1, by=0.1, dtin = dat_input, k=20000, seed=2000, plvl = 1);
proc export data=result
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example1\SAS\result1-20000.csv"
    dbms=csv replace;
run;


