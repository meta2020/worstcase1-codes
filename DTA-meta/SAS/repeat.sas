
%let kk = 1000; %let rr = 10; %let pl = 0;
ods html file="IVDresult2-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\";

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
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\IVDresult2-K&kk.-R&rr..csv"
    dbms=csv replace;
run;



%let kk = 2000; 
ods html file="IVDresult2-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\";

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
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\IVDresult2-K&kk.-R&rr..csv"
    dbms=csv replace;
run;



%let kk = 5000; 
ods html file="IVDresult2-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\";

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
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\IVDresult2-K&kk.-R&rr..csv"
    dbms=csv replace;
run;



%let kk = 10000; 
ods html file="IVDresult2-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\";

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
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\IVDresult2-K&kk.-R&rr..csv"
    dbms=csv replace;
run;



%let kk = 20000; 
ods html file="IVDresult2-K&kk.-R&rr..html" path="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\";

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
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SAS\IVDresult2-K&kk.-R&rr..csv"
    dbms=csv replace;
run;

