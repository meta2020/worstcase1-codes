
/*WORST-CASE BOUND OF META-ANALYSIS OF PROGNOSIS STUDIES*/
/*EXAMPLE 3*/

/*CLOSE ALL NOTES*/
ods _all_ close;
options NONOTES;


/*IMPORT THE CALCULATED DATA: 3 YEAR AND 5 YEAR*/
proc import  out=raw_data 
datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example3\data-35y.csv"
dbms = csv replace;
getnames = yes;
run;
proc import  out=par_all_wide 
datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example3\ml-par3.csv"
dbms = csv replace;
getnames = yes;
run;

%let t = 3;
data prog_meta;
	set raw_data;
	where ty=&t;
	keep  ty u_sen u_spe  v_sen v_spe v_senspe;
	rename ty=t v_sen=v11  v_spe=v22 v_senspe=v12;
run;

/*ASSUMPTION 3*/
proc sort data = prog_meta;
/*by descending v11 descending v22 descending v12;*/
/*by descending v22 descending v11 descending v12;*/
/*by descending v12;*/
run;


data prog_meta;
set prog_meta;
study = _N_;
run;

data cov_mat;
set prog_meta;
keep t study v11 v12 v22;
run;


/*sigma_wide: create sigma matrix*/
data tau_wide;
set par_all_wide;
keep t tau11 tau12 tau22;
run;
proc sort data = tau_wide; by t; run;
data sigma_wide;
merge cov_mat tau_wide;
by t;
sigma11 = v11+tau11;
sigma12 = v12+tau12;
sigma22 = v22+tau22;
/*keep study sens spec sigma1 sigma2;*/
run;
proc sort data=sigma_wide (keep = t study sigma11 sigma12 sigma22); by t study; run;


data sigma_mat;
set sigma_wide;
sens = 1; spec = 0; output;
sens = 0; spec = 1; output;
run;

data sigma_mat;
set sigma_mat;
if sens = 1 & spec = 0 then sigma22=sigma12;
if sens = 0 & spec = 1 then sigma11=sigma12;
drop sigma12;
run;


/*Data Z, S, C********************************/

/*dat_z: generate (z1, z2)~N(0, 1)*/
%let K = 10000;
%let seed=100;

/* create a TYPE=COV data set */
data dat_z_cov(type=COV);
input _TYPE_ $ 1-8 _NAME_ $ 9-16 z1 z2;
datalines;
COV     z1      1 0
COV     z2      0 1
MEAN            0 0
;
run;
proc simnormal data=dat_z_cov outsim=dat_z
nr = &K        			 /* size of sample */
seed = &seed;	 /* random number seed */
var z1-z2;
run;

/*dat_k: K=1000*/
data dat_k;
k = &K;
run;
/*dat_s: # of observed studies*/
data dat_s;
set prog_meta end=eof;
if eof then  output;
keep study;
rename study=s;
run;


data dat_c;
set par_all_wide;
c1=1;
c2 = -tau12/tau22;
run;
/***********************************************/

proc iml;
/*  Read data into IML ;*/
use dat_s;
read all;
close dat_s;

use dat_c;
read all;
close dat_c;

c = c1 || c2;

sum_sigma_inv = {0 0 ,0 0};
do i=1 to s;

  use sigma_mat;
  read all var {sigma11 sigma22} where(study=i);  
  close sigma_mat;

  sigma = sigma11 || sigma22;   
  sum_sigma_inv = sum_sigma_inv+ inv(sigma);  

end;

cleft =  c * inv(sum_sigma_inv);  
create dat_cleft var {sigma_c1 sigma_c2};
append from cleft;
close dat_cleft;

quit;

/*sigma_invhalfz_ik: sigma power -1/2*/

proc iml;

use dat_s;
read all;
close dat_s;

/*calculate sigma_invhalf*/
free X;

do i=1 to s;

use sigma_mat;
read all var {sigma11 sigma22} where(study=i);  
close sigma_mat;

sigma = sigma11 || sigma22;    /*    print sigma;*/

/*eigen vectors %*% (values^power * t(vectors))*/
eigval2 = eigval(sigma)##(-0.5) || eigval(sigma)##(-0.5);
sigma_invhalf =  eigvec(sigma) * (  eigval2 # t(eigvec(sigma)) );
/*sigma_invhalf_z_ik = sigma_invhalf_i * t(z_k);*/
X = X // sigma_invhalf;

end;

/*calculate sigma_invhalf*/

create sigma_invhalf from X[colname={sigma_invhalf1 sigma_invhalf2}];
append from X;
close;

quit;

/*sigma_invhalf: add study label*/
data sigma_invhalf;
merge sigma_mat(keep = t study sens spec) sigma_invhalf;
run;

/*Sigma1234*/
/*sig1i sig2i sig3i sig4i from sigma_invhalf*/
data sigma1234;
   merge sigma_invhalf (where=(sens=1) rename=(SIGMA_INVHALF1 = sig1 SIGMA_INVHALF2 = sig2))
         sigma_invhalf (where=(spec=1) rename=(SIGMA_INVHALF1 = sig3 SIGMA_INVHALF2 = sig4));
   by study;
   drop sens spec;
run;

%let p=0.9;

/*sauc: calculate sauc*/

%macro tsauc_bound(start= , end= , by= );

%local i;
%do i = 1 %to %eval(%sysfunc(Ceil(%sysevalf((&end - &start ) / &by ) ) ) +1) ;
 %let p=%sysevalf(( &start - &by ) + ( &by * &i )) ;
   %if &p <=&end %then %do;
      %put &p;

proc optmodel Printlevel=0;
	number s;
	number k;
	number cs1;
	number cs2;

	var psk{1..s,1..k};
	var ps{1..s};

	num z1{1..k};
	num z2{1..k};
	num sig1{1..s};
	num sig2{1..s};
	num sig3{1..s};
	num sig4{1..s};

	read data dat_k into k=k;
	read data dat_s into s=s;
	read data sigma1234 into [_n_] sig1=sig1 sig2=sig2 sig3=sig3 sig4=sig4;

  	read data dat_cleft into cs1=sigma_c1 cs2=sigma_c2;
 	read data dat_z into [_n_] z1=z1 z2=z2;

   con cons1{i in 1..s, j in 1..k}: 	0<=psk[i,j]<=1;
   con cons2{i in 1..s}: 				0<=ps[i]<=1;
   con cons3{i in 1..s}: 				ps[i]=(1/k)*sum{j in 1..k} psk[i,j];
   con cons4: 							&p.<=1/s*sum{i in 1..s} ps[i];

/* monotone assumption for sgm */
   con cons5{i in 1..(s-1)}: 			ps[i]<=ps[i+1]; 


   max maxb=(1/k)*cs1*sum{i in 1..s, j in 1..k} ((sig1[i]*z1[j]+sig2[i]*z2[j])*psk[i,j]/ps[i])
				+(1/k)*cs2*sum{i in 1..s, j in 1..k} ((sig3[i]*z1[j]+sig4[i]*z2[j])*psk[i,j]/ps[i]);

   solve;

  create data maxb from p=&p maxb=maxb k=k;

  min minb=(1/k)*cs1*sum{i in 1..s, j in 1..k} ((sig1[i]*z1[j]+sig2[i]*z2[j])*psk[i,j]/ps[i])
				+(1/k)*cs2*sum{i in 1..s, j in 1..k} ((sig3[i]*z1[j]+sig4[i]*z2[j])*psk[i,j]/ps[i]);

   solve;

  create data minb from p=&p minb=minb k=k;

quit;

data out_bound;
merge maxb minb;
by p k;
run;

/*calculate sauc*/
proc iml;

use par_all_wide;
read all;
close par_all_wide;

use out_bound;
read all;
close out_bound;

theta1 = theta1;
theta2 = theta2;
tau12 = tau12;
tau22 = tau22;
p=p;
k=k;
maxb=maxb;
minb=minb;

/*sauc*/
start SROC(x) global(theta1, theta2, tau12, tau22);
	f = logistic( theta1- (tau12/tau22)*(log(x/(1-x)) + theta2) );
   return(f);
finish;

/*Upper bound for sauc*/
start SROC_ub(x) global(theta1, theta2, tau12, tau22, maxb);
	f = logistic( theta1- (tau12/tau22)*(log(x/(1-x)) + theta2) + maxb );
   return(f);
finish;

/*Lower bound for sauc*/
start SROC_lb(x) global(theta1, theta2, tau12, tau22, minb);
	f = logistic( theta1- (tau12/tau22)*(log(x/(1-x)) + theta2) + minb );
   return(f);
finish;

/*integral on a*/
a ={0 1};    
call quad(sauc, "SROC", a);
call quad(sauc_ub, "SROC_ub", a);
call quad(sauc_lb, "SROC_lb", a);

if p=1 then sauc=p||sauc||sauc||maxb;
else sauc=p||sauc_lb||sauc_ub||maxb;

create sauc var { p sauc_lb sauc_ub maxb};
append from sauc;
close sauc;

quit;

data result;
	set result sauc; 
run;

 %end;
%end ;

%mend tsauc_bound;

data result;
run;
%tsauc_bound(start = 0.1 , end = 1, by = 0.1);


data result2;
  set result;
  if cmiss(of p)<1 ;
run;

/*EXPORT RESULTS*/
proc export data=result2
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example3\MCbound3.csv"
    dbms=csv replace;
run;

