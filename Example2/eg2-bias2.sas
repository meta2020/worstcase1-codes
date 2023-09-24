
/*WORST-CASE BOUND OF META-ANALYSIS OF DIAGNOSTIC TEST ACCURACY (SE & SP)*/
/*EXAMPLE 2*/

/*CLOSE ALL NOTES*/
proc datasets library=work kill noprint;
quit;
ods _all_ close;
options NONOTES;

proc import out=data_input 
datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\eg2-IVD-cc.csv"
dbms = csv replace;
getnames = yes;
run;
proc import out=par_all_wide 
datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\ml-par.csv"
dbms = csv replace;
getnames = yes;
run;

data dat_c1;
input c c1 c2;
datalines;
	1   1 0
	2   0 1
;
run;

data dat_c2;
set par_all_wide;
c  = 3;
c1 = 1;
c2 = -tau12/tau22;
keep c c1 c2;
run;

data dat_c;
set dat_c1 dat_c2;
run;



data dat_z_cov(type=COV);
input _TYPE_ $ 1-8 _NAME_ $ 9-16 z1 z2;
datalines;
COV     z1      1 0
COV     z2      0 1
MEAN            0 0
;
run;


/*VARIANCE OF y1 y2*/
data dat_v;
set data_input (keep=TP FN TN FP);
v11 = 1/TP+1/FN;
v22 = 1/TN+1/FP;
v3  = (v11+v22);
run;

%macro bias2(start=0.1, end=1, by=0.1);

%local i;
%do i = 1 %to %eval(%sysfunc(Ceil(%sysevalf((&end - &start ) / &by ) ) ) +1) ;
%let p=%sysevalf(( &start - &by ) + ( &by * &i )) ;
%if &p <=&end %then %do;
%put &p;

proc optmodel Printlevel=0;
	number n;
	number k;
	number cs1;
	number cs2;

	var psk{1..n,1..k};
	var ps{1..n};

	num z1{1..k};
	num z2{1..k};
	num sig1{1..n};
	num sig2{1..n};
	num sig3{1..n};
	num sig4{1..n};

	read data dat_k 		into k=k;
	read data dat_n 		into n=n;
	read data sigma1234 	into [_n_] sig1=sig1 sig2=sig2 sig3=sig3 sig4=sig4;
	read data dat_cleft 	into cs1=sigma_c1 cs2=sigma_c2;
	read data dat_z 		into [_n_] z1=z1 z2=z2;

	/* CONSTAINTS (1)-(3) */
	con cons1{i in 1..n, j in 1..k}: 		0<=psk[i,j]<=1;
	con cons2{i in 1..n}: 					0<=ps[i]<=1;
	con cons3{i in 1..n}: 					ps[i]=(1/k)*sum{j in 1..k} psk[i,j];
	con cons4: 								&p.<=1/n*sum{i in 1..n} ps[i];
	/* ASSUMPTION 2 */
	con cons5{i in 1..(n-1)}: 				ps[i]<=ps[i+1]; 

	max maxb=(1/k)*cs1*sum{i in 1..n, j in 1..k} ((sig1[i]*z1[j]+sig2[i]*z2[j])*psk[i,j]/ps[i])+(1/k)*cs2*sum{i in 1..n, j in 1..k} ((sig3[i]*z1[j]+sig4[i]*z2[j])*psk[i,j]/ps[i]);
	solve;
	create data maxb from p=&p maxb=maxb k=k;

	min minb=(1/k)*cs1*sum{i in 1..n, j in 1..k} ((sig1[i]*z1[j]+sig2[i]*z2[j])*psk[i,j]/ps[i])+(1/k)*cs2*sum{i in 1..n, j in 1..k} ((sig3[i]*z1[j]+sig4[i]*z2[j])*psk[i,j]/ps[i]);
	solve;
	create data minb from p=&p minb=minb k=k;

	data outdata; merge maxb minb; by p k; run;


quit;

data result; set result outdata; run;

%end;
%end ;

%mend bias2;

%macro repeat2(n, k, c, as);

	/*ASSUMPTION 2*/
	%if &as = 1 %then %do;
		proc sort data = dat_v; by descending v11; run;
	%end;
	%if &as = 2 %then %do;
		proc sort data = dat_v; by descending v22; run;
	%end;
	%if &as = 3 %then %do;
		proc sort data = dat_v; by descending v3; run;
	%end;

	data dat_v;
		set dat_v;
		study = _N_;
	run;

	/*EXTRACT TAU*/
	data tau_wide;
		set par_all_wide;
		keep tau11 tau12 tau22;
	run;

	/*CALCULATE SIGMA*/
	data sigma_wide;
		set dat_v;
		if _n_ = 1 then do;
		  set tau_wide;
		end;
		sigma11 = v11+tau11;
		sigma12 = tau12;
		sigma22 = v22+tau22;
		keep study sigma11 sigma12 sigma22;
	run;

	data sigma_mat;
		set sigma_wide; sens = 1; spec = 0; output;
		set sigma_wide; sens = 0; spec = 1; output;
	run;
	data sigma_mat;
		set sigma_mat;
		if sens = 1 & spec = 0 then sigma22=sigma12;
		if sens = 0 & spec = 1 then sigma11=sigma12;
		drop sigma12;
	run;

	data dat_n;
		set dat_v end=eof;
		if eof then  output;
		/*keep id;*/
		rename study=n;
	run;
	proc iml;
		/*  Read data into IML ;*/
		use dat_n; read all; close dat_n;
		use dat_c; read all where (c=&c); close dat_c;

		c = c1 || c2;

		sum_sigma_inv = {0 0 ,0 0};
		do i=1 to n;
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

	%local i;
   	%do i=1 %to &n;
	/* GENERATE THE BIVARIATE RANDOM VARIABLE Z */
	%if &n = 1 %then %do;
	proc simnormal data=dat_z_cov out=dat_z nr = &K seed = 2023; 
		var z1-z2;
	run;
	%end;
	%else %do;
	proc simnormal data=dat_z_cov out=dat_z nr = &K; 
		var z1-z2;
	run;
	%end;
	data dat_k; k = &K; run;

	/*Sigma_invhalfz_ik*/
	/*CALCULATE SIGMA^(-1/2) USING EIGENVALUE DECOMPOSITION*/
	proc iml;

		use dat_n; read all; close dat_n;

		free X;

		do i=1 to n;

			use sigma_mat;
			read all var {sigma11 sigma22} where(study=i);  
			close sigma_mat;

			sigma = sigma11 || sigma22;   

			use dat_z;
			read all;
			close dat_z;

			z_k = z1||z2;

			eigval2 = eigval(sigma)##(-0.5) || eigval(sigma)##(-0.5);
			sigma_invhalf =  eigvec(sigma) * (  eigval2 # t(eigvec(sigma)) );
			X = X // sigma_invhalf;

		end;

		create sigma_invhalf from X[colname={sigma_invhalf1 sigma_invhalf2}];
		append from X;
		close;

	quit;
	/*sigma_invhalf: add study label*/
	data sigma_invhalf;
		merge sigma_mat(keep = study sens spec) sigma_invhalf;
	run;
	/*sigma1234:  sig1i sig2i sig3i sig4i*/
	data sigma12;
		set sigma_invhalf;
		where sens=1;
		rename SIGMA_INVHALF1 = sig1 SIGMA_INVHALF2 = sig2;
	run;
	data sigma34;
		set sigma_invhalf;
		where spec=1;
		rename SIGMA_INVHALF1 = sig3 SIGMA_INVHALF2 = sig4;
	run;
	data sigma1234;
		merge sigma12 sigma34;
	run;

	/*CALCULATE BIAS*/
	data result; run;
	%bias2;

	data compare;
	set result;
	group = &i;
	run;

	proc sort data=compare; by descending p; quit;

	data repeat;
		set repeat compare;
	run;

%end;
%mend repeat2;

/*STORE RESULTS*/
/*==================================*/
/*	Parameter n: repeat times
	Parameter c: which b to maximize                            
 		c = 1: max b for logit-Se                   
		c = 2: max b for logit-Sp    
		c = 3: max b for SROC
	Parameter K: the length of random variable Z
	Parameter as: assumption
		as = 1: sort by v11
		as = 2: sort by v22
		as = 3: sort by v3=v11+v22
/*==================================*/

%let R_value  = 1;
%let K_value  = 2000;
%let c_value  = 3;
%let sort_value = 3;
data repeat;
run;
%repeat2(n = &R_value, k = &K_value, c = &c_value, as = &sort_value);

data repeat;
	set repeat;
	if cmiss(of p)<1;
run;

data b0;
	set repeat;
	by group;
	if first.group then output;
	keep p group maxb minb;
	rename maxb = maxb0 minb = minb0;
run;
data repeat2;
	merge repeat b0;
	by group;
	maxb2 = maxb - maxb0;
	minb2 = minb - minb0;
run;

/*EXPORT RESULTS*/
%let date = %SYSFUNC(PUTN(%sysevalf(%SYSFUNC(TODAY())),DATE9.));
%let time = %sysevalf(%SYSFUNC(TIME()), ceil);
proc export data=repeat2
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example2\SASresult\result2-R&R_value-K&K_value-c&&c_value-as&&sort_value-&date&time..csv"
    dbms=csv replace;
run;
