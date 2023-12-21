
/*WORST-CASE BOUND OF META-ANALYSIS OF INTERVENTION*/
/*EXAMPLE 1*/
proc datasets library=work kill noprint;quit;
/*CLOSE ALL NOTES*/
ods _all_ close;
options nonotes;

/*IMPORT RAW DATASET*/
/*VARIABLES: LOGOR AND PRECISION (=1/SIGMA)*/
proc import  out=dat_input
	datafile = "C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example1\eg1-raw.csv"
	dbms = csv replace;
	getnames = yes;
run;

/*CALCULTAE SUM OF INVERSE SIGMA^2 AND THE NUMBER OF PUBLISHED STUDIES*/
data dat2;
	set dat_input;
	isgm1=precision;
	isgm2=precision**2;
run;
proc univariate data=dat2 noprint;
	var isgm1 isgm2;
	output out=dat_sum n=n1 n2 sum=sisgm1 sisgm2;
run;

/* ANALYTICAL COPAS-JACKSON BOUND (tau=0)*/
data dat_cj;
	set dat_sum;
	do p=0.1 to 1 by 0.1;
	group=0;
	k=1;
	n=n1;
	m=n*(1-p)/p;
	u=probit(p);
	bound=1/p/sqrt(6.28)*exp(-0.5*u*u)*sisgm1/sisgm2;
	output;
	end;
	keep p bound group k;
run;


/* MC BASED BOUND */

%macro bias1(start=0.1, end=1, by=0.1);

%local i;
%do  i = 1 %to %eval(%sysfunc(Ceil(%sysevalf((&end - &start ) / &by ) ) ) +1) ;
%let p=%sysevalf(( &start - &by ) + ( &by * &i )) ;
%if  &p <=&end %then %do;
/*%put &p;*/

/*SORT DATA BY INCREASING VARIANCE*/
proc sort data=dat_input; by precision; run;

/*MC BASED BOUND*/
proc optmodel Printlevel=0;

	number n;
	number k;

	var psk{1..n,1..k};
	var ps{1..n};

	num zk{1..k};
	num precision{1..n};

	num si1;
	num si2;

	read data dat_k 		into k=k;
	read data dat_sum 		into si1=sisgm1 si2=sisgm2 n=n1;
	read data dat_input 	into [_n_] precision=precision;
	read data dat_z 		into [_n_] zk=z;

	/* CONSTAINTS (1)-(3) */
	con cons1 {i in 1..n, j in 1..k}: 	0<=psk[i,j]<=1;
	con cons2 {i in 1..n}: 				0<=ps[i]<=1;
	con cons3 {i in 1..n}: 				ps[i]=(1/k)*sum{j in 1..k} psk[i,j];
	con cons4: 							&p.<=1/n*sum{i in 1..n} ps[i];

	/* ASSUMPTION1 */
	con cons5 {i in 1..(n-1)}: 			ps[i]<=ps[i+1]; 

	/* MAXIMUM OF b*/
	max maxb = 1/si2*(1/k)*sum{i in 1..n, j in 1..k} precision[i]*zk[j]*psk[i,j]/ps[i];
	solve;
	create data maxb from p=&p maxb=maxb k=k;

	min minb =1/si2*(1/k)*sum{i in 1..n, j in 1..k} precision[i]*zk[j]*psk[i,j]/ps[i];
	solve;
	create data minb from p=&p minb=minb k=k;

	data outdata; merge maxb minb; by p k; run;

quit;

data result; set result outdata; run;

%end;
%end;

%mend bias1;


%macro repeat1(n, k);

	/* SET THE LENGTH OF Z VECTOR*/
	data dat_k; k=&k; run;

	%local i;
   	%do i=1 %to &n;
	/* GENREATE RANDOM VARIABLE Z*/
	%if &n = 1 %then %do;
	data dat_z;
		set dat_k;
		do j=1 to k;
			z=RAND('NORMal');
			z=rannor(2000);
			output;
		end;
	run;
	%end;
	%else %do;
	data dat_z;
		set dat_k;
		do j=1 to k;
		call streaminit(0);
/*		seed = ceil( (2**31 - 1)*rand("uniform") ); */
		seed = %sysevalf(2000 + &i);
		z=RAND('NORMal');
		z=rannor(seed);
		output;
		end;
	run;
	%end;

	data result;
	run;
	%bias1;

	proc sort data=dat_cj; by descending p; quit;
	proc sort data=result; by descending p; quit;

	data compare;
	merge dat_cj result;
/*	bias = bound-maxb;*/
	group = &i;
	run;

	proc sort data=compare; by descending p; quit;

	data repeat;
		set repeat compare;
	run;

%end;
%mend repeat1;

/*STORE RESULTS*/
%let K_value=2000;
%let R_value=10; 
data repeat;
run;
%repeat1(n=&R_value, k=&K_value);


data repeat;
	set repeat;
	if cmiss(of p)<1;
run;


/*EXPORT RESULTS*/
%let date = %SYSFUNC(PUTN(%sysevalf(%SYSFUNC(TODAY())),DATE9.));
%let time = %sysevalf(%SYSFUNC(TIME()), ceil);
proc export data=repeat
    outfile="C:\Users\zhouy\Documents\GitHub-Bios\worstcase1-codes\Example1\SASresult\result1-R&R_value-K&K_value-&date&time..csv"
    dbms=csv replace;
run;

