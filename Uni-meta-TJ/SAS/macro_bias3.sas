/*---------------------------------------------------------------------------------*/
/* MACRO FOR OPTIMIZING MAXIMUM/MINIMUM BIAS GIVE MARGINAL SELECTION PROBABILITIES */
/*---------------------------------------------------------------------------------*/

%macro bias1opt(start=0.1, end=1, by=0.1, plvl=0);
%local i;
%do  i = 1 %to %eval(%sysfunc(Ceil(%sysevalf((&end - &start ) / &by ) ) ) +1) ;
%let p=%sysevalf(( &start - &by ) + ( &by * &i )) ;
%if  &p <=&end %then %do;
%put &p;

/*SORT DATA BY INCREASING 1/VARIANCE*/
proc sort data=&dtin; by precision; run;

/*MC BASED BOUND*/
proc optmodel Printlevel= &plvl;

	number n;
	number k;
	str solstatus;
	str status;

	var pskk{1..n,1..k,1..k};
	var psk{1..n,1..k};
	var ps{1..n};
/*	var d{1..k} BINARY;*/

	num z1{1..k};
	num z2{1..k};
	num precision{1..n};

	num si1;
	num si2;

	read data dat_k 		into k=k;
	read data dat_sum 		into si1=sisgm1 si2=sisgm2 n=n1;
	read data &dtin 		into [_n_] precision=precision;
	read data dat_z 		into [_n_] z1=z1 z2=z2;

	/* CONSTAINTS (1)-(3) */
	con cons1 {i in 1..n, j in 1..k}: 				0<=psk[i,j]<=1;
	con cons2 {i in 1..n, j in 1..k, l in 1..k}: 	0<=pskk[i,j,l]<=1;
	con cons3 {i in 1..n}: 							0<=ps[i]<=1;
	con cons4 {i in 1..n, j in 1..k}: 				psk[i,j]=(1/k)*sum{l in 1..k} pskk[i,j,l];
	con cons5 {i in 1..n}: 							ps[i]=(1/k)*sum{j in 1..k} psk[i,j];
	con cons6: 										&p<=1/n*sum{i in 1..n} ps[i];


	/* ASSUMPTION1 */
	con cons7 {i in 1..(n-1), j in 1..k}: 			psk[i,j]>=psk[i+1,j]; 
/*	con cons8 {i in 1..(n-1), j in 1..k}: 			psk[i,j]<=psk[i+1,j]; */

	/* MAXIMUM OF b*/
	max maxb = 1/si2*(1/(k*k))*(sum{i in 1..n, j in 1..k, l in 1..k} precision[i]*z1[j]*pskk[i,j,l]/ps[i]);
	solve;
	status = _STATUS_;
	solstatus = _SOLUTION_STATUS_;
	create data maxb from p=&p k=k maxb=maxb status1=status solstatus1=solstatus s=(floor(n/&p.));

	data outdata; set maxb; run;

quit;

data result; set result outdata; run;

%end;
%end;
%mend bias1opt;

/*---------------------------------------*/
/* MAIN MACRO FOR STANDARD META-ANALYSIS */
/*---------------------------------------*/
data dat_z_cov(type=COV);
input _TYPE_ $ 1-8 _NAME_ $ 9-16 z1 z2;
datalines;
COV     z1      1 0
COV     z2      0 1
MEAN            0 0
;
run;

%macro bias1(
/*SET A SEQUENCE OF MARGINAL SELECTION PROBABILITY, NUMBER OF REPEAT*/
start=0.1, end=1, by=0.1, repeatn=1,
/*INPUT DATA, REQUIRED VAIRABLES: study y precision (=1/srqt(var+tau2))*/
dtin = dat_input, 
/*NUMBER OF RANDOM NUMBERS IN MONTE-CARLO SIMULATION */
k=200, 
/*SEED OF THE RANDOM NUMBER*/
seed=2000, 
/*PRINTLEVEL IN PROC OPTMODEL, */
plvl = 0);

data dat2;
	set &dtin;
	isgm1=precision;
	isgm2=precision**2;
run;
proc univariate data=dat2 noprint;
	var isgm1 isgm2;
	output out=dat_sum n=n1 n2 sum=sisgm1 sisgm2;
run;


data repeat; run;
/* SET THE LENGTH OF Z VECTOR*/
data dat_k; k=&k; run;
%local i;
%do i=1 %to &repeatn;
	/* GENREATE RANDOM VARIABLE Z*/
	%if &repeatn = 1 %then %do;
		proc simnormal data=dat_z_cov out=dat_z nr = &K seed = &seed; 
			var z1-z2;
		run;
		%end;

		%else %do;
		proc simnormal data=dat_z_cov out=dat_z nr = &K seed = %sysevalf(&seed + &i); 
			var z1-z2;
		run;
		%end;
		data dat_k; k = &k; run;
	%end;

	/*PREPARE RESULT DATA*/
	data result; run;
	%bias1opt(start=&start, end=&end, by=&by, plvl=&plvl);

	proc sort data=result; by descending p; quit;
	data result; set result; group = &i; if cmiss(of p)<1; run;
	data repeat; set repeat result; if cmiss(of p)<1; run;
%end;

%mend bias1;
