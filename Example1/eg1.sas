**************************************************;
* MC BOUND FOR EXAMPLE 1
*****************************************************;


ods _all_ close;
options nonotes;


data dat1;
input study y precision;
datalines;
	1   -1.55 1.57
	2   -1.49 1.07
	3   -1.33 1.71
	4   -0.35 2.72
	5   -0.19 1.98
	6   -0.43 1.88
	7   -0.61 2.11
	8   -0.97 1.48
	9   -1.64 1.10
	10  -1.19 0.60
	11  -0.28 1.47
	12   0.03 2.00
	13  -0.06 3.90
	14  -0.54 4.56
	;
	run;
proc univariate data=dat1 noprint;
	var study;
	output out=n_study n=s;
run;


/* ANALYTICAL COPAS-JACKSON BOUND (tau=0)*/
data dat2;
	set dat1;
	isgm1=precision;
	isgm2=precision**2;
run;

proc univariate data=dat2 noprint;
	var isgm1 isgm2;
	output out=out1 n=n1 n2 sum=sisgm1 sisgm2;
run;

data CJbound;
	set out1;
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
data dat_k;
	kkk=10000;
run;

%macro bias1(start= , end= , by= );

data dat12;
	set dat_k;
	do j=1 to kkk;
/*	RANDOM SEED*/
/*	call streaminit(0);*/
/*	seed = ceil( (2**31 - 1)*rand("uniform") ); */
/*	z=RAND('NORMal');*/
/*	z=rannor(seed);*/
	z=rannor(2023);
	output;
	end;
/*	drop j;*/
run;

proc sort data=dat12;by z;run;

%local i;
%do  i = 1 %to %eval(%sysfunc(Ceil(%sysevalf((&end - &start ) / &by ) ) ) +1) ;
%let p=%sysevalf(( &start - &by ) + ( &by * &i )) ;
%if  &p <=&end %then %do;
%put &p;
proc sort data=dat1;
by precision;run;

proc optmodel Printlevel=0;

	number s;
	number k;

	var psk{1..s,1..k};
	var ps{1..s};

	num zk{1..k};
	num precision{1..s};

	num si1;
	num si2;

	read data dat_k into k=kkk;
	read data n_study into s=s;
	read data dat1 into [_n_] precision=precision;
  	read data out1 into si1=sisgm1 si2=sisgm2;

 	read data dat12 into [_n_] zk=z;

	/* CONSTAINTS (1)-(4) */
   con cons1 {i in 1..s, j in 1..k}: 	0<=psk[i,j]<=1;
   con cons2 {i in 1..s}: 				0<=ps[i]<=1;
   
   con cons3 {i in 1..s}: 				ps[i]=(1/k)*sum{j in 1..k} psk[i,j];
   con cons4: 							&p.<=1/s*sum{i in 1..s} ps[i];

	/* ASSUMPTION1 */
   con cons5 {i in 1..(s-1)}: 			ps[i]<=ps[i+1]; 
	/* MAXIMUM OF b*/
   max maxb = 1/si2*(1/k)*sum{i in 1..s, j in 1..k} precision[i]*zk[j]*psk[i,j]/ps[i];

   solve;

   create data outdata from p=&p maxb = maxb k=k;

quit;

data result;
	set result outdata;
run;

%end;
%end;

%mend bias1;

/*STORE RESULTS*/
data result;
run;
%bias1(start = 0.1 , end = 1 , by = 0.1);




