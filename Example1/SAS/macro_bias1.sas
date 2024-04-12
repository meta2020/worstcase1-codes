/**/
/*MACRO FOR DERIVING MAXIMUM/MINIMUM BIAS GIVE MARGINAL SELECTION PROBABILITIES*/
/**/
%macro bias1(
/*SET A SEQUENCE OF MARGINAL SELECTION PROBABILITY*/
start=0.1, end=1, by=0.1, 
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

/* SET THE LENGTH OF Z VECTOR*/
data dat_k; k=&k; run;
/* GENREATE RANDOM VARIABLE Z*/
data dat_z;
	set dat_k;
	do j=1 to k;
		z=RAND('NORMal');
		z=rannor(&seed);
		output;
	end;
run;
/*PREPARE RESULT DATA*/
data result; run;

%local i;
%do  i = 1 %to %eval(%sysfunc(Ceil(%sysevalf((&end - &start ) / &by ) ) ) +1) ;
%let p=%sysevalf(( &start - &by ) + ( &by * &i )) ;
%if  &p <=&end %then %do;
/*%put &p;*/

/*SORT DATA BY INCREASING 1/VARIANCE*/
proc sort data=&dtin; by precision; run;

/*MC BASED BOUND*/
proc optmodel Printlevel= &plvl;

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
	read data &dtin 		into [_n_] precision=precision;
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
	solve with nlp;
	create data maxb from p=&p maxb=maxb k=k;

	min minb =1/si2*(1/k)*sum{i in 1..n, j in 1..k} precision[i]*zk[j]*psk[i,j]/ps[i];
	solve with nlp;

	create data minb from p=&p minb=minb k=k;

	data outdata; merge maxb minb; by p k; run;

quit;

data result; set result outdata; run;

%end;
%end;
proc sort data=result; by descending p; quit;
data result;
	set result;
	if cmiss(of p)<1;
run;

%mend bias1;
