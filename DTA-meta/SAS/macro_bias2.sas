/*------------------------------------------*/
/* GENERATING COVARIANCE OF RANDOM VARIABLE */
/*------------------------------------------*/
data dat_z_cov(type=COV);
input _TYPE_ $ 1-8 _NAME_ $ 9-16 z1 z2;
datalines;
COV     z1      1 0
COV     z2      0 1
MEAN            0 0
;
run;

/*---------------------------------------------------------------------------------*/
/* MACRO FOR OPTIMIZING MAXIMUM/MINIMUM BIAS GIVE MARGINAL SELECTION PROBABILITIES */
/*---------------------------------------------------------------------------------*/

%macro bias2opt(start=0.1, end=1, by=0.1, plvl=0);

%local i;
%do i = 1 %to %eval(%sysfunc(Ceil(%sysevalf((&end - &start ) / &by ) ) ) +1) ;
%let p=%sysevalf(( &start - &by ) + ( &by * &i )) ;
%if &p <=&end %then %do;
/*%put &p;*/

proc optmodel Printlevel=&plvl;
	number n;
	number k;
	number cs1;
	number cs2;
	str solstatus;
	str status;

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
	con cons4: 								&p<=1/n*sum{i in 1..n} ps[i];
	/* ASSUMPTION 2 */
	con cons5{i in 1..(n-1)}: 				ps[i]<=ps[i+1]; 

	max maxb=(1/k)*cs1*sum{i in 1..n, j in 1..k} ((sig1[i]*z1[j]+sig2[i]*z2[j])*psk[i,j]/ps[i])+(1/k)*cs2*sum{i in 1..n, j in 1..k} ((sig3[i]*z1[j]+sig4[i]*z2[j])*psk[i,j]/ps[i]);
	solve;
	status = _STATUS_;
	solstatus = _SOLUTION_STATUS_;
	create data maxb from p=&p k=k maxb=maxb status1=status solstatus1=solstatus;

	min minb=(1/k)*cs1*sum{i in 1..n, j in 1..k} ((sig1[i]*z1[j]+sig2[i]*z2[j])*psk[i,j]/ps[i])+(1/k)*cs2*sum{i in 1..n, j in 1..k} ((sig3[i]*z1[j]+sig4[i]*z2[j])*psk[i,j]/ps[i]);
	solve;
	status = _STATUS_;
	solstatus = _SOLUTION_STATUS_;

	create data minb from p=&p k=k minb=minb status2=status solstatus2=solstatus;

	data outdata; merge maxb minb; by p k; run;


quit;

data result; set result outdata; run;

%end;
%end;

%mend bias2opt;

/*---------------------------------------*/
/* MAIN MACRO FOR STANDARD META-ANALYSIS */
/*---------------------------------------*/

%macro bias2(/*SET A SEQUENCE OF MARGINAL SELECTION PROBABILITY*/
start=0.1, end=1, by=0.1, repeatn=1,
/*INPUT DATA, REQUIRED VAIRABLES: TP FN TN FP*/
dtin = dat_input, 
/*INPUT DATA2, THE ML ESTIMATES OF TAU11 TAU22 TAU12*/
dpar = par_all_wide,
/*NUMBER OF RANDOM NUMBERS IN MONTE-CARLO SIMULATION */
k=200, 
/*SEED OF THE RANDOM NUMBER*/
seed=2000, 
dat_z_cov = dat_z_cov,
/*PRINTLEVEL IN PROC OPTMODEL, */
plvl = 0,
assumption = 1);

	/*VARIANCE OF y1 y2*/
	data dat_v;
	set &dtin (keep=TP FN TN FP);
	v11 = 1/TP+1/FN;
	v22 = (1/TN+1/FP);
	v3  = 1/TP+1/FN+1/TN+1/FP;
	run;

	/*ASSUMPTION 2*/
	%if &assumption = 1 %then %do;
		proc sort data = dat_v; by descending v11 v22; run;
	%end;
	%if &assumption = 2 %then %do;
		proc sort data = dat_v; by descending v22 v11; run;
	%end;
	%if &assumption = 3 %then %do;
		proc sort data = dat_v; by descending v3 v11 v22; run;
	%end;

	data dat_v;
		set dat_v;
		study = _N_;
	run;

	/*EXTRACT TAU*/
	data tau_wide;
		set &dpar;
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

	data dat_c;
		set tau_wide;
		c1 = 1;
		c2 = -tau12/tau22;
	run;

	proc iml;
		/*  Read data into IML ;*/
		use dat_n; read all; close dat_n;
		use dat_c; read all; close dat_c;

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

	data repeat; run;
	%local i;
   	%do i=1 %to &repeatn;
		/* GENERATE THE BIVARIATE RANDOM VARIABLE Z */
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
		%bias2opt(start=&start, end=&end, by=&by, plvl=&plvl);

		proc sort data=result; by descending p; quit;
		data result; set result; group = &i; if cmiss(of p)<1; run;
		data repeat; set repeat result; if cmiss(of p)<1; run;

	%end;
%mend bias2;
