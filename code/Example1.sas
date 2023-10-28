/*Lettuce data from van Eeuwijk 1992 TAG 85:92-100*/

data b;
input env @@;
do gen="Pa", "DM", "Pi", "GT", "RW", "Wi", "Tr", "Ls";
  input yield @; output;
end;
datalines;
1 3.113 2.835 2.629 1.988 2.199 2.414 1.248 2.380
2 3.379 3.222 2.848 2.823 3.002 2.950 2.176 3.196
3 3.067 2.326 2.511 2.120 2.692 2.598 1.032 2.355
4 3.202 2.663 2.230 1.638 2.187 2.171 1.062 1.599
5 3.921 3.365 3.028 2.653 2.935 2.931 2.007 2.942
6 4.153 3.970 3.444 2.813 2.865 3.232 2.341 3.289
7 4.851 4.512 4.010 3.504 3.135 3.624 3.080 3.612
8 4.547 4.203 3.429 2.944 2.616 3.052 2.817 3.070
9 3.721 3.505 3.337 2.425 2.177 2.525 1.917 2.830
10 3.581 3.298 3.287 2.389 2.159 2.681 1.744 2.726
11 3.312 3.130 2.959 2.280 1.797 2.152 1.365 2.178
12 3.439 3.329 3.254 2.561 2.843 3.035 1.927 3.058
13 3.195 3.047 2.948 2.696 2.610 2.902 1.914 3.138
14 2.890 2.297 2.295 2.237 1.930 2.414 1.462 2.274
15 2.700 2.430 2.172 2.004 2.194 2.392 1.374 2.144
16 3.143 2.710 2.429 2.260 2.406 2.438 1.536 2.464
17 2.746 2.470 2.226 2.126 2.332 2.185 1.287 2.621
18 3.273 2.384 2.555 2.167 2.545 2.386 1.616 2.813
;
data wide;  /*for PLS*/
input env gen1-gen8;
datalines;
1 3.113 2.835 2.629 1.988 2.199 2.414 1.248 2.380
2 3.379 3.222 2.848 2.823 3.002 2.950 2.176 3.196
3 3.067 2.326 2.511 2.120 2.692 2.598 1.032 2.355
4 3.202 2.663 2.230 1.638 2.187 2.171 1.062 1.599
5 3.921 3.365 3.028 2.653 2.935 2.931 2.007 2.942
6 4.153 3.970 3.444 2.813 2.865 3.232 2.341 3.289
7 4.851 4.512 4.010 3.504 3.135 3.624 3.080 3.612
8 4.547 4.203 3.429 2.944 2.616 3.052 2.817 3.070
9 3.721 3.505 3.337 2.425 2.177 2.525 1.917 2.830
10 3.581 3.298 3.287 2.389 2.159 2.681 1.744 2.726
11 3.312 3.130 2.959 2.280 1.797 2.152 1.365 2.178
12 3.439 3.329 3.254 2.561 2.843 3.035 1.927 3.058
13 3.195 3.047 2.948 2.696 2.610 2.902 1.914 3.138
14 2.890 2.297 2.295 2.237 1.930 2.414 1.462 2.274
15 2.700 2.430 2.172 2.004 2.194 2.392 1.374 2.144
16 3.143 2.710 2.429 2.260 2.406 2.438 1.536 2.464
17 2.746 2.470 2.226 2.126 2.332 2.185 1.287 2.621
18 3.273 2.384 2.555 2.167 2.545 2.386 1.616 2.813
;

data e;
input env x1-x8;
datalines;
1 2.1  1136  993  911  881 14.75 11.23 13.15
2 2.1  1345 1277 1250 1815 11.78 13.28 14.93
3 2.1  1700 2191 2586 2556 16.14 16.48 16.37
4 2.1  1076 1090 1323 1065 14.81 13.61 12.74
5 2.2   960  779  539  457 13.61 12.46 10.89
6 2.0   316  482  421  556 12.81 11.23  9.04
7 2.0   145  117  102   42 11.91 10.29  8.05
8 1.9   109   93  127   42 10.96  8.71  7.73
9 2.2   555  504  415  383  9.64  8.50  9.90
10 2.0  641  663  596  780  8.45  7.98 11.78
11 2.0  676  666  541  546  7.70  9.04 12.60
12 1.6 1951 2427 2413 2286  9.22 12.05 14.39
13 1.5 1651 1789 1276 1518 11.78 13.01 15.21
14 1.5 2281 2359 2376 2514 13.15 14.45 15.62
15 1.5 1244 1456 1604 1398 14.07 15.37 16.23
16 2.3 1398 1852 2719 2975 14.51 15.96 16.45
17 1.5 2041 1515 1350  988 14.93 16.23 16.48
18 1.5 1326 1416 1779 1580 15.26 16.39 16.41
;
proc standard data=e out=f mean=0 std=1;
var x1-x8;
run;

%let p=0.0;  /*p=0 for complete data; use 0 < p < 1 to simulate missing data.*/

data a;
merge b f;
by env;
miss=ranuni(1); /*simulate MCAR*/
if miss<&p then yield=.;
run;

/*BASELINE*/
proc mixed data=g covtest method=REML;
class gen env;
model yield=gen/noint solution;
random env;
repeated env/ sub=gen type=ar(1);
run;

/*FA*/

proc mixed data=a lognote;
class gen env;
model yield=gen env;
run;

proc mixed data=a lognote;
class gen env;
model yield=gen env;
repeated env/sub=gen type=ar(1);
run;

     /*FA0(1)*/

proc mixed data=a lognote covtest method=ML;
ods output covparms=cp0;
class gen env;
model yield=gen env;
random x1-x8/sub=gen type=FA0(1);
parms / pdata=cp;
parms (.1)(.1)(.1)(.1)(.1)(.1)(.1)(.1)
      (.1)(.1);
repeated env/sub=gen type=ar(1);
run;

proc mixed data=a lognote covtest method=REML;
ods output covparms=cp;
class gen env;
model yield=gen env;
random x1-x8/sub=gen type=FA0(1);
parms / pdata=cp0;
repeated env/sub=gen type=ar(1);
run;

data cp2;
set cp;
output;
if _N_<9 and _N_>1 then do;
  estimate=0.1; output;
end;
run;

proc mixed data=a lognote;
ods output covparms=cp3;
class gen env;
model yield=gen env;
random x1-x8/sub=gen type=FA0(2);
parms / pdata=cp2;
repeated env/sub=gen type=ar(1);
run;

data cp;
set cp;
if _N_<9;
array slope slope1-slope8;
do i=1 to 8;
  slope[i]=0;
end;
slope[_N_]=estimate;
run;

proc means noprint data=cp;
var slope1-slope8;
output out=slopes sum=;
run;

data g;
set a;
if _N_=1 then set slopes;
array slope slope1-slope8;
array x x1-x8;
z=0;
do i=1 to 8;
  z=z+slope[i]*x[i];
end;
run;

proc mixed data=g method=ML;
class gen env;
model yield=gen z gen*z /ddfm=KR;
random env;
repeated env/ sub=gen type=ar(1);
run;

proc mixed data=g method=REML covtest;
class gen env;
model yield=gen gen*z/noint solution;
random env;
repeated env/ sub=gen type=ar(1);
run;

        /*FA0(2)*/

data cp3;
set cp3;
array s1 s11-s18;
array s2 s21-s28;
do i=1 to 8;
  s1[i]=0;
  s2[i]=0;
end;
if _N_=1 then s11=estimate;
if _N_=2 then s12=estimate;
if _N_=3 then s22=estimate; s21=0;
if _N_=4 then s13=estimate;
if _N_=5 then s23=estimate;
if _N_=6 then s14=estimate;
if _N_=7 then s24=estimate;
if _N_=8 then s15=estimate;
if _N_=9 then s25=estimate;
if _N_=10 then s16=estimate;
if _N_=11 then s26=estimate;
if _N_=12 then s17=estimate;
if _N_=13 then s27=estimate;
if _N_=14 then s18=estimate;
if _N_=15 then s28=estimate;

run;

proc means noprint data=cp3;
var s11-s18 s21-s28;
output out=slopes2 sum=;
run;

data gg;
set a;
if _N_=1 then set slopes2;
array s1 s11-s18;
array s2 s21-s28;
array x x1-x8;
z1=0;
z2=0;
do i=1 to 8;
  z1=z1+s1[i]*x[i];
  z2=z2+s2[i]*x[i];
end;
run;

proc mixed data=gg method=ML;
class gen env;
model yield=gen z1 z2 gen*z1 gen*z2 /ddfm=KR;
random env;
repeated env/ sub=gen type=ar(1);
run;

proc mixed data=gg covtest;
class gen env;
model yield=gen gen*z1 gen*z2/noint solution;
random env;
repeated env/ sub=gen type=ar(1);
run;


/*PLS*/

proc standard data=e out=f2 mean=0 std=1;
var x1-x8;
run;

data h;
merge f2 wide;
by env;
run;

proc pls data=h method=pls noscale nocenter;
model gen1-gen8=x1-x8;
output out=z XSCORE=z;
run;

proc reg data=z;
model z1=x1-x8;
run;

data k;
merge z b;
by env;
run;

proc mixed data=k;
class gen env;
model yield=gen z1 gen*z1 /ddfm=KR;
random env;
repeated env/ sub=gen type=ar(1);
run;

proc mixed data=k covtest method=REML;
class gen env;
model yield=gen gen*z1/noint solution;
random env;
repeated env/ sub=gen type=ar(1);
run;

proc mixed data=k covtest method=ML;
class gen env;
model yield=gen gen*z1/noint;
random env;
repeated env/ sub=gen type=ar(1);
run;

proc mixed data=k covtest method=REML;
class gen env;
model yield=gen gen*z1 gen*z2/noint solution;
random env;
repeated env/ sub=gen type=ar(1);
run;

proc mixed data=k covtest method=ML;
class gen env;
model yield=gen gen*z1 gen*z2/noint;
random env;
repeated env/ sub=gen type=ar(1);
run;


/*FACTORIAL REGRESSION*/

proc mixed data=a covtest method=REML;
ods output solutionF=svd0;
class gen env;
model yield=gen x1-x8 gen*x1 gen*x2 gen*x3 gen*x4 gen*x5 gen*x6 gen*x7 gen*x8/noint solution ddfm=KR outpm=SVD;
random env;
repeated env/ sub=gen type=ar(1);
parms (.1)(.1)(.1)/lowerb=0,0,0;
run;

proc mixed data=a covtest method=ML;
ods output solutionF=svd0;
class gen env;
model yield=gen x1-x8 gen*x1 gen*x2 gen*x3 gen*x4 gen*x5 gen*x6 gen*x7 gen*x8/noint solution ddfm=KR outpm=SVD;
random env;
repeated env/ sub=gen type=ar(1);
parms (.1)(.1)(.1)/lowerb=0,0,0;
run;


/*CRISS-CROSS REGRESSION*/

       /*adjust variable names to code for Macholdt et al. 2022*/
data total;
set a;
year=env;
Treatment=gen;
Total_yield=yield;
run;

       /*get environmental means to initialize*/
proc glimmix data=total;
ods output lsmeans=yearmeans;
class year treatment;
model Total_yield=year treatment;
lsmeans year;
run;

data yearmeans;
set yearmeans;
yearmean=estimate;
keep yearmean year;
run;

proc sort data=total;
by year;
run;

data temp;
merge total yearmeans;
by year;

proc sort data=temp;
by treatment;
run;

     /*regress on yield means to initialize treatment intercepts and slopes
       (this is classical FW regression)*/

proc hpmixed data=temp lognote noprofile;
ods output parameterestimates=PE_init CovParms=cp;
class year treatment;
model Total_yield=treatment treatment*yearmean/noint solution;
random year;
repeated year/subject=treatment type=ar(1);
parms (.2)(.2)(.2)/
lowerb=0,0,0;
nloptions maxiter=500;
run;

data PE_init;
set PE_init;
trt_int=0;
trt_slope=0;
if effect="Treatment"          then trt_int  =estimate;
if effect="yearmean*Treatment" then trt_slope=estimate;
if Treatment ne '';
run;

proc sort data=PE_init;
by treatment;
run;

proc means data=PE_init noprint;
by treatment;
var trt_int trt_slope;
output out=trt_int_slope sum=;
run;

%macro iter;

%let it_max=15;
%let tiny=1e-4;

data fit_old;
value_old=-1e12;
run;

%do iter=1 %to &it_max;

     /*now can start iterating*/

     /*criss*/
proc sort data=total;
by treatment;
run;

data criss;
merge total trt_int_slope;
by treatment;
total_yield_minus_offset=total_yield-trt_int;
run;

proc hpmixed data=criss lognote noprofile;
ods output parameterestimates=PE_year;
class year treatment;
model total_yield_minus_offset
                 =trt_slope
                  trt_slope*x1
                  trt_slope*x2
                  trt_slope*x3
                  trt_slope*x4
                  trt_slope*x5
                  trt_slope*x6
                  trt_slope*x7
                  trt_slope*x8
/noint solution;
random year;
repeated year/subject=treatment type=ar(1);
parms / pdata=cp noiter;
run;

data PE_year;
set PE_year;
intercept=0;
x1_slope=0;
x2_slope=0;
x3_slope=0;
x4_slope=0;
x5_slope=0;
x6_slope=0;
x7_slope=0;
x8_slope=0;
if effect="trt_slope"    then intercept=estimate;
if effect="trt_slope*x1" then x1_slope=estimate;
if effect="trt_slope*x2" then x2_slope=estimate;
if effect="trt_slope*x3" then x3_slope=estimate;
if effect="trt_slope*x4" then x4_slope=estimate;
if effect="trt_slope*x5" then x5_slope=estimate;
if effect="trt_slope*x6" then x6_slope=estimate;
if effect="trt_slope*x7" then x7_slope=estimate;
if effect="trt_slope*x8" then x8_slope=estimate;
run;

proc means data=PE_year noprint;
var intercept x1_slope x2_slope x3_slope x4_slope x5_slope x6_slope x7_slope x8_slope;
output out=yr_int_slope sum=;
run;

data cross;
set total;
if _N_=1 then set yr_int_slope;
yearmean=intercept + x1_slope*x1 + x2_slope*x2 + x3_slope*x3 + x4_slope*x4 + x5_slope*x5 + x6_slope*x6 + x7_slope*x7 + x8_slope*x8;
run;

      /*cross*/

proc hpmixed data=cross lognote noprofile;
ods output parameterestimates=PE_trt FitStatistics=fit_new CovParms=cp2;
class year treatment;
model Total_yield=treatment treatment*yearmean/noint solution;
random year;
repeated year/subject=treatment type=ar(1);
parms / pdata=cp;
run;

data cp;
set cp2;

data fit_new;
set fit_new;
value_new=value;
if _N_=1;
run;

data PE_trt;
set PE_trt;
trt_int=0;
trt_slope=0;
if effect="Treatment"          then trt_int  =estimate;
if effect="yearmean*Treatment" then trt_slope=estimate;
if Treatment ne '';
run;

proc sort data=PE_trt;
by treatment;
run;

proc means data=PE_trt noprint;
by treatment;
var trt_int trt_slope;
output out=trt_int_slope sum=;
run;

       /*standardize slopes to mean=1*/
proc means data=trt_int_slope noprint;
var trt_slope;
output out=mean_slope mean=mean_slope;
run;

data trt_int_slope;
set trt_int_slope;
if _N_=1 then set mean_slope;
trt_slope=trt_slope/mean_slope;
drop mean_slope;
run;

data check;
merge fit_old fit_new;
if abs(value_old-value_new)<&tiny then done='yes'; else done='no ';
iter=&iter;
if done='yes' then do;
  put "iter=" iter;
  dummy=&it_max;
  call symput('iter', dummy);
  put done;
end;
run;

data fit_old;
set fit_new;
value_old=value_new;
keep value_old;
run;

proc print data=check;
run;

%end;

proc print data=cp ; run;
proc print data=cp2; run;

%mend;

%iter;

proc mixed data=cross lognote covtest;
class year treatment;
model Total_yield=treatment treatment*yearmean/ddfm=KR noint solution;
random year;
repeated year/subject=treatment type=ar(1);
parms / pdata=cp2;
run;

proc mixed data=cross lognote covtest;
class year treatment;
model Total_yield=treatment yearmean treatment*yearmean/ddfm=KR;
random year;
repeated year/subject=treatment type=ar(1);
parms / pdata=cp2;
run;

proc mixed data=cross lognote noprofile method=ML;
class year treatment;
model Total_yield=treatment yearmean treatment*yearmean;
random year;
repeated year/subject=treatment type=ar(1);
parms / pdata=cp2;
run;


/*SVD*/

data svd;
set svd;
keep gen env pred;
run;

proc iml;
use svd; 
read all var _num_ into pred;
pred=pred[,2];
close svd;
use svd0; 
read all var _num_ into int;
int=int[1:8,1];
print int;
close svd0;
p=shape(pred,18,8); p=t(p);
p_centered=p-int*j(1,18,1);
print pred p p_centered;
call svd(u,q,v,p_centered);
print u q v;
a=u*diag(q)*t(v);
print a;
create z_svd from v;
append from v;
close z_svd;
run;
quit;

data z_svd;
set z_svd;
z1=col1;
z2=col2;
z3=col3;
env=_N_;
run;

       /*check*/
data w;
merge e z_svd;
by env;
run;

proc reg;
model col1=x1-x8;
run;

data v;
merge a z_svd;
by env;
run;

proc mixed data=v method=REML covtest;
class gen env;
model yield=gen gen*z1/noint solution;
random env;
repeated env/sub=gen type=ar(1);
run;

proc mixed data=v method=ML;
class gen env;
model yield=gen gen*z1;
random env;
repeated env/sub=gen type=ar(1);
run;

proc mixed data=v method=REML;
class gen env;
model yield=gen gen*z1 gen*z2;
random env;
repeated env/sub=gen type=ar(1);
run;

proc mixed data=v method=ML;
class gen env;
model yield=gen gen*z1 gen*z2;
random env;
repeated env/sub=gen type=ar(1);
run;
