# delimit ; 
clear; 
set mem 300m;	
set more 1 ;  
drop _all;
program drop _all;
log using "C:\research\hrs\wealth\death.log", replace; 
use "c:\research\hrs\wealth\dataprep2";
gen FullModel=1;


sort wave;

*keep PI and lagged health around;
gen PIholder=PI;
gen lhealage=age1*lheal;
gen maleage=age1*male;
gen surv=1-died;

***********************GENERATE SOME PROGRAMS HERE****************;
program define genvals;
  local i = 70;
  while `i' <=100 { ;
    quietly sum _UNCOND if age == `i';
    scalar UNCONDITIONAL`i' = r(mean);
    quietly sum _20th_percentile if age == `i';
    scalar _20p`i' = r(mean);
    quietly sum _40th_percentile if age == `i';
    scalar _40p`i' = r(mean);
    quietly sum _60th_percentile if age == `i';
    scalar _60p`i' = r(mean);
    quietly sum _80th_percentile if age == `i';
    scalar _80p`i' = r(mean);
    local i = `i'+1;
  } ;
end;

program define genexpect;
  gen i = 70;
  gen i2 = 71;
  gen probUNCONDITIONAL=1;
  gen prob20 = 1;
  gen prob40 = 1;
  gen prob60 = 1;
  gen prob80 = 1;
  gen probUNCONDITIONALh=1; 
  gen prob20h = 1;
  gen prob40h = 1;
  gen prob60h = 1;
  gen prob80h = 1;
  gen cprobUNCONDITIONAL=0;
  gen cprob20 = 0;
  gen cprob40 = 0;
  gen cprob60 = 0;
  gen cprob80 = 0;

  local i = 70;
  while `i' <=100 { ;
	replace probUNCONDITIONAL=probUNCONDITIONALh*(1-UNCONDITIONAL`i');
	replace cprobUNCONDITIONAL=cprobUNCONDITIONAL+probUNCONDITIONAL;
	replace probUNCONDITIONALh=probUNCONDITIONAL;

	replace prob20=prob20h*(1-_20p`i');
	replace cprob20=cprob20+prob20;
	replace prob20h=prob20;

	replace prob40=prob40h*(1-_40p`i');
	replace cprob40=cprob40+prob40;
	replace prob40h=prob40;

	replace prob60=prob60h*(1-_60p`i');
	replace cprob60=cprob60+prob60;
	replace prob60h=prob60;

	replace prob80=prob80h*(1-_80p`i');
	replace cprob80=cprob80+prob80;
	replace prob80h=prob80;

    local i = `i' + 1;
  } ;
* procedure above kills people at the start of the time period -- i want to kill tham at the end;
* this means add another year to their life;
replace cprob20=cprob20+1;
replace cprob40=cprob40+1;
replace cprob60=cprob60+1;
replace cprob80=cprob80+1;
end;


*this program generates a polynomial in adist;
program define makedist;
        version 3.1;
drop PI2 PI3 PIage;
gen PI2=PI*PI;
gen PI3= PI2*PI;
gen PIage=PI*age1;
end;
set trace off ;
makedist;


program define genPIshift;
	version 3.1;
* now generate percentiles of the distribution;
	replace PI=Pip1;
	makedist;
	gen bp1=ageshift+(PI*bPI1)+(PI2*bPI2)+(PIage*bPIage);
	gen ebp1=exp(bp1);
	gen survp1=ebp1/(1+ebp1);
	replace PI=Pip2;
	makedist;
	gen bp2=ageshift+(PI*bPI1)+(PI2*bPI2)+(PIage*bPIage);
	gen ebp2=exp(bp2);
	gen survp2=ebp2/(1+ebp2);
	replace PI=Pip3;
	makedist;
	gen bp3=ageshift+(PI*bPI1)+(PI2*bPI2)+(PIage*bPIage);
	gen ebp3=exp(bp3);
	gen survp3=ebp3/(1+ebp3);
	replace PI=Pip4;
	makedist;
	gen bp4=ageshift+(PI*bPI1)+(PI2*bPI2)+(PIage*bPIage);
	gen ebp4=exp(bp4);
	gen survp4=ebp4/(1+ebp4);

gen _20th_percentile=1-sqrt(survp1);
gen _40th_percentile=1-sqrt(survp2);
gen _60th_percentile=1-sqrt(survp3);
gen _80th_percentile=1-sqrt(survp4);

drop bp1 bp2 bp3 bp4 ebp1 ebp2 ebp3 ebp4 survp1 survp2 survp3 survp4; * drops intermediate output from genPIshift;

end;


*****************;
* UNCONDITIONAL ESTIMATIONS;
***not conditioning on anything (other than age);
*****************;
* ALL WOMEN FIRST;

logit surv age1 age2 age3 age4 if male==0;
predict survpf;

*****************;
* NOW MEN        ;

*reg surv age1 age2 age3 age4 if male==1;
logit surv age1 age2 age3 age4 if male==1;
predict survpm;

gen survpfo=sqrt(survpf);
gen survpmo=sqrt(survpm);
gen diedpmo=1-survpmo;
gen diedpfo=1-survpfo;





******************************************;
* NOW THE FULL MODEL                      ;
******************************************;
replace PI=PIholder; makedist;

if FullModel==1{;
logit surv age1 age2 age3 PI PI2 PIage lheal lhealage male maleage;

gen ageshift=_b[_cons]+(_b[age1]*age1)+(_b[age2]*age2) +(_b[age3]*age3); 
gen ageshiftFINAL=ageshift; * output ageshifter;
gen healshift=_b[lheal]+_b[lhealage]*age1;
gen maleshift=_b[male]+_b[maleage]*age1;
gen PIshift=_b[PI]+_b[PIage]*age1;
gen bPI1 =_b[PI];
gen bPI2=_b[PI2];
gen bPIage=_b[PIage];
};

if FullModel==0{;
logit surv age1 age2 age3 male maleage;

gen ageshift=_b[_cons]+(_b[age1]*age1)+(_b[age2]*age2) +(_b[age3]*age3); 
gen ageshiftFINAL=ageshift; * output ageshifter;
gen healshift=0;
gen maleshift=_b[male]+_b[maleage]*age1;
gen PIshift=0;
gen bPI1 =0;
gen bPI2=0;
gen bPIage=0;
};


drop if age1<70|age1>102;
sort age;
drop if age==age[_n-1];
save c:\research\hrs\wealth\death, replace;
use c:\research\hrs\wealth\death, replace;
drop if age>100;

***********************************;
* GRAPHS                           ;
***********************************;

set scheme s1mono;

* check out mortality rates, by age;
* women versus men;

twoway connected diedpmo diedpfo age,ysc(r(0 .5)) ylabel(0(.1).5) msymbol(T O) lpattern (solid dash) ti("Probability of Death, Women and Men", size(medsmall)) ytitle("") saving(menandwomendeath, replace);

/*Command msymbol(T O) determines the shape of the marker, where T stands for a big triange and O stands for a big circle, (small t and o would
stand for small triangle and circle, respectively);*/

graph use menandwomendeath.gph;
graph display, ysize(6) xsize(8.25);*this command changes the size of the graph, delete it in order to get a default size;

graph export c:\research\papers\wealth\menandwomendeath.eps, replace;


**********************************************;
* generate predictions for healthy women;
**********************************************;
genPIshift;
gen _UNCOND=1-sqrt(survpf);

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, ysc(r(0 .5)) ylabel(0(.1).5) msymbol(T D S O) lpattern(solid longdash shortdash longdash_dot) ti("Probability of Death, by Permanent Income Percentile, Women in Good Health", size(medsmall)) ytitle("") saving(womendeathg, replace);
graph use womendeathg.gph;
graph display, ysize(6) xsize(8.25);

graph export c:\research\papers\wealth\womendeathg.eps, replace;
graph export c:\research\papers\wealth\womendeathg.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;

genvals;
genexpect;
sum cprob*;
drop cprob* prob* i* ;
drop _*;
**********************************************;
* now generate predictions for unhealthy women;
**********************************************;
replace ageshift=ageshift+healshift; * use this to generate profiles for bad health;
genPIshift;
replace ageshift=ageshift-healshift; * use this to generate profiles for bad health;
gen _UNCOND=1-sqrt(survpf);

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ysc(r(0 .5)) ylabel(0(.1).5) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Probability of Death, by Permanent Income Percentile, Women in Bad Health", size(medsmall)) ytitle("") saving(womendeathb, replace);
graph use womendeathb.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\womendeathb.eps, replace;
graph export c:\research\papers\wealth\womendeathb.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;

genvals;
genexpect;
sum cprob*;
drop cprob* prob* i* ;
drop _*;
**********************************************;
* generate predictions for healthy men;
**********************************************;
replace ageshift=ageshift+maleshift; * use this to generate profiles for bad health;
genPIshift;
replace ageshift=ageshift-maleshift; * use this to generate profiles for bad health;
gen _UNCOND=1-sqrt(survpm);

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ysc(r(0 .5)) ylabel(0(.1).5) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Probability of Death, by Permanent Income Percentile, Men, Good Health", size(medsmall)) ytitle("") saving(mendeathg, replace);
graph use mendeathg.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\mendeathg.eps, replace;
graph export c:\research\papers\wealth\mendeathg.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;

genvals;
genexpect;
sum cprob*;
drop cprob* prob* i* ;
drop _*;
**********************************************;
* now generate predictions for unhealthy men;
**********************************************;
replace ageshift=ageshift+maleshift+healshift; * use this to generate profiles for bad health;
genPIshift;
replace ageshift=ageshift-maleshift-maleshift; * use this to generate profiles for bad health;
gen _UNCOND=1-sqrt(survpm);

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ysc(r(0 .5)) ylabel(0(.1).5) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Probability of Death, by Permanent Income Percentile, Men, Bad Health", size(medsmall)) ytitle("") saving(mendeathb, replace);
graph use mendeathb.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\mendeathb.eps, replace;
graph export c:\research\papers\wealth\mendeathb.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;

genvals;
genexpect;
sum cprob*;
drop cprob* prob* i* ;
drop _*;

use c:\research\hrs\wealth\death, replace;

* output for DP model;
keep age ageshiftFINAL healshift maleshift PIshift bPI2;
order age ageshiftFINAL healshift maleshift PIshift bPI2;
outsheet using c:\research\hrs\wealth\deathprof,nonames replace;

*save c:\research\hrs\wealth\death, replace;
program drop _all;
log close;
