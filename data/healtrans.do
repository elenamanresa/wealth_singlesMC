# delimit ; 
clear; 
set mem 300m;	
set more 1 ;  
drop _all;
program drop _all;
capture log close;
log using c:\research\hrs\wealth\healtrans.log, replace ; 

use c:\research\hrs\wealth\dataprep2;

* this program generates both a               ;
* 2 year and a 1 year health transition matrix;

* separate program for the dead people (death.do); 
drop if dead==1;
*keep PI and lagged health around;
gen PIholder=PI;
gen lhealage=lheal*age1;
gen maleage=male*age1;

*IMPORTANT:_20pg stands for gg (NOT 1-gg!!!!) for 20th percentile and _20pb stands for bb for 20th percentile;

program define genvals;
  local i = 70;
  while `i' <=100 { ;
    * quietly sum g_20 if age == `i';
 *   scalar _20pg`i' = r(mean);

	gen _20pg`i'=.;
	replace _20pg`i'= (1-g_20) if age==`i';
	gen _20pb`i'=.;
	replace _20pb`i'=b_20 if age==`i';
	
	gen _40pg`i'=.;
	replace _40pg`i'= (1-g_40) if age==`i';
	gen _40pb`i'=.;
	replace _40pb`i'=b_40 if age==`i';



	gen _60pg`i'=.;
	replace _60pg`i'= (1-g_60) if age==`i';
	gen _60pb`i'=.;
	replace _60pb`i'=b_60 if age==`i';	


	gen _80pg`i'=.;
	replace _80pg`i'= (1-g_80) if age==`i';
	gen _80pb`i'=.;
	replace _80pb`i'=b_80 if age==`i';

    local i = `i'+1;
*sum _*;
  } ;
end;


program define ONEyear;
sum _20pb*;
*  gen i = 70;
*  gen i2 = 71;

 gen prob20g = 1;
  gen prob40g = 1;
  gen prob60g = 1;
  gen prob80g = 1;

  gen prob20b = 1;
  gen prob40b = 1;
  gen prob60b = 1;
  gen prob80b = 1;

gen A=0;
gen B=0;
gen C=0;

  local i = 70;
  while `i' <=100 { ;

	replace A=(2-_20pb`i' - _20pg`i');
	replace B=(2*_20pg`i'-2);
	replace C=1+(_20pb`i'*_20pb`i') - (_20pb`i' + _20pg`i' );	

	replace prob20b=(-B+sqrt( (B*B) - 4*(A*C) ))/(2*A) if age==`i';
	replace prob20g=((prob20b*prob20b)+1-prob20b-_20pb`i')/(1-prob20b) if age==`i';


	replace A=(2-_40pb`i' - _40pg`i');
	replace B=(2*_40pg`i'-2);
	replace C=1+(_40pb`i'*_40pb`i') - (_40pb`i' + _40pg`i' );	
	replace prob40b=(-B+sqrt( (B*B) - 4*(A*C) ))/(2*A) if age==`i';
	replace prob40g=((prob40b*prob40b)+1-prob40b-_40pb`i')/(1-prob40b) if age==`i';


	replace A=(2-_60pb`i' - _60pg`i');
	replace B=(2*_60pg`i'-2);
	replace C=1+(_60pb`i'*_60pb`i') - (_60pb`i' + _60pg`i' );	
	replace prob60b=(-B+sqrt( (B*B) - 4*(A*C) ))/(2*A) if age==`i';
	replace prob60g=((prob60b*prob60b)+1-prob60b-_60pb`i')/(1-prob60b) if age==`i';

	replace A=(2-_80pb`i' - _80pg`i');
	replace B=(2*_80pg`i'-2);
	replace C=1+(_80pb`i'*_80pb`i') - (_80pb`i' + _80pg`i' );	
	replace prob80b=(-B+sqrt( (B*B) - 4*(A*C) ))/(2*A) if age==`i';
	replace prob80g=((prob80b*prob80b)+1-prob80b-_80pb`i')/(1-prob80b) if age==`i';


    local i = `i' + 1;
  } ;
drop A B C;
end;


 *this program generates a polynomial in adist;
program define makedist;
        version 3.1;
	drop PI2 PI3 PIage ;
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
	gen heal1=ebp1/(1+ebp1);
	replace PI=Pip2;
	makedist;
	gen bp2=ageshift+(PI*bPI1)+(PI2*bPI2)+(PIage*bPIage);
	gen ebp2=exp(bp2);
	gen heal2=ebp2/(1+ebp2);
	replace PI=Pip3;
	makedist;
	gen bp3=ageshift+(PI*bPI1)+(PI2*bPI2)+(PIage*bPIage);
	gen ebp3=exp(bp3);
	gen heal3=ebp3/(1+ebp3);
	replace PI=Pip4;
	makedist;
	gen bp4=ageshift+(PI*bPI1)+(PI2*bPI2)+(PIage*bPIage);
	gen ebp4=exp(bp4);
	gen heal4=ebp4/(1+ebp4);

	gen _20th_percentile=heal1;
	gen _40th_percentile=heal2;
	gen _60th_percentile=heal3;
	gen _80th_percentile=heal4;
	drop bp1 bp2 bp3 bp4 ebp1 ebp2 ebp3 ebp4 heal1 heal2 heal3 heal4; * drops intermediate output from genPIshift;
end;

**************ESTIMATION HERE ******************************;

logit heal age1 age2 age3 lheal lhealage PI PI2 PIage male maleage;

gen ageshift=_b[_cons]+(_b[age1]*age1)+(_b[age2]*age2) +(_b[age3]*age3);
gen ageshiftFINAL=_b[_cons]+(_b[age1]*age1)+(_b[age2]*age2) +(_b[age3]*age3);
gen healshift=_b[lheal]+_b[lhealage]*age1;
gen PIshift=_b[PI]+_b[PIage]*age1;
gen bPI1=_b[PI];
gen bPI2=_b[PI2];
gen bPIage=_b[PIage];
gen bmale=_b[male];
gen maleshift=_b[male]+_b[maleage]*age1;

*rename age1 age;
sort age HHID;
drop if age>102|age<70;
drop if age==age[_n-1];
save c:\research\hrs\wealth\healtrans, replace;
drop if age>100;

***********************************;
* GRAPHS                           ;
***********************************;

set scheme s1mono;

********** generate predicted transitions, WOMEN IN GOOD HEALTH;
genPIshift;
gen g_20=_20th_percentile;
gen g_40=_40th_percentile;
gen g_60=_60th_percentile;
gen g_80=_80th_percentile;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, ysc(r(0 .3)) ylabel(0(.1).3) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Probability of Being in Bad Health, by Permanent Income Percentile,", size(medsmall) span  color(black)) subtitle("In Good Health 2 Years Ago, Women", size(medsmall) span color(black)) l1(Prob(health=bad)) ytitle("") saving(healg, replace);

graph use healg.gph;
graph display, ysize(6) xsize(8.25);
*graph export c:\research\papers\wealth\healgw.eps, replace;
*graph export c:\research\papers\wealth\healgw.emf, replace;
drop _*;

********** generate predicted transitions, WOMEN IN BAD HEALTH;
replace ageshift=ageshift+healshift;
genPIshift;
replace ageshift=ageshift-healshift;
gen b_20=_20th_percentile;
gen b_40=_40th_percentile;
gen b_60=_60th_percentile;
gen b_80=_80th_percentile;

genvals;
ONEyear;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ysc(r(0.7 .95)) ylabel(0.7(.1).95) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Probability of Being in Bad Health, by Permanent Income Percentile,", size(medsmall) span color(black)) subtitle("In Bad Health 2 Years Ago, Women", siz(medsmall) span color(black)) l1(Prob(health=bad)) ytitle("") saving(healb, replace);
graph use healb.gph;
graph display, ysize(6) xsize(8.25);
*graph export c:\research\papers\wealth\healbw.eps, replace;
*graph export c:\research\papers\wealth\healbw.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;


sort age;
by age: sum prob*;

*ONE YEAR GRAPHS;

drop _*;

gen _20th_percentile = 1-prob20g;
gen _40th_percentile = 1-prob40g;
gen _60th_percentile = 1-prob60g;
gen _80th_percentile = 1-prob80g;
twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, ysc(r(0 .3)) ylabel(0(.1).3) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Probability of Being in Bad Health, by Permanent Income Percentile,", size(medsmall) span color(black)) subtitle("In Good Health 1 Year Ago, Women", siz(medsmall) span color(black)) l1(Prob(health=bad)) ytitle("") saving(healg1, replace);
graph use healg1.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\healgw1.eps, replace;
graph export c:\research\papers\wealth\healgw1.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;
drop _*;


gen _20th_percentile = prob20b;
gen _40th_percentile = prob40b;
gen _60th_percentile = prob60b;
gen _80th_percentile = prob80b;
twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ysc(r(0.7 .95)) ylabel(0.7(.1).95) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot)  ti("Probability of Being in Bad Health, by Permanent Income Percentile,", size(medsmall) span color(black)) subtitle("In Bad Health 1 Year Ago, Women", siz(medsmall) span color(black)) l1(Prob(health=bad)) ytitle("") saving(healb1, replace);
graph use healb1.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\healbw1.eps, replace;
graph export c:\research\papers\wealth\healbw1.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;


drop _*;

drop g_* b_* prob*;

*end;

********** generate predicted transitions, MEN IN GOOD HEALTH;


replace ageshift=ageshift+maleshift;
genPIshift;
replace ageshift=ageshift-maleshift;
gen g_20=_20th_percentile;
gen g_40=_40th_percentile;
gen g_60=_60th_percentile;
gen g_80=_80th_percentile;


twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ysc(r(0 .3)) ylabel(0(.1).3) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot)  ti("Probability of Being in Bad Health, by Permanent Income Percentile,", size(medsmall) span color(black)) subtitle("In Good Health 2 Years Ago, Men", size(medsmall) span color(black)) l1(Prob(health=bad)) ytitle("") saving(healgm, replace);
graph use healgm.gph;
graph display, ysize(6) xsize(8.25);
*graph export c:\research\papers\wealth\healgm.eps, replace; 
*graph export c:\research\papers\wealth\healgm.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;
drop _*;



********** generate predicted transitions, MEN IN BAD HEALTH;



replace ageshift=ageshift+maleshift+healshift;
genPIshift;
replace ageshift=ageshift-maleshift-healshift;

gen b_20=_20th_percentile;
gen b_40=_40th_percentile;
gen b_60=_60th_percentile;
gen b_80=_80th_percentile;

genvals;
ONEyear;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ysc(r(0.7 .95)) ylabel(0.7(.1).95) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Probability of Being in Bad Health, by Permanent Income Percentile,", size(medsmall) span color(black)) subtitle("In Bad Health 2 Years Ago, Men", size(medsmall) span color(black)) l1(Prob(health=bad)) ytitle("") saving(healbm, replace);
graph use healbm.gph;
graph display, ysize(6) xsize(8.25);
*graph export c:\research\papers\wealth\healbm.eps, replace;
*graph export c:\research\papers\wealth\healbm.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;
sort age;
by age: sum prob*;

*ONE YEAR GRAPHS;

drop _*;

gen _20th_percentile = 1-prob20g;
gen _40th_percentile = 1-prob40g;
gen _60th_percentile = 1-prob60g;
gen _80th_percentile = 1-prob80g;
twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ysc(r(0 .3)) ylabel(0(.1).3) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Probability of Being in Bad Health, by Permanent Income Percentile,", size(medsmall) span color(black)) subtitle("In Good Health 1 Year Ago, Men", siz(medsmall) span color(black)) l1(Prob(health=bad)) ytitle("") saving(healgm1, replace);
graph use healgm1.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\healgm1.eps, replace;
graph export c:\research\papers\wealth\healgm1.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;
drop _*;


gen _20th_percentile = prob20b;
gen _40th_percentile = prob40b;
gen _60th_percentile = prob60b;
gen _80th_percentile = prob80b;
twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ysc(r(0.7 .95)) ylabel(0.7(.1).95) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Probability of Being in Bad Health, by Permanent Income Percentile,", size(medsmall) span color(black)) subtitle("In Bad Health 1 Year Ago, Men", siz(medsmall) span color(black)) l1(Prob(health=bad)) ytitle("") saving(healbm1, replace);
graph use healbm1.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\healbm1.eps, replace;
graph export c:\research\papers\wealth\healbm1.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;
drop _*;

drop g_* b_* prob*;

use c:\research\hrs\wealth\healtrans, replace;


keep age1 ageshiftFINAL healshift maleshift PIshift bPI2;
order age1 ageshiftFINAL healshift maleshift PIshift bPI2;

outsheet using c:\research\hrs\wealth\healthprof,nonames replace;

drop _all;
program drop _all;
log close;
