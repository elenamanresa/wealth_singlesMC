# delimit ; 
clear; 
set mem 300m;	
set more 1 ;  
drop _all;
program drop _all;
capture log close;
log using c:\research\hrs\wealth\incprof.log, replace ; 
use c:\research\hrs\wealth\dataprep2;
gen useheal=0 ;
********************************************;
* IMPORTANT NOTE ABOUT THE DATA             ;
* WAVE 1 OF THE AHEAD IS PRETTY SCREWED UP  ;
* THE QUESTIONS ASKED ARE LESS COMPREHENSIVE;
* THAN IN OTHER WAVES, AND SO THE MEANS OF  ;
* DIFFERENT VARIABLES ARE LOWER             ;
********************************************;
*drop if (wave==2);
*drop if (wave==3);

replace lhhinc=. if (wave==2|wave==3);
replace lanninc=. if (wave==2|wave==3);
replace lsocy=. if (wave==2|wave==3);
*replace linc=. if (wave==2|wave==3);




sort HHID  realyear female;
drop if HHID==HHID[_n-1] & realyear==realyear[_n-1];

gen healage=age1*heal;
gen maleage=age1*male;


******************PROCEDURES HERE ***********************;
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


*this program generates average income at different percentiles of the income distribution;
program define genPIshift;
        version 3.1;
	replace PI=Pip1;
	makedist;
	gen lannincp1=ageshiftvarc+PI*bPI+PI2*bPI2+(PIage*bPIage);
	replace PI=Pip2;
	makedist;
	gen lannincp2=ageshiftvarc+PI*bPI+PI2*bPI2+(PIage*bPIage);
	replace PI=Pip3;
	makedist;
	gen lannincp3=ageshiftvarc+PI*bPI+PI2*bPI2+(PIage*bPIage);
	replace PI=Pip4;
	makedist;
	gen lannincp4=ageshiftvarc+PI*bPI+PI2*bPI2+(PIage*bPIage);

	gen _20th_percentile=exp(lannincp1);
	gen _40th_percentile=exp(lannincp2);
	gen _60th_percentile=exp(lannincp3);
	gen _80th_percentile=exp(lannincp4);
	drop lannincp*;
end;

if FE==0{;
* OLS;
reg lhhinc age1 age2 age3 age4 age5;
predict lhhincp;

reg lanninc age1 age2 age3 age4 age5;
predict lannincp;

reg lsocy age1 age2 age3 age4 age5;
predict lsocyp;

reg linc age1 age2 age3 age4 if lfpr==0;
replace lfpr=0;
predict lincpl;
};

if FE==1{;
* fixed effects ;

* do fixed-effects regression and get predictions;
/*
xtreg lanninc age1 age2 age3 heal healage PIage maleage, fe i(HHID);
predict lannincp;
gen lannincres=lanninc-lannincp;
gen ageshift=_b[_cons]+(_b[age1]*age1)+(_b[age2]*age2) +(_b[age3]*age3);
gen bheal=_b[heal];
gen healshift=_b[heal]+_b[healage]*age1;
gen maleshift=_b[maleage]*age1;
gen PIshift=_b[PIage]*age1;
gen bPIage=_b[PIage];
*/
if useheal==1{;
xtreg lanninc age1 age2 age3 heal healage , fe i(HHID);
predict lannincp;
gen lannincres=lanninc-lannincp;
gen ageshift=_b[_cons]+(_b[age1]*age1)+(_b[age2]*age2) +(_b[age3]*age3);
gen bheal=_b[heal];
gen healshift=_b[heal]+_b[healage]*age1;
gen maleshift=0;
gen PIshift=0;
gen bPIage=0;
};
if useheal==0{;
xtreg lanninc age1 age2 age3, fe i(HHID);
predict lannincp;
gen lannincres=lanninc-lannincp;
gen ageshift=_b[_cons]+(_b[age1]*age1)+(_b[age2]*age2) +(_b[age3]*age3);
gen bheal=0;
gen healshift=0;
gen maleshift=0;
gen PIshift=0;
gen bPIage=0;
};

*sort age;
*by age: sum lanninc* age*;
*exit;

* now run a regression on the residuals;
reg lannincres PI PI2 male;
predict respred;
replace ageshift=ageshift+_b[_cons];
gen ageshiftFINAL=ageshift;
replace maleshift=maleshift+_b[male];
replace PIshift=PIshift+_b[PI];
gen bPI =_b[PI];
gen bPI2=_b[PI2];
gen bmale=_b[male];

* see if cohort effects are important to worry about;
reg lannincres PI PI2 male cohort2 cohort3 cohort4 cohort5 cohort6 cohort7; * now generate percentiles of the distribution;

* figure out variance of residual;
gen resres=lannincres-respred;
egen sdres=sd(resres);
gen varc=(sdres*sdres)/2;
sum varc;
gen ageshiftvarc=ageshift+varc;

* regress some other income measures on age;
xtreg lsocy age1 age2 age3 age4 age5, fe i(HHID);
predict lsocyp;

xtreg lhhinc age1 age2 age3 age4 age5, fe i(HHID);
predict lhhincp;

xtreg lhhinc age1 age2 age3 age4 age5 if lfpr==0, fe i(HHID);
replace lfpr=0;
predict lincpl;
};

gen hhincp=exp(lhhincp);
gen incp=exp(lincp);
gen annincp=exp(lannincp);
gen socyp=exp(lsocyp);
gen incpl=exp(lincpl);



drop if age1>102|age1<70;
sort age1;
sort age1;
drop if age1==age1[_n-1];
save c:\research\hrs\wealth\incprof, replace;
drop if age>100;

***********************************;
* GRAPHS                           ;
***********************************;

set scheme s1mono;

* GRAPHS;
* by type of income;
twoway connected hhincp annincp incp socyp age, msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) c(ll) ti("Income, by type") saving(income, replace);
graph use income.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\incomebytype.eps, replace;

*for anninc, by asset quantile;
* by type of income;
if FE==1{;

****************** healthy women;
genPIshift;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, ysc(r(7000 20000)) ylabel(7000 10000 15000 20000) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Income, by Permanent Income Percentile, Healthy Women",size(medlarge)) l1(1998 dollars) ytitle("") saving(income, replace);
graph use income.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\incomegw.eps, replace;
graph export c:\research\papers\wealth\incomegw.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;

*************** unhealthy women;
drop _*;
replace ageshiftvarc=ageshiftvarc+healshift;
genPIshift;
replace ageshiftvarc=ageshiftvarc-healshift;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ysc(r(7000 20000)) ylabel(7000 10000 15000 20000) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Income, by Permanent Income Percentile, Unhealthy Women", size(medlarge)) l1(1998 dollars) ytitle("") saving(income, replace);
graph use income.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\incomebw.eps, replace;
graph export c:\research\papers\wealth\incomebw.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;

******************** healthy men;
drop _*;
replace ageshiftvarc=ageshiftvarc+maleshift;
genPIshift;
replace ageshiftvarc=ageshiftvarc-maleshift;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ysc(r(7000 20000)) ylabel(7000 10000 15000 20000) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Income, by Permanent Income Percentile, Healthy Men", size(medlarge)) l1(1998 dollars) ytitle("") saving(income, replace);
graph use income.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\incomegm.eps, replace;
graph export c:\research\papers\wealth\incomegm.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;

******************* unhealthy men;
drop _*;
replace ageshiftvarc=ageshiftvarc+healshift+maleshift;
genPIshift;
replace ageshiftvarc=ageshiftvarc-healshift-maleshift;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ysc(r(7000 20000)) ylabel(7000 10000 15000 20000) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Income, by Permanent Income Percentile, Unhealthy Men", size(medlarge)) l1(1998 dollars) ytitle("") saving(income, replace);
graph use income.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\incomebm.eps, replace;
graph export c:\research\papers\wealth\incomebm.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;

};

*replace ageshift=ageshift;
*keep age ageshift bPI bPI2;
*order age ageshift bPI bPI2;

use c:\research\hrs\wealth\incprof, replace;


keep age ageshiftFINAL healshift maleshift PIshift bPI2;
order age ageshiftFINAL healshift maleshift PIshift bPI2;

sum;


outsheet using c:\research\hrs\wealth\incprof,nonames replace;
drop _all;
program drop _all;
log close;


