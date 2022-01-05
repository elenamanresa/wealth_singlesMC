# delimit ; 
clear; 
set mem 300m;	
set more 1 ;  
drop _all;
program drop _all;
log using c:\research\hrs\wealth\medex.log, replace ; 


use c:\research\hrs\wealth\dataprep2;

* log medcost;
gen lmedcost=ln(medcost);
gen healage=heal*age1;
*gen PIage=heal*age1;
gen maleage=male*age1;


* drop married women;
sort HHID female realyear;
drop if HHID==HHID[_n-1] & realyear==realyear[_n-1];

sort wave;
by wave: sum age1 dead medcost;

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


sum P*;



*this program generates a polynomial in adist;
program define genPIshift;
        version 3.1;
	replace PI=Pip1;
	makedist;
	gen lmedcostp1=ageshiftvarc+PI*bPI+PI2*bPI2+(PIage*bPIage);
	replace PI=Pip2;
	makedist;
	gen lmedcostp2=ageshiftvarc+PI*bPI+PI2*bPI2+(PIage*bPIage);
	replace PI=Pip3;
	makedist;
	gen lmedcostp3=ageshiftvarc+PI*bPI+PI2*bPI2+(PIage*bPIage);
	replace PI=Pip4;
	makedist;
	gen lmedcostp4=ageshiftvarc+PI*bPI+PI2*bPI2+(PIage*bPIage);


	gen _20th_percentile=exp(lmedcostp1);
	gen _40th_percentile=exp(lmedcostp2);
	gen _60th_percentile=exp(lmedcostp3);
	gen _80th_percentile=exp(lmedcostp4);
drop lmedcostp*;
end;



* OLS VERSUS FIXED EFFECTS;
* OLS is not consistent with what JOHN, Cristina and I have decided on--it is only descriptive;
if FE==0{;
reg medcost age1 age2 age3 age4 if male==0 & aa1==1;
predict medpf1;

reg medcost age1 age2 age3 age4 if male==0 & aa2==1;
predict medpf2;

reg medcost age1 age2 age3 age4 if male==0 & aa3==1;
predict medpf3;

reg medcost age1 age2 age3 age4 if male==0 & aa4==1;
predict medpf4;
};

if FE==1{;
* fixed effects ;
*xtreg lmedcost age1 age2 age3 age4 age5, fe i(HHID);
*xtreg lmedcost age1 age2 age3 age4 age5 heal, fe i(HHID);

xtreg lmedcost age1 age2 age3  heal healage maleage PIage, fe i(HHID);
predict lmedcostp;
* get predictions;
gen lmedcostres=lmedcost-lmedcostp;
gen ageshift=_b[_cons]+(_b[age1]*age1)+(_b[age2]*age2) +(_b[age3]*age3);
gen healshift=_b[heal]+_b[healage]*age1;
gen maleshift=_b[maleage]*age1;
gen PIshift=_b[PIage]*age1;
gen bPIage=_b[PIage];
* now run an OLS regression on the residuals;
reg lmedcostres PI PI2 male;
predict respred;
gen b=_b[_cons];
replace ageshift=ageshift+_b[_cons];
gen ageshiftFINAL=ageshift;
gen bPI =_b[PI];
gen bPI2=_b[PI2];
gen bmale=_b[male];
replace maleshift=maleshift+_b[male];
replace PIshift=PIshift+_b[PI]; * do not add in quadratic term, as we keep track of that separately;
* see if cohort effects are important to worry about (as it turns out, they are);
reg lmedcostres PI PI2 cohort2 cohort3 cohort4 cohort5 cohort6 cohort7; * now generate percentiles of the distribution;

* figure out variance of residual;
gen resres=lmedcostres-respred;
egen sdres=sd(resres);
gen varc=(sdres*sdres);

/*
gen varc1=resres*resres;
reg varc1 age1 ;
reg varc1 age1 age2 male heal PI;
sum varc*;
*/


gen ageshiftvarc=ageshift+(varc/2);

* gen medcostp=exp(lmedcostp); * this term is meaningless, in the aggregate;
};


*rename age1 age;
sort age HHID;

drop if age>100|age<70;
drop if age==age[_n-1];

* WOMEN, GOOD HEALTH;

*In order to take off a blue or pink background delete "plotregion(ifcolor(pink*0.08))" in a line below;

genPIshift;
twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, ylabel(0(5000)25000) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) plotregion(ifcolor(pink*0.08)) ti("Medical Expenses, by Permanent Income Percentile, Women in Good Health", size(medsmall)) l1(1998 dollars) ytitle("") saving(medex, replace);
graph use medex.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\medexgw.eps, replace;
graph export c:\research\papers\wealth\medexgw.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;
drop _*;

* WOMEN, BAD HEALTH;

replace ageshiftvarc =ageshiftvarc+healshift;
genPIshift;
replace ageshiftvarc =ageshiftvarc-healshift;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, ylabel(0(5000)25000) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) plotregion(ifcolor(pink*0.08)) ti("Medical Expenses, by Permanent Income Percentile, Women in Bad Health",size(medsmall)) l1(1998 dollars) ytitle("") saving(medex, replace);
graph use medex.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\medexbw.eps, replace;
graph export c:\research\papers\wealth\medexbw.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;
drop _*;

* MEN, GOOD HEALTH;

replace ageshiftvarc =ageshiftvarc+maleshift;
genPIshift;
replace ageshiftvarc =ageshiftvarc-maleshift;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, ylabel(0(5000)25000) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) plotregion(ifcolor(blue*0.08)) ti("Medical Expenses, by Permanent Income Percentile, Men in Good Health", size(medsmall)) l1(1998 dollars) ytitle("") saving(medex, replace);
graph use medex.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\medexgm.eps, replace;
graph export c:\research\papers\wealth\medexgm.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;
drop _*;

* MEN, BAD HEALTH;

replace ageshiftvarc =ageshiftvarc+maleshift+healshift;
genPIshift;
replace ageshiftvarc =ageshiftvarc-maleshift-healshift;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, ylabel(0(5000)25000) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) plotregion(ifcolor(blue*0.08)) ti("Medical Expenses, by Permanent Income Percentile, Men in Bad Health", size(medsmall)) l1(1998 dollars) ytitle("") saving(medex, replace);
graph use medex.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\medexbm.eps, replace;
graph export c:\research\papers\wealth\medexbm.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;


keep age ageshiftFINAL healshift maleshift PIshift bPI2;
order age ageshiftFINAL healshift maleshift PIshift bPI2;

*keep age ageshift healshift bPI bPI2 bmale;
*order age ageshift healshift bPI bPI2 bmale;
sum;

outsheet using c:\research\hrs\wealth\medexprof,nonames replace;
drop _all;
program drop _all;
log close;
