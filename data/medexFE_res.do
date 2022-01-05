# delimit ; 

************************************ FE*****************************;

clear; 
set mem 300m;	
set more 1 ;  
drop _all;
program drop _all;
capture log close;
log using c:\research\hrs\wealth\medex.log, replace ; 

*----------------programs here----------------------------;
*---------------------------------------------------------;
program define fixresid;
			reg res b22-b40 ag11-ag65
ab11-ab66;

predict predres;

fixcoh;
summ predres;
predict meanres; summ predres meanres res;
end;

program define fixcoh;
	version 3.1;
		local i = 22;
		while `i' <=40{;
		replace b`i'=0;
		local i=`i'+1 };

replace b40=1;

fixabag;
end;

*this program generates a polynomial in adist;
program define genPIshift;
        version 3.1;
	replace PI=Pip1;
	makedist;
	gen lmedcostp1=ageshiftvarc+PI*bPI+PI2*bPI2+(PIage*bPIage);
	replace lmedcostp1=lmedcostp1+(vage+ PI*bvPI+ PI2*bvPI2  +PIage*bvPIage )/2;
	replace PI=Pip2;
	makedist;
	gen lmedcostp2=ageshiftvarc+PI*bPI+PI2*bPI2+(PIage*bPIage);
	replace lmedcostp2=lmedcostp2+(vage+ PI*bvPI+ PI2*bvPI2  +PIage*bvPIage )/2;
	replace PI=Pip3;
	makedist;
	gen lmedcostp3=ageshiftvarc+PI*bPI+PI2*bPI2+(PIage*bPIage);
	replace lmedcostp3=lmedcostp3+(vage+ PI*bvPI+ PI2*bvPI2  +PIage*bvPIage )/2;
	replace PI=Pip4;
	makedist;
	gen lmedcostp4=ageshiftvarc+PI*bPI+PI2*bPI2+(PIage*bPIage);
	replace lmedcostp4=lmedcostp4+(vage+ PI*bvPI+ PI2*bvPI2  +PIage*bvPIage )/2;

	gen _20th_percentile=exp(lmedcostp1);
	gen _40th_percentile=exp(lmedcostp2);
	gen _60th_percentile=exp(lmedcostp3);
	gen _80th_percentile=exp(lmedcostp4);
drop lmedcostp*;
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

*------------------------------------------------------------;


use c:\research\hrs\wealth\dataprep2;

gen cohort8=0;
replace cohort8=1 if cohort==.;
replace cohort=8 if cohort8==1;

replace medcost=rmedcost;
replace medcost=. if wave<3;
*drop if wave<3;
*drop if wave>7;
replace medcost=250 if medcost<250;

* log medcost;
gen lmedcost=ln(medcost);
*drop if lmedcost==.;
replace heal=1 if dead==1; * in dataprep2 i do this if died==1.  changing it to "if dead==1" gets at medical expenses of those who died in previous waves;
drop if heal==.; * drops 3 obs;
sum *medc* PI male heal age1;


sort age1;
by age1: sum medcost lmedcost;

gen healage=heal*age1;
*gen PIage=heal*age1;
gen maleage=male*age1;


* drop married women;
sort HHID female realyear;
drop if HHID==HHID[_n-1] & realyear==realyear[_n-1];

sort wave;
by wave: sum age1 dead medcost;

makedist;


* OLS VERSUS FIXED EFFECTS;
* OLS is not consistent with what JOHN, Cristina and I have decided on--it is only descriptive;
if FE==0{;
reg lmedcost age1 age2 age3  heal healage maleage PIage;

reg lmedcost age1 age2 age3  heal healage maleage PIage PI PI2 male;
gen  ageshiftOLS=_b[_cons]+(_b[age1]*age1)+(_b[age2]*age2) +(_b[age3]*age3);
sort age1;
by age1: sum ageshift*;
};

if FE==1{;
* fixed effects ;
*xtreg lmedcost age1 age2 age3 age4 age5, fe i(HHID);
*xtreg lmedcost age1 age2 age3 age4 age5 heal, fe i(HHID);

xtreg lmedcost age1 age2 age3 age4 heal healage maleage PIage , fe i(HHID);

predict lmedcostp;
* get predictions;

gen ageshift=_b[_cons]+(_b[age1]*age1)+(_b[age2]*age2) +(_b[age3]*age3) +(_b[age4]*age4) ;
gen healshift=_b[heal]+_b[healage]*age1;
gen maleshift=_b[maleage]*age1;
gen PIshift=_b[PIage]*age1;
gen bPIage=_b[PIage];

drop lmedcostp;
gen lmedcostp=ageshift+(healshift*heal)+(maleshift*male)+(PIshift*PI);
gen lmedcostres=lmedcost-lmedcostp;



*sum lmedcost*;
* now run an OLS regression on the residuals;
reg lmedcostres PI PI2 male cohort2 cohort3 cohort4 cohort6 cohort7 cohort8;
predict respred;
* set cohort dummies to 0;
replace cohort2=0;
replace cohort3=0;
replace cohort4=0;
replace cohort6=0;
replace cohort7=0;
replace cohort8=0;
predict meanres;
gen lmedcadj=lmedcost+meanres-respred;
gen medcadj=exp(lmedcadj);

save C:\research\hrs\wealth\medresid, replace;

gen b=_b[_cons];
replace ageshift=ageshift+_b[_cons];
gen ageshiftFINAL=ageshift;
gen bPI =_b[PI];
gen bPI2=_b[PI2];
gen bmale=_b[male];
replace maleshift=maleshift+_b[male];
replace PIshift=PIshift+_b[PI]; * do not add in quadratic term, as we keep track of that separately;
* see if cohort effects are important to worry about (as it turns out, they are);


* figure out variance of residual;
gen resres=lmedcostres-respred;
*drop if resres==.;
egen sdres=sd(resres);
gen varc=(sdres*sdres);

* figure out the unconditional variance;
egen sdmedc=sd(lmedcost);
gen varcmed=sdmedc*sdmedc;

* R2 measure;
gen R2meas=1-varc/varcmed;

gen varc2=lmedcostres*lmedcostres;
gen varc1=resres*resres;
reg varc1;
reg varc1 age1 ;
reg varc1 age1 age2 male heal PI;



*xtreg varc1 age1 age2 age3 heal healage PIage, fe i(HHID);
*gen vageFE=_b[_cons]+(_b[age1]*age1)+(_b[age2]*age2) +(_b[age3]*age3);




reg varc1 age1 age2 age3 age4 male maleage heal healage PI PI2 PIage ;
* reg varc1 age1 age2 age3 male maleage heal healage PI PI2 PIage cohort2 cohort3 cohort4  cohort6 cohort7 cohort8;
predict varc1p;


*exit;

gen vage=_b[_cons]+(_b[age1]*age1)+(_b[age2]*age2) +(_b[age3]*age3) +(_b[age4]*age4);





* gen vage=_b[_cons]+(_b[age1]*age1);
gen vmale=_b[male]+(_b[maleage]*age1);
*gen vmale=0;
gen vheal=_b[heal]+(_b[healage]*age1);
gen vPI=_b[PI]+(_b[PIage]*age1);
gen vPI2=_b[PI2];

gen bvPI=_b[PI];
gen bvPIage=_b[PIage];
gen bvPI2=_b[PI2];


* sort age1;
* by age1: sum vage*;
*two conn vage vageFE age1;

sum varc* R2meas;



gen ageshiftvarc=ageshift;

* gen medcostp=exp(lmedcostp); * this term is meaningless, in the aggregate;
};


*rename age1 age;
sort age HHID;

drop if age>102|age<70;
drop if age==age[_n-1];
save c:\research\hrs\wealth\medexprof, replace;
drop if age>100;

***********************************;
* GRAPHS                           ;
***********************************;

set scheme s1mono;

* WOMEN, GOOD HEALTH;

*In order to take off a blue or pink background delete "plotregion(ifcolor(pink*0.08))" in a line below;

genPIshift;
twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, ylabel(0(10000)50000) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Medical Expenses, by Permanent Income Percentile, Women in Good Health", size(medsmall)) l1(1998 dollars) ytitle("") saving(medex, replace);
graph use medex.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\medexgw.eps, replace;
graph export c:\research\papers\wealth\medexgw.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;
drop _*;

* WOMEN, BAD HEALTH;

replace ageshiftvarc =ageshiftvarc+healshift+vheal/2;
genPIshift;
replace ageshiftvarc =ageshiftvarc-healshift-vheal/2;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ylabel(0(10000)50000) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Medical Expenses, by Permanent Income Percentile, Women in Bad Health",size(medsmall)) l1(1998 dollars) ytitle("") saving(medex, replace);
graph use medex.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\medexbw.eps, replace;
graph export c:\research\papers\wealth\medexbw.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;
drop _*;

* MEN, GOOD HEALTH;

replace ageshiftvarc =ageshiftvarc+maleshift+vmale/2;
genPIshift;
replace ageshiftvarc =ageshiftvarc-maleshift-vmale/2;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ylabel(0(10000)50000) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Medical Expenses, by Permanent Income Percentile, Men in Good Health", size(medsmall)) l1(1998 dollars) ytitle("") saving(medex, replace);
graph use medex.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\medexgm.eps, replace;
graph export c:\research\papers\wealth\medexgm.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;
drop _*;

* MEN, BAD HEALTH;

replace ageshiftvarc =ageshiftvarc+maleshift+healshift+vheal/2+vmale/2;
genPIshift;
replace ageshiftvarc =ageshiftvarc-maleshift-healshift-vheal/2-vmale/2;

twoway connected _20th_percentile _40th_percentile _60th_percentile _80th_percentile age, legend(off) ylabel(0(10000)50000) msymbol(T S D O) lpattern(solid longdash shortdash longdash_dot) ti("Medical Expenses, by Permanent Income Percentile, Men in Bad Health", size(medsmall)) l1(1998 dollars) ytitle("") saving(medex, replace);
graph use medex.gph;
graph display, ysize(6) xsize(8.25);
graph export c:\research\papers\wealth\medexbm.eps, replace;
graph export c:\research\papers\wealth\medexbm.emf, replace;
**********numbers behind the graphs****************;
sort age;
by age: sum _20th_percentile _40th_percentile _60th_percentile _80th_percentile age;

use c:\research\hrs\wealth\medexprof, replace;


keep age ageshiftFINAL healshift maleshift PIshift bPI2 vage vheal vmale vPI vPI2;
order age ageshiftFINAL healshift maleshift PIshift bPI2 vage vheal vmale vPI vPI2;

sum;

outsheet using c:\research\hrs\wealth\medexprof,nonames replace;
drop _all;
program drop _all;
log close;
