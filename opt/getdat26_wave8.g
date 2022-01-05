/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

fn logitsr(x) = sqrt(exp(x)./(1+exp(x)));
fn logit(x)   = exp(x)./(1+exp(x));
fn logitrv(p) = ln(p./(1-p));

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
GET1YR:  Converts two-year transition probs into 1-year probs
         Using formula by O. Nartova
         gg = pr(h_t+2=good|h_t=good); bb = pr(h_t+2=bad|h_t=bad)
         _g = pr(h_t+1=good|h_t=good); _b = pr(h_t+1=bad|h_t=bad)
*/
proc(2) = get1yr(gg,bb); 

    local big_A, big_B, big_C, _g, _b;

    big_A = 2-gg-bb;
    big_B = 2*(gg-1);
    big_C = 1-bb-gg+(bb^2);

    _b    = -big_B + sqrt(big_B^2 - 4*big_A.*big_C);
    _b    = _b./(2*big_A);
    _g    = ( (1-bb) - _b.*(1-_b) )./(1-_b);

retp(_g,_b); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
GETAGES:  gets ages right
*/
proc(3) = getages(bornage,dieage,dat);

    local j,k,bage2,ib;

    k = indexcat(dat[.,1],bornage);
    if scalmiss(k) == 1; k = 1; endif;
    bage2 = dat[k,1];
    ib = 1 +maxc((bage2-bornage)|0);
    j = indexcat(dat[.,1],dieage);
    if scalmiss(j) == 1;
        j = rows(dat);
    endif;
retp(k,j,ib); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
FIXOBS:  Adjusts data to have a common set of years
         Adds filler rows if data doesn't cover [bornage,dieage]
*/
proc(1) = fixobs(dat,bornage,dieage,k,j,lowmiss,himiss);

    local cn, rn, rn2, bage2, dage2;

    cn    = cols(dat);
    bage2 = dat[k,1];
    dage2 = dat[j,1];
    dat   = dat[k:j,.];

    if bage2 > bornage;
        rn  = bage2-bornage;
        dat = (lowmiss*ones(rn,cn))|dat;
        dat[1:rn,1] = seqa(bornage,1,rn);
    endif;
    rn2 = j-k+bage2-bornage+1;                      /* Number of rows in dat */
    if dage2 < dieage;
        rn  = dieage-dage2;
        dat = dat|(himiss*ones(rn,cn));
        dat[rn2+1:_tr,1] = seqa(dage2+1,1,rn);
    endif;

retp(dat); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
proc(7) = getprofs();

    local data, datstr, datah, dataw, k, j, ib, mortshg, mortshb, mortswg, 
          mortswb, morts, mort_pi, pctle, pieffcts, survc0, survc, alivec, 
          i, healshg, healshb, healswg, healswb, healdats, heal_pi, 
          ggbb, gg, bb, _g, _b, yprofsh, yprofsw, yprof, y_pi, 
          rorshk, wlthshk, onevec;

    format /ro 10,4;
    onevec = ones(_tr-1,1);
    pctle  = 0.4;

   /*----------------------Mortality and Survival Rates----------------------*/
    if mortdif==0;
        datstr = datapath$+"deathprof_nodiff.out";
    elseif mortdif==1;
        datstr = datapath$+"deathprof.out";
    endif;
    loadm data[] = ^datstr;
    data     = reshape(data,rows(data)/6,6);
    data[.,1]= data[.,1]-1;/* Adjust for Eric's Indexing and 2-year interval */

   /* Order: age, b_const, b_badhealth b_male, b_PI, b_PI2,                  */
   /*        there are separate values of b_* for each age                   */
   /* Pr(alive_t+1|..) = sqrt( logit(b_const[age] + b_badhealh[age] + ...) ) */

    datah    = data[.,1]~(data[.,2]+data[.,4])~data[.,3]~data[.,5:6];
    {k,j,ib} = getages(bornage,dieage-1,datah);
    datah    = fixobs(datah,bornage,dieage,k,j,0,-1e10);
    mort_pi  = datah[.,4:5];
    dataw    = data[.,1]~data[.,2]~data[.,3]~data[.,5:6];
    dataw    = fixobs(dataw,bornage,dieage,k,j,0,-1e10);
    mort_pi  = mort_pi~dataw[.,4:5];

    pieffcts = (mort_pi[.,1]*pctle+mort_pi[.,2]*(pctle^2))*ones(1,2);
    pieffcts = pieffcts~(mort_pi[.,3]*pctle+mort_pi[.,4]*(pctle^2))*ones(1,2);

    "Mortality rates for single individuals at the ";; pctle;; " PI Percentile";
    "husband age      hbad     hgood       wbad      wgood";;
    mortshg  = datah[.,2];
    mortshb  = datah[.,2]+datah[.,3];
    mortswg  = dataw[.,2];
    mortswb  = dataw[.,2]+dataw[.,3];
    morts    = mortshb~mortshg~mortswb~mortswg;
    ageseq~(1-logitsr(morts+pieffcts));?;  

   /*---------------------Health Transition Probabilities---------------------*/
    datstr   = datapath$+"healthprof.out";
    loadm data[] = ^datstr;
    data     = reshape(data,rows(data)/6,6);
    data[.,1]= data[.,1]-1;/* Adjust for Eric's Indexing and 2-year interval */

   /* Order: age, b_const, b_badhealth b_male, b_PI, b_PI2,                  */
   /*        there are separate values of b_* for each age                   */
   /* Pr(h_t+1=bad|..) = logit(b_const[age] + b_badhealh[age] + ...)         */

    datah    = data[.,1]~(data[.,2]+data[.,4])~data[.,3]~data[.,5:6];
    {k,j,ib} = getages(bornage,dieage-1,datah);
    datah    = fixobs(datah,bornage,dieage,k,j,0,0);
    heal_pi  = datah[.,4:5];
    dataw    = data[.,1]~data[.,2]~data[.,3]~data[.,5:6];
    dataw    = fixobs(dataw,bornage,dieage,k,j,0,0);
    heal_pi  = heal_pi~dataw[.,4:5];

    pieffcts = (heal_pi[.,1]*pctle+heal_pi[.,2]*(pctle^2))*ones(1,2);
    pieffcts = pieffcts~(heal_pi[.,3]*pctle+heal_pi[.,4]*(pctle^2))*ones(1,2);

    "Good health probabilities for single individuals at the ";; pctle;; " PI Percentile";
    "husband age      hbad     hgood       wbad      wgood";;
    healshg  = datah[.,2];
    healshb  = datah[.,2]+datah[.,3];
    healswg  = dataw[.,2];
    healswb  = dataw[.,2]+dataw[.,3];
    healdats = healshb~healshg~healswb~healswg;
    ggbb     = logit(healdats+pieffcts);
    gg       = 1 - ggbb[.,2|4];
    bb       = ggbb[.,1|3];
    {_g,_b}  = get1yr(gg,bb);
    _b       = 1 - _b; 
    ageseq~_b[.,1]~_g[.,1]~_b[.,2]~_g[.,2];?;

  /*-----------------------------Income Profiles-----------------------------*/
    datstr   = datapath$+"incprof.out";
    loadm data[] = ^datstr;
    data     = reshape(data,rows(data)/6,6);

   /* Order: age, b_const, b_badhealth b_male, b_PI, b_PI2,                  */
   /*        there are separate values of b_* for each age                   */
   /* y = exp(b_const[age] + b_badhealh[age] + ... )                         */

    datah    = data[.,1]~(data[.,2]+data[.,4])~data[.,3]~data[.,5:6];
    {k,j,ib} = getages(bornage,dieage,datah);
    datah    = fixobs(datah,bornage,dieage,k,j,0,-1e10);
    y_pi     = datah[.,4:5];
    dataw    = data[.,1]~data[.,2]~data[.,3]~data[.,5:6];
    dataw    = fixobs(dataw,bornage,dieage,k,j,0,-1e10);
    y_pi     = y_pi~dataw[.,4:5];

    "Income Profiles at the ";; pctle;; " PI Percentile";
    "husband age     males   females";;
    yprofsh  = datah[.,2];
    yprofsw  = dataw[.,2];
    yprof    = yprofsh~yprofsw;
    pieffcts = y_pi[.,1]*pctle+y_pi[.,2]*(pctle^2);
    pieffcts = pieffcts~(y_pi[.,3]*pctle+y_pi[.,4]*(pctle^2));
    ageseq~exp(yprof+pieffcts);?;  

  /*------------------Get observed date-specific ROR shocks-------------------*/
    datstr   = datapath$+"wlthshk8.txt";

    loadm wlthshk[] = ^datstr;
    wlthshk  = reshape(wlthshk,rows(wlthshk)/5,5);
    {k,j,ib} = getages(momyr1,momyr2+1,wlthshk);
    wlthshk  = fixobs(wlthshk,momyr1,momyr2+1,k,j,0,0);
    rorshk   = wlthshk[.,1+rshktype];

    if rshktype==0;
        rorshk = 0*rorshk;
    endif;
    if rshktype>2;
        rorshk = rorshk - mu_r;
    endif;

    "Rate of return shocks:";;
    (yearseq|(momyr2+1))~rorshk;?; 

retp(morts, healdats, yprof, mort_pi, heal_pi, y_pi, rorshk); 
endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(4) = getmxtab(mxcoef,pctle,smplyrs);

    local datstr, data, k, j , ib, m_const, m_badh, m_male, m_PI, m_PI2, 
          v_const, v_badh, v_male, v_PI, v_PI2, mnmxs, mnmx_pi, stdmxs, 
          varmx_pi, mxmean, mxstd, mx995, title1;

    if smplyrs /= 0;
        smplyrs = smplyrs-bornage+1;
    endif;

 /* First for the means */

    if findmedex==0;

        datstr   = datapath$+"medexprofX.out";
        loadm data[] = ^datstr;
        data     = reshape(data,rows(data)/11,11@6,6@);
        {k,j,ib} = getages(bornage,dieage,data);
        data     = fixobs(data,bornage,dieage,k,j,0,-1e10);

     /* Order: age, m_const, m_badhealth m_male, m_PI, m_PI2,                */
     /*        v_const, v_badhealth v_male, v_PI, v_PI2                      */
     /*        there are separate values of m_* for each age                 */
     /* mx = exp(m_const[age] + m_badhealh[age] + ... )                      */

        m_const = data[.,2];
        m_badh  = data[.,3];
        m_male  = data[.,4];
        m_PI    = data[.,5];
        m_PI2   = data[.,6];
        v_const = data[.,7];
        v_badh  = data[.,8];
        v_male  = data[.,9];
        v_PI    = data[.,10];
        v_PI2   = data[.,11];

    else;

        m_const = mxcoef[1,1] + mxcoef[2,1]*ageseq + mxcoef[3,1]*(ageseq^2/100)
                  + mxcoef[4,1]*(ageseq^3/10000);
        m_badh  = mxcoef[5,1] + mxcoef[6,1]*ageseq;
        m_male  = mxcoef[7,1] + mxcoef[8,1]*ageseq;
        m_PI    = mxcoef[9,1] + mxcoef[10,1]*ageseq;
        m_PI2   = mxcoef[11,1]*ones(_TR,1);
        v_const = mxcoef[1,2] + mxcoef[2,2]*ageseq + mxcoef[3,2]*(ageseq^2/100)
                  + mxcoef[4,2]*(ageseq^3/10000);
        v_badh  = mxcoef[5,2] + mxcoef[6,2]*ageseq;
        v_male  = mxcoef[7,2] + mxcoef[8,2]*ageseq;
        v_PI    = mxcoef[9,2] + mxcoef[10,2]*ageseq;
        v_PI2   = mxcoef[11,2]*ones(_TR,1);

    endif;

    mnmxs    = (m_const[1:_TR]+m_badh[1:_TR])~m_const[1:_TR];
    mnmxs    = (m_const+m_badh+m_male)~(m_const+m_male)~mnmxs;
    mnmx_pi  = m_PI[1:_TR]~m_PI2[1:_TR];
    mnmx_pi  = m_PI~m_PI2~mnmx_pi;
    stdmxs   = (v_const[1:_TR]+v_badh[1:_TR])~v_const[1:_TR];
    stdmxs   = (v_const+v_badh+v_male)~(v_const+v_male)~stdmxs;
    varmx_pi = v_PI[1:_TR]~v_PI2[1:_TR];
    varmx_pi = v_PI~v_PI2~varmx_pi;

    if minc(minc(stdmxs))<0; /* penalize negative variances in GMM criterion */
        negpen = 1;
    else;
        negpen = 0;
    endif;

    mnmxs    = mnmxs + (mxvarscl.==0)*discrtzn; /* Ad Hoc discretization adjustment */
    stdmxs   = abs(stdmxs);
    mnmxs    = mnmxs + (1-mxvarscl)*(stdmxs)/2; /* Mean-preserving adjustment */
    stdmxs   = sqrt(mxvarscl*stdmxs);
    mnmx_pi  = mnmx_pi + (1-mxvarscl)*(varmx_pi)/2;
    mnmx_pi  = mnmx_pi*exp((mxvarscl.==0)*discrtzn);/* Ad Hoc discretization adjustment */
    varmx_pi = mxvarscl*varmx_pi;

    if pctle /= -1;

        {mxmean,mxstd,mx995} = mxprofs(mnmxs,stdmxs,mnmx_pi,varmx_pi,pctle);
    
        title1   = "         age        hbad         hgood         wbad       wgood";
        ?;" Analytical Mean Health Care Costs at the ";; pctle;; " PI Percentile";
        title1;;
        ageseq[smplyrs]~mxmean[smplyrs,.];?;

        " Analytical Std. deviation Health Care Costs at the ";; pctle;; " PI Percentile";
        title1;;
        ageseq[smplyrs]~mxstd[smplyrs,.];?;

        " 99.5-th Percentile Annual Health Care Costs at the ";; pctle;; " PI Percentile";
        title1;;
        ageseq[smplyrs]~mx995[smplyrs,.];?;

    endif;

retp(mnmxs,stdmxs,mnmx_pi,varmx_pi); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(3) = mxprofs(mnlnmxs,stdmxs,mnmx_pi,varmx_pi,pctle);
 
    local pieffcts, mnmxs2, stdmxs2, pct995, mxmean, mxstd, mx995;

    pieffcts = (mnmx_pi[.,1]*pctle+mnmx_pi[.,2]*(pctle^2))*ones(1,2);
    pieffcts = pieffcts~(mnmx_pi[.,3]*pctle+mnmx_pi[.,4]*(pctle^2))*ones(1,2);
    mnmxs2   = mnlnmxs + pieffcts;
    pieffcts = (varmx_pi[.,1]*pctle+varmx_pi[.,2]*(pctle^2))*ones(1,2);
    pieffcts = pieffcts~(varmx_pi[.,3]*pctle+varmx_pi[.,4]*(pctle^2))*ones(1,2);
    stdmxs2  = sqrt(stdmxs^2 + pieffcts);
    
    mxmean   = exp(mnmxs2 + (stdmxs2.^2)/2);
    mxstd    = sqrt( exp(2*mnmxs2 + stdmxs2.^2).*(exp(stdmxs2.^2)-1) );
    pct995   = 2.5758293;
    mx995    = exp(mnmxs2 + pct995*stdmxs2);

retp(mxmean,mxstd,mx995); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* 
GETDATA:  load and clean up data
*/

proc(15) = getdata();

    local datstr, obsdat, wlthdat, totobs, agedat, PIdat, asstdat, MStatdat, 
          beqdat, mxdat, mxobsdat, hsdat, hsobsdat, agedat96, datawgts, 
          chrtcnts, chrttype, avgage96, i, ashft2, toobig, isim94;

    datstr = datapath$+"wlthmat13b.out";
    loadm wlthdat[2688,89] = ^datstr;

    if (ageshft<2);
        ashft2=0;
    elseif (ageshft==2);
        ashft2=1;
    elseif (ageshft==4);
        ashft2=2;
    endif;

    obsdat   = (wlthdat[.,76+ashft2].>dieage) + (wlthdat[.,76+ashft2].<bornage);
@   obsdat   = (obsdat+wlthdat[.,27+ashft2]).>0;  @ @ Drop if initial assets missing:  not necessary here @
    obsdat   = 1-obsdat;
    wlthdat  = selif(wlthdat,obsdat);
    totobs   = rows(wlthdat);"Number of observations";;totobs;?;

    HHIDdat  = wlthdat[.,1];
    MStatdat = wlthdat[.,6+ashft2:12];
    agedat   = wlthdat[.,76+ashft2:82];
    MStatdat = MStatdat.*(agedat.<=dieage); /* Kill if too old for model */
    asstdat  = wlthdat[.,20+ashft2:26];
    toobig   = asstdat.>asstmax;
    asstdat  = asstdat.*(1-toobig) + toobig.*asstmax;
    obsdat   = (wlthdat[.,27+ashft2:33].==0);
    obsdat   = obsdat.*(MStatdat.>0);
    beqdat   = obsdat*0;
/*
    mxdat    = wlthdat[.,48+ashft2:54];
    mxobsdat = (wlthdat[.,55+ashft2:61].==0);
*/
    mxdat    = wlthdat[.,34+ashft2:40];
    mxobsdat = (wlthdat[.,41+ashft2:47].==0);

    hsdat    = wlthdat[.,62+ashft2:68];
    hsobsdat = (wlthdat[.,69+ashft2:75].==0);
    hsobsdat = hsobsdat.*(Mstatdat.>0);

    if wgtddata==0;
        PIdat    = wlthdat[.,5];
        datawgts = ones(totobs,1);
    elseif wgtddata==1;
        PIdat    = wlthdat[.,4];
        datawgts = wlthdat[.,2];
        datawgts = datawgts*totobs/sumc(datawgts);
    endif;

    MStatdat = MStatdat[.,1:mmtyrs];
    asstdat  = asstdat[.,1:mmtyrs];
    agedat   = agedat[.,1:mmtyrs];
    obsdat   = obsdat[.,1:mmtyrs];
    beqdat   = beqdat[.,1:mmtyrs];
    mxdat    = mxdat[.,1:mmtyrs];
    mxobsdat = mxobsdat[.,1:mmtyrs];
    hsdat    = hsdat[.,1:mmtyrs];
    hsobsdat = hsobsdat[.,1:mmtyrs];
    agedat96 = agedat[.,1];

    {chrtcnts,chrttype} = getchrt(agedat[.,1],cohorts);
    avgage96 = zeros(chrtnum,1);
    i=1; do until i>chrtnum;
         avgage96[i] = ((chrttype.==i)'*agedat[.,1])/chrtcnts[i,4];
    i=i+1; endo;
    avgage96 = round(avgage96);
    "Average Cohort Age"; seqa(1,1,chrtnum)~avgage96;
    isim94 = wlthdat[.,76|(5-wgtddata)|6|62|62|69];   /* Use for life exp. calculations */
    isim94 = selif(isim94[.,1:5],isim94[.,6].==0);
    save path =^shkpath isim94;

retp(agedat,PIdat,asstdat,beqdat,MStatdat,obsdat,mxdat,mxobsdat,hsdat,
     hsobsdat,totobs,datawgts,agedat96,avgage96,HHIDdat); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(1) = getninc(assets, labinc);

    local totror, totinc, netinc, taxBrk, taxMar, taxdim, i, tottax, 
          inBrk, blwBrk, abvBrk, inctax, Brk_low, Brk_hi, stuff;

    totror = mu_r + swchROR*rorshk[1];
    totinc = labinc+assets*totror;

    taxBrk = 0|6250|40200|68400|93950|148250|284700|(1e15); 
    taxMar = 0.0765|0.2616|0.4119|0.3499|0.3834|0.4360|0.4761; 
    taxdim = rows(taxMar);

    tottax = 0*totinc;
    blwBrk = tottax;

    i=1; do until i>taxdim;
        Brk_low = taxBrk[i];
        Brk_hi  = taxBrk[i+1];
        blwBrk  = totinc.<Brk_low;
        inBrk   = (1-blwBrk).*(totinc.<=Brk_hi);
        abvBrk  = totinc.>Brk_hi;
        inctax  = inBrk.*(totinc-Brk_low) + abvBrk*(Brk_hi-Brk_low);
        inctax  = inctax*taxMar[i];
        tottax  = tottax+inctax;
    i=i+1; endo; 
   
    netinc = totinc - swchTax*tottax;

/*
    stuff = assets~labinc~totinc~tottax~netinc; stuff[1:200,.]; 
    meanc(stuff); wait;
*/

retp(netinc); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* 
INITDIST:  This gives the initial distribution of assets, wages, health,
           permanent income, and persistent health costs
*/
proc(0) = initdist(isimage);

    local ageindx, gotobs, gotobs2, num0, rn, num, agesim96, pisim96, mstat96,  
          asim96, mxsim96, hsimh96, hsimw96, goth96, gotw96, simwgt96, HHIDsim,
          meanhsh, hstemp, meanhsw, hhdat, hwdat, msdat, hmissdat, ai, asstsign, 
          oldasst, lnabasst, varasst, bta, asstnoiz, projerr, ln_inc, yPI_96, 
          incsim96, cohsim96, asimx, incsim96x, mxsim96x, cohsimx, agesimx, 
          pisimx, simwgtx, HHIDx, astnoizx, hsimhx, hsimwx, mstatx, mstatdx, 
          hhdx, hwdx, hmissdx;

    ageindx  = ones(totobs,1);
    if simtype==2;
        ageindx  = (agedat[.,1].>=(isimage[1]+ageshft)).*(agedat[.,1].<=(isimage[2]+ageshft));
    endif;
    gotobs   = (obsdat[.,1].*mxobsdat[.,1].*hsobsdat[.,1]).*ageindx; 
    num0     = seqa(1,1,totobs);
    num0     = selif(num0,gotobs);  /* restrict ourselves to non-missing obs */
    rn       = rows(num0);

    "Number of observations in initial distribution =";; rn;?;
    num      = rndu(_nn,1)*rn;            /* random draws of indices */
    num      = floor(num)+onesim;
    num      = num0[num];

  /*--------Draw assets, wages and health from empirical distribution--------*/
    agesim96 = agedat[num,1];
    pisim96  = PIdat[num];
    HHIDsim  = HHIDdat[num];
    mstat96  = Mstatdat[num,1];
    gotobs2  = hsobsdat[num,1].*(mstat96.>0);
    goth96   = gotobs2.*((mstat96.==1)+(mstat96.==3));
    gotw96   = gotobs2.*((mstat96.==2)+(mstat96.==3));
    asim96   = asstdat[num,1];
    asim96   = asim96 + (asim96.>asstmax).*(asstmax-asim96);
    mxsim96  = mxdat[num,1];
    hsimh96  = 1-hsdat[num,1];
    hsimw96  = 1-hsdat[num,1];
    simwgt96 = datawgts[num,1];

    meanhsh  = meanc(hsimh96.*goth96)/meanc(goth96);
    hstemp   = rndu(_nn,1).<meanhsh;
    hsimh96  = hsimh96.*goth96 + hstemp.*(1-goth96); @ should be redundant @

    meanhsw  = meanc(hsimw96.*gotw96)/meanc(gotw96);    
    hstemp   = rndu(_nn,1).<meanhsw;
    hsimw96  = hsimw96.*gotw96 + hstemp.*(1-gotw96); @ should be redundant @

    msdat    = mstatdat[num,1:mmtyrs];
    hhdat    = 1-hsdat[num,1:mmtyrs];
    hwdat    = hhdat;
    hmissdat = msdat.<1; /* Identify observations with missing data or MS = 0 */ 
    hmissdat = ((hsobsdat[num,1:mmtyrs].==0)+hmissdat).>0;
    hmissdat = cumsumc(hmissdat')';
    hmissdat = hmissdat.>0;

    "Fraction of HH with men   ";; meanc(goth96);
    "Fraction of HH with women ";; meanc(gotw96);
    "Fraction in good health (men and women), 1996 = ";; meanhsh;; meanhsw;

  /*---Now alter asset and coh distribution to allow for measurement error---*/ 
  /*--------Model:  coh = a + v, v = zero mean error, orthogonal to a--------*/
  /*-----------------a^ = a + u, u orthogonal to a, v------------------------*/
    asstsign = -1+2*(asim96.>=0);
    oldasst  = abs(asim96).<1;                         /* recode zero assets */
    oldasst  = asim96.*(1-oldasst) + oldasst.*asstsign;
    lnabasst = ln(abs(oldasst));
    varasst  = vcx(lnabasst) - assterr;
    bta      = varasst/(varasst+assterr);
    asim96   = lnabasst*bta + (1-bta)*meanc(lnabasst);
    projerr  = sqrt(bta*assterr)*rndn(_nn,1); /* Now add in omitted "true" variation */
    asim96   = asim96 + projerr;
    asim96   = asstsign.*exp(asim96);
    ai       = agesim96-bornage+1;
    yPI_96   = y_PI[ai,1:2].*(mstat96.==1) + y_PI[ai,3:4].*(mstat96.==2);
    ln_inc   = yprof[ai,1].*(mstat96.==1) + yprof[ai,2].*(mstat96.==2)
               + yPI_96[.,1].*pisim96+ yPI_96[.,2].*(pisim96^2);
    incsim96 = exp(ln_inc);
    incsim96 = getninc(asim96,incsim96);
    cohsim96 = asim96 +incsim96 - mxsim96;

    asstnoiz = exp(sqrt(assterr)*rndn(_nn,simyrs));
    asstnoiz[.,1] = oldasst./asim96;
    format /ro 12,4;
    "Mean assets and std dev of measured and 'true' assets, 1996 = ";
    meanc(oldasst)~stdc(oldasst)~meanc(asim96)~stdc(asim96);?;
    "Std. Deviation of Measurement error";; 
    if simtype==1;
        yearseq~stdc(asstnoiz);
    elseif simtype==2;
        ageseqx~stdc(asstnoiz);
    endif;

    if simtype==1;
        save path= ^shkpath asim96, incsim96, mxsim96, cohsim96, mstat96, 
                            pisim96, agesim96, asstnoiz, hsimh96, hsimw96, 
                            msdat, hhdat, hwdat, hmissdat, simwgt96, 
                            HHIDsim;
    elseif simtype==2;
        asimx   = asim96;    incsim96x = incsim96;  mxsim96x = mxsim96;
        cohsimx = cohsim96;  mstatx    = mstat96;   pisimx   = pisim96;
        agesimx = agesim96;  astnoizx  = asstnoiz;  hsimhx   = hsimh96;   
        hsimwx  = hsimw96;   mstatdx   = mstatdat;  hhdx     = hhdat;
        hwdx    = hwdat;     hmissdx   = hmissdat;  simwgtx  = simwgt96; 
        HHIDx   = HHIDsim;
  
        save path= ^shkpath asimx, incsim96x, mxsim96x, cohsimx, mstatx, 
                            pisimx, agesimx, astnoizx, hsimhx, hsimwx, 
                            mstatdx, hhdx, hwdx, hmissdx, simwgtx, HHIDx;
    endif;

retp(); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
proc(3) = gethms(ai2, mstatsim, healsimh, healsimw, pisim96, heal_pi, mort_pi, 
                 morts, healdats, _nn);
                          
    local healshk, too_old, pieffct, healshb, healshg, healswb, healswg, 
          mortshg, mortshb, mortswg, mortswb, bhshb, bhshg, bhswg, bhswb, 
          mstatshk, survshg, survshb, survswg, survswb,hsimpoh, hsimpow, 
          mstpo;

    mortshb = morts[.,1];     mortshg = morts[.,2]; 
    mortswb = morts[.,3];     mortswg = morts[.,4];

    healshb = healdats[.,1];  healshg = healdats[.,2]; 
    healswb = healdats[.,3];  healswg = healdats[.,4];

    healshk = rndu(_nn,2);                           /* for appropriate age */
    too_old = ai2.>_TR;
    ai2     = ai2.*(1-too_old)+too_old*_TR;

/*  All simulation vectors are for a single year                             */
/*  Note:  1 denotes good health, while transition probs are for bad health  */

    pieffct = heal_pi[ai2,1].*pisim96 + heal_pi[ai2,2].*(pisim96^2); /* Men */
    bhshb   = logit(healshb[ai2] + pieffct);
    bhshg   = logit(healshg[ai2] + pieffct);
    {bhshg, bhshb} = get1yr(1-bhshg,bhshb);
    bhshg   = 1 - bhshg;

    pieffct = heal_pi[ai2,3].*pisim96 + heal_pi[ai2,4].*(pisim96^2); /* Women */
    bhswb   = logit(healswb[ai2] + pieffct);
    bhswg   = logit(healswg[ai2] + pieffct);
    {bhswg, bhswb} = get1yr(1-bhswg,bhswb);
    bhswg   = 1 - bhswg;

/*  Note:  Having a > inequality is essential */

    hsimpoh = (healsimh.==0).*bhshb + (healsimh.==1).*bhshg;
    hsimpoh = healshk[.,1].>hsimpoh;

    hsimpow = (healsimw.==0).*bhswb + (healsimw.==1).*bhswg;
    hsimpow = healshk[.,2].>hsimpow;

    mstatshk= rndu(_nn,2);
    pieffct = mort_pi[ai2,1].*pisim96 + mort_pi[ai2,2].*(pisim96^2); /* Men */
    survshg = logitsr( mortshg[ai2] + pieffct );
    survshb = logitsr( mortshb[ai2] + pieffct );

    pieffct = mort_pi[ai2,3].*pisim96 + mort_pi[ai2,4].*(pisim96^2); /* Women */
    survswg = logitsr( mortswg[ai2] + pieffct );
    survswb = logitsr( mortswb[ai2] + pieffct );

    mstpo   = 1*(mstatsim.==1).*             /* husband at time t */
                ( healsimh.*(mstatshk[.,1].<survshg)  +
                 (1-healsimh).*(mstatshk[.,1].<survshb)   ) +
              2*(mstatsim.==2).*                /* wife at time t */
                ( healsimw.*(mstatshk[.,2].<survswg)  +
                 (1-healsimw).*(mstatshk[.,2].<survswb)   ); 
    
    too_old = (ai2+1).>_TR;
    mstpo   = mstpo.*(1-too_old);

retp(hsimpoh, hsimpow, mstpo); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
hsimpute:  Impute health status probabilities for years between AHEAD waves
           Here data are sorted by year, with different ages alive each year
           Note:  if health status probabilities vary by marital status, this
                  imputation procedure will NOT work.
*/
proc(0) = hsimpute(healdats, heal_pi, _nn);

    local healshg, healshb, healswg, healswb, mstat96, pisim96, agesim96, pieffct, 
          ii, ageindx, ai2, too_old, gg_t, g_sh, g_sw, g_t, bb_t,b_sh, b_sw, b_t, 
          B_gg_sh, B_gg_sw, B_gb_sh, B_gb_sw, B_bg_sh, B_bg_sw, B_bb_sh, B_bb_sw,
          gg_sh, gg_sw, gb_sh, gb_sw, bg_sh, bg_sw, bb_sh, bb_sw, g1, g2, g3, g4, cn;

    healshb  = healdats[.,1];  healshg = healdats[.,2]; 
    healswb  = healdats[.,3];  healswg = healdats[.,4];

    load path=^shkpath pisim96, agesim96, mstat96;  

    ageindx  = agesim96-bornage; /* Offset used to get age-appropriate probs */ 
/*
    First Step:  get one- and two-year-ahead transition probs                
        g_ij  = Pr(h_t+1 = good | h_t = good, type ij)
        gg_ij = Pr(h_t+2 = good | h_t = good, type ij)
        b_ij  = Pr(b_t+1 = bad | h_t = bad, type ij)
        bb_ij = Pr(h_t+2 = bad | h_t = bad, type ij)
*/
    g_sh = {};  g_sw = {}; b_sh = {};  b_sw = {};
/*
    Second Step:  Do one-year-imputations for all years 
        Use formulae in handout to get
        B_ef_ij = Pr(h_t+1 = bad | h_t = e, h_t+2 = f, type ij)  
        Add dummy column to convert results to
        B_ef_ij = Pr(h_t = bad | h_t-1 = e, h_t+1 = f, type ij)
*/
    B_gg_sh = -ones(_nn,1);  B_gb_sh = -ones(_nn,1);  B_bb_sh = -ones(_nn,1);  B_bg_sh = -ones(_nn,1);
    B_gg_sw = -ones(_nn,1);  B_gb_sw = -ones(_nn,1);  B_bb_sw = -ones(_nn,1);  B_bg_sw = -ones(_nn,1);

    save path =^shkpath B_gg_sh, B_gb_sh, B_bg_sh, B_bb_sh, B_gg_sw, B_gb_sw, B_bg_sw, B_bb_sw;
    clear B_gg_sh, B_gb_sh, B_bg_sh, B_bb_sh, B_gg_sw, B_gb_sw, B_bg_sw, B_bb_sw;

    ii=1; do until ii > simyrs;
        ai2     = ageindx+ii;          /* Use age-appropriate probabilities */
        too_old = ai2.>_TR;
        ai2     = ai2.*(1-too_old)+too_old*_TR;

    /*  Note:  Transition probs are for bad health */
    
        pieffct = heal_pi[ai2,1].*pisim96 + heal_pi[ai2,2].*(pisim96^2); /* Men */
        bb_t    = logit(healshb[ai2] + pieffct);
        gg_t    = 1 - logit(healshg[ai2] + pieffct);
        {g_t, b_t} = get1yr(gg_t,bb_t);

        if ii<simyrs;
            load path =^shkpath B_gg_sh, B_gb_sh, B_bg_sh, B_bb_sh;
            {g1, g2, g3, g4} = onestep(g_t, g_t, gg_t, b_t, bb_t);
            B_gg_sh = B_gg_sh~(1-g1);
            B_gb_sh = B_gb_sh~(1-g2);
            B_bg_sh = B_bg_sh~(1-g3);
            B_bb_sh = B_bb_sh~(1-g4);
            save path =^shkpath B_gg_sh, B_gb_sh, B_bg_sh, B_bb_sh;
            clear B_gg_sh, B_gb_sh, B_bg_sh, B_bb_sh;
        endif;

        if ii>1;  load path =^shkpath g_sh, b_sh;  endif;
        g_sh     = g_sh~g_t;
        b_sh     = b_sh~b_t;
        save path =^shkpath g_sh, b_sh;  clear g_sh, b_sh;

        pieffct  = heal_pi[ai2,3].*pisim96 + heal_pi[ai2,4].*(pisim96^2); /* Women */

        bb_t     = logit(healswb[ai2] + pieffct);
        gg_t     = 1 - logit(healswg[ai2] + pieffct);
        {g_t, b_t} = get1yr(gg_t,bb_t);
        if ii<simyrs;
            load path =^shkpath B_gg_sw, B_gb_sw, B_bg_sw, B_bb_sw;
            {g1, g2, g3, g4} = onestep(g_t, g_t, gg_t, b_t, bb_t);
            B_gg_sw = B_gg_sw~(1-g1);
            B_gb_sw = B_gb_sw~(1-g2);
            B_bg_sw = B_bg_sw~(1-g3);
            B_bb_sw = B_bb_sw~(1-g4);
            save path =^shkpath B_gg_sw, B_gb_sw, B_bg_sw, B_bb_sw;
            clear B_gg_sw, B_gb_sw, B_bg_sw, B_bb_sw;
        endif;

        if ii>1;  load path =^shkpath g_sw, b_sw;  endif;
        b_sw     = b_sw~b_t;
        g_sw     = g_sw~g_t; 
        save path =^shkpath g_sw, b_sw;  clear g_sw, b_sw;

    ii=ii+1; endo;

retp; endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
proc(4) = onestep(g_t, g_tp1, gg_t, b_t, bb_t);

    local g1, g2, g3, g4;

    g1 = g_t.*g_tp1./gg_t;       /* Pr(h_t+1 = good | h_t = good, h_t+2 = good) */
    g2 = g_t.*(1-g_tp1)./(1-gg_t);/* Pr(h_t+1 = good | h_t = good, h_t+2 = bad) */
    g3 = (1-b_t).*g_tp1./(1-bb_t);/* Pr(h_t+1 = good | h_t = bad, h_t+2 = good) */
    g4 = (1-b_t).*(1-g_tp1)./bb_t; /* Pr(h_t+1 = good | h_t = bad, h_t+2 = bad) */
retp(g1, g2, g3, g4); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
proc(4) = twostep(_g, _b, cn);

    local g_t, g_tp1, g_tp2, ggg_t, b_t, b_tp1, b_tp2, bbb_t, 
          G_gg, G_gb, G_bg, G_bb;

    g_t   = _g[.,cn];
    g_tp1 = _g[.,cn+1];
    g_tp2 = _g[.,cn+2];
    b_t   = _b[.,cn];
    b_tp1 = _b[.,cn+1];
    b_tp2 = _b[.,cn+2];

    ggg_t = g_t.*g_tp1.*g_tp2 + g_t.*(1-g_tp1).*(1-b_tp2) +
            (1-g_t).*(1-b_tp1).*g_tp2 + (1-g_t).*b_tp1.*(1-b_tp2);

    bbb_t = b_t.*b_tp1.*b_tp2 + b_t.*(1-b_tp1).*(1-g_tp2) +
            (1-b_t).*(1-g_tp1).*b_tp2 + (1-b_t).*g_tp1.*(1-g_tp2);

    G_gg  = ( g_t.*g_tp1.*g_tp2 )~( g_t.*(1-g_tp1).*(1-b_tp2) );
    G_gg  = G_gg~( (1-g_t).*(1-b_tp1).*g_tp2 )~( (1-g_t).*b_tp1.*(1-b_tp2) );
    G_gg  = G_gg./ggg_t;

    G_gb  = ( g_t.*g_tp1.*(1-g_tp2) )~( g_t.*(1-g_tp1).*b_tp2 );
    G_gb  = G_gb~( (1-g_t).*(1-b_tp1).*(1-g_tp2) )~( (1-g_t).*b_tp1.*b_tp2 );
    G_gb  = G_gb./(1-ggg_t);

    G_bg  = ( (1-b_t).*g_tp1.*g_tp2 )~( (1-b_t).*(1-g_tp1).*(1-b_tp2) );
    G_bg  = G_bg~( b_t.*(1-b_tp1).*g_tp2 )~( b_t.*b_tp1.*(1-b_tp2) );
    G_bg  = G_bg./(1-bbb_t);

    G_bb  = ( (1-b_t).*g_tp1.*(1-g_tp2) )~( (1-b_t).*(1-g_tp1).*b_tp2 );
    G_bb  = G_bb~( b_t.*(1-b_tp1).*(1-g_tp2) )~( b_t.*b_tp1.*b_tp2 );
    G_bb  = G_bb./bbb_t;
retp(G_gg, G_gb, G_bg, G_bb); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
INITMX:  Generate medical expenditure shocks.  Here, "medex shocks" are  
         uniformly-distributed variables applied to Markov chain probabilities
*/

proc(0) = initmx(_ns);

    local ztacdfsim96, epscdfsim, xicdfsim, ztacdfx, epscdfx, xicdfx;

    ztacdfsim96 = rndu(_ns,1);          /*  Initial draw of AR(1)          */
    epscdfsim   = rndu(_ns,simyrs+1);   /*  simulated innovation on AR(1)  */
    xicdfsim    = rndu(_ns,simyrs+1);   /*  simulated white noise shock    */

    if simtype==1;
        save path=^shkpath ztacdfsim96, xicdfsim, epscdfsim;
    elseif simtype==2;
        ztacdfx = ztacdfsim96;  xicdfx = xicdfsim;  epscdfx = epscdfsim;
        save path=^shkpath ztacdfx, xicdfx, epscdfx;
    endif;           
	 
retp; endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
INITSIM:  Simulate sequences of health, health cost and demographic shocks
          Here data are sorted by year, with different ages alive each year
*/
proc(0) = initsim(morts, mort_pi, healdats, heal_pi);

    local hsimh96, hsimw96, mxsim96, mstat96, pisim96, agesim96, ageindx,
          hsimhx, hsimwx, mstatx, mxsim96x, pisimx, agesimx, ii, ai2,  
          healsimh, healsimw, hsimpoh, hsimpow, mstatsim, mstpo, agesim, 
          hlsimhx, hlsimwx;

    if simtype==1;
        load path=^shkpath hsimh96, hsimw96, mxsim96, mstat96, pisim96, 
                           agesim96;  
    elseif simtype==2;
        load path=^shkpath hsimhx, hsimwx, mxsim96x, mstatx, pisimx, agesimx;           
        hsimh96 = hsimhx;    hsimw96 = hsimwx;  mstat96  = mstatx[.,1];
        mxsim96 = mxsim96x;  pisim96 = pisimx;  agesim96 = agesimx[.,1]; 
        clear hsimhx, hsimwx, mstatx, mxsim96x, pisimx, agesimx;
    endif;    

    mstatsim = mstat96;        /* 0=>defunct, 1=>husband, 2=>wife, 3=>couple */
    agesim   = agesim96;
 
    healsimh = hsimh96;                     /* 0=>bad health, 1=>good health */
    healsimw = hsimw96;

    ageindx  = agesim96-bornage; /* Offset used to get age-appropriate probs */ 

    ii=1; do until ii > simyrs;      /* Update with transition probabilities */
        ai2      = ageindx+ii;                        /* for appropriate age */
        {hsimpoh, hsimpow, mstpo} = gethms(ai2, mstatsim[.,ii], healsimh[.,ii], 
                                           healsimw[.,ii], pisim96, heal_pi, 
                                           mort_pi, morts, healdats,_nn);
        healsimh = healsimh~hsimpoh;
        healsimw = healsimw~hsimpow;
        mstatsim = mstatsim~mstpo;
        agesim   = agesim~(agesim[.,ii]+(mstpo.>0));
    ii=ii+1; endo;
    agesim = agesim.*(mstatsim.>0);

    ?;"Number of simulated observations";; _nn;

    if simtype==1;
        save path=^shkpath agesim, healsimh, healsimw, mstatsim;
    elseif simtype==2;
        agesimx = agesim;    mstatx  = mstatsim;  
        hlsimhx = healsimh;  hlsimwx = healsimw;
        save path=^shkpath agesimx, hlsimhx, hlsimwx, mstatx;
    endif;

    initmx(_nn);                           /* Generate and save medex shocks */

retp; endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
INITSIM2:  Simulate sequences of health, health cost and demographic shocks
           Here data are sorted by year, with different ages alive each year
           This version uses health and mortality shocks observed in the data,
           with missing (no-wave) values imputed using Baye's Rule.
*/
proc(0) = initsim2(morts, mort_pi, healdats, heal_pi);

    local hsimh96, hsimw96, mxsim96, mstat96, pisim96, agesim96, msdat, 
          hhdat, hwdat, hmissdat, healsimh, healsimw, hsimhx, hsimwx, 
          mstatx, mxsim96x, pisimx, agesimx, mstatdx, hhdx, hwdx, hmissdx,
          hsimpoh, hsimpow, mstatsim, mstpo, agesim, ageindx, jj, ii, 
          ii_L, ii_H, clndryr, hmiss, msmiss, ai2, too_old, healshk,
          obs_ggh, obs_gbh, obs_bgh, obs_bbh, obs_ggw, obs_gbw, obs_bgw, obs_bbw, 
          B_gg_sh, B_gg_sw, B_gb_sh, B_gb_sw, B_bg_sh, B_bg_sw, B_bb_sh, B_bb_sw,
          gg_sh, gg_sw, gb_sh, gb_sw, bg_sh, bg_sw, bb_sh, bb_sw,
          iprob, ip2_L, ip2_H, ip3, prtrgt, hlsimhx, hlsimwx;

  /*---Construct health status imputation probabilities using Baye's Rule----*/

    hsimpute(healdats, heal_pi, _nn);

    load path=^shkpath hsimh96, hsimw96, mxsim96, mstat96, pisim96, agesim96, 
                       msdat, hhdat, hwdat, hmissdat; 

    if simtype==1;
        load path=^shkpath hsimh96, hsimw96, mstat96, mxsim96, pisim96, 
                           agesim96, msdat, hhdat, hwdat, hmissdat;
    elseif simtype==2;
        load path=^shkpath hsimhx, hsimwx, mstatx, mxsim96x, pisimx, agesimx, 
                           mstatdx, hhdx, hwdx, hmissdx;   
        hsimh96  = hsimhx;    hsimw96 = hsimwx;  mstat96  = mstatx[.,1];  
        mxsim96  = mxsim96x;  pisim96 = pisimx;  agesim96 = agesimx[.,1];  
        msdat    = mstatdx;   hhdat   = hhdx;    hwdat    = hwdx; 
        hmissdat = hmissdx;
        clear hsimhx, hsimwx, mstatx, mxsim96x, pisimx, agesimx, mstatdx, 
              hhdx, hwdx, hmissdx;
    endif;     

    mstatsim = mstat96;                   /* 0=>defunct, 1=>husband, 2=>wife */
    agesim   = agesim96;
 
    healsimh = hsimh96;                     /* 0=>bad health, 1=>good health */
    healsimw = hsimw96;

    ageindx  = agesim96-bornage; /* Offset used to get age-appropriate probs */ 

    jj = mmtyrs; do while jj>1;
        msmiss = (msdat[.,jj-1].==0).*(msdat[.,jj].>0);
        msdat[.,jj-1] = msdat[.,jj-1].*(1-msmiss) + msmiss.*msdat[.,jj];
    jj=jj-1; endo;

    jj=1; do until jj > mmtyrs;

      /*------------First, do simulations to find missing values-------------*/

        if jj < mmtyrs;
            ii_L = mmtcols[jj];
            ii_H = mmtcols[jj+1]-1;
        else;              /* Tack on an extra year beyond the sample period */
            ii_L = cols(healsimh);
            ii_H = ii_L;
        endif;
 
        ii=ii_L; do until ii > ii_H; /* Update with transition probabilities */
            ai2      = ageindx+ii;                    /* for appropriate age */
            too_old  = ai2.>_TR;
            ai2      = ai2.*(1-too_old)+too_old*_TR;
            {hsimpoh, hsimpow, mstpo} = gethms(ai2, mstatsim[.,ii], healsimh[.,ii], 
                                               healsimw[.,ii], pisim96, heal_pi, 
                                               mort_pi, morts, healdats,_nn);
            healsimh = healsimh~hsimpoh;
            healsimw = healsimw~hsimpow;
            mstatsim = mstatsim~mstpo;
        ii=ii+1; endo;

      /*-------Replace simulations with observed values, when possible--------*/
        if jj >= mmtyrs;
            jj=jj+1; continue;    /* At this point, no more data is available */
        endif;
        ii_H     = ii_H+1;
        hmiss    = hmissdat[.,jj+1];
        healsimh[.,ii_H] = hhdat[.,jj+1].*(1-hmiss) + healsimh[.,ii_H].*hmiss;
        healsimw[.,ii_H] = hwdat[.,jj+1].*(1-hmiss) + healsimw[.,ii_H].*hmiss;
        mstatsim[.,ii_H] = msdat[.,jj+1];

      /*----------Now fill in non-wave years, using imputation probs----------*/
        obs_ggh = (healsimh[.,ii_L].>0.5).*(healsimh[.,ii_H].>0.5);
        obs_gbh = (healsimh[.,ii_L].>0.5).*(healsimh[.,ii_H].<0.5);
        obs_bgh = (healsimh[.,ii_L].<0.5).*(healsimh[.,ii_H].>0.5);
        obs_bbh = (healsimh[.,ii_L].<0.5).*(healsimh[.,ii_H].<0.5);
        obs_ggw = (healsimw[.,ii_L].>0.5).*(healsimw[.,ii_H].>0.5);
        obs_gbw = (healsimw[.,ii_L].>0.5).*(healsimw[.,ii_H].<0.5);
        obs_bgw = (healsimw[.,ii_L].<0.5).*(healsimw[.,ii_H].>0.5);
        obs_bbw = (healsimw[.,ii_L].<0.5).*(healsimw[.,ii_H].<0.5);
        healshk = rndu(_nn,2);

        mstatsim[.,ii_L+1] = mstatsim[.,ii_H];
        load path=^shkpath B_gg_sh, B_gg_sw, B_gb_sh, B_gb_sw,
                           B_bg_sh, B_bg_sw, B_bb_sh, B_bb_sw;

        iprob  =  obs_ggh.*B_gg_sh[.,ii_L+1] + obs_gbh.*B_gb_sh[.,ii_L+1] +
                  obs_bgh.*B_bg_sh[.,ii_L+1] + obs_bbh.*B_bb_sh[.,ii_L+1];
        healsimh[.,ii_L+1] = (healshk[.,1].>iprob);
        iprob  = obs_ggw.*B_gg_sw[.,ii_L+1] + obs_gbw.*B_gb_sw[.,ii_L+1] +
                 obs_bgw.*B_bg_sw[.,ii_L+1] + obs_bbw.*B_bb_sw[.,ii_L+1];
        healsimw[.,ii_L+1] = (healshk[.,2].>iprob);

        clear B_gg_sh, B_gg_sw, B_gb_sh, B_gb_sw, B_bg_sh, B_bg_sw, B_bb_sh, B_bb_sw;

    jj=jj+1; endo;
        
    ii=1; do until ii > (simyrs);
        agesim   = agesim~( agesim[.,ii]+(mstatsim[.,ii+1].>0) );
        too_old  = (agesim[.,ii+1].>dieage);
        agesim[.,ii+1] = agesim[.,ii+1].*(1-too_old) + dieage*too_old;
        mstatsim[.,ii+1] = mstatsim[.,ii+1].*(1-too_old);
    ii=ii+1; endo;

    agesim = agesim.*(mstatsim.>0);

    if simtype==1;
        save path=^shkpath agesim, healsimh, healsimw, mstatsim;
    elseif simtype==2;
        agesimx = agesim;    mstatx  = mstatsim;  
        hlsimhx = healsimh;  hlsimwx = healsimw;
        save path=^shkpath agesimx, hlsimhx, hlsimwx, mstatx;
    endif;

    initmx(_nn);                          /* Generate and save medex shocks */
		 
retp; endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(0) = sumshks();

    local agesim96, healsimh, healsimw, mstatsim, agesim, frh, frw, frc, 
          frdead, yearseq2, hsm, hsw, cn, cp1, cp2, cph, cpw;

    yearseq2 = yearseq|(maxc(yearseq)+1);
    load path= ^shkpath healsimh, healsimw, mstatsim, agesim;

    frdead = (mstatsim.==0);
    frh    = (mstatsim.==1);
    frw    = (mstatsim.==2);
    frc    = (mstatsim.==3);

    "Fraction of Households that are husbands, wives, couples or dead";;
    yearseq2~meanc(frh)~meanc(frw)~meanc(frc)~meanc(frdead);?;

    hsm = (mstatsim.==1)+(mstatsim.==3);/*ones(_nn,1); */
    hsw = (mstatsim.==2)+(mstatsim.==3);/*ones(_nn,1); */

    "Average health status (1=>good) for men and women that are still alive";;
    yearseq2~(meanc(healsimh.*hsm)./meanc(hsm))~(meanc(healsimw.*hsw)./meanc(hsw))~maxc(agesim);?;

    cn  = cols(mstatsim);
    cp1 = healsimh[.,1:cn-1];
    cp2 = healsimh[.,2:cn];
    cph = meanc(cp2.*cp1)./meanc(cp1);
    cp1 = 1-cp1;
    cph = cph~( meanc(cp2.*cp1)./meanc(cp1) );

    cp1 = healsimw[.,1:cn-1];
    cp2 = healsimw[.,2:cn];
    cpw = meanc(cp2.*cp1)./meanc(cp1);
    cp1 = 1-cp1;
    cpw = cpw~( meanc(cp2.*cp1)./meanc(cp1) );

    "Time-t Probability of Good Health at time-t+1";
    "                   hgood         hbad         wgood        wbad";;
    yearseq2[1:cn-1]~cph~cpw;

retp; endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
proc(0) = sumshks2();

    local hlsimhx, hlsimwx, mstatx, frh, frw, frc, frdead, hsm, hsw, ageseq2;

    load path= ^shkpath hlsimhx, hlsimwx, mstatx;

    frdead = (mstatx.==0);
    frh    = (mstatx.==1);
    frw    = (mstatx.==2);
    frc    = (mstatx.==3);
    ageseq2 = ageseqx|(maxc(ageseqx)+1);

    "Fraction of Households that are husbands, wives, couples or dead";;
    ageseq2~meanc(frh)~meanc(frw)~meanc(frc)~meanc(frdead);?;

    hsm = (mstatx.==1)+(mstatx.==3);
    hsw = (mstatx.==2)+(mstatx.==3);

    "Average health status (1=>good) for men and women that are still alive";;
    ageseq2~(meanc(hlsimhx.*hsm)./meanc(hsm))~(meanc(hlsimwx.*hsw)./meanc(hsw));?;

retp; endp;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(3) = compvinv(optwgts);

    local mmtmtx, mnvals, totobs, vcv, vinv, vdiag, rn, _W;

    load path =^datapath mmtmtx;
    mnvals  = 0;
    mnvals  = meanc(mmtmtx);      /* This data should already be zero-mean */
    mmtmtx  = mmtmtx - mnvals';

    totobs  = rows(mmtmtx);
    vcv     = mmtmtx'*mmtmtx;
    clear mmtmtx;
    vcv     = vcv/totobs;

   @-----Take the principal diagonal and form a diagonal weighting matrix-----@

    vinv    = inv(vcv);  
    rn      = rows(vcv);
    vdiag   = diag(vcv);
    save path=^datapath vcv;
    vdiag   = diagrv(eye(rn), vdiag);
    vdiag   = inv(vdiag);

    if optwgts==0;
        _W  = eye(rn);
    elseif optwgts==1;
        _W = vdiag;
    elseif optwgts==2;
        _W = vinv;
    endif;    

    format /ro 12,4;
    "Principal diagonals of V, Vinv (V^{-1}), V_diag^{-1}, and _W = ";;
    seqa(1,1,rn)~diag(vcv)~diag(vinv)~diag(vdiag)~diag(_W);?;

retp(VCV,vinv,_W); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(10) = getWmtx(agedat96,PIdat,asstdat,MStatdat,obsdat,mxdat,mxobsdat,
                   datawgts,optwgts,xtrasst);

    local aqntdat, aqntcnts, mxqntdat, mxqntcnts, mxmnsdat, mxmnscnts, 
          lmxmnsdat, lmxmnscnts, lmxstddat, mxcrldat1, mxcrlcnt1, mxcrldat2, 
          mxcrlcnt2, mmtvec, obsvec, qntvec, pdfvec, mmttype, VCV, vinv, _W,
          Wdiag;

    {aqntdat,aqntcnts} = simqunts(PIdat, agedat96, asstdat, MStatdat, obsdat,
                                  pistate_a, cohorts_a, quants_a, mmtyrs, 0, 
                                  datawgts);

    {mxqntdat,mxqntcnts} = simqunts(PIdat, agedat96, mxdat, MStatdat, mxobsdat,
                                    pistate_m, cohorts_m, quants_m, mmtyrs, 0, 
                                    datawgts);

    {mxmnsdat,mxmnscnts} = simqunts(PIdat, agedat96, mxdat, MStatdat, mxobsdat,
                                    pistate_m, cohorts_m, 0, mmtyrs, 0, datawgts);

    {lmxmnsdat, lmxmnscnts, lmxstddat, mxcrldat1, mxcrlcnt1, mxcrldat2, mxcrlcnt2}
                         = simcrrl(PIdat, agedat96, ln(mxdat+(1-mxobsdat)), MStatdat, 
                                   mxobsdat, pistate_m, cohorts_m, mmtyrs, 0, datawgts);

    ?;"Asset quantiles from Data";; print aqntdat;?; 
    ?;"Medex quantiles from Data";; print mxqntdat;?; 
    ?;"Medex means from Data";; print mxmnsdat;?; print mxmnscnts;?; 

    ?;"Log Medex means from Data ";; print lmxmnsdat;?; 
    ?;"Log Medex standard deviations from Data";; print lmxstddat;?; 
    ?;"Log Medex 1st autocorrelations from Data";; print mxcrldat1;?; 
    ?;"Log Medex 2nd autocorrelations from Data";; print mxcrldat2;?; 

    grphmtx(aqntdat,1,0,0,qnum_a,pinum_a,chrtnum_a);
    grphmtx(mxqntdat,3,0,0,qnum_m,pinum_m,chrtnum_m);
    grphmtx(mxmnsdat,3,0,0,0,pinum_m,chrtnum_m);

    savemmts = 1;                                           /* Modify Global */
    {mmtvec,obsvec,qntvec,pdfvec,mmttype} 
        = makemmts(PIdat,agedat[.,1],asstdat,MStatdat,obsdat,datawgts,aqntdat,
                   mxdat,mxobsdat,mxqntdat,mxmnsdat,lmxmnsdat,lmxstddat,mxcrldat1,
                   mxcrldat2,0);
    savemmts = 0;
    {VCV,vinv,_W} = compvinv(optwgts);

    if optwgts<2;
        Wdiag = diag(_W);
        Wdiag = Wdiag + (xtrasst-1).*(mmttype.==1).*Wdiag;
        _W    = diagrv(_W,Wdiag);
    endif;    

    {aqntdat,aqntcnts} = simqunts(PIdat, agedat96, asstdat, MStatdat, obsdat, 
                                  pistate_a, cohorts_a, quants_a, mmtyrs, 1, 
                                  datawgts); 
    {mxqntdat,mxqntcnts} = simqunts(PIdat, agedat96, mxdat, MStatdat, mxobsdat,
                                    pistate_m, cohorts_m, quants_m, mmtyrs, 1, 
                                    datawgts);
    {mxmnsdat,mxmnscnts} = simqunts(PIdat, agedat96, mxdat, MStatdat, mxobsdat,
                                    pistate_m, cohorts_m, 0, mmtyrs, 1, datawgts);

    grphmtx(aqntdat,1,1,0,qnum_a,pinum_a,chrtnum_a);   /* Adjusts for composition bias */
    grphmtx(mxqntdat,3,1,0,qnum_m,pinum_m,chrtnum_m);
    grphmtx(mxmnsdat,3,1,0,0,pinum_m,chrtnum_m);

retp(VCV,vinv,_W,lmxmnsdat,lmxstddat,mmtvec,obsvec,qntvec,pdfvec,mmttype); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(0) = grphmtx(dataprfs,vartype,comptype,datatype,qnum_j,pinum_j,chrtnum_j);

    local iMStat, MSIndex, gmat, cn, iChrt, iPI, iQunt, tempprf, skipem, 
          _tr2, rn, rn2, gotsome, ageseq2, ageseq3, name1, name2, name3, 
          name4, name5, iYear, jYear, fnamestr;

    if vartype==1;                                 /* consumption vs. assets */
        name1 = "a";
    elseif vartype==2;
        name1 = "c";
    elseif vartype==3;
        name1 = "m";
    endif;    

    if comptype==0;                        /* All observations vs. survivors */
        name2 = "all";
    elseif comptype==1;
        name2 = "srv";
    endif;    

    if datatype==0;                                  /* Data vs. simulations */
        name3 = "dt";        
    elseif datatype==1;
        name3 = "sm";        
    endif;

    if qnum_j==0;                                        /* looking at means */
        iQunt=0;
    else;
        iQunt=1;
    endif;

    _tr2 = maxc(avgage96)+simyrs - bornage;
    _tr2 = maxc(_tr|_tr2);
    ageseq2 = seqa(bornage,1,_tr2);

    do until iQunt>qnum_j;

        name5 = ftos(iQunt,"%*.*lf",1,0);

        iMStat=1; do until iMStat>MSnum;

            if MSsplit==0;                 /* No Marital Status Distinctions */
                MSIndex = 1;
            elseif MSsplit==1;          /* Single M vs. Single F vs. Couples */
                MSIndex = iMStat;
            endif;

            gmat = ones(_tr2,chrtnum_j*pinum_j)*miss(0,0);
            gmat = ageseq2~gmat;
            gotsome = zeros(_tr2,1);        /* Records ages with observations */
            cn   = 1;

            iChrt=1; do until iChrt>chrtnum_j;
                iPI=1; do until iPI>pinum_j;
                    if iQunt==0; /* Means */
                        tempprf = getmatrix(dataprfs,iChrt|iPI|MSIndex|1); /* a row vector */
                    else;
                        tempprf = getmatrix(dataprfs,iChrt|iPI|MSIndex|iQunt); /* a row vector */
                    endif;
                    skipem  = (tempprf.==mvcode);
                    skipem  = skipem*(skipem');
                    cn      = cn+1;
                    rn      = mmtcols+avgage96[iChrt]-bornage;
                    rn2     = mmtyrs - skipem;

                    if rn2>0;
                        rn = rn[1:rn2];
                        rn  = selif(rn,rn.<=_tr2);
                        rn2 = rows(rn);
                        gmat[rn,cn] = tempprf[1:rn2]';
                        gotsome[rn] = ones(rn2,1); 
                    endif;
                iPI=iPI+1; endo;
            iChrt=iChrt+1; endo;

            iYear= 1; do until iYear>rn-1;
                if gotsome[iYear]==1;
                    jYear=_tr2; do until jYear==iYear;
                        if gotsome[jYear]==1; 
                            gotsome[iYear:jYear]=ones(jYear-iYear+1,1); 
                            jYear=iYear+1;
                        endif;
                    jYear=jYear-1; endo;
                iYear=rn-1;
                endif;
            iYear=iYear+1; endo;

            gmat    = selif(gmat,gotsome); /* Drop ages with no observations */
            ageseq3 = selif(ageseq2,gotsome);
            cn = 2; do until cn>cols(gmat);  /* Interpolation for "interior" */
                gmat[.,cn] = fillin(gmat[.,cn],ageseq3);   /* missing values */
            cn=cn+1; endo;

            if MSsplit==0;
               name4 = "";
            elseif MSsplit==1;
                if iMStat==1;
                    name4 = "m";
                elseif IMStat==2;
                    name4 = "f";
                endif;
            endif;
  
            fnamestr = grphpath $+ name1 $+ name2 $+ name3 $+ name4 $+ name5;
            save ^fnamestr = gmat; 
            if (datatype*basecase) == 1; /* Save for comparison graphs */
                fnamestr = grphpath $+ name1 $+ name2 $+ "bn" $+ name4 $+ name5;
                save ^fnamestr = gmat; 
            endif;

        iMStat=iMStat+1; endo;
    iQunt=iQunt+1; endo;

retp; endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
grphmtx2:  Set up files for making graphs.  Here we are doing experiments
           where we track a single cohort through time
*/
proc(0) = grphmtx2(dataprfs,vartype,comptype,datatype,qnum_j,pinum_j);

    local iMStat, MSIndex, tempprf, gmat, cn, iChrt, iPI, iQunt, name1, name2, 
          name3, name4, name5, name6, fnamestr;

    if vartype==1;                                 /* consumption vs. assets */
        name1 = "a";
    elseif vartype==2;
        name1 = "c";
    elseif vartype==3;
        name1 = "m";
    elseif vartype==4;
        name1 = "i";
    endif;  

    if comptype==0;                        /* All observations vs. survivors */
        name2 = "all";
    elseif comptype==1;
        name2 = "srv";
    endif;    

    if datatype==0;                                  /* Data vs. simulations */
        name3 = "dt";        
    elseif datatype==1;
        name3 = "sm";        
    endif;

    name6 = "x";

    if qnum_j==0;                                        /* looking at means */
        iQunt=0;
    else;
        iQunt=1;
    endif;

    do until iQunt>qnum_j;

        name5 = ftos(iQunt,"%*.*lf",1,0)$+name6;

        iMStat=1; do until iMStat>MSnum;

            if MSsplit==0;                 /* No Marital Status Distinctions */
                MSIndex = 1;
            elseif MSsplit==1;          /* Single M vs. Single F vs. Couples */
                MSIndex = iMStat;
            endif;

            gmat = ones(_trexpr,pinum_j)*miss(0,0);
            gmat = ageseqx~gmat;

            iPI=1; do until iPI>pinum_j;
                if iQunt==0; /* Means */
                    tempprf = getmatrix(dataprfs,1|iPI|MSIndex|1); /* a row vector */
                else;
                    tempprf = getmatrix(dataprfs,1|iPI|MSIndex|iQunt); /* a row vector */
                endif;
                gmat[.,iPI+1] = tempprf';
            iPI=iPI+1; endo;

            if MSsplit==0;
               name4 = "";
            elseif MSsplit==1;
                if iMStat==1;
                    name4 = "m";
                elseif IMStat==2;
                    name4 = "f";
                endif;
            endif;
  
            fnamestr = grphpath $+ name1 $+ name2 $+ name3 $+ name4 $+ name5;
            save ^fnamestr = gmat; 
            if (datatype*basecase) == 1; /* Save for comparison graphs */
                fnamestr = grphpath $+ name1 $+ name2 $+ "bn" $+ name4 $+ name5;
                save ^fnamestr = gmat; 
            endif;

        iMStat=iMStat+1; endo;
    iQunt=iQunt+1; endo;

retp; endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
proc(1)=fillin(x,ageseq2);

    local i, j, frac, mv, rn;
                                             
    mv = miss(0,0);
    rn = rows(x);
    i=2; do until i>rn;
        if (x[i]==mv)and(x[i-1]/=mv);
            j=i+1; do until j>rn;
                if (x[j]/=mv); break; endif;
            j=j+1; endo;
            if j<=rn;
                frac = (ageseq2[i]-ageseq2[i-1])/(ageseq2[j]-ageseq2[i-1]);
                x[i] = x[i-1] + frac*(x[j]-x[i-1]);
            endif;
        endif;
    i=i+1; endo;

retp(x); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
MAKENAME.G:  Returns a filename string with quantile index
             Note:  We are assuming quantile index is between 1 and 9.
*/
proc(1) = makename(path,rootname,qunum);

    local quantstr, fnamestr;

    quantstr = ftos(qunum,"%*.*lf",1,0);
    fnamestr = path $+ rootname $+ quantstr;

    retp(fnamestr);
endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
proc(0) = makgrph(quants_j,pinum_j,chrtnum_j,grphtype);

    local fs0, iQunt, fnamestr, quntstr, quntstr2, typestr, typestr2, titlstr, 
          tinum, PIstr, alldat, srvdat, allsim, allboth, allbench, srvsim, cn, 
          rn, z, chrtcol1, chrtcol2, chrtcol3, chrtcol4, plc0, plc1, plt0, 
          plt1, plt2, plt3, plt4, qnum_j, gotdata, cn2;

  /*--------------------------Adjust graph settings--------------------------*/

    fs0 = "-C=1 -CF=" $+ grphpath;
    z   = "%*.*lf";

    if grphtype==1; 
        typestr  = "Assets";
        typestr2 = "a";
    elseif grphtype==2; 
        typestr  = "Consumption";
        typestr2 = "c";
    elseif grphtype==3; 
        typestr  = "Medex";
        typestr2 = "m";
    endif;

    titlstr = ":  ";
    tinum   = 0;
    if pinum_j>1;
        titlstr = "Income"$+titlstr;
        tinum   = 1;
    endif;
    if chrtnum_j>1;
        if tinum==0;
            titlstr = "Cohort"$+titlstr;
        else;
            titlstr = "Cohort and "$+titlstr;
        endif;
        tinum = 1;
    endif;
    if tinum>0;
        titlstr = " by "$+titlstr;
    endif;

    titlstr  = typestr$+titlstr;
    PIstr    = "_"$+ftos(pinum_j,z,1,0)$+"PI_";

    _pdate   = "";
    _pmcolor = zeros(8,1)|15;             /* Black letters, white background */
    fs0      = "-C=1 -CF=" $+ grphpath $+typestr2;
    z         = "%*.*lf";

    chrtcol1 = seqa(2,1,pinum_j)|seqa(2*pinum_j+2,1,pinum_j);
    if (pinum_j<3);
         chrtcol1 = chrtcol1|seqa(4*pinum_j+2,1,pinum_j);
    endif;
    chrtcol2 = seqa(pinum_j+2,1,pinum_j)|seqa(3*pinum_j+2,1,pinum_j);
    chrtcol3 = seqa(2,1,chrtnum_j) @2|4|6@;
    chrtcol4 = seqa(2,1,chrtnum_j-1);

    plc0     = (2|0|4|1|2|0|4|1);

    if pinum_j==1;   
      plc1 = ones(3,1).*.plc0[1:chrtnum_j];        /* Common colors for Cohorts */
   /* plc1 = ones(3,1).*.plc0[1|3|4]; */      /* Common colors for Cohorts */

    else;
        plc1 = ones(chrtnum_j,1).*.plc0[1:pinum_j]; /* Common colors for Quantiles */
    endif;
    
    _pcolor  = plc1;
    plt0     = (6|3|6|3|6|3|6|3);                                   /* line types */     
    plt1     = plt0.*.ones(pinum_j,1);                          
    plt2     = (6|3).*.ones(2*pinum_j,1);
    if pinum_j==1;
        plt2 = (6|3).*.ones(chrtnum_j,1); 
     /* plt2 = (6|3).*.ones(3,1); */
 
    elseif pinum_j<3;  
        plt2 = (6|3).*.ones(3*pinum_j,1); 
    endif;
    plt3     = (6|3).*.ones(2*pinum_j,1);
    plt4     = (6|3).*.ones(pinum_j,1);
    _plwidth = 8*ones(2*chrtnum_j*pinum_j,1);                    /* line widths */
 
    if quants_j==0;                                         /* looking at means */
        qnum_j=0;
        iQunt=0;
    else;
        qnum_j=rows(quants_j);
        iQunt=1;
    endif;

    do until iQunt>qnum_j;

        if iQunt==0;
            quntstr = "Mean ";
            quntstr2 = "avg";        
        elseif quants_j[iQunt]==0.5;
            quntstr = "Median ";
            quntstr2 = ftos(50,z,2,0);
        else;
            quntstr2 = ftos(100*quants_j[iQunt],z,2,0);
            quntstr  = quntstr2$+"th Percentile ";
        endif;

        ylabel(typestr $+ " (000s of 1998 dollars)");
        xlabel("Age");

        if grphtype==1;                                      /* Set up Y axis */ 
            if iQunt==0;
                if pinum_j==1;
                    ytics(0,160,20,miss(0,0));                      
                elseif pinum_j<6;
                    ytics(0,300,50,miss(0,0));                  
                else;
                    ytics(0,500,50,miss(0,0)); 
                endif;
            elseif quants_j[iQunt]==0.5;
                if pinum_j==1; 
                    ytics(0,80,20,miss(0,0));
                elseif pinum_j<6;
                    ytics(0,200,50,miss(0,0));                
                else;
                    ytics(0,250,50,miss(0,0)); 
                endif;
            elseif quants_j[iQunt]==0.9;
                if pinum_j==1;
                    ytics(0,600,100,miss(0,0));
                else;
                    ytics(0,800,200,miss(0,0));
                endif;
            endif;
        elseif grphtype==2;
            if iQunt==0;
                if pinum_j==1;
                    ytics(0,30,5,miss(0,0));
                else;
                    ytics(0,40,5,miss(0,0)); 
                endif;
            elseif quants_j[iQunt]==0.5;
                if pinum_j==1;
                    ytics(0,20,5,miss(0,0));
                else;
                    ytics(0,25,5,miss(0,0));
                 endif;
            elseif quants_j[iQunt]==0.9;
                if pinum_j==1;
                    ytics(0,60,10,miss(0,0));
                else;
                    ytics(0,100,10,miss(0,0));
                endif;
            endif;
        elseif grphtype==3;
            if iQunt==0;
                if pinum_j==1;
                    ytics(0,30,5,miss(0,0)); 
                else;
                    ytics(0,60,5,miss(0,0));
                endif;
            elseif quants_j[iQunt]==0.5;
                if pinum_j==1;
                    ytics(0,20,5,miss(0,0));
                else;
                    ytics(0,50,10,miss(0,0));
                endif;
            elseif quants_j[iQunt]==0.9;
                if pinum_j==1;
                    ytics(0,60,10,miss(0,0));
                else;
                    ytics(0,120,10,miss(0,0));
                endif;
            endif;
        endif;

        if grphtype==2;
            fnamestr = makename(grphpath,typestr2$+"allsm",iQunt);
            load allsim = ^fnamestr;
            _ptitlht = 0.18;
            _pltype  = plt1;
 
            dograph(allsim,quntstr$+titlstr$+"Model",fs0$+"sim"$+quntstr2$+PIstr$+"1.eps");

            _ptitlht = 0.14;
            _pltype  = plt2;

            fnamestr = makename(grphpath,typestr2$+"srvsm",iQunt);
            load srvsim = ^fnamestr;
            rn       = minc(rows(allsim)|rows(srvsim));
            if pinum_j==1;
                cn   = chrtcol3;
            else;
                cn   = chrtcol1;
            endif;
            allboth  = allsim[1:rn,1|cn]~srvsim[1:rn,cn];
            dograph(allboth,quntstr$+typestr$+":  Everyone in Sims (Solid) vs. Survivors (Dashed)",
                    fs0$+"sim"$+quntstr2$+PIstr$+"4.eps");

            iQunt=iQunt+1; 
            continue;
        endif;

        _ptitlht = 0.18;
        _pltype  = plt1;

        fnamestr = makename(grphpath,typestr2$+"alldt",iQunt);
        load alldat = ^fnamestr;

        cn2 = cols(alldat);

        if pinum_j>3;
           cn2 = cn2 - pinum_j; /* drop oldest cohort */
        endif;

        dograph(alldat[.,1:cn2],quntstr$+titlstr$+"Data",fs0$+"data"$+quntstr2$+PIstr$+"1.eps");

        fnamestr = makename(grphpath,typestr2$+"allsm",iQunt);
        load allsim = ^fnamestr;
  
        dograph(allsim,quntstr$+titlstr$+"Model",fs0$+"sim"$+quntstr2$+PIstr$+"1.eps");

        if pinum_j>1;
            _pltype = plt4;
            cn      = chrtcol1;
            dograph(alldat[.,1|cn],quntstr$+titlstr$+"Data",
                    fs0$+"data"$+quntstr2$+PIstr$+"2.eps");
            cn      = chrtcol2;
            dograph(alldat[.,1|cn],quntstr$+titlstr$+"Data",
                    fs0$+"data"$+quntstr2$+PIstr$+"3.eps");
            cn      = chrtcol1;
        else;
            cn = chrtcol3;
            _pltype  = plt2;
        endif;

        _pltype = plt2;
        rn      = minc(rows(alldat)|rows(allsim));
        gotdata = alldat[1:rn,cn]+0.001;
        gotdata = gotdata./gotdata;
        allboth = alldat[1:rn,1|cn]~(allsim[1:rn,cn].*gotdata);
 
        dograph(allboth,quntstr$+typestr$+":  Data (Solid) vs. Model (Dashed)",
                fs0$+"both"$+quntstr2$+PIstr$+"2.eps");

        if pinum_j>1;
            cn      = chrtcol2;
            _pltype = plt3;
            gotdata = alldat[1:rn,cn]+0.001;
            gotdata = gotdata./gotdata;
            allboth = alldat[1:rn,1|cn]~(allsim[1:rn,cn].*gotdata);
            dograph(allboth,quntstr$+typestr$+":  Data (Solid) vs. Model (Dashed)",
                    fs0$+"both"$+quntstr2$+PIstr$+"3.eps");
        endif;

        if basecase==0;
            _ptitlht = 0.16;
            fnamestr = makename(grphpath,typestr2$+"allbn",iQunt);
            load allbench = ^fnamestr;
            if pinum_j==1;
                cn   = chrtcol3;
            else;
                cn   = chrtcol1;
            endif;
            _pltype  = plt1;
            rn       = minc(rows(allsim)|rows(allbench));
            allboth = allsim[1:rn,1|cn]~allbench[1:rn,cn];
            dograph(allboth,quntstr$+typestr$+":  Experiment (Solid) vs. Baseline (Dashed)",
                    fs0$+"comp"$+quntstr2$+PIstr$+"2.eps");
        endif;

        _ptitlht = 0.15;

        fnamestr = makename(grphpath,typestr2$+"srvdt",iQunt);
        load srvdat = ^fnamestr;

        if pinum_j==1;
/*
            cn      = chrtcol4;                              /* drop oldest cohort */
            _pcolor = ones(3,1).*.plc0[1:chrtnum_j-1];       /* Common colors for Cohorts */
            _pltype = (6|3).*.ones(chrtnum_j-1,1); 
*/
            cn      = 2|4;                              /* drop oldest cohort */
            _pcolor = ones(3,1).*.plc0[3|5];       /* Common colors for Cohorts */
            _pltype = (6|3).*.ones(2,1);

        else;
            cn      = chrtcol1;
            _pltype = plt2;
        endif;

        if grphtype==1;                                      /* Set up Y axis */ 
            if quants_j[iQunt]==0.5;
                if pinum_j==1; 
                    ytics(0,100,20,miss(0,0));
                elseif pinum_j<6;
                    ytics(0,250,50,miss(0,0));                
                else;
                    ytics(0,300,50,miss(0,0)); 
                endif;
            endif;
        endif;

        rn       = minc(rows(alldat)|rows(srvdat));

        allboth  = alldat[1:rn,1|cn]~srvdat[1:rn,cn];
        dograph(allboth,quntstr$+typestr$+":  Everyone in Data (Solid) vs. Survivors (Dashed)",
                fs0$+"data"$+quntstr2$+PIstr$+"4.eps");

        fnamestr = makename(grphpath,typestr2$+"srvsm",iQunt);
        load srvsim = ^fnamestr;
        rn       = minc(rows(allsim)|rows(srvsim));
        allboth  = allsim[1:rn,1|cn]~srvsim[1:rn,cn];
        dograph(allboth,quntstr$+typestr$+":  Everyone in Sims (Solid) vs. Survivors (Dashed)",
                fs0$+"sim"$+quntstr2$+PIstr$+"4.eps");

    iQunt=iQunt+1; endo;

retp; endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(0) = makgrphx(quants_j,pinum_j,grphtype);

    local fs0, fs1, iQunt, fnamestr, quntstr, quntstr2, typestr, typestr2, 
          titlstr, PIstr, allsim, srvsim, allbench, allboth, cn0, cn, rn, z, 
          plc0, plc1, plt0, plt1, qnum_j, rn2, pls0, pls1, psc0, psc1;

  /*--------------------------Adjust graph settings--------------------------*/

    fs0 = "-C=1 -CF=" $+ grphpath;
    z   = "%*.*lf";

    if grphtype==1; 
        typestr  = "Assets";
        typestr2 = "a";
        _plegctl = {2,5,6.5,4.1};
        fs1      = ":  Model";
    elseif grphtype==2; 
        typestr  = "Consumption";
        typestr2 = "c";
        _plegctl = {2,5,6.6,4.2};
        fs1      = ":  Model";
    elseif grphtype==3; 
        typestr  = "Medical Expenses";
        typestr2 = "m";
        _plegctl = {2,5,1.2,4};
        fs1      = "";
    elseif grphtype==4; 
        typestr  = "Income";
        typestr2 = "i";
        _plegctl = {2,5,1.1,4.2};
        fs1      = "";
    endif;

    _plegstr = "top \000second \000third \000fourth \000bottom ";

    titlstr = "";
    if pinum_j>1;
        titlstr = " by Income Quintile";
    endif;

    titlstr  = typestr$+titlstr;
    PIstr    = "_"$+ftos(pinum_j,z,1,0)$+"PI_x";

    _pdate   = "";
    _pmcolor = zeros(8,1)|15;             /* Black letters, white background */
    fs0      = "-C=1 -CF=" $+ grphpath $+typestr2;
    z         = "%*.*lf";

    if pinum_j==1;   
        plc0 = 3;
    else;
        plc0 = (2|0|4|1|2|0|4|1);
        plc0 = plc0[1:pinum_j];
    endif;

    _pcolor  = ones(2,1).*.plc0;
    pls0     = (0|0|11|8|10|11|8|10);                           /* symbol types */  
    pls0     = pls0[1:pinum_j]; 
    pls1     = (0|0).*.ones(pinum_j,1);                          
    _pstype  = pls0;
    psc0     = (0|0|1|1|1|1|1|1);                          /* lines with symbol */
    psc0     = psc0[1:pinum_j]; 
    psc1     = (0|0).*.ones(pinum_j,1); 
    _plctrl  = psc0;
    plt0     = (6|3|6|3|6|3|6|3);                                 /* line types */  
    plt0     = plt0[1:pinum_j];   
    plt1     = (6|3).*.ones(pinum_j,1);                          
    _pltype  = plt0;
    _plwidth = 8*ones(2*pinum_j,1);                              /* line widths */
 
    if quants_j==0;                                         /* looking at means */
        qnum_j=0;
        iQunt=0;
    else;
        qnum_j=rows(quants_j);
        iQunt=1;
    endif;

    do until iQunt>qnum_j;

        if iQunt==0;
            quntstr = "Mean ";
            quntstr2 = "avg";        
        elseif quants_j[iQunt]==0.5;
            quntstr = "Median ";
            quntstr2 = ftos(50,z,2,0);
        else;
            quntstr2 = ftos(100*quants_j[iQunt],z,2,0);
            quntstr  = quntstr2$+"th Percentile ";
        endif;

        ylabel(typestr $+ " (000s of 1998 dollars)");
        xlabel("Age");

        if grphtype==1;                                      /* Set up Y axis */ 
            if iQunt==0;
                if pinum_j==1;
                    ytics(0,100,20,miss(0,0));                      
                elseif pinum_j<6;
                    ytics(0,300,50,miss(0,0));                  
                else;
                    ytics(0,500,50,miss(0,0)); 
                endif;
            elseif quants_j[iQunt]==0.5;
                if pinum_j==1; 
                    ytics(0,80,10,miss(0,0));
                elseif pinum_j<6;
                    ytics(0,200,50,miss(0,0));                
                else;
                    ytics(0,250,50,miss(0,0)); 
                endif;
            elseif quants_j[iQunt]==0.9;
                if pinum_j==1;
                    ytics(0,600,100,miss(0,0));
                else;
                    ytics(0,800,200,miss(0,0));
                endif;
            endif;
        elseif grphtype==2;
            if iQunt==0;
                if pinum_j==1;
                    ytics(0,30,5,miss(0,0));
                else;
                    ytics(0,40,5,miss(0,0)); 
                endif;
            elseif quants_j[iQunt]==0.5;
                if pinum_j==1;
                    ytics(0,20,5,miss(0,0));
                else;
                    ytics(0,25,5,miss(0,0));
                 endif;
            elseif quants_j[iQunt]==0.9;
                if pinum_j==1;
                    ytics(0,60,10,miss(0,0));
                else;
                    ytics(0,100,10,miss(0,0));
                endif;
            endif;
        elseif grphtype==3;
            if iQunt==0;
                if pinum_j==1;
                    ytics(0,30,5,miss(0,0)); 
                else;
                    ytics(0,45,5,miss(0,0));
                endif;
            elseif quants_j[iQunt]==0.5;
                if pinum_j==1;
                    ytics(0,20,5,miss(0,0));
                else;
                    ytics(0,50,10,miss(0,0));
                endif;
            elseif quants_j[iQunt]==0.9;
                if pinum_j==1;
                    ytics(0,60,10,miss(0,0));
                else;
                    ytics(0,70,10,miss(0,0));
               endif;
            endif;
        elseif grphtype==4;
            if iQunt==0;
                if pinum_j==1;
                    ytics(0,15,2.5,miss(0,0));
                else;
                    ytics(0,25,5,miss(0,0)); 
                endif;
            elseif quants_j[iQunt]==0.5;
                if pinum_j==1;
                    ytics(0,15,5,miss(0,0));
                else;
                    ytics(0,25,5,miss(0,0));
                 endif;
            elseif quants_j[iQunt]==0.9;
                if pinum_j==1;
                    ytics(0,25,5,miss(0,0));
                else;
                    ytics(0,50,10,miss(0,0));
                endif;
            endif;
        endif;

        fnamestr = makename(grphpath,typestr2$+"allsm",iQunt) $+ "x";
        load allsim = ^fnamestr;
        _ptitlht = 0.18;
        _pstype  = pls0;
        _plctrl  = psc0;
        _pltype  = plt0;
        
        rn2 = 1 + (grphtype==3);
        rn  = rows(allsim);
        cn0 = seqa(2,1,cols(allsim)-1);
        cn0 = rev(cn0);                       @ reverse column order @
        cn  = 1|cn0;

        dograph(allsim[rn2:rn,cn],quntstr$+titlstr$+fs1,fs0$+"sim"$+quntstr2$+PIstr$+"1.eps");

        _ptitlht = 0.14;
        _pltype  = plt1;
        _pstype  = pls1;
        _plctrl  = psc1;
        _plegctl = {0,0,0,0};

        fnamestr = makename(grphpath,typestr2$+"srvsm",iQunt)$+ "x";
        load srvsim = ^fnamestr;
        rn       = minc(rows(allsim)|rows(srvsim));
        allboth  = allsim[1:rn,cn]~srvsim[1:rn,cn0];

        dograph(allboth[rn2:rn,.],quntstr$+typestr$+":  Everyone in Sims (Solid) vs. Survivors (Dashed)",
                fs0$+"sim"$+quntstr2$+PIstr$+"4.eps");

        if (grphtype==1) and (basecase==0);
            _ptitlht = 0.16;
            fnamestr = makename(grphpath,typestr2$+"allbn",iQunt)$+ "x";
            load allbench = ^fnamestr;
            rn       = minc(rows(allsim)|rows(allbench));
            allboth = allsim[1:rn,cn]~allbench[1:rn,cn0];
            dograph(allboth[rn2:rn,.],quntstr$+typestr$+":  Experiment (Solid) vs. Baseline (Dashed)",
                    fs0$+"comp"$+quntstr2$+PIstr$+"2.eps");
        endif;

    iQunt=iQunt+1; endo;

retp; endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(0) = dograph(dat,figtitle,figname);

    local x,y;

    title(figtitle);
    figtitle;;
    dat;?;
    if useMPI==0;
        x  = dat[.,1];
        y  = dat[.,2:cols(dat)]/1000;
        graphprt(figname);
        xy(x,y); 
    endif;
    pause(4);
retp; endp;
