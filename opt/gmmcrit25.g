/*---------------------------------------------------------------------------*/
/*  Procedures for MSM Estimation:  Uses decision rules to find simulated    */
/*      moments and compares to data moments.                                */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(1) = fixcoh(cohsim96);

    local toopoor;

    toopoor  = cohsim96.<cfloor;               /* Bound at consumption floor */
    cohsim96 = cohsim96.*(1-toopoor) + toopoor*cfloor;
    format /ro 12,4;
    ?;"Mean and std dev of cash-on-hand, year 1996 = ";; 
    meanc(cohsim96)~stdc(cohsim96);?;

retp(cohsim96); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(0) = savevecs();

    local agevec, simvec, swchvec, medexvec, prefvec, asstvec, mxpicoef,
          cohsim96, ztacdfsim96, pisim96, agesim96, healsimh, healsimw, 
          mstatsim, xicdfsim, epscdfsim, cohsimx, ztacdfx, pisimx, agesimx, 
          hlsimhx, hlsimwx, mstatx, xicdfx, epscdfx;

    agevec   = bornage|dieage|_tr;
    simvec   = _nn|simyrs|momyr1|momyr2;
    swchvec  = swchMort|swchBeta|swchY|swchmxst|swchTax|swchZeta|swchXi;
    swchvec  = swchvec|swchROR|swchBeq|swchGdif;
    medexvec = rhomx|fracar1|fracar1i|fracwn|minmedex; 
    prefvec  = _delta|_beta|nu|phi0|K0;
    asstvec  = cfloor|tauBeq|exBeq|mu_r|sigma_r; 

    save path = ^iopath agevec;
    save path = ^iopath simvec;
    save path = ^iopath swchvec;
    save path = ^iopath medexvec;
    save path = ^iopath prefvec;
    save path = ^iopath asstvec;
    save path = ^iopath rorshk;

    mortprfs = vecr(mortprfs);  save path = ^iopath mortprfs;
    mort_pi  = vecr(mort_pi);   save path = ^iopath mort_pi;
    hsprobs  = vecr(hsprobs);   save path = ^iopath hsprobs;
    heal_pi  = vecr(heal_pi);   save path = ^iopath heal_pi;
    yprof    = vecr(yprof);     save path = ^iopath yprof;
    y_pi     = vecr(y_pi);      save path = ^iopath y_pi;

    mnlnmxs  = vecr(mnlnmxs);   save path = ^iopath mnlnmxs;
    stdlnmxs = vecr(stdlnmxs);  save path = ^iopath stdlnmxs;
    mxpicoef = vecr(mnmx_pi)|vecr(varmx_pi);  save path = ^iopath mxpicoef;

    if simtype==1;
        load path=^shkpath cohsim96, ztacdfsim96, pisim96, agesim96, 
                           healsimh, healsimw, mstatsim, xicdfsim, 
                           epscdfsim;
    elseif simtype==2;
        load path=^shkpath cohsimx, ztacdfx, pisimx, agesimx, hlsimhx, 
                           hlsimwx, mstatx, xicdfx, epscdfx;
        cohsim96 = cohsimx[.,1];  ztacdfsim96 = ztacdfx;  pisim96   = pisimx;
        agesim96 = agesimx[.,1];  healsimh    = hlsimhx;  healsimw  = hlsimwx;
        mstatsim = mstatx;        xicdfsim    = xicdfx;   epscdfsim = epscdfx;
    endif;

    if allalive==1;
        mstatsim = mstatsim[.,1]*ones(1,cols(mstatsim));
    endif;

    healsimh  = vec(healsimh);         /*  simulated husband's health status */ 
    healsimw  = vec(healsimw);         /*  simulated wife's health status    */ 
    mstatsim  = vec(mstatsim);         /*  simulated marital status          */
    xicdfsim  = vec(xicdfsim);         /*  simulated innovation on AR(1)     */
    epscdfsim = vec(epscdfsim);        /*  simulated white noise shock       */
    cohsim96  = fixcoh(cohsim96);      /*  Bound below by consumption floor  */

    save path=^iopath cohsim96, ztacdfsim96, pisim96, agesim96, healsimh, 
                      healsimw, mstatsim, xicdfsim, epscdfsim;

retp(); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
PUNSCALE_P:  Converts a vector of transformed parameters into the parameters
             used in the life-cycle model
*/
proc(6) = punscale_p(parmvec);

    local bigR;

    _delta = parmvec[1];
    _beta  = parmvec[2];
    nu     = parmvec[3];
    cfloor = exp(parmvec[4]);
    phi0   = parmvec[5];
    phi0   = maxc(minc(phi0|1)|0.000001);
    K0     = parmvec[6]*1000;
    bigR   = 1+mu_r;
    phi0   = (( bigR*(1-phi0)/phi0 )^nu)/(bigR*_beta); 

retp(_delta,_beta,nu,cfloor,phi0,K0); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
PUNSCALE_M:  Converts a vector of transformed parameters into the parameters
             used ARMA model of medex shocks
*/
proc(2) = punscale_m(shkparms);

    local rhomx, fracar1;

    rhomx   = shkparms[1];
    fracar1 = shkparms[2];

    if pscaled_m == 1;
        rhomx   = logit(rhomx);
        fracar1 = logit(fracar1);
    endif;

retp(rhomx, fracar1); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
GETPARMS:   Converts scaled parameter vector into set of parameters
            Note that individual parameters are globals
*/
proc(0) = getparms(parmvec);

    local uparams,rn;

    if pscaled_p==1; 
        {_delta, _beta, nu, cfloor, phi0, K0}
          = punscale_p(parmvec);
    elseif pscaled_p==0;
        _delta = parmvec[1];
        _beta  = parmvec[2];
        nu     = parmvec[3];
        cfloor = parmvec[4];
        phi0   = parmvec[5];
        K0     = parmvec[6];
    endif;

    if swchBeq==0;
        phi0 = 0;  K0 = 1;
    endif;

    if prnres>0;
        uparams = _delta|_beta|nu|cfloor|phi0|K0;
        ?;"Transformed parameters, and parameters actually used:";
        rn=1; do until rn>rows(parmvec);
            plabel[rn];; parmvec[rn];; uparams[rn];
        rn=rn+1; endo;
        ?; "Other key parameters:";
        "taubeq (estate tax rate)     ";; taubeq;
        "exbeq (estate tax exemption) ";; exbeq;
        "mu_r (mean ROR in dec rules) ";; mu_r;
        "sigma_r (std dev ...)        ";; sigma_r; ?;
    endif;

retp; endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(0) = loadsims();

local cohsim, asstsim, ztasim, ztaindexsim, xisim, xiindexsim, Medicaidsim, 
      medexsim, conssim, beqsim, mssim2, alivesim, MStatsim, transfersim, 
      aliveavg, simavg, netIncomesim, medex1yrsim, toosmall, asim96, mxsim96, 
      incsim96, cohsim96, asimx, incsim96x, mxsim96x, mxsimx, incsimx, 
      cohsimx, asstsimx, ztasimx, ztaindx, xisimx, xiindx, Medicaidx, 
      transferx, mx1simx, conssimx, beqsimx, mssim2x, alivex;

    load path=^iopath cohsim, ztasim, ztaindexsim, xisim, xiindexsim, 
               Medicaidsim, medexsim, conssim, beqsim, mssim2, asstsim,
               netIncomesim, transfersim;

    cohsim   = reshape(cohsim,simyrs+1,_nn)';
    asstsim  = reshape(asstsim,simyrs+1,_nn)';
    ztasim   = reshape(ztasim,simyrs+1,_nn)';
    xisim    = reshape(xisim,simyrs+1,_nn)';      
    medexsim = reshape(medexsim,simyrs+1,_nn)';
    ztaindexsim = reshape(ztaindexsim,simyrs+1,_nn)';
    xiindexsim  = reshape(xiindexsim,simyrs+1,_nn)';
    Medicaidsim = reshape(Medicaidsim,simyrs+1,_nn)';
    transfersim = reshape(transfersim,simyrs+1,_nn)';
    medex1yrsim = medexsim;
    conssim  = reshape(conssim,simyrs+1,_nn)';
    beqsim   = reshape(beqsim,simyrs+1,_nn)';
    netIncomesim = reshape(netIncomesim,simyrs+1,_nn)';
    mssim2   = reshape(mssim2,simyrs+1,_nn)';
    alivesim = mssim2.>0;
    aliveavg = meanc(alivesim);

    if simtype==1;
        load path=^shkpath asim96, mxsim96, incsim96, cohsim96;
    elseif simtype==2;
        load path=^shkpath asimx, incsim96x, mxsim96x, cohsimx;
        asim96 = asimx[.,1];  mxsim96 = mxsim96x[.,1];  incsim96 = incsim96x[.,1];
        cohsim96 = cohsimx[.,1];
    endif;

    asstsim[.,1]  = asim96;
    medexsim[.,1] = mxsim96;
    cohsim[.,1]   = cohsim96;
    netIncomesim[.,1] = incsim96;
    medexsim = (medexsim[.,1:simyrs]+medexsim[.,2:simyrs+1])/2; @ 2-year averages @
    medexsim = mxsim96~medexsim;
    toosmall = medexsim.<minmedex;   @ bottom coding @
    medexsim = medexsim.*(1-toosmall) + toosmall*minmedex;

    if prnres>1;
        "     Year        Survival    cash-o-h       zeta+xi      Medicaid    h. costs   consumption    bequests        assets   net income   tot 1yr mx";;
        simavg = aliveavg~(meanc(cohsim.*alivesim)./aliveavg)~(meanc((ztasim+xisim).*alivesim)./aliveavg);
        simavg = simavg~(meanc(Medicaidsim.*alivesim)./aliveavg);
        simavg = simavg~(meanc(medexsim.*alivesim)./aliveavg)~(meanc(conssim.*alivesim)./aliveavg);
        simavg = simavg~(meanc(beqsim.*alivesim)./aliveavg)~(meanc(asstsim.*alivesim)./aliveavg);
        simavg = simavg~(meanc(netIncomesim.*alivesim)./aliveavg)~(meanc((Medicaidsim+medex1yrsim).*alivesim)./aliveavg);
        seqa(momyr1,1,simyrs+1)~simavg;
    endif;

    cohsim   = cohsim[.,1:simyrs];
    asstsim  = asstsim[.,1:simyrs];
    ztasim   = ztasim[.,1:simyrs];
    medexsim = medexsim[.,1:simyrs];
    conssim  = conssim[.,1:simyrs];
    mssim2   = mssim2[.,1:simyrs];
    beqsim   = beqsim[.,1:simyrs];
    alivesim = alivesim[.,1:simyrs];
    netIncomesim = netIncomesim[.,1:simyrs];
    medex1yrsim  = medex1yrsim[.,1:simyrs];
    Medicaidsim  = Medicaidsim[.,1:simyrs];
    transfersim  = transfersim[.,1:simyrs];

    save path=^shkpath cohsim, asstsim, ztasim, ztaindexsim, xisim, xiindexsim, 
                       Medicaidsim, medexsim, medex1yrsim, conssim, beqsim, 
                       mssim2, alivesim, netIncomesim, transfersim;
    if simtype==1;
        save path=^shkpath cohsim, asstsim, ztasim, ztaindexsim, xisim, 
                           xiindexsim, Medicaidsim, medexsim, medex1yrsim, 
                           conssim, beqsim, mssim2, alivesim, netIncomesim;
    elseif simtype==2;
        cohsimx   = cohsim;       asstsimx = asstsim;   ztasimx = ztasim;
        ztaindx   = ztaindexsim;  xisimx   = xisim;     xiindx  = xiindexsim;
        Medicaidx = Medicaidsim;  mxsimx   = medexsim;  mx1simx = medex1yrsim;
        transferx = transfersim;  conssimx = conssim;   beqsimx = beqsim;    
        mssim2x   = mssim2;       alivex   = alivesim;  incsimx = netIncomesim;
        save path=^shkpath cohsimx, asstsimx, ztasimx, ztaindx, xisimx, xiindx,
                           Medicaidx, transferx, mxsimx, mx1simx, conssimx, beqsimx, 
                           mssim2x, alivex, incsimx;
    endif;


retp(); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
PIQuant:  Finds Permanent Income quantiles and Identifies PI interval
*/
proc(3) = PIquant(data,wgts,quants,chktie);

    local qnum, rn, cndmnum, qnts, srtddata, srtd2, k, rn1, cdf, qn, qno, 
          qntype, i, j, qnt, rn2;

    qnum     = rows(quants);
    rn       = rows(data);
    if quants[qnum]==0;
        qnts    = 1~maxc(data);
        cndmnum = 1~rn;
        qntype  = ones(rn,1);
        goto pidone;
    endif;

    qnts     = zeros(qnum,2);
    qnts[.,1]= quants;
    cndmnum  = zeros(qnum+1,2);
    cndmnum[.,1] = (quants|1);

  /*-----------------------Remove Missing Observations-----------------------*/
    data     = missrv(data,mvcode)~wgts~seqa(1,1,rn);
    srtddata = sortc(data,1);
    k        = (srtddata[.,1].==mvcode);
    k        = k'k; 
    if k > 0;
        srtd2    = srtddata[1:k,.];
    elseif k==0;
        srtd2 = {};
    endif;
    srtddata = srtddata[k+1:rn,.];

    rn1      = rn-k;                  /* Number of non-missing observations */
    cdf      = srtddata[.,2];                /* Base distribution on Weights */
    cdf      = cumsumc(cdf)/sumc(cdf);
    qno      = 0;
    qntype   = ones(rn1,1);
        
    i=1; do until i>qnum;
        qn   = cdf.< quants[i]; 
        qn   = qn'qn+1;
        qnt  = srtddata[qn,1];
        qnts[i,2] = qnt;

      /*--------Adjust for ties:  shouldn't happen w/ continuous dist--------*/
        if chktie == 1;
            j  = (srtddata[qn:rn1,1].==qnt);
            j  = j'j -1; 
            qn = qn+j;
        endif;

        if qn==qno;  
            i=i+1;
            continue;
        endif;

        rn2 = qn-qno;                     /* Counts for interval (q_i-1,q_i] */
        cndmnum[i,2] = rn2;
        qntype[qno+1:qn] = i*ones(rn2,1);     /* observations in (q_i-1,q_i] */

        qno = qn;
    i=i+1; endo;                               /* End loop through quantiles */

  /*----------------Count observations above highest quantile----------------*/
    if qn < rn1;
        qntype[qn+1:rn1] = (qnum+1)*ones(rn1-qn,1);
        cndmnum[qnum+1,2] = rn1 - qn;        
    endif;

  /*---------------------Add in Markers for Missing Data---------------------*/

    if k>0;  qntype = zeros(k,1)|qntype;  endif;
    data   = srtd2|srtddata;
    data   = data~qntype;
    data   = sortc(data,3);
    qntype = data[.,4];

pidone:

retp(qnts,cndmnum,qntype); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
PIQuant2:  Finds Permanent Income quantiles and Identifies PI interval
           Assumes data have been normalized to [0,1] interval
*/
proc(3) = PIquant2(data,quants,chktie);

    local qnum, rn, cndmnum, qnts, srtddata, srtd2, k, rn1, cdf, qn, qno, 
          qntype, i, j, qnt, rn2;

    qnum     = rows(quants);
    rn       = rows(data);
    if quants[qnum]==0;
        qnts    = 1~maxc(data);
        cndmnum = 1~rn;
        qntype  = ones(rn,1);
        goto pidone2; 
    endif;

    qnts     = zeros(qnum,2);
    qnts[.,1]= quants;
    cndmnum  = zeros(qnum+1,2);
    cndmnum[.,1] = (quants|1);

  /*-----------------------Remove Missing Observations-----------------------*/
    data     = missrv(data,mvcode)~seqa(1,1,rn);
    srtddata = sortc(data,1);
    k        = (srtddata[.,1].==mvcode);
    k        = k'k; 
    if k > 0;
        srtd2    = srtddata[1:k,.];
    elseif k==0;
        srtd2 = {};
    endif;
    srtddata = srtddata[k+1:rn,.];

    rn1      = rn-k;                  /* Number of non-missing observations */
    cdf      = srtddata[.,1];
    qno      = 0;
    qntype   = ones(rn1,1);
        
    i=1; do until i>qnum;
        qn   = cdf.<= quants[i]; 
        qn   = qn'qn;
        qnt  = srtddata[qn,1];
        qnts[i,2] = qnt;

      /*--------Adjust for ties:  shouldn't happen w/ continuous dist--------*/
        if chktie == 1;
            j  = (srtddata[qn:rn1,1].==qnt);
            j  = j'j -1; 
            qn = qn+j;
        endif;

        if qn==qno;  
            i=i+1;
            continue;
        endif;

        rn2 = qn-qno;                     /* Counts for interval (q_i-1,q_i] */
        cndmnum[i,2] = rn2;
        qntype[qno+1:qn] = i*ones(rn2,1);     /* observations in (q_i-1,q_i] */

        qno = qn;
    i=i+1; endo;                               /* End loop through quantiles */

  /*----------------Count observations above highest quantile----------------*/
    if qn < rn1;
        qntype[qn+1:rn1] = (qnum+1)*ones(rn1-qn,1);
        cndmnum[qnum+1,2] = rn1 - qn;        
    endif;

  /*---------------------Add in Markers for Missing Data---------------------*/

    if k>0;  qntype = zeros(k,1)|qntype;  endif;
    data   = srtd2|srtddata;
    data   = data~qntype;
    data   = sortc(data,2);
    qntype = data[.,3];

pidone2:

retp(qnts,cndmnum,qntype); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
Getquant:  Finds quantiles, after removing missing values
           Can handle weighted data
*/
proc(2) = getquant(data,wgts,quants,chktie);

    local qnum, rn, cndmnum, qnts, srtddata, srtd2, k, rn1, cdf, qn, qno, 
          i, j, qnt, rn2;

    qnum     = rows(quants);
    rn       = rows(data);
    qnts     = zeros(qnum,2);
    qnts[.,1]= quants;
    cndmnum  = zeros(qnum+1,2);
    cndmnum[.,1] = (quants|1);

    data     = packr(data~wgts);              /* Remove Missing Observations */
    srtddata = sortc(data,1);
    rn1      = rows(srtddata);         /* Number of non-missing observations */
    cdf      = srtddata[.,2];                /* Base distribution on Weights */
    cdf      = cumsumc(cdf)/sumc(cdf);
    qno      = 0;
        
    i=1; do until i>qnum;
        qn   = cdf.< quants[i]; 
        qn   = qn'qn+1;
        qnt  = srtddata[qn,1];
        qnts[i,2] = qnt;

      /*--------Adjust for ties:  shouldn't happen w/ continuous dist--------*/
        if chktie == 1;
            j  = (srtddata[qn:rn1,1].==qnt);
            j  = j'j -1; 
            qn = qn+j;
        endif;

        if qn==qno;  
            i=i+1;
            continue;
        endif;

        rn2 = qn-qno;                     /* Counts for interval (q_i-1,q_i] */
        cndmnum[i,2] = rn2;
        qno = qn;
    i=i+1; endo;                               /* End loop through quantiles */

  /*----------------Count observations above highest quantile----------------*/
    if qn < rn1;  cndmnum[qnum+1,2] = rn1 - qn;  endif;

retp(qnts,cndmnum); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
Getqunt2:  Finds quantiles, after removing missing values
           Uses GAUSS procedures
*/
proc(2) = getqunt2(data,quants);

    local qnum, rn, cndmnum, qnts, k, qn, qno, i, qnt, rn2;

    qnum      = rows(quants);
    qnts      = zeros(qnum,2);
    qnts[.,1] = quants;
    cndmnum   = zeros(qnum+1,2);
    cndmnum[.,1] = (quants|1);

    data      = packr(data);                  /* Remove Missing Observations */
    rn        = rows(data);
    if rn>1;
        qnts[.,2] = quantile(data,quants);
    else;
        qnts[.,2] = -1e6;
    endif;
           
    qno      = 0;
    i=1; do until i>qnum;
        qn   = data.<=qnts[i,2]; 
        qn   = qn'qn;
        if qn==qno;  
            i=i+1;
            continue;
        endif;

        rn2 = qn-qno;                     /* Counts for interval (q_i-1,q_i] */
        cndmnum[i,2] = rn2;
        qno = qn;
    i=i+1; endo;                               /* End loop through quantiles */

  /*----------------Count observations above highest quantile----------------*/
    if qn < rn;  cndmnum[qnum+1,2] = rn - qn;  endif;

retp(qnts,cndmnum); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
Getmean:  Finds means, after removing missing values
          Uses GAUSS procedures
*/
proc(2) = getmean(data,wgts);

    local mns, rn, cndmnum;

    mns     = zeros(1,2);
    cndmnum = zeros(1,2);
    data    = packr(data~wgts);              /* Remove Missing Observations */
    if ismiss(data) == 1;
        mns[2] = -1e6;
    else; 
        wgts   = data[.,2];
        wgts   = wgts/sumc(wgts);
        data   = data[.,1].*wgts;                     /* Use weighted data */
        mns[2] = sumc(data);
        cndmnum[2] = rows(data);
    endif;

retp(mns,cndmnum); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(2) = getchrt(data,cohorts_j);

    local rn, chrtnum_j, chrtcnts, srtddata, srtd2, k, rn1, cn, cno, chrttype, 
          h, j, rn2, rn3, agevec;

    rn            = rows(data);
    chrtnum_j     = rows(cohorts_j)-1;
    chrtcnts      = zeros(chrtnum_j,4);
    chrtcnts[.,1] = seqa(1,1,chrtnum_j);
    chrtcnts[.,2] = cohorts_j[1:chrtnum_j]+1;
    chrtcnts[.,3] = cohorts_j[2:chrtnum_j+1];

  /*-----------------------Remove Missing Observations-----------------------*/
    data     = missrv(data,mvcode)~seqa(1,1,rn);
    srtddata = sortc(data,1);
    k        = (srtddata[.,1]==mvcode);
    k        = k'k; 
    if k > 0;
        srtd2    = srtddata[1:k,.];
    elseif k==0;
        srtd2 = {};
    endif;
    srtddata = srtddata[k+1:rn,.];

    rn1      = rn-k;                  /* Number of non-missing observations */
    chrttype = ones(rn1,1);
    cno      = 0;

    agevec   = srtddata[.,1] - ageshft;            /* Use age in first wave */
         
    h=chrtnum_j; do until h<1; 
        cn  = agevec.>cohorts_j[h]; /* Note that cohorts has chrtnum+1 elements */
        cn  = cn'cn;

        if cn==cno;  
            h=h-1;
            continue;
        endif;

        rn2 = cn-cno;                /* Counts for interval (age_h,age_h+1) */
        rn3 = rn1 - cn;
        chrttype[rn3+1:rn3+rn2] = h*ones(rn2,1); /* observations in (age_h,age_h+1) */
        cno = cn;
    h=h-1; endo;                              /* End loop through quantiles */

  /*---------------------------Zero Out Outliers----------------------------*/ 
    chrttype = chrttype.*(1-(agevec.>cohorts_j[chrtnum_j+1]));
    chrttype = chrttype.*(1-(agevec.<=cohorts_j[1]));

    h=1; do until h>chrtnum_j;
        cn = (chrttype.==h);
        chrtcnts[h,4] = cn'cn;
    h=h+1; endo;

  /*--------------------Add in Markers for Missing Data---------------------*/
    if k>0;  chrttype = zeros(k,1)|chrttype;  endif;
    data     = srtd2|srtddata;
    data     = data~chrttype;
    data     = sortc(data,2);
    chrttype = data[.,3];

retp(chrtcnts,chrttype); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(2) = simqunts(pisim96, agesim96, datasim, mssim2, alivesim, pistate_j, 
                   cohorts_j, quants_j, mmtyrs_j, comptype, wgts);

    local pinum_j, PIqnts, PIcnts, PItype, chrtnum_j, qnum_j, chrtcnts, chrttype, 
          indicat0, indicat, iChrt, iMStat, iPI, iQunt, iYear, dataprfs, datacnts, 
          data, tempprf, tempnum, recsize, cn, nummtx;

    pinum_j   = rows(pistate_j) + (pistate_j /= 0);
    chrtnum_j = rows(cohorts_j)-1;
    qnum_j    = rows(quants_j); 

    {PIqnts,PIcnts,PItype} = PIquant2(pisim96,pistate_j,chktie);

    {chrtcnts,chrttype} = getchrt(agesim96,cohorts_j);
    if prnres>1;
        ?;"PI counts";;PIcnts;"Cohort counts";;chrtcnts;?;
    endif;

    cn = cols(alivesim);
    if comptype==1;
        ?;"Looking at Survivors Only!!";
        alivesim = minc(alivesim[.,1|cn]')*ones(1,cn);
    endif;

    nummtx   = zeros(chrtnum_j,pinum_j); 
    recsize  = chrtnum_j|pinum_j|MSnum|qnum_j|mmtyrs_j;
    dataprfs = arrayinit(recsize,mvcode);
    datacnts = arrayinit(recsize,0);

    iChrt=1; do until iChrt>chrtnum_j;
        iPI=1; do until iPI>pinum_j;
            indicat0 = (chrttype.==iChrt).*(PItype.==iPI);
            indicat  = indicat0.*alivesim[.,1];
            nummtx[iChrt,iPI] = rows(indicat)*meanc(indicat);

            iMStat=1; do until iMStat>MSnum;
                if MSsplit==0;             /* No Marital Status Distinctions */
                    indicat = (indicat0*ones(1,mmtyrs_j));
                elseif MSsplit==1;                 /* Single M vs. Single Fs */
                    indicat = (indicat0*ones(1,mmtyrs_j)).*(mssim2.==iMStat);
                endif;

                indicat = indicat.*alivesim;
/*  "Cohort = ";;iCHrt;; "PI Quintile = ";;iPI;; meanc(indicat)';  */
                data    = datasim.*miss(indicat,0);
                iYear=1; do until iYear>mmtyrs_j;
                    if sumc(indicat[.,iYear])<quantmin; /* quantile() won't work */
                        break;
                    endif;
                    if (quants_j==0); /* using means, rather than quantiles */
                        {tempprf,tempnum} = getmean(data[.,iYear],wgts);
                    else;
                        if wgtddata==0;
                            {tempprf,tempnum} = getqunt2(data[.,iYear],quants_j);
                        else;
                            {tempprf,tempnum} = getquant(data[.,iYear],wgts,quants_j,chktie);
                        endif;
                    endif;
                    iQunt=1; do until iQunt>qnum_j;
                        dataprfs[iChrt,iPI,iMStat,iQunt,iYear] = tempprf[iQunt,2];
                        datacnts[iChrt,iPI,iMStat,iQunt,iYear] = tempnum[iQunt,2];
                    iQunt=iQunt+1; endo;
                iYear=iYear+1; endo;
            iMStat=iMStat+1; endo;
        iPI=iPI+1; endo;
    iChrt=iChrt+1; endo;

    if prnres>1;
        "Cohort-PI counts (from first wave in moment criterion)";;nummtx;
    endif;

retp(dataprfs,datacnts); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(7) = simcrrl(pisim96, agesim96, datasim, mssim2, alivesim, pistate_j, 
                  cohorts_j, mmtyrs_j, comptype, wgts);

    local pinum_j, PIqnts, PIcnts, PItype, chrtnum_j, chrtcnts, chrttype, 
          indicat0, indicat, iChrt, iMStat, iPI, iYear, icY0, icY1, icY2,
          meanprfs, meancnts, stdprfs, crrlprf1, crrlcnt1, crrlprf2, crrlcnt2,

          data0, dY0, dY1, dY2, tempprf, tempnum, recsize, cn;

    pinum_j   = rows(pistate_j) + (pistate_j /= 0);
    chrtnum_j = rows(cohorts_j)-1;

    {PIqnts,PIcnts,PItype} = PIquant2(pisim96,pistate_j,chktie);

    {chrtcnts,chrttype} = getchrt(agesim96,cohorts_j);
    if prnres>1;
        ?;"PI counts";;PIcnts;"Cohort counts";;chrtcnts;?;
    endif;

    cn = cols(alivesim);
    if comptype==1;
        ?;"Looking at Survivors Only!!";
        alivesim = minc(alivesim[.,1|cn]')*ones(1,cn);
    endif;

    recsize  = chrtnum_j|pinum_j|MSnum|mmtyrs_j;     @ means @ 
    meanprfs = arrayinit(recsize,mvcode);
    meancnts = arrayinit(recsize,mvcode);
    stdprfs  = arrayinit(recsize,mvcode);
    recsize  = chrtnum_j|pinum_j|MSnum|(mmtyrs_j-1); @ First autocorrelations @
    crrlprf1 = arrayinit(recsize,mvcode);
    crrlcnt1 = arrayinit(recsize,0);
    recsize  = chrtnum_j|pinum_j|MSnum|(mmtyrs_j-2); @ Second autocorrelations @
    crrlprf2 = arrayinit(recsize,mvcode);
    crrlcnt2 = arrayinit(recsize,0);


    iChrt=1; do until iChrt>chrtnum_j;
        iPI=1; do until iPI>pinum_j;
            indicat0 = (chrttype.==iChrt).*(PItype.==iPI);
            indicat  = indicat0.*alivesim[.,1];

            iMStat=1; do until iMStat>MSnum;
                if MSsplit==0;             /* No Marital Status Distinctions */
                    indicat = (indicat0*ones(1,mmtyrs_j));
                elseif MSsplit==1;                 /* Single M vs. Single Fs */
                    indicat = (indicat0*ones(1,mmtyrs_j)).*(mssim2.==iMStat);
                endif;

                indicat = indicat.*alivesim;
            /*  "Cohort = ";;iCHrt;; "PI Quintile = ";;iPI;; meanc(indicat)';  */
                data0   = datasim.*miss(indicat,0);
                dY0     = datasim[.,1]*0;
                dY1     = dY0;
                dY2     = dY0;

                iYear=1; do until iYear>mmtyrs_j;

                    if sumc(indicat[.,iYear])<quantmin; /* quantile() won't work */
                        break;
                    endif;

                    dY0 = data0[.,iYear];
                    {tempprf,tempnum} = getmean(dY0,wgts);
                    meanprfs[iChrt,iPI,iMStat,iYear] = tempprf[1,2];
                    meancnts[iChrt,iPI,iMStat,iYear] = tempnum[1,2];

                    dY0 = dY0 - tempprf[1,2];     /* Convert into residual */
                    {tempprf,tempnum} = getmean(dY0^2,wgts);
                    stdprfs[iChrt,iPI,iMStat,iYear] = sqrt(tempprf[1,2]);
                    dY0 = dY0/sqrt(tempprf[1,2]);    /* Normalize residual */

                    if iYear>1;
                        {tempprf,tempnum} = getmean(dY0.*dY1,wgts);
                        crrlprf1[iChrt,iPI,iMStat,iYear-1] = tempprf[1,2];
                        crrlcnt1[iChrt,iPI,iMStat,iYear-1] = tempnum[1,2];
                    endif;

                    if iYear>2;
                        {tempprf,tempnum} = getmean(dY0.*dY2,wgts);
                        crrlprf2[iChrt,iPI,iMStat,iYear-2] = tempprf[1,2];
                        crrlcnt2[iChrt,iPI,iMStat,iYear-2] = tempnum[1,2];
                    endif;

                    dY2 = dY1;
                    dY1 = dY0;                    
 
                iYear=iYear+1; endo;
            iMStat=iMStat+1; endo;
        iPI=iPI+1; endo;
    iChrt=iChrt+1; endo;

retp(meanprfs, meancnts, stdprfs, crrlprf1, crrlcnt1, crrlprf2, crrlcnt2); 
endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
getpdf:  Finds kernel densities, using non-missing data
         Kernel density estimator written by Ruud Koenig
*/

proc(1) = getpdf(data,indicat,xvals);

    local rn, k, iWave, datac, pdfs, kdns, dkdns, bw;

    indicat = indicat.>0;
    indicat = miss(indicat,0);
    data    = data.*indicat;
    pdfs    = {};
    rn      = rows(data);

    iWave=1; do until iWave>cols(data);
        kdns  = 0*xvals;
        datac = data[.,iWave];
        datac = packr(datac);                 /* Remove Missing Observations */
        if ismiss(datac) /= 1;            
            {kdns, dkdns, bw} = ukernel(xvals[iWave],datac,0,1,&k_gauss); 
        endif;
        pdfs  = pdfs|kdns;
    iWave=iWave+1; endo;

retp(pdfs); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(6) = makemmts_q(PIdat, agedat, data, MStatdat, obsdat, wgts, dataprfs,
                     pistate_j, cohorts_j, quants_j, mmtyrs_j, keepinit, addpdfs);

    local nobs, PIqnts, PIcnts, PItype, pinum_j, chrtcnts, chrttype, chrtnum_j, 
          alivedat, indicat0, indicat, iChrt, iPI, iMStat, iQunt, iYear, cnum,
          qnum_j, tempprf, temppdf, mmtdata, mmtmtx, obsmtx, mmtvec, iM, mmtskip,
          mmtskeep, zn, obsvec, qntvec, pdfvec;

    nobs      = rows(PIdat);
    pinum_j   = rows(pistate_j) + (pistate_j /= 0);
    chrtnum_j = rows(cohorts_j)-1;
    qnum_j    = rows(quants_j);
    if quants_j==0;  addpdfs=0; endif;  /* Looking at means */

    {PIqnts,PIcnts,PItype} = PIquant2(PIdat,pistate_j,chktie);
    {chrtcnts,chrttype}    = getchrt(agedat,cohorts_j);

    mmtmtx   = {};                       /* rows->obs, cols->moment function */
    obsmtx   = {};
    qntvec   = {};
    pdfvec   = {};
    mmtskip  = {};
    mmtskeep = {};
    iM       = 0;

    iChrt=1; do until iChrt>chrtnum_j;

        iPI=1; do until iPI>pinum_j;

            indicat0 = (chrttype.==iChrt).*(PItype.==iPI);

            iMStat=1; do until iMStat>MSnum;
                if MSsplit==0;             /* No Marital Status Distinctions */
                    indicat = indicat0*ones(1,mmtyrs_j);
                elseif MSsplit==1;                  /* Single M vs. Single F */
                    indicat = (indicat0*ones(1,mmtyrs_j)).*(MStatdat.==iMStat);
                endif;
                indicat = indicat.*obsdat;
                indicat = indicat.*wgts;
                iQunt=1; do until iQunt>qnum_j; 
                    tempprf = getmatrix(dataprfs,iChrt|iPI|iMStat|iQunt); /* a row vector */
                    qntvec  = qntvec|(tempprf');
                    if quants_j==0;   /* looking at means */
                        mmtdata = (data - tempprf);
                    else;
                        mmtdata = (data.<=tempprf) - quants_j[iQunt];
                    endif;
                    mmtdata = mmtdata.*indicat;
                    mmtmtx  = mmtmtx~mmtdata;
                    obsmtx  = obsmtx~indicat;
                    temppdf = -tempprf';
                    if addpdfs==1;
                        temppdf = getpdf(data,indicat,tempprf');
                    endif;
                    pdfvec  = pdfvec|temppdf;
                    cnum    = sumc(indicat);
                 /* Identify moments to drop:  These include moments using   */
                 /* the initial asset distribution, or moments with too few  */
                 /* observations */
                    iYear=1; do until iYear>mmtyrs_j; 
                        iM = iM+1;
                        mmtskeep = mmtskeep|iM;
                        if cnum[iYear] < cellmin;
                            mmtskip = mmtskip|iM;
                            cnum[iYear:mmtyrs_j] = zeros(mmtyrs_j-iYear+1,1); /* Skip all future years */
                        elseif (iYear==1) and (keepinit==0);
                            mmtskip = mmtskip|iM;
                        endif;
                    iYear = iYear+1; endo;
                iQunt=iQunt+1; endo;
            iMStat=iMStat+1; endo;
        iPI=iPI+1; endo;
    iChrt=iChrt+1; endo;

    if rows(mmtskip)>0;
        mmtskeep[mmtskip] = zeros(rows(mmtskip),1);
    endif;
    mmtskeep = sortc(mmtskeep,1);
    zn       = sumc(mmtskeep.< 0.99);
    mmtskeep = mmtskeep[zn+1:iM];

    mmtmtx = mmtmtx[.,mmtskeep];     /* Delete moments derived from initial */
    obsmtx = obsmtx[.,mmtskeep];                  /* distribution of assets */
    qntvec = qntvec[mmtskeep];
    pdfvec = pdfvec[mmtskeep];

    mmtvec = meanc(mmtmtx);
    obsvec = meanc(obsmtx);

retp(qntvec,pdfvec,mmtskeep,iM,mmtmtx,obsmtx); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(5) = makemmts_c(PIdat, agedat, data0, MStatdat, obsdat, wgts, meanprfs, 
                     stdprfs, crrlprf1, crrlprf2, pistate_j, cohorts_j, mmtyrs_j);

    local nobs, PIqnts, PIcnts, PItype, pinum_j, chrtcnts, chrttype, chrtnum_j, 
          alivedat, indicat0, indicat, iChrt, iPI, iMStat, iQunt, iYear, cnum,
          meanvec, data1, stdvec, dY0, dY1, dY2, icY0, icY1, icY2, indicati,
          tempprf, mmtdata, mmtmtx, obsmtx, mmtvec, iM, mmtskip,
          mmtskeep, zn, obsvec, crrlvec;

    nobs      = rows(PIdat);
    pinum_j   = rows(pistate_j) + (pistate_j /= 0);
    chrtnum_j = rows(cohorts_j)-1;

    {PIqnts,PIcnts,PItype} = PIquant2(PIdat,pistate_j,chktie);
    {chrtcnts,chrttype}    = getchrt(agedat,cohorts_j);

    mmtmtx   = {};                       /* rows->obs, cols->moment function */
    obsmtx   = {};
    crrlvec  = {};
    mmtskip  = {};
    mmtskeep = {};
    iM       = 0;

    iChrt=1; do until iChrt>chrtnum_j;

        iPI=1; do until iPI>pinum_j;

            indicat0 = (chrttype.==iChrt).*(PItype.==iPI);

            iMStat=1; do until iMStat>MSnum;
                if MSsplit==0;             /* No Marital Status Distinctions */
                    indicat = indicat0*ones(1,mmtyrs_j);
                elseif MSsplit==1;                  /* Single M vs. Single F */
                    indicat = (indicat0*ones(1,mmtyrs_j)).*(MStatdat.==iMStat);
                endif;
                indicat = indicat.*obsdat;
                indicat = indicat.*sqrt(wgts);      @ terms are multiplied below @

             /* meanvec and stdvec should be from data as the model is matching  */
             /* level means, not log means or deviations.                        */

                meanvec = getmatrix(meanprfs,iChrt|iPI|iMStat);  /* a row vector */
                data1   = data0 - meanvec;             /* Convert into residuals */               
                stdvec  = getmatrix(stdprfs,iChrt|iPI|iMStat);  
                data1   = data1./stdvec;                    /* Normalize residual */
                dY0     = data1[.,1]*0;
                dY1     = dY0;
                dY2     = dY0;
                icY0    = dY0;
                icY1    = dY0;
                icY2    = dY0;

                iYear=1; do until iYear>mmtyrs_j;
                    dY0  = data1[.,iYear];
                    icY0 = indicat[.,iYear];

                    if iYear>2; @ We lack 2-year averages for year 1 @
                        tempprf  = getscalar4d(crrlprf1,iChrt,iPI,iMStat,iYear-1);
                        crrlvec  = crrlvec|tempprf;
                        mmtdata  = dY0.*dY1 - tempprf;
                        indicati = icY0.*icY1;
                        mmtdata  = mmtdata.*indicati;
                        mmtmtx   = mmtmtx~mmtdata;
                        obsmtx   = obsmtx~indicati;
                        cnum     = sumc(indicati.>0);
                        iM       = iM+1;
                        mmtskeep = mmtskeep|iM;
                        if cnum < cellmin; /* Drop moments with too few observations */
                            mmtskip = mmtskip|iM;
                        endif;
                    endif;

                    if iYear>3;
                        tempprf  = getscalar4d(crrlprf2,iChrt,iPI,iMStat,iYear-2);
                        crrlvec  = crrlvec|tempprf;
                        mmtdata  = dY0.*dY2 - tempprf;
                        indicati = icY0.*icY2;
                        mmtdata  = mmtdata.*indicati;
                        mmtmtx   = mmtmtx~mmtdata;
                        obsmtx   = obsmtx~indicati;
                        cnum     = sumc(indicati.>0);
                        iM       = iM+1;
                        mmtskeep = mmtskeep|iM;
                        if cnum < cellmin;
                            mmtskip = mmtskip|iM;
                        endif;
                    endif;

                    dY2  = dY1;
                    dY1  = dY0; 
                    icY2 = icY1;
                    icY1 = icY0;                   
 
                iYear=iYear+1; endo;
            iMStat=iMStat+1; endo;
        iPI=iPI+1; endo;
    iChrt=iChrt+1; endo;

    if rows(mmtskip)>0;
        mmtskeep[mmtskip] = zeros(rows(mmtskip),1);
    endif;
    mmtskeep = sortc(mmtskeep,1);
    zn       = sumc(mmtskeep.< 0.99);
    mmtskeep = mmtskeep[zn+1:iM];

    mmtmtx  = mmtmtx[.,mmtskeep];  
    obsmtx  = obsmtx[.,mmtskeep];  
    crrlvec = crrlvec[mmtskeep];

    mmtvec  = meanc(mmtmtx);
    obsvec  = meanc(obsmtx);

retp(crrlvec,mmtskeep,iM,mmtmtx,obsmtx); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(5) = makemmts(PIdat, agedat, asstdat, MStatdat, obsdat, wgts, asstqnts,
                   mxdat, mxobsdat, mxquants, mxmeans, lmxmnsdat, lmxstddat, 
                   mxcrrls1, mxcrrls2, addpdfs);

    local mmtvec_a, obsvec_a, qntvec_a, pdfvec_a, mmtskeep_a, iM_a, mmtmtx_a, obsmtx_a,
          mmtvec_m, obsvec_m, qntvec_m, pdfvec_m, mmtskeep_m, iM_m, mmtmtx_m, obsmtx_m,
          meanvec_m, crrlvec_m, mmtvec, obsvec, qntvec, pdfvec, mmtskeep, iM, mmtmtx, 
          obsmtx, mmttype, nmom, MSMval, i;


/*  First, evaluate quantiles for assets */

    {qntvec_a,pdfvec_a,mmtskeep_a,iM_a,mmtmtx_a,obsmtx_a}
        = makemmts_q(PIdat,agedat[.,1],asstdat,MStatdat,obsdat,datawgts,
                     asstqnts,pistate_a,cohorts_a,quants_a,mmtyrs,0,addpdfs);

    mmtmtx     = mmtmtx_a;
    obsmtx     = obsmtx_a;
    qntvec     = qntvec_a;
    pdfvec     = pdfvec_a;
    mmtskeep   = mmtskeep_a;
    iM         = iM_a;
    mmttype    = ones(iM_a,1);

    if medexmmts==1;

    /*  Next, evaluate quantiles for Medex */

        {qntvec_m,pdfvec_m,mmtskeep_m,iM_m,mmtmtx_m,obsmtx_m}
            = makemmts_q(PIdat,agedat[.,1],mxdat,MStatdat,mxobsdat,datawgts,
                         mxquants,pistate_m,cohorts_m,quants_m,mmtyrs,0,addpdfs);

        mmtmtx     = mmtmtx~mmtmtx_m;
        obsmtx     = obsmtx~obsmtx_m;
        qntvec     = qntvec|qntvec_m;
        pdfvec     = pdfvec|pdfvec_m;
        mmtskeep_m = mmtskeep_m + iM_a;
        mmtskeep   = mmtskeep|mmtskeep_m;
        iM         = iM+iM_m;
        mmttype    = mmttype|(2*ones(iM_m,1));

    /*  Next, evaluate Medex means */

        {meanvec_m,pdfvec_m,mmtskeep_m,iM_m,mmtmtx_m,obsmtx_m}
            = makemmts_q(PIdat,agedat[.,1],mxdat,MStatdat,mxobsdat,datawgts,
                         mxmeans,pistate_m,cohorts_m,0,mmtyrs,0,0);

        mmtmtx     = mmtmtx~mmtmtx_m;
        obsmtx     = obsmtx~obsmtx_m;
        qntvec     = qntvec|meanvec_m;
        mmtskeep_m = mmtskeep_m + iM;
        mmtskeep   = mmtskeep|mmtskeep_m;
        iM         = iM+iM_m;
        mmttype    = mmttype|(3*ones(iM_m,1));

    /*  Next, look at log Medex autocorrelations */

        {crrlvec_m,mmtskeep_m,iM_m,mmtmtx_m,obsmtx_m}
            = makemmts_c(PIdat,agedat[.,1],ln(mxdat+(1-mxobsdat)),MStatdat,
                         mxobsdat,datawgts,lmxmnsdat,lmxstddat,mxcrrls1,mxcrrls2,
                         pistate_m,cohorts_m,mmtyrs);

        mmtmtx     = mmtmtx~mmtmtx_m;
        obsmtx     = obsmtx~obsmtx_m;
        qntvec     = qntvec|crrlvec_m;
        mmtskeep_m = mmtskeep_m + iM;
        mmtskeep   = mmtskeep|mmtskeep_m;
        mmttype    = mmttype|(4*ones(iM_m,1));

    endif;
    
    mmtvec     = meanc(mmtmtx);
    obsvec     = meanc(obsmtx);
    nmom       = rows(mmtvec);
    mmttype    = mmttype[mmtskeep];
    MSMval     = (totobs*(mmtvec^2).*diag(_W));


    if savemmts==1;                                /* N.B. This is a global */
        save path=^datapath mmtmtx, obsmtx;
    endif;
    if prnres>1;
        ?;"Check of Moment Conditions";
        "                               Average   % Observed     Adjusted    MSM Value         Type";;
        seqa(1,1,nmom)~mmtskeep~mmtvec~obsvec~(mmtvec./obsvec)~MSMval~mmttype;?;
        "Summed squared avg. differences = ";; mmtvec'mmtvec;?;
    endif;
    if prnres>0;
        i=1; do until i>maxc(mmttype);
            "Contribution from moment category";; i;; 
            if i==1;
                "Asset Quantiles       ";;
            elseif i==2;
                "Medex Quantiles       ";;
            elseif i==3;
                "Medex Means           ";;
            elseif i==4;
                "Medex Autocorrelations";;
            endif;
            sumc((mmttype.==i).*MSMval);
        i=i+1; endo;
        ?;
    endif;
    
retp(mmtvec,obsvec,qntvec,pdfvec,mmttype); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(1) = getcrit(parmvec);

    local PIsim96, agesim96, simwgt96, asim96, asstsim, conssim, medexsim, 
          mssim2, alivesim, asstnoiz,  execret, aqntsim, aqntcnts, cprfsim, 
          cprfcnts, mxprfsim, mxprfcnts, mxmnssim, mxmnscnts, lmedexsim, 
          lmxmnssim, lmxmnscnts, mxstdsim, mxcrlsim1, mxcrlcnt1, mxcrlsim2, 
          mxcrlcnt2, mmtvec, obsvec, qntvec, pdfvec, mmttype, criter;

    getparms(parmvec);
    savevecs();                           /* Save input vector for C program */
    if prnout==1; output off; endif;

 /*   if useMPI < 2;
        execret = exec(rulecall," ");
    elseif useMPI == 2;
        execret = exec("mpirun", rulecall);
    endif;
        */

    if prnout==1; output on; endif;
    loadsims();              /* Load simulation results into local directory */

 /* Get model-predicted summary stats */
    load path=^shkpath pisim96, agesim96, simwgt96, asim96, asstsim, mssim2, 
                       conssim, alivesim, medexsim, asstnoiz;
    asstsim  = asstsim.*asstnoiz;
    asstsim  = asstsim[.,mmtcols];                   /* AHEAD Survey Years */
    conssim  = conssim[.,mmtcols];
    mssim2   = mssim2[.,mmtcols];
    alivesim = alivesim[.,mmtcols];
    medexsim = medexsim[.,mmtcols];

    {aqntsim,aqntcnts}   = simqunts(pisim96,agesim96,asstsim,mssim2,alivesim, 
                                    pistate_a,cohorts_a,quants_a,mmtyrs,0,simwgt96);

    {cprfsim,cprfcnts}   = simqunts(pisim96,agesim96,conssim,mssim2,alivesim, 
                                    pistate_a,cohorts_a,quants_a,mmtyrs,0,simwgt96);

    {mxprfsim,mxprfcnts} = simqunts(pisim96,agesim96,medexsim,mssim2,alivesim, 
                                    pistate_m,cohorts_m,quants_m,mmtyrs,0,simwgt96);

    {mxmnssim,mxmnscnts} = simqunts(pisim96,agesim96,medexsim,mssim2,alivesim, 
                                    pistate_m,cohorts_m,0,mmtyrs,0,simwgt96);

    lmedexsim = ln(medexsim.*alivesim + (1-alivesim)); @ missing values are assigned large negative numbers @
    {lmxmnssim, lmxmnscnts, mxstdsim, mxcrlsim1, mxcrlcnt1, mxcrlsim2, mxcrlcnt2}
                         = simcrrl(pisim96,agesim96,lmedexsim,mssim2,alivesim,
                                   pistate_m,cohorts_m,mmtyrs,0,simwgt96);
    if prnres>1;
        ?;"Asset Profiles from the Model";; print aqntsim;
        ?;"Medex Profiles from the Model";; print mxprfsim;
        ?;"Medex Means from the Model";; print mxmnssim;
        ?;"Log Medex Means from the Model";; print lmxmnssim;
        ?;"Log Medex 1st autocorrelations from Model";; print mxcrlsim1;?; 
        ?;"Log Medex 2nd autocorrelations from Model";; print mxcrlsim2;?; 
    endif;
    if prngrph==1; 
        grphmtx(aqntsim,1,0,1,qnum_a,pinum_a,chrtnum_a);
        grphmtx(cprfsim,2,0,1,qnum_a,pinum_a,chrtnum_a);
        grphmtx(mxprfsim,3,0,1,qnum_m,pinum_m,chrtnum_m);
        grphmtx(mxmnssim,3,0,1,0,pinum_m,chrtnum_m);
    endif;

 /* We use lmxmnsdat and lmxstddat because we are not matching logged means or std deviations. */

    {mmtvec,obsvec,qntvec,pdfvec,mmttype} 
        = makemmts(PIdat,agedat[.,1],asstdat,MStatdat,obsdat,datawgts,aqntsim,
                   mxdat,mxobsdat,mxprfsim,mxmnssim,lmxmnsdat,lmxstddat,
                   mxcrlsim1,mxcrlsim2,0);

    criter = (1+negpen)*totobs*mmtvec'*_W*mmtvec;
    "GMM criterion value";; criter;

    if prngrph==1;

        {aqntsim,aqntcnts}   = simqunts(pisim96,agesim96,asstsim,mssim2,alivesim,
                                        pistate_a,cohorts_a,quants_a,mmtyrs,1,simwgt96);
        {cprfsim,cprfcnts}   = simqunts(pisim96,agesim96,conssim,mssim2,alivesim,
                                        pistate_a,cohorts_a,quants_a,mmtyrs,1,simwgt96);
        {mxprfsim,mxprfcnts} = simqunts(pisim96,agesim96,medexsim,mssim2,alivesim, 
                                        pistate_m,cohorts_m,quants_m,mmtyrs,1,simwgt96);
        {mxmnssim,mxmnscnts} = simqunts(pisim96,agesim96,medexsim,mssim2,alivesim, 
                                        pistate_m,cohorts_m,0,mmtyrs,1,simwgt96);

        grphmtx(aqntsim,1,1,1,qnum_a,pinum_a,chrtnum_a); /* Survivors Only */
        grphmtx(cprfsim,2,1,1,qnum_a,pinum_a,chrtnum_a);
        grphmtx(mxprfsim,3,1,1,qnum_m,pinum_m,chrtnum_m);
        grphmtx(mxmnssim,3,1,1,0,pinum_m,chrtnum_m);
    endif;

    sttime = timerec(sttime,"Doing one function evaluation");
    "We are at the end of function evaluation";; feval;?;
    feval = feval+1;

retp(criter); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
proc(0) = prnmedex(rhomx,fracar1,mxcoef);

    local j;

    ?;"Medex Coefficients";
    "rho_mx          ";;rhomx;
    "V(AR1)/V(AR+WN) ";;fracar1;?;

    "                       Mean      Variance";
    j=1; do until j>rows(mxcoef);
        mxlabel[j];;mxcoef[j,.];
    j=j+1; endo;
    ?;

retp; endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(1) = outrloop(mxvec);

    local mxmcoef, mxvcoef, mxcoef, parmvec, criter, func, iter, tim;

    ?;"We are at the beginning of function evaluation";; feval;

    mxvec    = mxvec.*zerovec[1:24] + fixvals[1:24].*(1-zerovec[1:24]);
    {rhomx, fracar1} = punscale_m(mxvec[1:2]);
    fracar1i = fracar1*(1-rhomx^2);    /* sigma of AR(1) innovations: global */           
    fracwn   = 1-fracar1;                    /* sigma of white noise: global */
    mxmcoef  = mxvec[3:13];
    mxvcoef  = mxvec[14:24];
    mxcoef   = mxmcoef~mxvcoef;
    prnmedex(rhomx,fracar1,mxcoef);

    {mnlnmxs, stdlnmxs, mnmx_pi, varmx_pi} = getmxtab(mxcoef,0.4,smplyrs);  /* Outputs are globals */

    parmvec  = fixvals[25:30];
    criter   = -100;

    if job==2;
        prnres  = 1;
        {psimplx_p,func,iter,tim} = AMOEBA(psimplx_p,ftol,maxsec,maxiter,&getcrit,prnum);
        parmvec = psimplx_p[minindc(func),.]';
        prnres  = 2;
        "Parameter estimates:";; getparms(parmvec);
        "MODEL EVALUATED AT INITIAL PARAMETER VALUES";
        criter  = getcrit(parmvec);

    elseif job==4;
        prnres  = 1;
        criter  = getcrit(parmvec);

    else;
        "WRONG JOB!!!";?;

    endif;

retp(criter); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(1) = oneloop(allparms);

    local mxmcoef, mxvcoef, mxcoef, parmvec, criter;

    ?;"We are at the beginning of function evaluation";; feval;

    allparms = allparms.*zerovec + fixvals.*(1-zerovec);
    {rhomx, fracar1} = punscale_m(allparms[1:2]);
    fracar1i = fracar1*(1-rhomx^2);    /* sigma of AR(1) innovations: global */           
    fracwn   = 1-fracar1;                    /* sigma of white noise: global */
    mxmcoef  = allparms[3:13];
    mxvcoef  = allparms[14:24];
    mxcoef   = mxmcoef~mxvcoef;
    prnmedex(rhomx,fracar1,mxcoef);

    {mnlnmxs, stdlnmxs, mnmx_pi, varmx_pi} = getmxtab(mxcoef,0.4,smplyrs);  /* Outputs are globals */

    parmvec  = allparms[25:30]; 
    criter   = -100;

    if job /= 4;
        criter  = getcrit(parmvec);
    else;
        "WRONG JOB!!!";?;
    endif;

retp(criter); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(0) = getexpr(allparms);
    
    local mxmcoef, mxvcoef, mxcoef, parmvec, execret, pisimx, agesimx, asstsimx, 
          mssim2x, conssimx, incsimx, alivex, mxsimx, astnoizx, simwgtx, aqntsim, 
          aqntcnts, cprfsim, cprfcnts, mxprfsim, mxprfcnts, mxmnssim, mxmnscnts, 
          imnssim, imnscnts, lmedexsim, lmxmnssim, lmxmnscnts, mxstdsim, 
          mxcrlsim1, mxcrlcnt1, mxcrlsim2, mxcrlcnt2;

    ?;"We are at the beginning of function evaluation";; feval;

    allparms = allparms.*zerovec + fixvals.*(1-zerovec);
    {rhomx, fracar1} = punscale_m(allparms[1:2]);
    fracar1i = fracar1*(1-rhomx^2);    /* sigma of AR(1) innovations: global */           
    fracwn   = 1-fracar1;                    /* sigma of white noise: global */
    mxmcoef  = allparms[3:13];
    mxvcoef  = allparms[14:24];
    mxcoef   = mxmcoef~mxvcoef;
    prnmedex(rhomx,fracar1,mxcoef);

    {mnlnmxs, stdlnmxs, mnmx_pi, varmx_pi} = getmxtab(mxcoef,0.4,smplyrs);  /* Outputs are globals */

    parmvec  = allparms[25:30]; 
    getparms(parmvec);
    savevecs();                           /* Save input vector for C program */

    if useMPI < 2;
        execret = exec(rulecall," ");
    elseif useMPI == 2;
        execret = exec("mpirun", rulecall);
    endif;

    loadsims();              /* Load simulation results into local directory */

 /* Get model-predicted summary stats */
    load path=^shkpath pisimx, agesimx, simwgtx, asstsimx, mssim2x, conssimx, 
                       incsimx, alivex, mxsimx, astnoizx;
    asstsimx = asstsimx.*astnoizx;

/*  Calculate consumption variance over the life-cycle */
/*
    local conssimx2, cgrowth, cg2, alivex2, lspanx, cg1m, cg2m, cgstd, 
          i, pii, avgstd, PIqnts, PIcnts, PItype;
    conssimx2 = alivex.*conssimx + (1-alivex);
    cgrowth = conssimx2[.,2:simyrs-1]./conssimx2[.,1:simyrs-2];
    cg2 = cgrowth^2;
    alivex2 = alivex[.,2:simyrs-1];
    lspanx = meanc(alivex2'+1e-9);
    cg1m = meanc((cgrowth.*alivex2)')./lspanx;
    cg2m = meanc((cg2.*alivex2)')./lspanx;
    cgstd = cg2m-cg1m^2; @ person by person @
    cgstd = sqrt(cgstd);
    {PIqnts,PIcnts,PItype} = PIquant2(pisimx,pistate_a,chktie);
    i=1; do until i>pinum_a;
        pii    = PItype.==i;
        avgstd = meanc(cgstd.*pii)/meanc(pii);
        "std of cons in quintile";; i;; avgstd;; meanc(lspanx.*pii)/meanc(pii);
    i=i+1; endo;

*/
    {aqntsim,aqntcnts}   = simqunts(pisimx,agesimx[.,1],asstsimx,mssim2x,alivex, 
                                    pistate_a,cohortsx,quants_a,simyrsx,0,simwgtx);
    {cprfsim,cprfcnts}   = simqunts(pisimx,agesimx[.,1],conssimx,mssim2x,alivex,
                                    pistate_a,cohortsx,quants_a,simyrsx,0,simwgtx);
    {mxprfsim,mxprfcnts} = simqunts(pisimx,agesimx[.,1],mxsimx,mssim2x,alivex,
                                    pistate_m,cohortsx,quants_m,simyrsx,0,simwgtx);
    {mxmnssim,mxmnscnts} = simqunts(pisimx,agesimx[.,1],mxsimx,mssim2x,alivex, 
                                    pistate_m,cohortsx,0,simyrsx,0,simwgtx);
    {imnssim,imnscnts}   = simqunts(pisimx,agesimx[.,1],incsimx,mssim2x,alivex, 
                                    pistate_a,cohortsx,0,simyrsx,0,simwgtx);

    lmedexsim = ln(mxsimx.*alivex + (1-alivex)); @ missing values are assigned large negative numbers @
    {lmxmnssim, lmxmnscnts, mxstdsim, mxcrlsim1, mxcrlcnt1, mxcrlsim2, mxcrlcnt2}
                         = simcrrl(pisimx,agesimx[.,1],lmedexsim,mssim2x,alivex, 
                                   pistate_m,cohortsx,simyrsx,0,simwgtx);
    if prnres>1;
        ?;"Asset Profiles from the Model";; print aqntsim;
        ?;"Medex Profiles from the Model";; print mxprfsim;
        ?;"Medex Means from the Model";; print mxmnssim;
        ?;"Log Medex Means from the Model";; print lmxmnssim;
        ?;"Log Medex 1st autocorrelations from Model";; print mxcrlsim1;?; 
        ?;"Log Medex 2nd autocorrelations from Model";; print mxcrlsim2;?; 
    endif;
    if prngrph==1; 
        grphmtx2(aqntsim,1,0,1,qnum_a,pinum_a);
        grphmtx2(cprfsim,2,0,1,qnum_a,pinum_a);
        grphmtx2(mxprfsim,3,0,1,qnum_m,pinum_m);
        grphmtx2(mxmnssim,3,0,1,0,pinum_m);
        grphmtx2(imnssim,4,0,1,0,pinum_a);
    endif;

    {aqntsim,aqntcnts}   = simqunts(pisimx,agesimx[.,1],asstsimx,mssim2x,alivex,
                                    pistate_a,cohortsx,quants_a,simyrsx,1,simwgtx);
    {cprfsim,cprfcnts}   = simqunts(pisimx,agesimx[.,1],conssimx,mssim2x,alivex,
                                    pistate_a,cohortsx,quants_a,simyrsx,1,simwgtx);
    {mxprfsim,mxprfcnts} = simqunts(pisimx,agesimx[.,1],mxsimx,mssim2x,alivex,
                                    pistate_m,cohortsx,quants_m,simyrsx,1,simwgtx);
    {mxmnssim,mxmnscnts} = simqunts(pisimx,agesimx[.,1],mxsimx,mssim2x,alivex, 
                                    pistate_m,cohortsx,0,simyrsx,1,simwgtx);
    {imnssim,imnscnts}   = simqunts(pisimx,agesimx[.,1],incsimx,mssim2x,alivex, 
                                    pistate_a,cohortsx,0,simyrsx,1,simwgtx);

    if prngrph==1; 
        grphmtx2(aqntsim,1,1,1,qnum_a,pinum_a); /* Survivors Only */
        grphmtx2(cprfsim,2,1,1,qnum_a,pinum_a);
        grphmtx2(mxprfsim,3,1,1,qnum_m,pinum_m);
        grphmtx2(mxmnssim,3,1,1,0,pinum_m);
        grphmtx2(imnssim,4,1,1,0,pinum_a);
    endif;

retp(); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
elife:  Simulate sequences of health, health cost and demographic shocks
        Here data are sorted by age, with each year having a single age
*/
proc(2) = elife(morts, mort_pi, healdats, heal_pi, calctype, initages, 
                initPI, initMS, initHS, survages, nn_3);

    local initage, isim94, goodage, agesimx, hlsimhx, hlsimwx, mssimx, pisimx, gotPI, 
          hsimpoh, hsimpow, mstpo, ii, ai2, frh, frw, frdead, agedist, 
          agedisth, agedistw, elspan, ageseqx2, _trx2, num, fralive, rs, 
          svec, survprob;

    initage  = initages[1];    
    _trx2    = dieage-initage+1;
    ageseqx2 = seqa(initage,1,_trx2+1);

    if calctype == 0;
        load path=^shkpath isim94;
        goodage = (isim94[.,1].>=initages[1]).*(isim94[.,1].<=initages[2]);
        isim94   = selif(isim94,goodage);
        if initMS > -1;
            isim94   = selif(isim94, isim94[.,3].==initMS);
        endif;
        if initHS[1] > -1;
            isim94   = selif(isim94, isim94[.,4].==(1-initHS[1]));
        endif;
        if initHS[2] > -1;
            isim94   = selif(isim94, isim94[.,5].==(1-initHS[2]));
        endif;
        if initPI[1] > -1;
            if initPI[2] > initPI[1];
                gotPI  = (isim94[.,2].>=initPI[1]).*(isim94[.,2].<initPI[2]);
                isim94 = selif(isim94, gotPI);
            else;
                isim94[.,2] = initPI[1]*ones(rows(isim94),1);
            endif;
        endif;

        "Number of observations";; rows(isim94);

        if rows(isim94)>4;
            num      = rndu(nn_3,1)*rows(isim94); /* random draws of indices */
            num      = floor(num)+ones(nn_3,1);
            pisimx   = isim94[num,2];
            agesimx  = isim94[num,1];
            mssimx   = isim94[num,3];     /* 0=>defunct, 1=>husband, 2=>wife */
            hlsimhx  = 1-isim94[num,4];     /* 0=>bad health, 1=>good health */
            hlsimwx  = 1-isim94[num,5];
        else;
            agesimx  = initage*ones(nn_3,1); 
            hlsimhx  = initHS[1]*ones(nn_3,1);  /* 0=>bad health, 1=>good health */
            hlsimwx  = initHS[2]*ones(nn_3,1);
            pisimx   = meanc(initPI)*ones(nn_3,1);
            mssimx   = initMS*ones(nn_3,1);   /* 0=>defunct, 1=>husband, 2=>wife */
        endif;

    elseif calctype == 1;
        agesimx  = initage*ones(nn_3,1); 
        hlsimhx  = initHS[1]*ones(nn_3,1);  /* 0=>bad health, 1=>good health */
        hlsimwx  = initHS[2]*ones(nn_3,1);
        pisimx   = meanc(initPI)*ones(nn_3,1);
        mssimx   = initMS*ones(nn_3,1);   /* 0=>defunct, 1=>husband, 2=>wife */
    endif;

    ii=1; do until ii >  _trx2;      /* Update with transition probabilities */
        ai2     = initage-bornage+ii;                 /* for appropriate age */
        {hsimpoh, hsimpow, mstpo} = gethms(ai2, mssimx[.,ii], hlsimhx[.,ii], 
                                           hlsimwx[.,ii], pisimx, heal_pi, 
                                           mort_pi, morts, healdats, nn_3);
        hlsimhx = hlsimhx~hsimpoh;
        hlsimwx = hlsimwx~hsimpow;
        mssimx  = mssimx~mstpo;
    ii=ii+1; endo;

    frdead  = (mssimx.==0);
    frh     = (mssimx.==1);
    frw     = (mssimx.==2);

    agedist  = meanc(frdead[.,2:_trx2+1])-meanc(frdead[.,1:_trx2]);
    agedist  = 0|agedist;
    agedisth = meanc(frh[.,1:_trx2])-meanc(frh[.,2:_trx2+1]);
    agedisth = agedisth/meanc(frh[.,1]);
    agedisth = 0|agedisth;
    agedistw = meanc(frw[.,1:_trx2])-meanc(frw[.,2:_trx2+1]);
    agedistw = agedistw/meanc(frw[.,1]);
    agedistw = 0|agedistw;
    agedist  = agedist~agedisth~agedistw;

    elspan   = (ageseqx2-initage)'agedist;

/*  Find Survival probabilities for various ages  */

    fralive  = 1 - meanc(frdead);
    rs       = rows(survages);
    survprob = zeros(rs,1);
 
    ii=1; do until ii>rs;
        svec = ageseqx2.==survages[ii];
        survprob[ii] = fralive'svec;
    ii=ii+1; endo;

    format /ro 12,4;
    "Number of simulated observations  ";; nn_3;
    "Calculation type (0 => data dist) ";; calctype;
    "Initial age                       ";; initage;; initages';
    "Initial Health Status (M,F)       ";; initHS';
    "Initial Marital Status            ";; initMS;
    "Permanent Income Rank             ";; initPI';
    "Expected lifespan                 ";; elspan;
    "Survival probabilities by age     ";; survages~survprob;?;

    "Fraction of Households that are men, women or dead";;
    ageseqx2~meanc(frh)~meanc(frw)~meanc(frdead)~agedist;?;

    "Average health status (1=>good) for men and women that are still alive";;
    ageseqx2~(meanc(hlsimhx.*frh)./meanc(frh))~(meanc(hlsimwx.*frw)./meanc(frw));?;
		 
retp(elspan,survprob); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(0) = getelife(_nn3,mortprfs,mort_pi,hsprobs,heal_pi,initages,survages); 

    local elspan, lifetab, lifetab2, rn, rs, recsize, recsize2, lifetab3, 
          lifetab4, lifetab5, lifetab6, survprob, tablab, iPI, iPI2, PIstep,
          iS, iR, iC;

    lifetab  = zeros(8,6); @ PI percentiles @
    lifetab2 = zeros(8,6); @ PI quintiles   @

    rs       = rows(survages);
    recsize  = rows(lifetab)|cols(lifetab)|rs;
    lifetab3 = arrayinit(recsize,0);
    recsize  = recsize[3|1|2];
    lifetab4 = arrayinit(recsize,0);
    recsize2 = rows(lifetab2)|cols(lifetab2)|rs;
    lifetab5 = arrayinit(recsize2,0);
    recsize2 = recsize2[3|1|2];
    lifetab6 = arrayinit(recsize2,0);

/*  onesim = ones(_nn3,1); initdist(exprage);  */
    {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                              -1, -1, -ones(2,1), survages, _nn3);
    lifetab[6,6] = elspan[1];
    lifetab[7,6] = elspan[1];
    lifetab[8,6] = elspan[1];

    lifetab2[6,6] = elspan[1];
    lifetab2[7,6] = elspan[1];
    lifetab2[8,6] = elspan[1];

    lifetab3 = putarray(lifetab3,6|6,survprob');
    lifetab3 = putarray(lifetab3,7|6,survprob');
    lifetab3 = putarray(lifetab3,8|6,survprob');

    lifetab5 = putarray(lifetab5,6|6,survprob');
    lifetab5 = putarray(lifetab5,7|6,survprob');
    lifetab5 = putarray(lifetab5,8|6,survprob');

/*  Males */
    {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                              -1, 1, -ones(2,1), survages, _nn3);
    lifetab[7,2:3]  = elspan[1]*ones(1,2);
    lifetab2[7,2:3] = elspan[1]*ones(1,2);
    lifetab3 = putarray(lifetab3,7|2,survprob');
    lifetab3 = putarray(lifetab3,7|3,survprob');
    lifetab5 = putarray(lifetab5,7|2,survprob');
    lifetab5 = putarray(lifetab5,7|3,survprob');

    {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                              -1, 1, zeros(2,1), survages, _nn3);
    lifetab[6,3]  = elspan[1];
    lifetab2[6,3] = elspan[1];
    lifetab3 = putarray(lifetab3,6|3,survprob');
    lifetab5 = putarray(lifetab5,6|3,survprob');

    {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                              -1, 1, ones(2,1), survages, _nn3);
    lifetab[6,2]  = elspan[1];
    lifetab2[6,2] = elspan[1];
    lifetab3 = putarray(lifetab3,6|2,survprob');
    lifetab5 = putarray(lifetab5,6|2,survprob');

/*  Females */
    {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                              -1, 2, -ones(2,1), survages, _nn3);
    lifetab[7,4:5]  = elspan[1]*ones(1,2);
    lifetab2[7,4:5] = elspan[1]*ones(1,2);
    lifetab3 = putarray(lifetab3,7|4,survprob');
    lifetab3 = putarray(lifetab3,7|5,survprob');
    lifetab5 = putarray(lifetab5,7|4,survprob');
    lifetab5 = putarray(lifetab5,7|5,survprob');

    {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                              -1, 2, zeros(2,1), survages, _nn3);
    lifetab[6,5]  = elspan[1];
    lifetab2[6,5] = elspan[1];
    lifetab3 = putarray(lifetab3,6|5,survprob');
    lifetab5 = putarray(lifetab5,6|5,survprob');

    {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                              -1, 2, ones(2,1), survages, _nn3);
    lifetab[6,4]  = elspan[1];
    lifetab2[6,4] = elspan[1];
    lifetab3 = putarray(lifetab3,6|4,survprob');
    lifetab5 = putarray(lifetab5,6|4,survprob');

    {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                              -1, -1, zeros(2,1), survages, _nn3);
    lifetab[8,3|5]  = elspan[1]*ones(1,2);
    lifetab2[8,3|5] = elspan[1]*ones(1,2);
    lifetab3 = putarray(lifetab3,8|3,survprob');
    lifetab3 = putarray(lifetab3,8|5,survprob');
    lifetab5 = putarray(lifetab5,8|3,survprob');
    lifetab5 = putarray(lifetab5,8|5,survprob');

    {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                              -1, -1, ones(2,1), survages, _nn3);
    lifetab[8,2|4]  = elspan[1]*ones(1,2);
    lifetab2[8,2|4] = elspan[1]*ones(1,2);
    lifetab3 = putarray(lifetab3,8|2,survprob');
    lifetab3 = putarray(lifetab3,8|4,survprob');
    lifetab5 = putarray(lifetab5,8|2,survprob');
    lifetab5 = putarray(lifetab5,8|4,survprob');

    rn = 0;
    PIstep = 0.2;
    iPI=PIstep; do until iPI > 1;
        rn     = rn+1;
        iPI2   = iPI-PIstep;
     /* First, consider observations where PI falls in an interval */

        lifetab2[rn,1] = iPI/PIstep;
        lifetab5 = putarray(lifetab5,rn|1,iPI*ones(1,rs)/PIstep);

        {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                                  (iPI2|iPI), -1, -ones(2,1), survages, _nn3);
        lifetab2[rn,6] = elspan[1];
        lifetab5 = putarray(lifetab5,rn|6,survprob');

        {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                                  (iPI2|iPI), 1, zeros(2,1), survages, _nn3);
        lifetab2[rn,3] = elspan[1];
        lifetab5 = putarray(lifetab5,rn|3,survprob');

        {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                                  (iPI2|iPI), 1, ones(2,1), survages, _nn3);
        lifetab2[rn,2] = elspan[1];
        lifetab5 = putarray(lifetab5,rn|2,survprob');

        {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                                  (iPI2|iPI), 2, zeros(2,1), survages, _nn3);
        lifetab2[rn,5] = elspan[1];
        lifetab5 = putarray(lifetab5,rn|5,survprob');

        {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                                  (iPI2|iPI), 2, ones(2,1), survages, _nn3);
        lifetab2[rn,4] = elspan[1];
        lifetab5 = putarray(lifetab5,rn|4,survprob');

     /* Next, use all observations, but set PI to pre-determined value */

        iPI2 = iPI - PIstep/2;

        lifetab[rn,1] = iPI2;
        lifetab3 = putarray(lifetab3,rn|1,iPI2*ones(1,rs));

        {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 0, initages, 
                                  iPI2|iPI2, -1, -ones(2,1), survages, _nn3);
        lifetab[rn,6] = elspan[1];
        lifetab3 = putarray(lifetab3,rn|6,survprob');

        {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 1, initages, 
                                  iPI2, 1, zeros(2,1), survages, _nn3);
        lifetab[rn,3] = elspan[1];
        lifetab3 = putarray(lifetab3,rn|3,survprob');

        {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 1, initages, 
                                  iPI2, 1, ones(2,1), survages, _nn3);
        lifetab[rn,2] = elspan[1];
        lifetab3 = putarray(lifetab3,rn|2,survprob');

        {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 1, initages, 
                                  iPI2, 2, zeros(2,1), survages, _nn3);
        lifetab[rn,5] = elspan[1];
        lifetab3 = putarray(lifetab3,rn|5,survprob');

        {elspan,survprob} = elife(mortprfs, mort_pi, hsprobs, heal_pi, 1, initages, 
                                  iPI2, 2, ones(2,1), survages, _nn3);
        lifetab[rn,4] = elspan[1];
        lifetab3 = putarray(lifetab3,rn|4,survprob');

    iPI=iPI+PIstep; endo;

    tablab = "     PI rank        MGood         MBad        FGood         FBad          All";
    
    "Life Expectancy Tables ";
    tablab;;
    lifetab;?;

    tablab;;
    lifetab2;?;

    save path=^shkpath lifetab;

    iS=1; do until iS>rs;
        iR=1; do until iR>recsize[2];
            iC=1; do until iC>recsize[3];
                lifetab4[iS,iR,iC] = lifetab3[iR,iC,iS];
            iC=iC+1; endo;
        iR=iR+1; endo;

        "Probability of being alive at age";; survages[iS];
        tablab;;
        getmatrix(lifetab4,iS);?;
    iS=iS+1; endo;

    iS=1; do until iS>rs;
        iR=1; do until iR>recsize2[2];
            iC=1; do until iC>recsize2[3];
                lifetab6[iS,iR,iC] = lifetab5[iR,iC,iS];
            iC=iC+1; endo;
        iR=iR+1; endo;

        "Probability of being alive at age";; survages[iS];
        tablab;;
        getmatrix(lifetab6,iS);?;
    iS=iS+1; endo;

retp; endp;

