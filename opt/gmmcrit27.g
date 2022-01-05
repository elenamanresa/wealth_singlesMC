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

proc(1) = getcrit(parmvec);

    local PIsim96, agesim96, simwgt96, asim96, asstsim, conssim, medexsim, 
          mssim2, alivesim, asstnoiz,  execret, aqntsim, aqntcnts, cprfsim, 
          cprfcnts, mxprfsim, mxprfcnts, mxmnssim, mxmnscnts, lmedexsim, 
          lmxmnssim, lmxmnscnts, mxstdsim, mxcrlsim1, mxcrlcnt1, mxcrlsim2, 
          mxcrlcnt2, mmtvec, obsvec, qntvec, pdfvec, mmttype, criter;

    getparms(parmvec);
    savevecs();                           /* Save input vector for C program */
    if prnout==1; output off; endif;

    if useMPI < 2;
        execret = exec(rulecall," ");
    elseif useMPI == 2;
        execret = exec("mpirun", rulecall);
    endif;

    if prnout==1; output on; endif;
    loadsims();              /* Load simulation results into local directory */

    sttime = timerec(sttime,"Doing one function evaluation");

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

    {mnlnmxs, stdlnmxs, mnmx_pi, varmx_pi} = getmxtab(mxcoef,0.1,smplyrs);  /* Outputs are globals */

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

    {mnlnmxs, stdlnmxs, mnmx_pi, varmx_pi} = getmxtab(mxcoef,0.1,smplyrs);  /* Outputs are globals */

    parmvec  = allparms[25:30]; 
    getparms(parmvec);
    savevecs();                           /* Save input vector for C program */

    if useMPI < 2;
        execret = exec(rulecall," ");
    elseif useMPI == 2;
        execret = exec("mpirun", rulecall);
    endif;

    loadsims();              /* Load simulation results into local directory */


retp(); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
