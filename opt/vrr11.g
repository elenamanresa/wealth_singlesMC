/*---------------------------------------------------------------------------*/
/*  Procedures for imputing medical expense shock process:                   */
/*      moments and compares to data moments.                                */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(4) = infershk(maxrat);

local asstsim, ztasim, medexsim, Medicaidsim, healsimh, healsimw, healsim, 
      conssim, alivesim, agesim, mstatsim, malesim, pisim96, posvalue, lamsim, 
      i, ismale, alivei, livefrac, rnrec, rn, cn, agei, malesim2, PIsim2, 
      agesim2, healsim2, asstsim2, lamsim2, mxsim2, conssim2, X, ck, rscl, 
      im, yvec, xTx, xTy, prfmcoef, ypred, errvec, stderr, Rsq, resids, 
      prfvcoef, prfcoef, mnprfs, stdprfs, mnprf_pi, varprf_pi, stuff;
      
    load path=^shkpath healsimh, healsimw, mstatsim;
    malesim = (mstatsim.==1);
    healsim = healsimh.*malesim + healsimw.*(1-malesim);

    clear mstatsim, healsimh, healsimw; 
    load path=^shkpath asstsim, medexsim, Medicaidsim, conssim, alivesim, agesim, pisim96;

    medexsim = medexsim+Medicaidsim; /* Add in expenditures picked up by Medicaid */
    posvalue = medexsim.>0;
    medexsim = medexsim.*posvalue + (1-posvalue);
    posvalue = conssim.>0;
    conssim  = conssim.*posvalue + (1-posvalue);    
    lamsim   = medexsim./conssim;        /* Infer from FOC */
    posvalue = lamsim.>maxrat;
    lamsim   = lamsim.*(1-posvalue) + posvalue*maxrat; /* Chop off outliers */

    clear Medicaidsim, posvalue;

    malesim2 = {}; 
    PIsim2   = {};
    agesim2  = {};
    healsim2 = {};
    asstsim2 = {};
    mxsim2   = {};
    conssim2 = {};
    lamsim2  = {};
    rnrec    = ones(simyrs,2);
    livefrac = zeros(simyrs,1);

    i=1; do until i>simyrs;
        alivei   = alivesim[.,1];
        livefrac[i] = meanc(alivei);
        alivei   = selif(seqa(1,1,rows(alivei)),alivei); 
        rn       = rows(alivei);
        if i==1;
            rnrec[i,2] = rn;
        else;
            rnrec[i,1] = rnrec[i-1,2]+1;
            rnrec[i,2] = rnrec[i-1,2]+rn;
        endif;
        malesim2 = malesim2|malesim[alivei,1];
        PIsim2   = PIsim2|pisim96[alivei,1];
        agesim2  = agesim2|agesim[alivei,1];
        healsim2 = healsim2|(healsim[alivei,1].==0);  /* Good Health is the omitted category */
        conssim2 = conssim2|conssim[alivei,1];
        mxsim2   = mxsim2|medexsim[alivei,1];
        asstsim2 = asstsim2|(asstsim[alivei,1]/1000);
        lamsim2  = lamsim2|lamsim[alivei,1];
        
        cn = simyrs-i+1;
                                                             
        if i<simyrs;
            malesim  = malesim[.,2:cn];
            agesim   = agesim[.,2:cn];
            healsim  = healsim[.,2:cn];
            asstsim  = asstsim[.,2:cn];
            conssim  = conssim[.,2:cn];
            lamsim   = lamsim[.,2:cn]; 
            alivesim = alivesim[.,2:cn];  
            medexsim = medexsim[.,2:cn];         
        endif;
    i=i+1; endo;

    stuff = conssim2~mxsim2~lamsim2;
    meanc(stuff)~stdc(stuff);?;
    quantile(stuff,0.25|0.5|0.9|0.99|0.999);?;

    clear conssim2, asstsim2, mxsim2;      


   /*-----Construct the X matrix, containing PIA & Assets interacted w/ age-----*/

    rn   = rows(healsim2);
    x    = ones(rn,1)~agesim2~(agesim2^2/100)~(agesim2^3/10000);
    x    = x~healsim2~(healsim2.*agesim2)~malesim2~(malesim2.*agesim2);
    x    = x~PIsim2~(PIsim2.*agesim2)~(PIsim2^2);
    cn   = cols(x);

    clear agesim2, malesim2, healsim2, PIsim2;

    yvec = ln(lamsim2);
    
   /*----------Compute X'X and X'y in stages, for memory constraints---------*/
    ck   = 2000;
    rscl = 1;
    rscl = rscl*ck;
    im   = floor(rn/ck);
    xTx  = 0;
    xTy  = 0;
    i=1; do until i > im;
        xTx = xTx + x[(i-1)*ck+1:i*ck,.]'x[(i-1)*ck+1:i*ck,.]/rscl;
        xTy = xTy + x[(i-1)*ck+1:i*ck,.]'yvec[(i-1)*ck+1:i*ck]/rscl;
    i=i+1; endo;
    if (im*ck) < rn;
        xTx = xTx + x[im*ck+1:rn,.]'x[im*ck+1:rn,.]/rscl;
        xTy = xTy + x[im*ck+1:rn,.]'yvec[im*ck+1:rn]/rscl;
    endif;

    prfmcoef = inv(xTx)*(xTy);    
    ypred    = x*prfmcoef;

    errvec = yvec-ypred;
    yvec   = yvec-meanc(yvec);
    Rsq    = 1 - vcx(errvec)/vcx(yvec);
    stderr = invpd(xTx)*vcx(errvec);
    stderr = sqrt(diag(stderr));
    ?;"Preference Shock Coefficients";
    i=1; do until i>cn;
        mxlabel[i];;prfmcoef[i]~stderr[i];
    i=i+1; endo;
    ?;"Residual Variance =            ";; vcx(errvec);
    "R-squared for the regression = ";; Rsq;?;
    
   /*------------Next, do a regression for the residual variance-------------*/

    yvec = errvec^2;

    xTy  = 0;
    i=1; do until i > im;
        xTy = xTy + x[(i-1)*ck+1:i*ck,.]'yvec[(i-1)*ck+1:i*ck]/rscl;
    i=i+1; endo;
    rn = rows(x);
    if (im*ck) < rows(x);
        xTy = xTy + x[im*ck+1:rn,.]'yvec[im*ck+1:rn]/rscl;
    endif;

    prfvcoef = inv(xTx)*(xTy);    
    ypred    = x*prfvcoef;
   
    errvec = yvec-ypred;
    yvec   = yvec-meanc(yvec);
    Rsq    = 1 - (errvec'errvec)/(yvec'yvec);
    stderr = invpd(xTx)*vcx(errvec);
    stderr = sqrt(diag(stderr));
    ?;"Residual Variance Coefficients"; 
    i=1; do until i>cn;
        mxlabel[i];;prfvcoef[i]~stderr[i];
    i=i+1; endo;
    ?;"Residual Variance =            ";; vcx(errvec);
    "R-squared for the regression = ";; Rsq;?;

    prfcoef = prfmcoef~prfvcoef;
    {mnprfs,stdprfs,mnprf_pi,varprf_pi} = getmxtab(prfcoef,-1,0);

    "Mean Coefficients";
     "         age        hbad         hgood         wbad       wgood          PI        PI^2";;
     ageseq~mnprfs~mnprf_pi;?;

    "Variance (Square Root) Coefficients";
     "         age        hbad         hgood         wbad       wgood          PI        PI^2";;
     ageseq~stdprfs~varprf_pi;?;


retp(mnprfs,stdprfs,mnprf_pi,varprf_pi); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
proc(2) = simmxd(mnlnmxs,stdmxs,mnmx_pi,varmx_pi,twoyrs,details);

    local ztasim96, pisim96, agesim96, healsimh, healsimw, mstatsim, xisim, 
          alivesim, mxinnovp, ztasim, lnmxsim, ageindx, i, ai, too_old, 
          mnlnmxh, mnlnmxw, mnlnmx, stdlnmxh, stdlnmxw, stdlnmx, PIeffect,  
          lnmxsimh, lnmxsimw, mxsimh, mxsimw, mxsim, mxsim1, tempprf, tempnum, 
          mxqnts, qnttab, fracaliv, meanmx, avgage, malesim, healsim, agesim, 
          malesim2, pisim2, agesim2, healsim2, mxsim2, yearsim2, alivei, rn, cn, 
          mxdata, fmlsim, rn2, Xmat, xystats, label, yvec, ck, rscl, im, xTx, 
          xTy, Qmat, _beta, ypred, errvec, Rsq, stderr, stderr_w, xTerr,
          mmtcols2, ser, simyrs2;
    
    load path=^shkpath pisim96, agesim96, healsimh, healsimw, mstatsim;


    simyrs2  = @simyrs @25;
    ztasim96 = rndn(1,_nn)*sqrt(fracar1);
    xisim    = rndn(_nn,simyrs2)*sqrt(fracwn);


/*
agesim96 = 70*ones(_nn,1);
mstatsim = mstatsim[.,1]*ones(1,simyrs2);
healsimh = healsimh[.,1]*ones(1,simyrs2);
healsimw = healsimw[.,1]*ones(1,simyrs2);
*/

    mstatsim = mstatsim[.,1:simyrs2];
    ageindx  = agesim96-bornage;
                         
    mxinnovp = rndn(simyrs2,_nn)*sqrt(fracar1i);
    ztasim   = recserar(mxinnovp, ztasim96, rhomx*ones(1,_nn));
    ztasim   = ztasim';
    lnmxsim  = ztasim+xisim;
    clear mxinnovp, ztasim, xisim;
    lnmxsimh = {};
    lnmxsimw = {};
    
    i=1; do until i> simyrs2;  
        ai       = ageindx+i; 
        too_old  = ai.>_TR;
        ai       = ai.*(1-too_old)+too_old*_TR;
       
        mnlnmxh  = mnlnmxs[ai,1].*(1-healsimh[.,i]) + mnlnmxs[ai,2].*healsimh[.,i];    
        PIeffect = mnmx_pi[ai,1].*pisim96+mnmx_pi[ai,2].*(pisim96^2);
        mnlnmxh  = mnlnmxh+pieffect;
      
        mnlnmxw  = mnlnmxs[ai,3].*(1-healsimw[.,i]) + mnlnmxs[ai,4].*healsimw[.,i]; 
        PIeffect = mnmx_pi[ai,3].*pisim96+mnmx_pi[ai,4].*(pisim96^2);
        mnlnmxw  = mnlnmxw+pieffect;
 
        stdlnmxh = stdmxs[ai,1].*(1-healsimh[.,i]) + stdmxs[ai,2].*healsimh[.,i]; 
        PIeffect = varmx_pi[ai,1].*pisim96+varmx_pi[ai,2].*(pisim96^2);
        stdlnmxh = sqrt(stdlnmxh^2+pieffect);
   
        stdlnmxw = stdmxs[ai,3].*(1-healsimw[.,i]) +  stdmxs[ai,4].*healsimw[.,i];
        PIeffect = varmx_pi[ai,3].*pisim96+varmx_pi[ai,4].*(pisim96^2);
        stdlnmxw = sqrt(stdlnmxw^2+pieffect);

        lnmxsimh = lnmxsimh~(mnlnmxh+stdlnmxh.*lnmxsim[.,i]);
        lnmxsimw = lnmxsimw~(mnlnmxw+stdlnmxw.*lnmxsim[.,i]);
    i=i+1; endo;       

    
    mxsimh   = exp(lnmxsimh);
    mxsimw   = exp(lnmxsimw);

    mxqnts   = (0.1|0.25|0.50|0.75|0.90|0.95|0.99);

    malesim  = (mstatsim.==1);
    fmlsim   = (mstatsim.==2);
    alivesim = (mstatsim.>0);

    agesim   = agesim96+seqa(0,1,simyrs2)';
    too_old  = agesim.>dieage;
    agesim   = agesim.*(1-too_old)+too_old*dieage;
    healsim  = healsimh[.,1:simyrs2].*malesim + healsimw[.,1:simyrs2].*fmlsim;
    mxsim    = mxsimh.*malesim + mxsimw.*fmlsim;
    clear mstatsim, healsimh, healsimw, mxsimh, mxsimw, lnmxsimh, lnmxsimw; 

    mxsim1   = zeros(_nn,1)~mxsim[.,2:simyrs2];@ Convert annual expenses into 2-year averages @
    i=2; do until i> simyrs2;
        mxsim1[.,i] = (mxsim[.,i-1]+mxsim[.,i])/2; 
    i=i+1; endo;
    mxsim1[.,1] = mxsim1[.,2];

    fracaliv = meanc(alivesim);
    avgage   = meanc(agesim.*alivesim)./fracaliv;
    meanmx   = meanc(mxsim.*alivesim)./fracaliv;
    mxsim    = mxsim.*miss(alivesim,0);
    mxqnts   = (0.1|0.25|0.50|0.75|0.90|0.95|0.99);

    qnttab = {};
    i=1; do until i> simyrs2;  
        {tempprf,tempnum} = getqunt2(mxsim[.,i],mxqnts);
        qnttab = qnttab~tempprf[.,2];
    i=i+1; endo;  

    if details==1;    
        ?;"Simulated one-year health cost quantiles by simulation year";;
        mxqnts~qnttab;
        "    Mean mx";; meanmx';
        "    %alive ";; fracaliv';
        "    Avg Age";; avgage';
    endif;

    malesim2 = {}; 
    PIsim2   = {};
    agesim2  = {};
    healsim2 = {};
    yearsim2 = {};
    rn2      = {};
    mxsim2   = {};

    if twoyrs==0; i=1;  elseif twoyrs==1; i=2; endif;


    i=1; do until i>simyrs2;
        alivei   = alivesim[.,i];
        alivei   = selif(seqa(1,1,_nn),alivei);
        rn2      = rn2|alivei;
        rn       = rows(alivei);
        yearsim2 = yearsim2|((momyr1+i-1)*ones(rn,1));
        malesim2 = malesim2|malesim[alivei,i];
        PIsim2   = PIsim2|pisim96[alivei];
        agesim2  = agesim2|agesim[alivei,i];
        healsim2 = healsim2|(healsim[alivei,i].==0); /* reverse sign */
        if twoyrs==0;
            mxsim2   = mxsim2|mxsim[alivei,i];
        elseif twoyrs==1;
            mxsim2   = mxsim2|mxsim1[alivei,i];
        endif;
    i=i+1; endo;

    mxdata = mxsim2~rn2~agesim2~PIsim2~malesim2~healsim2~yearsim2;
    save path=^shkpath mxdata;
    clear mxdata, rn2, yearsim2;

    fracaliv = meanc(alivesim);
    avgage   = meanc(agesim.*alivesim)./fracaliv;
    meanmx   = meanc(mxsim1.*alivesim)./fracaliv;
    mxsim1   = mxsim1.*miss(alivesim,0);
    mxqnts   = (0.1|0.25|0.50|0.75|0.90|0.95|0.99);

    qnttab = {};
    i=1; do until i>simyrs2;  
        {tempprf,tempnum} = getqunt2(mxsim1[.,i],mxqnts);
        qnttab = qnttab~tempprf[.,2];
    i=i+1; endo;  

    if details==1;    
        ?;"Simulated two-year health cost quantiles by AHEAD wave year";;
        mxqnts~qnttab;
        "    Mean mx";; meanmx';
        "    %alive ";; fracaliv';
        "    Avg Age";; avgage';
    endif;

    clear mxsim, mxsim1;

    Xmat = agesim2~(agesim2^2)~(agesim2^3)~(agesim2^4)~healsim2~(healsim2.*agesim2);
    Xmat = Xmat~(malesim2.*agesim2)~(PIsim2.*agesim2)~PIsim2~(PIsim2^2)~malesim2~ones(rows(mxsim2),1);

    yvec = ln(mxsim2);

    rn  = rows(Xmat);
    cn  = cols(Xmat);
    rn2 = rn - floor(rn/100)*100;
    ?;"Number of observations ";; rn;; rn2;?;

    let string label = {"Age        ","Age^2      ",
                        "Age^3      ","Age^4      ",
                        "Health     ","Health*Age ",
                        "Male*Age   ","PI*Age     ",
                        "PI         ","PI^2       ",
                        "Male       ","Constant   ",
                        "ln(medexs) "};

    xystats=(meanc(xmat)~stdc(xmat))|(meanc(yvec)~stdc(yvec));
    "Descriptive Statistics";
    "Variable         Mean    Std Dev";
    i=1; do until i> cn+1;  
        label[i];;xystats[i,.];
    i=i+1; endo;  

    save path=^shkpath Xmat;
    clear agesim2,PIsim2,malesim2,healsim2;
    
   /*----------Compute X'X and X'y in stages, for memory constraints---------*/
    ck   = 2000;
    rscl = 1;
    rscl = rscl*ck;
    im   = floor(rn/ck);
    xTx  = 0;
    xTy  = 0;
    i=1; do until i > im;
        xTx = xTx + Xmat[(i-1)*ck+1:i*ck,.]'Xmat[(i-1)*ck+1:i*ck,.]/rscl;
        xTy = xTy + Xmat[(i-1)*ck+1:i*ck,.]'yvec[(i-1)*ck+1:i*ck]/rscl;
    i=i+1; endo;
    if (im*ck) < rn;
        xTx = xTx + Xmat[im*ck+1:rn,.]'Xmat[im*ck+1:rn,.]/rscl;
        xTy = xTy + Xmat[im*ck+1:rn,.]'yvec[im*ck+1:rn]/rscl;
    endif;

    Qmat   = invpd(xTx);
    _beta  = Qmat*(xTy);
    
    ypred  = Xmat*_beta;
    errvec = yvec-ypred;
    yvec   = yvec-meanc(yvec);

    xTerr  = 0;
    Xmat   = Xmat.*errvec;
    i=1; do until i > im;
        xTerr = xTerr + Xmat[(i-1)*ck+1:i*ck,.]'Xmat[(i-1)*ck+1:i*ck,.]/rscl;
    i=i+1; endo;
    if (im*ck) < rn;
        xTerr = xTerr + Xmat[im*ck+1:rn,.]'Xmat[im*ck+1:rn,.]/rscl;
    endif;

    Rsq      = 1 - (errvec'errvec)/(yvec'yvec);
    stderr   = Qmat*vcx(errvec);
    stderr   = sqrt(diag(stderr)/rscl);
    stderr_w = Qmat*xTerr*Qmat;
    stderr_w = sqrt(diag(stderr_w)/rscl);
    ser      = stdc(errvec);
    ?;"Medex Coefficients";
    "Variable        Coeff    Std Err";
    i=1; do until i> cn;  
        label[i];;_beta[i];;stderr[i];;stderr_w[i];
    i=i+1; endo;  
    "Std error of regression =      ";; ser;
    "R-squared for the regression = ";; Rsq;?;

retp(_beta,ser); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
proc(0) = summxd(twoyrs,onlymyrs);

    local oopmedexsim, medexsim, Medicaidsim, transfersim, mssim2, pisim96, agesim96, 
          mxsim96, HHIDsim, healsimh, healsimw, alivesim, aliveavg, malesim, fmlsim, 
          healsim, agesim, too_old, malesim2, PIsim2, HHIDsim2, agesim2, healsim2, 
          yearsim, yearsim2, rn2, mxsim2, mdcdsim2, transim2, simyrs2, i0, i, rn, 
          alivei, mxdist;

    load path=^iopath oopmedexsim, Medicaidsim, transfersim, mssim2;
    load path=^shkpath pisim96, agesim96, healsimh, healsimw, mxsim96, HHIDsim;

    medexsim    = reshape(oopmedexsim,simyrs+1,_nn)';
    Medicaidsim = reshape(Medicaidsim,simyrs+1,_nn)';
    transfersim = reshape(transfersim,simyrs+1,_nn)';
    mssim2      = reshape(mssim2,simyrs+1,_nn)';
    alivesim    = mssim2.>0;
    aliveavg    = meanc(alivesim);
    malesim     = (mssim2[.,1].==1)*ones(1,simyrs+1);
    fmlsim      = (mssim2[.,1].==2)*ones(1,simyrs+1);
    healsim     = healsimh.*malesim + healsimw.*fmlsim;
    agesim      = agesim96+seqa(0,1,simyrs)';
    too_old     = agesim.>dieage;
    agesim      = agesim.*(1-too_old)+too_old*dieage;
    yearsim     = seqa(0,1,simyrs+1)+momyr1;

    medexsim[.,1] = mxsim96;

    if twoyrs==1;
        medexsim = (medexsim[.,1:simyrs]+medexsim[.,2:simyrs+1])/2; @ 2-year averages @
        medexsim = mxsim96~medexsim;
        Medicaidsim = (Medicaidsim[.,1:simyrs]+Medicaidsim[.,2:simyrs+1])/2; @ 2-year averages @
        Medicaidsim = zeros(_nn,1)~Medicaidsim;
        transfersim = (transfersim[.,1:simyrs]+transfersim[.,2:simyrs+1])/2; @ 2-year averages @
        transfersim = zeros(_nn,1)~transfersim;
    endif;

    if onlymyrs==1;                                             @ AHEAD Survey Years @
        alivesim    = alivesim[.,mmtcols];
        medexsim    = medexsim[.,mmtcols];
        Medicaidsim = Medicaidsim[.,mmtcols];
        transfersim = transfersim[.,mmtcols];
        malesim     = malesim[.,mmtcols];
        healsim     = healsim[.,mmtcols];
        agesim      = agesim[.,mmtcols];
        aliveavg    = aliveavg[mmtcols];
        yearsim     = yearsim[mmtcols];
    endif;

    mxsim2   = {};
    mdcdsim2 = {};
    transim2 = {};
    malesim2 = {};
    healsim2 = {}; 
    agesim2  = {};
    yearsim2 = {};
    PIsim2   = {};
    rn2      = {};
    HHIDsim2 = {};

/*  Note that there is no simulated medex for year 1.  Drop these observations  */

    simyrs2 = simyrs;
    i0 = 2;
    if onlymyrs==1;
        simyrs2 = rows(mmtcols);
    else;
        if twoyrs==1; i0=3; endif;
    endif;

    i=i0; do until i>simyrs2;
        alivei   = alivesim[.,i];
        alivei   = selif(seqa(1,1,_nn),alivei);
        rn       = rows(alivei);
        mxsim2   = mxsim2|medexsim[alivei,i];
        mdcdsim2 = mdcdsim2|Medicaidsim[alivei,i];
        transim2 = transim2|transfersim[alivei,i];
        malesim2 = malesim2|malesim[alivei,i];
        agesim2  = agesim2|agesim[alivei,i];
        yearsim2 = yearsim2|(yearsim[i]*ones(rn,1));
        healsim2 = healsim2|(healsim[alivei,i].==0); /* reverse sign */
        PIsim2   = PIsim2|pisim96[alivei];
        rn2      = rn2|alivei;
        HHIDsim2 = HHIDsim2|HHIDsim[alivei];
    i=i+1; endo;

    mxdist = mxsim2~mdcdsim2~transim2~agesim2~PIsim2~malesim2~healsim2~yearsim2~rn2~HHIDsim2;
    save path=^shkpath mxdist;

    "Mean values for simulated medex";; meanc(mxdist); 

retp; endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
proc(0) = sumbqd(onlymyrs);

    local asstsim, beqsim, mssim2, asim96, pisim96, agesim96, healsimh, 
          healsimw, HHIDsim, alivesim, aliveavg, diedsim, diedavg, malesim, 
          fmlsim, healsim, agesim, too_old, yearsim, cn, lasstsim, healsim2, 
          malesim2, agesim2, PIsim2, HHIDsim2, yearsim2, lasstsim2, beqsim2, 
          simyrs2, i, diedi, rn, bqdist;

    load path=^iopath asstsim, beqsim, mssim2;
    load path=^shkpath asim96, pisim96, agesim96, healsimh, healsimw, HHIDsim;

    asstsim  = reshape(asstsim,simyrs+1,_nn)';
    asstsim[.,1] = asim96;
    beqsim   = reshape(beqsim,simyrs+1,_nn)';
    mssim2   = reshape(mssim2,simyrs+1,_nn)';
    alivesim = mssim2.>0;
    aliveavg = meanc(alivesim);
    diedsim  = alivesim[.,1:simyrs] - alivesim[.,2:simyrs+1];
    diedsim  = zeros(_nn,1)~(diedsim.>0);
    diedavg  = meanc(diedsim);
    malesim  = (mssim2[.,1].==1)*ones(1,simyrs+1);
    fmlsim   = (mssim2[.,1].==2)*ones(1,simyrs+1);
    healsim  = healsimh.*malesim + healsimw.*fmlsim;
    agesim   = agesim96+seqa(0,1,simyrs)';
    too_old  = agesim.>dieage;
    agesim   = agesim.*(1-too_old)+too_old*dieage;
    yearsim  = seqa(0,1,simyrs+1)+momyr1;

    if onlymyrs==1;                                            
        cn = mmtcols-1; @ People are modeled as dying between AHEAD Survey Years @
    else;
        cn = seqa(1,1,simyrs);
    endif;

/*  There are no simulated bequests for year 1.  Drop these observations  */

    cn = cn[2:rows(cn)];

    diedsim  = diedsim[.,cn];
    diedavg  = diedavg[cn];
    beqsim   = beqsim[.,cn];
    lasstsim = asstsim[.,cn-1];
    healsim  = healsim[.,cn-1];    @ health of dead people is not observed! @
    malesim  = malesim[.,cn-1];
    agesim   = agesim[.,cn];
    yearsim  = yearsim[cn];   

    beqsim2   = {};
    lasstsim2 = {};
    healsim2  = {};
    malesim2  = {}; 
    agesim2   = {};
    yearsim2  = {};
    PIsim2    = {};
    HHIDsim2  = {};

    simyrs2 = rows(cn);

    i=1; do until i>simyrs2;
        diedi     = diedsim[.,i];
        diedi     = selif(seqa(1,1,_nn),diedi);
        rn        = rows(diedi);
        beqsim2   = beqsim2|beqsim[diedi,i];
        lasstsim2 = lasstsim2|lasstsim[diedi,i];
        healsim2  = healsim2|(healsim[diedi,i].==0); /* reverse sign */
        malesim2  = malesim2|malesim[diedi,i];
        agesim2   = agesim2|agesim[diedi,i];
        yearsim2  = yearsim2|(yearsim[i]*ones(rn,1));
        PIsim2    = PIsim2|pisim96[diedi];
        HHIDsim2  = HHIDsim2|HHIDsim[diedi];
    i=i+1; endo;

    bqdist = beqsim2~lasstsim2~agesim2~PIsim2~malesim2~healsim2~yearsim2~HHIDsim2;

    save path=^shkpath bqdist;

    "Mean values for simulated bequests";; meanc(bqdist); 

retp; endp;

