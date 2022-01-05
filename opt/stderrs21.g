/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
GETSE:   Compute standard errors
         Assumes PI is expressed as a rank
*/
proc(4)= getse(allparms,newmoms,xtrasst);

    local parmvec, _v, _w, vinv, mmtvec, obsvec, qntvec, pdfvec, mmttype, _m, 
          usepdfs, shortvec, gradstep, _D, DD, DDsc, DOmD, VCV, _P, _Q, Ovid, 
          stderr, tstat, tau, df, rn, rn2, rn3, zn, lmxmnsdat, lmxstddat,
          v0, vi0, W0;
    
    allparms  = allparms.*zerovec + fixvals.*(1-zerovec);
    {rhomx, fracar1} = punscale_m(allparms);/* Work with unscaled parameters */
    allparms[1:2] = rhomx|fracar1;
    fixvals[1:2]  = rhomx|fracar1;
    parmvec   = allparms[25:31];
    getparms(parmvec);                      /* Work with unscaled parameters */
    parmvec   = _delta|_beta|_nu|_omega|cfloor|phi0|K0;
    allparms[25:31] = parmvec;
    pvecse    = allparms;                /* Global Used for finding gradient */
    pscaled_p = 0;
    pscaled_m = 0;

    {_v, vinv, _W, lmxmnsdat,lmxstddat, mmtvec, obsvec, qntvec, pdfvec, mmttype}
    = getWmtx(agedat96,PIdat,asstdat,MStatdat,obsdat,mxdat,mxobsdat,
              datawgts,optwgts,xtrasst);

   /*--Update mmtvec, qntvec (and pdfvec, for consistency) to reflect model--*/
   /*---fits.  Recompute the vcv matrix and its inverse, vinv, using model---*/
   /*------------------means (but the same data as before)-------------------*/

    {v0, vi0, W0, mmtvec, obsvec, qntvec, pdfvec, mmttype} 
             = getWmtx2(allparms,xtrasst);  

    if newmoms==1;
        _v   = v0;
        vinv = vi0;
        _W   = W0;
    endif;

    _m      = rows(mmtvec);     /*  _V should be MxM, Parmvec should be px1  */
    usepdfs = mmttype.<3;                          /* quantile-based moments */
    if maxc(mmttype) > 2;
        pdfvec = pdfvec|ones(_m-sumc(usepdfs),1);
    endif;

 /*-----------Find Jacobian for model-predicted summary statistics-----------*/
 /*-----------------This should be Mxp, p = # of parameters------------------*/ 

    rn      = rows(allparms);
    colkeep = seqa(1,1,rn);                          /* Exclude fixed values */
    "Parameter Vector and Parameters Used";;
    colkeep~zerovec~allparms;?;
    colkeep = colkeep.*zerovec;                      /* This is a global     */
    colkeep = sortc(colkeep,1);
    zn      = sumc(colkeep.< 0.99);
    colkeep = colkeep[zn+1:rn];

    prnres   = 1;
    shortvec = allparms[colkeep];
    rn2      = rows(shortvec);
    rn3      = sumc(zerovec[1:24]);
    gradstep = ones(rn2-rn3,1);
    if rn3 > 0;
        gradstep = (0.1*ones(rn3,1))|gradstep;
    endif;
    gradstep = gradstep*0.025;

    "Parameters Used and Gradient Step Size";;
    colkeep~shortvec~gradstep;?;

    "Observation fraction and pdf vector ";;
    mmttype~obsvec~pdfvec;

    _D = gradp2(&modqnts,shortvec,gradstep); /* Rows are moments, columns are parameters */
    _D = (_D.*obsvec).*pdfvec;

    "D = " _D;?;    /* _D should be Mxp, where p is the number of parameters */
    if optwgts<2;
        DD   = _D'*_w*_D;                                /* DD should be pxp */
        DDsc = diagrv(eye(rn2),diag(DD)^-1); /* Improves numerical stability */
        DD   = DD*DDsc;
        DOMD = _D'*_w*_v*_w*_D*DDsc;
        VCV  = DDsc*inv(DD)*DOmD*inv(DD);
        "DD = ";; DD*inv(DDsc);
        "DOMD = ";; DOMD*inv(DDsc);                          
      /*-----Formula for Overidentification Test from Newey (JoE, 1995)------*/
        _P   = eye(_m) - _D*DDsc*inv(DD)*_D'_w;
        _Q   = _P*_v*_P';/* Use Moore-Penrose (generalized) inverse to invert */
        Ovid = mmtvec'pinv(_Q)*mmtvec;
    elseif optwgts==2; 
        DD   = _D'*vinv*_D;
        DDsc = diagrv(eye(rn2),diag(DD)^-1); /* Improves numerical stability */
        DD   = DD*DDsc;
        VCV  = DDsc*inv(DD);
        Ovid = mmtvec'*vinv*mmtvec;
    endif;

    tau  = totobs/_nn;
    VCV  = (1+tau)*VCV/totobs;
    stderr = zeros(rn,1);
    stderr[colkeep] = sqrt(diag(VCV));
    tstat= allparms./stderr;
    "standard errors and t-stats=" allparms~stderr~tstat;?;

    Ovid = totobs*Ovid/(1+tau);
    "Overidentification statistic = ";; Ovid;
    df = _m-rows(colkeep);
    "Degrees of Freedom =           ";; df;
    "p-value                        ";; cdfchic(ovid,df);

retp(VCV, stderr, tstat, Ovid); endp;

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
MODMOMS:  Calculates the vector of moments for GMM estimation
          N.B.  Parameters must be UNscaled
*/

proc(5) = modmoms(allparms);

    local mxmcoef, mxvcoef, mxcoef, parmvec, PIsim96, agesim96, asim96, 
          asstsim, conssim, medexsim, mssim2, alivesim, asstnoiz, execret, 
          simwgt96, aqntsim, aqntcnts, cprfsim, cprfcnts, mxprfsim, mxprfcnts, 
          mxmnssim, mxmnscnts, lmedexsim, lmxmnssim, lmxmnscnts, mxstdsim, 
          mxcrlsim1, mxcrlcnt1, mxcrlsim2, mxcrlcnt2, mmtvec, obsvec, qntvec, 
          pdfvec, mmttype;

    rhomx    = allparms[1];
    fracar1  = allparms[2];
    fracar1i = fracar1*(1-rhomx^2);    /* sigma of AR(1) innovations: global */           
    fracwn   = 1-fracar1;                    /* sigma of white noise: global */
    mxmcoef  = allparms[3:13];
    mxvcoef  = allparms[14:24];
    mxcoef   = mxmcoef~mxvcoef;
    ?;"Medex Coefficients";
    rhomx;;fracar1;
    mxcoef;?;

    {mnlnmxs, stdlnmxs, mnmx_pi, varmx_pi} = getmxtab(mxcoef,0.5,smplyrs);  /* Outputs are globals */

    parmvec  = allparms[25:31]; 
    getparms(parmvec);            /* Work with unscaled parameters */
    parmvec  = _delta|_beta|_nu|_omega|cfloor|phi0|K0;

    savevecs();                          /* Save input vector for C program */
    output off;

    if useMPI < 2;
        execret = exec(rulecall," ");
    elseif useMPI == 2;
        execret = exec("mpirun", rulecall);
    endif;

    output on;
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

 /* We use lmxmnsdat and lmxstddat because we are not matching logged means or std deviations. */

    {mmtvec,obsvec,qntvec,pdfvec,mmttype} 
        = makemmts(PIdat,agedat[.,1],asstdat,MStatdat,obsdat,datawgts,aqntsim,
                   mxdat,mxobsdat,mxprfsim,mxmnssim,lmxmnsdat,lmxstddat,
                   mxcrlsim1,mxcrlsim2,1);

    sttime = timerec(sttime,"Doing one function evaluation");

retp(mmtvec, obsvec, qntvec, pdfvec, mmttype); endp;

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
MODQNTS:  calculates the vector of moments for GMM estimation
*/

proc(1) = modqnts(shortvec);

    local allparm, mmtvec, obsvec, qntvec, pdfvec, mmttype;
 
    allparm = pvecse;
    allparm[colkeep] = shortvec;
    {mmtvec, obsvec, qntvec, pdfvec, mmttype} = modmoms(allparm);

retp(qntvec); endp;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

proc(8) = getWmtx2(allparms,xtrasst);

    local mmtvec, obsvec, qntvec, pdfvec, mmttype, mmtmtx, mnvals, totobs, 
          vcv, vinv, vdiag, rn, _W, Wdiag;

    savemmts = 1;
    {mmtvec, obsvec, qntvec, pdfvec, mmttype} = modmoms(allparms);
    savemmts = 0;

    load path =^datapath mmtmtx;
    mnvals   = meanc(mmtmtx);      /* This data should already be zero-mean */
    mmtmtx   = mmtmtx - mnvals';
    totobs   = rows(mmtmtx);
    vcv      = mmtmtx'*mmtmtx;
    clear mmtmtx;
    vcv      = vcv/totobs;

   /*----Take the principal diagonal and form a diagonal weighting matrix-----*/

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

    if optwgts<2;
        Wdiag = diag(_W);
        Wdiag = Wdiag + (xtrasst-1).*(mmttype.==1).*Wdiag;
        _W    = diagrv(_W,Wdiag);
    endif;   
 
/*
    format /ro 12,4;
    "Principal diagonals of V, Vinv (V^{-1}), V_diag^{-1}, and _W = ";;
    seqa(1,1,rn)~diag(vcv)~diag(vinv)~diag(vdiag)~diag(_W);?;
*/

retp(vcv,vinv,_W,mmtvec,obsvec,qntvec,pdfvec,mmttype); endp;

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/


/*
** gradp2.src
** modified:  03/21/07, JBJ
** (C) Copyright 1988-1998 by Aptech Systems, Inc.
** All Rights Reserved.
**
** This Software Product is PROPRIETARY SOURCE CODE OF APTECH
** SYSTEMS, INC.    This File Header must accompany all files using
** any portion, in whole or in part, of this Source Code.   In
** addition, the right to create such files is strictly limited by
** Section 2.A. of the GAUSS Applications License Agreement
** accompanying this Software Product.
**
** If you wish to distribute any portion of the proprietary Source
** Code, in whole or in part, you must first obtain written
** permission from Aptech Systems.
**
**> gradp
**
**  Purpose:    Computes the gradient vector or matrix (Jacobian) of a
**              vector-valued function that has been defined in a procedure.
**              Single-sided (forward difference) gradients are computed.
**
**  Format:     g = gradp(&f,x0,step);
**
**  Input:      f    scalar, procedure pointer to a vector-valued function:
**
**                                          f:Kx1 -> Nx1
**
**                   It is acceptable for f(x) to have been defined in terms of
**                   global arguments in addition to x, and thus f can return
**                   an Nx1 vector:
**
**                        proc f(x);
**                           retp( exp(x*b) );
**                        endp;
**
**              x0    Kx1 vector of points at which to compute gradient.
**
**              step  Kx1 vector of step size adjustments.  Gradp uses 1e-8.
**
**  Output:     g     NxK matrix containing the gradients of f with respect
**                    to the variable x at x0.
**
**  Remarks:    gradp will return a row for every row that is returned by f.
**              For instance, if f returns a 1x1 result, then gradp will
**              return a 1xK row vector. This allows the same function to be 
**              used where N is the number of rows in the result returned by f.
**              Thus, for instance, gradp can be used to compute the
**              Jacobian matrix of a set of equations.
**
**  Example:    proc myfunc(x);
**                 retp( x .* 2 .* exp( x .* x ./ 3 ));
**              endp;
**
**              x0 = { 2.5, 3.0, 3.5 };
**              y = gradp(&myfunc,x0);
**
**                           82.98901842    0.00000000    0.00000000
**                  y =       0.00000000  281.19752975    0.00000000
**                            0.00000000    0.00000000 1087.95414117
**   
**              It is a 3x3 matrix because we are passing it 3 arguments and
**              myfunc returns 3 results when we do that.  The off-diagonals
**              are zeros because the cross-derivatives of 3 arguments are 0.
**
**  Globals:    None
**
**  See Also:   hessp
*/

proc 1 = gradp2(f,x0,step);
    local f:proc;
    local n,k,grdd,dh,ax0,xdh,arg,dax0,i,f0;

    /* check for complex input */
    if iscplx(x0);
        if hasimag(x0);
            errorlog "ERROR: Not implemented for complex matrices.";
            end;
        else;
            x0 = real(x0);
        endif;
    endif;

    f0 = f(x0);
    n = rows(f0);
    k = rows(x0);
    grdd = zeros(n,k);

/* Computation of stepsize (dh) for gradient */

    ax0 = abs(x0);
    if x0 /= 0;
        dax0 = x0./ax0;
    else;
        dax0 = 1;
    endif;
    dh  = step.*(maxc((ax0~(1e-2)*ones(rows(x0),1))').*dax0);
    xdh = x0+dh;
    dh  = xdh-x0;    /* This increases precision slightly */
    arg = diagrv(reshape(x0,k,k)',xdh);

    i = 1;
    do until i > k;
        grdd[.,i] = f(arg[.,i]);
        i = i+1;
    endo;

    grdd = (grdd-f0)./(dh');

    retp(grdd);
endp;

