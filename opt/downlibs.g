/*  DOWNLIBS.G

  - Kernel density estimator code by Ruud Koning (1996)
  - Simplex algorithm (AMOEBA) by Bo Honore and Ekaterini Kyriazidou
*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/* KERNEL DENSITY LIBRARY ver 1.0
Author: Ruud H. Koning rhkoning@xs4all.nl.
Date: 13 July 1996
Provided without guarantee for public noncommercial use.
*/

/* kernel.src version 1.0 13 july 1996
Source file with functions of the kernel library:
    ukernel     procedure for univariate kernel estimation
    mkernel     procedure for multivariate kernel estimation
    bandw1      procedure for automatic bandwidth selection
    nw          procedure for nonparametric regression
    viewuknl    procedere for plotting the estimated density

*/

/* comment the following line out if you don't have the pgraph library
installed. You will not be able to use the procedure view_ukernel */
#include pgraph.ext

/* procedure to calculate a univariate kernel density estimate and its
derivative
usage:      {f, d, h} = ukernel(x, z, h, w, &kf);
input:      x:      T vector where density is to be estimated
            z:      n vector with observed data points
            h:      scalar bandwidth, if h<=0 bandwidth is determined by
                    procedure bandw1(z)
            w:      n vector with weights
            &kf:    pointer to weighting function
output:     f:      T vector with estimated density
            d:      T vector with estimated derivative
            h:      scaler bandwidth
*/

proc (3)=ukernel(x, z,h,w, &kf);
 local arg, i, n, f, d, k, kff, kfd, pkff,argl,argh,kffl,kffh,fl,fh,
  hl,hh;
 local kf:proc;

 /* error checking here */
 if (cols(x)>1);
  errorlog "ukernel.g: x has too many columns";
  retp(-1,-1,-1 );
 endif;
 if (cols(z)>1);
  errorlog "ukernel.g: z has too many columns";
  retp( -1,-1,-1 );
 endif;

 /* initialization */
 i = 1;
 n = rows(x);
 k = cols(x);

 /* determine bandwidth */
 if (h<=0);
  h=bandw1(z);
 endif;

 f=zeros(n,1);
 d=zeros(n,1);
 do while (i<=n);
  arg = (x[i]-z)/h;
  {kff, kfd} = kf(arg);
  f[i] = meanc(kff.*w)/h;
  d[i] = meanc(kfd.*w)/(h^2);
  i=i+1;
 endo;
 retp( f,d,h );
endp;


/* procedure to calculate a multivariate kernel density estimate and its
derivative
usage:      {f, d, h} = mkernel(x, z, h, w, &kf);
input:      x:      T x k matrix where density is to be estimated
            z:      n x k matrix with observed data points
            h:      scalar, bandwidth parameter (same bandwidth for all
                    components), if h<=0 it is determined by eq 4.14 of
                    Silverman
            w:      n vector with weights
            &kf:    pointer to weighting function
output:     f:      T vector with estimated density
            d:      T x k matrix with estimated derivatives
            h:      scalar, bandwidth used
*/

proc (3)=mkernel(x, z, h, w, &kf);
 local arg,i,n,f,d,k,kff,kfd,pkff,p,deneq0;
 local kf:proc;

 i = 1;
 n = rows(x);
 k = cols(x);
 if (h<=0);
   p=1/(k+4);
   h=(4/(k+2))^p*n^(-p);
 endif;
 f = zeros(n,1);
 d = zeros(n,k);
 do while (i<=n);
  arg = (x[i,.]-z)/h;
  {kff, kfd} = kf(arg);
  pkff = prodc(kff');
  f[i] = meanc(pkff.*w)/(h^k);
  deneq0=kff .==0;
  kff=kff+deneq0;
  d[i,.] = meanc((kfd.*pkff.*w)./kff)'/(h^(k+1));
  i = i+1;
 endo;
 retp( f,d,h );
endp;


/* computes Nadaraya-Watson kernel regression estimator

usage:  {g,d,h} = nw(x,z,y,h,&kf);
input:  x:          T-vector where regression function is evaluated
        z:          n-vector with data on independent variable
        y:          n-vector with data on dependent variable
        h:          scalar, bandwidth, if h<=0 the bandwidth is determined
                    automatically using procedure bandw1
        &kf:        pointer to univariate kernel density estimator
output: g:          T-vector with values of regression evaluated in x
        d:          T-vector with estimated derivative in x
        h:          scalar bandwidth used
*/

proc (3)=nw(x, z, y,h,&kf);
    local num, numder, denom, denomder, g, s, n, k, d,v;
    local kf: proc;

 n = rows(x);
 k = cols(x);
 g = zeros(n,1);
 d = zeros(n,1);

 if (k ne 1);
  print "error in nw.g: too many columns in x";
  retp( g, d );
 endif;
 if (cols(z) ne 1);
  print "error in nw.g: too many columns in z";
  retp( g, d );
 endif;

 {num,numder,h} = ukernel(x,z,h,y,&kf);
 {denom,denomder,h} = ukernel(x,z,h,1,&kf);
 s = (denom.==0);
 g = num./(denom+s);
 d = numder./(denom+s) - num.*denomder./(denom^2+s);
 retp( (1-s).*g, (1-s).*d,h );
endp;

/* bandw1
procedure to calculate the optimal bandwidth in kernel estimation of a density.
The optimal bandwidth is calculated according to eq. 3.31 of Silverman (1986)

usage:      h=bandw1(y);
input:      y:      n-vector whose density will be estimated;
output:     h:      scalar, optimal bandwidth choice;

*/

proc bandw1(y);
 local s, ys,n, a, iqr, qi1, qi3;

 if (cols(y)>1);
  errorlog "input error in bandw1.g: too many columns";
  retp( -1 );
 endif;
 s=sqrt(vcx(y));
 n=rows(y);
 ys=sortc(y,1);
 qi1=round(0.25*n);
 qi3=round(0.75*n);
 iqr=ys[qi3]-ys[qi1];
 retp( 0.9*minc(s|(iqr/1.34))/n^0.2 );
endp;

/* view_ukernel
procedure to plot kernel density estimate. Library pgraph must be activated.

usage: call view_ukernel(x,f,h);
input:      x:  n-vector with data points
            f:  n x k matrix (k=1 or k=3) with estimated densities for
                different bandwidths
            h:  k-vector (k=1 or k=3) with bandwidths
output:     none
globals:    all globals of the pgraph library

*/

proc (0)=view_ukernel(x,f,h);
 local data,k,xlow,xhigh,flow,fhigh,xlegend,ylegend;

 /* error checking */
 if (cols(f) ne rows(h));
  errorlog "error in viewkrnl.g: rows(f) unequal to cols(h)";
  retp();
 endif;
 /* set global variables pgraph */
 _pdate=0;

 k=cols(f);
 xlow=minc(x);
 xhigh=maxc(x);
 flow=minc(minc(f));
 fhigh=maxc(maxc(f));
 xlegend=xlow;
 ylegend=flow+0.8*(fhigh-flow);
 _plegctl=1|4|xlegend|ylegend;
 if (k==1);
  _plegstr="h=" $+ ftos(h,"%*.*lf",5,3);
  else;
  _plegstr="h=" $+ ftos(h[1],"%*.*lf",5,3) $+ "\0h="
   $+ftos(h[2],"%*.*lf",5,3)$+"\0h="$+ ftos(h[3],"%*.*lf",5,3);
 endif;
 data=sortc(x~f,1);
 xy(data[.,1],data[.,2:cols(data)]);
endp;

/* This file contains some kernel functions. All functions take an
n x k matrix u as their argument and return an n x k matrix with the
function evaluated in each point of x and an n x k matrix d with the
derivative of the function in each point of x.
    k_bw:     biweight kernel function
    k_epan:   Epanechnikov kernel
    k_gauss:  Gaussian kernel
    k_triang: triangular kernel
    k_rect:   rectangular kernel

usage:
    {f,d}=k_bw(u);

Author: Ruud H. Koning rhkoning@xs4all.nl.
Date: 13 July 1996

*/

proc (2)=k_bw(u);
 local select;
 select = abs(u).<=1;
 retp( 15/16*((1-u.^2).^2).*select, -15/4*u.*(1-u.^2).*select);
endp;

proc (2)=k_gauss(u);
    retp( pdfn(u), -u.*pdfn(u) );
endp;

proc (2)=k_epan(u);
    local c, s;
    s = abs(u).<sqrt(5);
    c = 0.75/sqrt(5);
    retp( c*(1-0.2*u.^2).*s, -0.4*c*u.*s );
endp;

proc (2)=k_rect(u);
    retp( 0.5*(abs(u).<1), 0*u );
endp;

proc (2)=k_trian(u);
    local s, a;
    a = abs(u);
    s = a.<1;
    retp( (1-a).*s, (-(u.>=0) + (u.<=0)).*s );
endp;


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/



/************* OPTIMIZATION ROUTINES  ******************************
The following are GAUSS versions of the optimization routines in
"Numerical Recipes" written by Bo Honore and Ekaterini Kyriazidou
of Northwestern University, with support from an NSF grant.
You are welcome to use and distribute them as long as you

include proper attribution to the authors.
You should know that you use the routines at your own risk.
THE ROUTINES FOLLOW:

***********************************************************/

PROC (4) = AMOEBA(p,ftol,maxsec,maxit,&fct,prnum);

/* This procedure minimizes a function called FCT using the simplex
   method. The procedure is tailored after AMOEBA in Numerical Recipes.
   We are minimizing over an NDIM-dimensional vector of parameters.
   Input is a matrix P, whose NDIM+1 rows are NDIM-dimensional vectors
   which are the vertices of the starting simplex.
   FTOL is the fractional convergence tolerance to be achieved in the function
   value.
   ALP, BET, GAM below, are parameters which define the expansions and
   contractions.

   INPUT

   p      ((ndim+1) x ndim)   starting simplex
   maxsec (1x1)               maximum number of seconds allowed
   maxit  (1x1)               maximum number of iteraions allowed
   ftol   (1x1)               fractional convergence tolerance
   fct    proc                function we want to minimize
   prnum (1x1)          print every prnnum iterations @ JBJ 5/11/98  @

   OUTPUT

   p      ((ndim+1) x ndim)   final simplex; first row is the minimizing
                  vector
   y      ((ndim+1) x 1)      vector of values of fct at final simplex;
                              first number is the value of the function at
                              the minimum
   iter   (1x1)               number of iterations taken
   tim    (1x1)               number of seconds of running time

   While running the following are printed on the screen:
   the number of the current iteration and of seconds of running time,
   the value of the function at the current simplex (transposed), and
   the current simplex (transposed).
   In the end the following are printed on the screen:
   the number of the final iteration, the number of seconds taken,
   the value of the function at the final simplex (transposed); the
   first number is the value of the function at the minimum
   the final simplex (transposed); first column is the minimizing vector.
*/

     local y,j,date1,tim,ndim,npts,pr,prr,pbar,ind,ihi,inhi,ilo,rtol,
       alp,bet,gam,ypr,yprr,i,iter,fct:proc   ;

     tim=0;
     date1=date;
     alp=1.; bet=0.5; gam=2.0;
     ndim=cols(p);
     npts=ndim+1;
     y=zeros(ndim+1,1);
     j=1;
       do while j<=npts;
       y[j,1]=fct(p[j,.]');
       j=j+1;
       endo;


     iter=0;
     begy:
      tim=ethsec(date1,date)/100;
      ind=sortind(y);
      ihi=ind[npts,1];
      inhi=ind[npts-1,1];
      ilo=ind[1,1];
      if (abs(y[ihi,1]+y[ilo,1]) > 1e-15);  /*  Added -8-4-96 (JBJ) */
        rtol=2.*abs(y[ihi,1]-y[ilo,1])/abs(y[ihi,1]+y[ilo,1]);
        else; rtol=2.*abs(y[ihi,1]-y[ilo,1])/(1e-15);
      endif;
      call monit(p,y,tim,iter,prnum);  /* Add this line to observe each iteration */
      if rtol<ftol;
      call monit(p,y,tim,iter,prnum);
      RETP(p,y,iter,tim);
      endif;
      if (iter==maxit);
      call monit(p,y,tim,iter,prnum);
      "Maximum number of iterations exceeded";
      RETP(p,y,iter,tim);
      endif;
      if (tim .ge maxsec);
      call monit(p,y,tim,iter,prnum);
      "Maximum number of seconds exceeded";
      RETP(p,y,iter,tim);
      endif;
      iter=iter+1;
      pbar=(sumc(p)-p[ihi,.]')/ndim;
      pr=(1.+alp)*pbar-alp*p[ihi,.]';
      ypr=fct(pr);
      if ypr <=y[ilo,1];
     prr=gam*pr+(1.-gam)*pbar;
     yprr=fct(prr);
     if yprr < y[ilo,1];
    p[ihi,.]=prr';
    y[ihi,1]=yprr;
     else;
    p[ihi,.]=pr';
    y[ihi,1]=ypr;
     endif;
      elseif ypr >= y[inhi,1];
     if ypr < y[ihi,1];
    p[ihi,.]=pr';
    y[ihi,1]=ypr;
     endif;
/*prr=bet*p[ihi,.]'+(1.-gam)*pbar; this is what I was given -- it seems screwed up, gam should be bet*/
     prr=bet*p[ihi,.]'+(1.-bet)*pbar;
     yprr=fct(prr);
     if yprr < y[ihi,1];
    p[ihi,.]=prr';
    y[ihi,1]=yprr;
     else;
    p=.5*(p+p[ilo,.]);
    i=1;
    do while i<=npts;
	 if i==ilo;
	  goto waydone;
	 endif;
       y[i,1]=fct(p[i,.]');
	  waydone:
       i=i+1;
    endo;
     endif;
      else;
     p[ihi,.]=pr';
     y[ihi,1]=ypr;
      endif;
      goto begy;

ENDP;

PROC (0) = monit(p,y,tim,iter,prnum);
     screen on;
     format /ro 9,7;
     if (iter/prnum == round(iter/prnum));
        "iteration: ";; iter; "hours elapsed: ";; tim/3600; ?;
        format /ro 5,6;
        "Value of function at current simplex: ";; y'; ?;
        "current simplex (transposed): ";; p'; ?;
     endif;
     format /ro 12,5;

ENDP;

/************************************************************************
Economics       is entirely responsible for the information
provided above. Please direct any comments to Alan G. Isaac at:
          E-Mail Address: aisaac@american.edu
          Phone Number:   (202) 885-3785
*/

