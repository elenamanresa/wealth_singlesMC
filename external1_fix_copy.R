######################   SUBPROGRAMS CODE FROM De NARDI, FRENCH, AND JONES #########################

#------------------------------Define Sub-programs--------------------------------#

logitsr <- function(x){
  sqrt(exp(x)/(1+exp(x)));
}

logit <- function(x){
  exp(x)/(1+exp(x));
}

logitrv <- function(p){
  log(p/(1-p));
}

#GET1YR:  Converts two-year transition probs into 1-year probs. Using formula by O. Nartova
#gg = pr(h_t+2=good|h_t=good); bb = pr(h_t+2=bad|h_t=bad)
#g = pr(h_t+1=good|h_t=good); b = pr(h_t+1=bad|h_t=bad)

get1yr <- function(gg,bb){
  big_A = 2-gg-bb;
  big_B = 2*(gg-1);
  big_C = 1-bb-gg+(bb^2);
  
  b    = -big_B + sqrt(big_B^2 - 4*big_A*big_C);
  b    = b/(2*big_A);
  g    = ( (1-bb) - b*(1-b) )/(1-b);
  output = list("b"=b, "g"=g)
}

#GETAGES:  gets ages right. Subfunction for "fixobs", get the position (row #) of bornage and dieage. 

getages <- function(bornage,dieage,dat){
  
  k = match(bornage, dat[,1]);
  if(is.na(k)){
    k=1;
  }
  bage2 = dat[k,1];
  ib = 1 + max((bage2-bornage),0);
  j = match(dieage, dat[,1]);
  if(is.na(j)){
    j = nrow(dat)
  }
  output = list("k"=k,"j"=j,"ib"=ib)
  return(output)
}

#FIXOBS:  Adjusts data to have a common set of years. Adds filler rows (with missing values) if data doesn't cover [bornage,dieage].
#Make the data contain years from bornage to dieage. 

fixobs <- function(dat,bornage,dieage,k,j,lowmiss,himiss){
  cn    = ncol(dat);
  bage2 = dat[k,1];
  dage2 = dat[j,1];
  dat   = dat[k:j,];
  
  if(bage2 > bornage){
    rn  = bage2-bornage;
    names(dat) = head(c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),cn);
    dat = rbind(lowmiss*matrix(1, rn, cn),dat);
    dat[1:rn,1] = seq(bornage,bornage+rn-1,by=1);
  }
  
  rn2 = j-k+bage2-bornage+1;                      # Number of rows in dat
  
  if(dage2 < dieage){
    rn  = dieage-dage2;
    names(dat) = head(c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),cn);
    dat = rbind(dat,himiss*matrix(1,rn,cn));
    dat[(rn2+1):TR,1] = seq(dage2+1,dage2+rn,by = 1);
  }
  return(dat)
}

#getprofs: Adjust and prepare the following datasets:
#a. Mortality and Surrival Rates: morts, morts_pi.
#b. Health Transition Probabilities: healdats, heal_pi.
#c. Income Profiles: y_pi, yprof.
#d. Observe date-specific ROR (rate of return) shocks: rorshk. 

getprofs <- function(){
  onevec = rep(1,TR-1);
  pctle  = 0.4;
  #----------------------Mortality and Survival Rates----------------------#
  if(mortdif==0){
    datstr = paste0(datapath,"deathprof_nodiff.out");
  }else if(mortdif==1){
    datstr = paste0(datapath,"deathprof.out");
  }
  data <- read.delim( datstr, header=FALSE)
  data[,1] <- data[,1]-1     # Adjust for Eric's Indexing and 2-year interval 
  
  # Order: age, b_const, b_badhealth b_male, b_PI, b_PI2,                
  #       there are separate values of b_* for each age                
  # Pr(alive_t+1|..) = sqrt( logit(b_const[age] + b_badhealh[age] + ...) ) 
  datah <- cbind(data[,1],data[,2]+data[,4],data[,3],data[,5:6])
  
  outcome = getages(bornage,dieage-1,datah);
  k = outcome$k;
  j = outcome$j;
  ib = outcome$ib;
  
  datah    = fixobs(datah,bornage,dieage,k,j,0,-1e10);
  mort_pi  = datah[,4:5];
  
  dataw    = cbind(data[,1],data[,2],data[,3],data[,5:6]);
  dataw    = fixobs(dataw,bornage,dieage,k,j,0,-1e10);
  mort_pi  = cbind(mort_pi,dataw[,4:5]);
  
  pieffcts = (mort_pi[,1]*pctle+mort_pi[,2]*(pctle^2))*matrix(1,nrow(mort_pi),2);
  pieffcts = cbind(pieffcts,(mort_pi[,3]*pctle+mort_pi[,4]*(pctle^2))*matrix(1,nrow(mort_pi),2));
  
  mortshg  = datah[,2];
  mortshb  = datah[,2]+datah[,3];
  mortswg  = dataw[,2];
  mortswb  = dataw[,2]+dataw[,3];
  morts    = cbind(mortshb,mortshg,mortswb,mortswg);
  
  
  #---------------------Health Transition Probabilities---------------------#
  datstr   = paste0(datapath,"healthprof.out");
  data <- read.delim( datstr, header=FALSE)
  
  data[,1]= data[,1]-1; # Adjust for Eric's Indexing and 2-year interval
  
  # Order: age, b_const, b_badhealth b_male, b_PI, b_PI2,                
  #        there are separate values of b_* for each age                  
  # Pr(h_t+1=bad|..) = logit(b_const[age] + b_badhealh[age] + ...)       
  
  datah    = cbind(data[,1],(data[,2]+data[,4]),data[,3],data[,5:6]);
  
  outcome = getages(bornage,dieage-1,datah);
  k = outcome$k;
  j = outcome$j;
  ib = outcome$ib;
  
  datah    = fixobs(datah,bornage,dieage,k,j,0,0);
  heal_pi  = datah[,4:5];
  dataw    = cbind(data[,1],data[,2],data[,3],data[,5:6]);
  dataw    = fixobs(dataw,bornage,dieage,k,j,0,0);
  heal_pi  = cbind(heal_pi,dataw[,4:5]);
  
  pieffcts = (heal_pi[,1]*pctle+heal_pi[,2]*(pctle^2))*matrix(1,nrow(heal_pi),2);
  pieffcts = cbind(pieffcts,(heal_pi[,3]*pctle+heal_pi[,4]*(pctle^2))*matrix(1,nrow(heal_pi),2));
  
  healshg  = datah[,2];
  healshb  = datah[,2]+datah[,3];
  healswg  = dataw[,2];
  healswb  = dataw[,2]+dataw[,3];
  healdats = cbind(healshb,healshg,healswb,healswg);
  ggbb     = logit(healdats+pieffcts);
  gg       = 1 - ggbb[,c(2,4)];
  bb       = ggbb[,c(1,3)];
  
  g = get1yr(gg,bb)$g;
  b = get1yr(gg,bb)$b;
  b= 1 - b; 
  
  
  #-----------------------------Income Profiles-----------------------------#
  datstr   = paste0(datapath,"incprof.out");
  data <- read.delim( datstr, header=FALSE)
  
  # Order: age, b_const, b_badhealth b_male, b_PI, b_PI2,                  
  #       there are separate values of b_* for each age                   
  # y = exp(b_const[age] + b_badhealh[age] + ... )                        
  
  datah    = cbind(data[,1],(data[,2]+data[,4]),data[,3],data[,5:6]);
  k = getages(bornage,dieage,datah)$k;
  j = getages(bornage,dieage,datah)$j;
  ib = getages(bornage,dieage,datah)$ib;
  
  datah    = fixobs(datah,bornage,dieage,k,j,0,-1e10);
  y_pi     = datah[,4:5];
  dataw    = cbind(data[,1],data[,2],data[,3],data[,5:6]);
  dataw    = fixobs(dataw,bornage,dieage,k,j,0,-1e10);
  y_pi     = cbind(y_pi,dataw[,4:5]);
  
  yprofsh  = datah[,2];
  yprofsw  = dataw[,2];
  yprof    = cbind(yprofsh,yprofsw);
  pieffcts = y_pi[,1]*pctle+y_pi[,2]*(pctle^2);
  pieffcts = cbind(pieffcts,(y_pi[,3]*pctle+y_pi[,4]*(pctle^2)));
  
  #------------------Get observed date-specific ROR shocks-------------------#
  
  datstr   = paste0(datapath,"wlthshk8.txt")  
  wlthshk <- read.table(datstr, quote="\"", comment.char="")
  
  outcome = getages(momyr1,momyr2+1,wlthshk);
  k = outcome$k;
  j = outcome$j;
  ib = outcome$ib;
  wlthshk  = fixobs(wlthshk,momyr1,momyr2+1,k,j,0,0);
  rorshk   = wlthshk[,1+rshktype];
  
  if(rshktype==0){
    rorshk = 0*rorshk;
  }else if(rshktype>2){
    rorshk = rorshk - mu_r;
  }
  
  output = list("morts"=morts, "healdats"=healdats, "yprof"=yprof, "mort_pi"=mort_pi, 
                "heal_pi"=heal_pi, "y_pi"=y_pi, "rorshk"=rorshk)
}


getchrt <- function(data,cohorts_j){
  rn            = length(data);
  chrtnum_j     = length(cohorts_j)-1;
  chrtcnts      = matrix(0,chrtnum_j,4);
  
  chrtcnts[,1] = seq(1,chrtnum_j, by=1);
  chrtcnts[,2] = cohorts_j[1:chrtnum_j]+1;
  chrtcnts[,3] = cohorts_j[2:(chrtnum_j+1)];
  
  
  #-----------------------Remove Missing Observations-----------------------#
  data[is.na(data)] <- mvcode
  data  = cbind(data,seq(1,rn,by=1));
  srtddata  = data[order(data[,1]), ];
  k  = (srtddata[,1]==mvcode);
  k  = t(k)%*%k; 
  
  if(k > 0){
    srtd2  = srtddata[1:k,];
  } else if(k==0){
    srtd2 = srtddata[FALSE,];
  }
  
  srtddata = srtddata[(k+1):rn,];
  rn1      = rn-k;                  # Number of non-missing observations 
  chrttype = rep(1,rn1);
  
  cno      = 0;
  
  agevec   = srtddata[,1] - ageshft;           # Use age in first wave
  
  h = chrtnum_j;
  
  while (h>=1) {
    cn  = (agevec>cohorts_j[h]); #Note that cohorts has chrtnum+1 elements 
    cn  = t(cn)%*%cn;
    if(cn==cno){
      h=h-1; 
      next;
    }
    
    rn2 = cn-cno;     # Counts for interval (age_h,age_h+1) 
    rn3 = rn1 - cn;
    chrttype[(rn3+1):(rn3+rn2)] = h*rep(1,rn2);  # observations in (age_h,age_h+1) 
    cno = cn;
    h = h-1;          # End loop through quantiles
  }
  
  #---------------------------Zero Out Outliers----------------------------#
  chrttype = chrttype*(1-(agevec>cohorts_j[chrtnum_j+1]));
  chrttype = chrttype*(1-(agevec<=cohorts_j[1]));
  
  h=1; 
  
  while(h <= chrtnum_j){
    cn = (chrttype == h);
    chrtcnts[h,4] = t(cn)%*%cn;
    h = h+1;
  }
  
  #--------------------Add in Markers for Missing Data---------------------#
  if(k>0){
    chrttype = c(rep(0,k),chrttype);
  }
  
  data     = rbind(srtd2,srtddata);
  data     = cbind(data,chrttype);
  data     = data[order(data[,2]), ]
  chrttype = data[,3];
  
  output = list("chrtcnts"=chrtcnts, "chrttype"=chrttype)
  
}

#getdata: Cut the large raw dataset into several small datasets including: agedat, PIdat, asstdat, beqdat, MStatdat, 
#obsdat, agedat96, avgage96 (average age for each cohort), HHIDdat. Also have these datasets prepared and cleaned. 

getdata <- function(){
  
  datstr   = paste0(datapath,"wlthmat13b.out");
  wlthdat <- read.delim( datstr, header=FALSE)
  
  if(ageshft<2){
    ashft2=0;
  }else if(ageshft==2){
    ashft2=1;
  }else if(ageshft==4){
    ashft2=2;
  }
  
  
  wlthdat <- subset(wlthdat, get(paste0("V",76+ashft2)) < dieage & get(paste0("V",76+ashft2)) > bornage)
  
  
  totobs   = nrow(wlthdat);   #Number of observations
  
  # household id
  HHIDdat  = wlthdat[,1]; 
  
  # age data
  agedat   = wlthdat[,(76+ashft2):82];
  # mortaility data
  MStatdat = wlthdat[,(6+ashft2):12];
  MStatdat = MStatdat*(agedat<=dieage); # Kill if too old for model
  
  # asset data
  asstdat  = wlthdat[,(20+ashft2):26];
  toobig   = asstdat>asstmax;
  asstdat  = asstdat*(1-toobig) + toobig*asstmax;
  
  #?#  # what is this??? observed data?
  obsdat   = (wlthdat[,(27+ashft2):33]==0);
  obsdat   = obsdat*(MStatdat>0);
  
  #?# # what is this?? bequests?? why is it all 0?
  beqdat   = obsdat*0;
  
  # medical expenditure
  mxdat    = wlthdat[,(34+ashft2):40];
  mxobsdat = 1-(1-(wlthdat[,(41+ashft2):47]==0));
  
  #?# # health data???
  hsdat    = wlthdat[,(62+ashft2):68];
  hsobsdat = (wlthdat[,(69+ashft2):75]==0);
  hsobsdat = hsobsdat*(MStatdat>0);
  
  # Permanent Income quantile. How was this computed?
  if(wgtddata==0){
    PIdat    = wlthdat[,5];
    datawgts = rep(1,totobs);
  }else if(wgtddata==1){
    PIdat    = wlthdat[,4];
    datawgts = wlthdat[,2];
    datawgts = datawgts*totobs/sum(datawgts);
  }
  
  MStatdat = MStatdat[,1:mmtyrs];
  asstdat  = asstdat[,1:mmtyrs];
  agedat   = agedat[,1:mmtyrs];
  obsdat   = obsdat[,1:mmtyrs];
  beqdat   = beqdat[,1:mmtyrs];
  mxdat    = mxdat[,1:mmtyrs];
  mxobsdat = mxobsdat[,1:mmtyrs];
  hsdat    = hsdat[,1:mmtyrs];
  hsobsdat = hsobsdat[,1:mmtyrs];
  agedat96 = agedat[,1];
  
  outcome = getchrt(agedat[,1],cohorts);
  chrtcnts = outcome$chrtcnts;
  chrttype = outcome$chrttype;
  
  avgage96 = rep(0,chrtnum);
  
  i = 1;
  while(i <= chrtnum){
    avgage96[i] = (t(chrttype==i)%*%agedat[,1])/chrtcnts[i,4];
    i = i+1
  }
  
  avgage96 = round(avgage96);
  
  isim94 = wlthdat[,c(76,(5-wgtddata),6,62,62,69)];   #Use for life exp. calculations 
  
  isim94 = isim94[isim94[,6]==0,1:5]
  
  write.table(isim94, file = paste0(shkpath,"isim94.csv"), sep=",", row.names = FALSE, col.names = FALSE)
  
  output = list("agedat"=agedat,"PIdat" = PIdat, "asstdat"= asstdat,"beqdat" = beqdat,"MStatdat" = MStatdat,"obsdat" = obsdat,
                "mxdat" = mxdat,"mxobsdat" = mxobsdat,"hsdat" = hsdat,"hsobsdat" = hsobsdat,"totobs" = totobs,
                "datawgts" = datawgts,"agedat96" = agedat96,"avgage96" = avgage96,"HHIDdat" = HHIDdat)
}



#---------------Generate the health and mortality shocks---------------#

#getninc: Adjust income considering tax. Obtain net after-tax income (labor income+assets)

getninc <- function(assets, labinc){
  
  totror = mu_r + swchROR*rorshk[1];
  totinc = labinc+assets*totror;
  
  taxBrk = c(0,6250,40200,68400,93950,148250,284700,1e15); 
  taxMar = c(0.0765,0.2616,0.4119,0.3499,0.3834,0.4360,0.4761); 
  taxdim = length(taxMar);
  
  tottax = 0*totinc;
  blwBrk = tottax;
  
  i = 1;
  while (i<=taxdim) {
    
    Brk_low = taxBrk[i];
    Brk_hi  = taxBrk[i+1];
    blwBrk  = totinc<Brk_low;
    inBrk   = (1-blwBrk)*(totinc<=Brk_hi);
    abvBrk  = totinc>Brk_hi;
    inctax  = inBrk*(totinc-Brk_low) + abvBrk*(Brk_hi-Brk_low);
    inctax  = inctax*taxMar[i];
    tottax  = tottax+inctax;
    
    i=i+1;
    
  }
  
  netinc = totinc - swchTax*tottax;
  
  return(netinc)
}


# INITDIST: This gives the initial distribution of assets, wages, health, permanent income, and persistent health costs.
# Generate simulated data (empirical draws) of assets, wages, health, permanent income, and persistent health costs. 
#(Add noises to assets and cohort distribution. Though noises are 0 here.)
#output datasets are: asim96, incsim96, mxsim96, cohsim96, mstat96, pisim96, agesim96, asstnoiz, hsimh96, hsimw96,
#msdat, hhdat, hwdat, hmissdat, simwgt96, HHIDsim.

initdist <- function(isimage){
  
  ageindx  = rep(1,totobs);
  
  if(simtype==2){
    ageindx  = (agedat[,1]>=(isimage[1]+ageshft))*(agedat[,1]<=(isimage[2]+ageshft));
  }
  
  gotobs  = (obsdat[,1]*mxobsdat[,1]*hsobsdat[,1])*ageindx;
  num0 = seq(1,totobs, by =1)
  num0 = as.vector(subset( cbind(num0, gotobs), gotobs == 1, select=c(num0))); # restrict ourselves to non-missing obs 
  rn  = length(num0);
  
  num      = runif(nn, min = 0, max = 1)*rn;            # random draws of indices 
  num      = floor(num)+onesim;
  num      = num0[num];
  
  #--------Draw assets, wages and health from empirical distribution--------#
  agesim96 = agedat[num,1];
  pisim96  = PIdat[num];
  HHIDsim  = HHIDdat[num];
  mstat96  = MStatdat[num,1];
  gotobs2  = hsobsdat[num,1]*(mstat96>0);
  goth96   = gotobs2*((mstat96==1)+(mstat96==3));
  gotw96   = gotobs2*((mstat96==2)+(mstat96==3));
  asim96   = asstdat[num,1];
  asim96   = asim96 + (asim96>asstmax)*(asstmax-asim96);
  mxsim96  = mxdat[num,1];
  hsimh96  = 1-hsdat[num,1];
  hsimw96  = 1-hsdat[num,1];
  simwgt96 = datawgts[num];
  
  meanhsh  = mean(hsimh96*goth96)/mean(goth96);
  hstemp   = runif(nn, min = 0, max = 1)<meanhsh;
  hsimh96  = hsimh96*goth96 + hstemp*(1-goth96); # should be redundant 
  hsimh96  = as.vector(hsimh96)
  meanhsw  = mean(hsimw96*gotw96)/mean(gotw96);     
  hstemp   = runif(nn, min = 0, max = 1)<meanhsw;
  hsimw96  = hsimw96*gotw96 + hstemp*(1-gotw96); # should be redundant 
  hsimw96  = as.vector(hsimw96)
  
  msdat    = MStatdat[num,1:mmtyrs];
  hhdat    = 1-hsdat[num,1:mmtyrs];
  hwdat    = hhdat;
  hmissdat = msdat<1; # Identify observations with missing data or MS = 0 
  hmissdat = ((1-hsobsdat[num,1:mmtyrs])+hmissdat)>0;
  #hmissdat = (hsobsdat[num,1:mmtyrs]+hmissdat)>0;
  hmissdat = t(apply(t(hmissdat), 2, cumsum))
  hmissdat = hmissdat>0;
  hmissdat = hmissdat*1
  
  
  #---Now alter asset and coh distribution to allow for measurement error---#
  #--------Model:  coh = a + v, v = zero mean error, orthogonal to a--------#
  #-----------------a^ = a + u, u orthogonal to a, v------------------------#
  
  asstsign = -1+2*(asim96>=0);
  oldasst  = abs(asim96)<1;                         # recode zero assets 
  oldasst  = asim96*(1-oldasst) + oldasst*asstsign;
  lnabasst = log(abs(oldasst));
  varasst  = var(lnabasst) - assterr;
  bta      = varasst/(varasst+assterr);
  asim96   = lnabasst*bta + (1-bta)*mean(lnabasst);
  projerr  = sqrt(bta*assterr)*rnorm(nn, mean = 0, sd = 1); # Now add in omitted "true" variation 
  asim96   = asim96 + projerr;
  asim96   = asstsign*exp(asim96);
  ai       = agesim96-bornage+1;
  yPI_96   = y_pi[ai,1:2]*(mstat96==1) + y_pi[ai,3:4]*(mstat96==2);
  ln_inc   = yprof[ai,1]*(mstat96==1) + yprof[ai,2]*(mstat96==2) + yPI_96[,1]*pisim96+ yPI_96[,2]*(pisim96^2);
  incsim96 = exp(ln_inc);
  
  incsim96 = getninc(asim96,incsim96);
  cohsim96 = asim96 +incsim96 - mxsim96;
  
  asstnoiz = exp(sqrt(assterr)*matrix(rnorm(nn*simyrs, mean = 0, sd = 1),nn,simyrs));
  asstnoiz[,1] = oldasst/asim96;
  
  if(simtype==1){
    
    datlist = c("asim96", "incsim96", "mxsim96", "cohsim96", "mstat96", 
                "pisim96", "agesim96", "asstnoiz", "hsimh96", "hsimw96", 
                "msdat", "hhdat", "hwdat", "hmissdat", "simwgt96", "HHIDsim");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      
      write.table(get(dat), file = paste0(shkpath,dat,".csv"), sep=",", row.names=FALSE, col.names = FALSE);
    }
    
  }else if(simtype==2){
    
    asimx   = asim96;    incsim96x = incsim96;  mxsim96x = mxsim96;
    cohsimx = cohsim96;  mstatx    = mstat96;   pisimx   = pisim96;
    agesimx = agesim96;  astnoizx  = asstnoiz;  hsimhx   = hsimh96;   
    hsimwx  = hsimw96;   mstatdx   = mstatdat;  hhdx     = hhdat;
    hwdx    = hwdat;     hmissdx   = hmissdat;  simwgtx  = simwgt96; 
    HHIDx   = HHIDsim;
    
    datlist = c("asimx", "incsim96x", "mxsim96x", "cohsimx", "mstatx", 
                "pisimx", "agesimx", "astnoizx", "hsimhx", "hsimwx", 
                "mstatdx", "hhdx", "hwdx", "hmissdx", "simwgtx", "HHIDx");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      
      write.table(get(dat), file = paste0(shkpath,dat,".csv"),sep=",",row.names=FALSE, col.names = FALSE);
    }
  }
}

#gethms: Get the 1-yr health transition probabilities for husbands and wifes. Specify the observation status: 
# mstat    0 => defunct   1 => husband  2 => wife  3 => couple 

gethms <- function(ai2, mstatsim, healsimh, healsimw, pisim96, heal_pi, mort_pi, 
                   morts, healdats, nn){
  
  mortshb = morts[,1];     mortshg = morts[,2]; 
  mortswb = morts[,3];     mortswg = morts[,4];
  
  healshb = healdats[,1];  healshg = healdats[,2]; 
  healswb = healdats[,3];  healswg = healdats[,4];
  
  healshk = matrix(runif(nn*2, min = 0, max = 1),nn,2);         # for appropriate age 
  too_old = ai2>TR;
  ai2     = ai2*(1-too_old)+too_old*TR;
  ai2     = as.matrix(ai2)
  
  # All simulation vectors are for a single year                           
  # Note:  1 denotes good health, while transition probs are for bad health  
  
  pieffct = heal_pi[ai2,1]*pisim96 + heal_pi[ai2,2]*(pisim96^2); # Men 
  bhshb   = logit(healshb[ai2] + pieffct);
  bhshg   = logit(healshg[ai2] + pieffct);
  b = get1yr(1-bhshg,bhshb)$b;
  g = get1yr(1-bhshg,bhshb)$g;
  bhshb = b;
  bhshg = g;
  bhshg   = 1 - bhshg;
  
  pieffct = heal_pi[ai2,3]*pisim96 + heal_pi[ai2,4]*(pisim96^2); # Women 
  bhswb   = logit(healswb[ai2] + pieffct);
  bhswg   = logit(healswg[ai2] + pieffct);
  b = get1yr(1-bhswg,bhswb)$b;
  g = get1yr(1-bhswg,bhswb)$g;
  bhswb = b;
  bhswg = g;
  bhswg   = 1 - bhswg;
  
  # Note:  Having a > inequality is essential 
  
  hsimpoh = (healsimh==0)*bhshb + (healsimh==1)*bhshg;
  hsimpoh = healshk[,1]>hsimpoh;
  
  hsimpow = (healsimw==0)*bhswb + (healsimw==1)*bhswg;
  hsimpow = healshk[,2]>hsimpow;
  
  mstatshk= matrix(runif(nn*2, min = 0, max = 1),nn,2);
  pieffct = mort_pi[ai2,1]*pisim96 + mort_pi[ai2,2]*(pisim96^2); # Men   
  survshg = logitsr( mortshg[ai2] + pieffct );
  survshb = logitsr( mortshb[ai2] + pieffct ); 
  
  pieffct = mort_pi[ai2,3]*pisim96 + mort_pi[ai2,4]*(pisim96^2); # Women 
  survswg = logitsr( mortswg[ai2] + pieffct );
  survswb = logitsr( mortswb[ai2] + pieffct );   
  
  mstpo   = 1*(mstatsim==1)*             # husband at time t
    ( healsimh*(mstatshk[,1]<survshg)  +
        (1-healsimh)*(mstatshk[,1]<survshb)   ) +
    2*(mstatsim==2)*               # wife at time t 
    ( healsimw*(mstatshk[,2]<survswg)  +
        (1-healsimw)*(mstatshk[,2]<survswb)   );       
  
  too_old = (ai2+1)>TR;
  mstpo   = mstpo*(1-too_old);   
  
  output = list("hsimpoh" = hsimpoh, "hsimpow" = hsimpow, "mstpo" = mstpo)
}


#INITSIM:  Simulate sequences of health, health cost and demographic shocks
#(get simulated age, health probabilities and observation status.)
# Here data are sorted by year, with different ages alive each year


initsim <- function(morts, mort_pi, healdats, heal_pi){
  
  if(simtype==1){
    
    datlist = c("hsimh96", "hsimw96", "mxsim96", "mstat96", "pisim96", "agesim96");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
      assign(dat,data) ;
    }
  } else if(simtype==2){
    
    datlist = c("hsimhx", "hsimwx", "mxsim96x", "mstatx", "pisimx", "agesimx");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
      assign(dat,data) ;
    }
    
    hsimh96 = hsimhx;    hsimw96 = hsimwx;  mstat96  = mstatx[,1];
    mxsim96 = mxsim96x;  pisim96 = pisimx;  agesim96 = agesimx[,1];
    rm(hsimhx, hsimwx, mstatx, mxsim96x, pisimx, agesimx);
  }
  
  mstatsim = mstat96;        # 0=>defunct, 1=>husband, 2=>wife, 3=>couple 
  agesim   = agesim96;
  
  healsimh = hsimh96;                     # 0=>bad health, 1=>good health 
  healsimw = hsimw96;
  
  ageindx  = agesim96-bornage; # Offset used to get age-appropriate probs
  
  ii = 1;
  while (ii <= simyrs) {
    
    ai2 = ageindx+ii;                        # for appropriate age 
    
    outcome = gethms(ai2, mstatsim[,ii], healsimh[,ii], healsimw[,ii], pisim96, heal_pi, mort_pi, morts, healdats,nn);
    hsimpoh = outcome$hsimpoh
    hsimpow = outcome$hsimpow
    mstpo = outcome$mstpo
    
    healsimh = cbind(healsimh,hsimpoh);
    healsimw = cbind(healsimw,hsimpow);
    mstatsim = cbind(mstatsim,mstpo);
    agesim   = cbind(agesim,(agesim[,ii]+(mstpo>0)));
    ii=ii+1;
  }
  
  agesim = agesim*(mstatsim>0);
  
  if(simtype==1){
    
    datlist = c( "agesim", "healsimh", "healsimw", "mstatsim");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      
      write.table(get(dat), file = paste0(shkpath,dat,".csv"),sep=",", row.names=FALSE, col.names = FALSE);
    }
  }else if(simtype==2){
    agesimx = agesim;    mstatx  = mstatsim;  
    hlsimhx = healsimh;  hlsimwx = healsimw;
    
    datlist = c( "agesimx", "hlsimhx", "hlsimwx", "mstatx");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      
      write.table(get(dat), file = paste0(shkpath,dat,".csv"), sep=",",row.names=FALSE,col.names = FALSE);
    }
  }
}


#Derive the health probability between two years. Given one year ahead and later probabilities.

onestep <- function(g_t, g_tp1, gg_t, b_t, bb_t){
  
  g1 = g_t*g_tp1/gg_t;           # Pr(h_t+1 = good | h_t = good, h_t+2 = good) 
  g2 = g_t*(1-g_tp1)/(1-gg_t);   #Pr(h_t+1 = good | h_t = good, h_t+2 = bad) 
  g3 = (1-b_t)*g_tp1/(1-bb_t);   # Pr(h_t+1 = good | h_t = bad, h_t+2 = good) 
  g4 = (1-b_t)*(1-g_tp1)/bb_t;   #Pr(h_t+1 = good | h_t = bad, h_t+2 = bad) 
  
  output = list("g1" = g1, "g2" = g2, "g3"= g3,"g4" = g4)
}


#hsimpute:  Impute health status probabilities for years between AHEAD waves
#Here data are sorted by year, with different ages alive each year
#Note:  if health status probabilities vary by marital status, this
#imputation procedure will NOT work.
#output datasets:  "B_gg_sh", "B_gb_sh", "B_bg_sh", "B_bb_sh", "B_gg_sw", "B_gb_sw", "B_bg_sw", "B_bb_sw",
#"g_sh", "b_sh", "g_sw", "b_sw".

hsimpute <- function(healdats, heal_pi, nn){
  
  healshb  = healdats[,1];  healshg = healdats[,2]; 
  healswb  = healdats[,3];  healswg = healdats[,4];
  
  datlist = c("pisim96", "agesim96", "mstat96");
  
  for (i in 1:length(datlist)) {
    
    dat = datlist[i];
    data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
    assign(dat,data) ;
  }
  
  ageindx  = agesim96-bornage; # Offset used to get age-appropriate probs 
  
  
  #First Step:  get one- and two-year-ahead transition probs                
  #g_ij  = Pr(h_t+1 = good | h_t = good, type ij)
  #gg_ij = Pr(h_t+2 = good | h_t = good, type ij)
  #b_ij  = Pr(b_t+1 = bad | h_t = bad, type ij)
  #bb_ij = Pr(h_t+2 = bad | h_t = bad, type ij)
  
  g_sh  <- vector();
  g_sw  <- vector();
  b_sh  <- vector();
  b_sw  <- vector();
  
  #Second Step:  Do one-year-imputations for all years 
  #Use formulae in handout to get
  #B_ef_ij = Pr(h_t+1 = bad | h_t = e, h_t+2 = f, type ij)  
  #Add dummy column to convert results to
  #B_ef_ij = Pr(h_t = bad | h_t-1 = e, h_t+1 = f, type ij)
  
  B_gg_sh = rep(-1,nn);  B_gb_sh = rep(-1,nn);  B_bb_sh = rep(-1,nn);  B_bg_sh = rep(-1,nn);
  B_gg_sw = rep(-1,nn);  B_gb_sw = rep(-1,nn);  B_bb_sw = rep(-1,nn);  B_bg_sw = rep(-1,nn);
  
  ii = 1
  while (ii<=simyrs) {
    
    ai2     = ageindx+ii;         # Use age-appropriate probabilities 
    too_old = ai2>TR;
    ai2     = ai2*(1-too_old)+too_old*TR;
    
    #  Note:  Transition probs are for bad health 
    
    pieffct = heal_pi[as.matrix(ai2),1]*pisim96 + heal_pi[as.matrix(ai2),2]*(pisim96^2); #Men 
    bb_t    = logit(healshb[as.matrix(ai2)] + pieffct);
    gg_t    = 1 - logit(healshg[as.matrix(ai2)] + pieffct);
    
    g_t = get1yr(gg_t,bb_t)$g;
    b_t = get1yr(gg_t,bb_t)$b;
    
    if(ii<simyrs){
      g1 = onestep(g_t, g_t, gg_t, b_t, bb_t)$g1;
      g2 = onestep(g_t, g_t, gg_t, b_t, bb_t)$g2;
      g3 = onestep(g_t, g_t, gg_t, b_t, bb_t)$g3;
      g4 = onestep(g_t, g_t, gg_t, b_t, bb_t)$g4;
      
      B_gg_sh = cbind(B_gg_sh,(1-g1));
      B_gb_sh = cbind(B_gb_sh,(1-g2));
      B_bg_sh = cbind(B_bg_sh,(1-g3));
      B_bb_sh = cbind(B_bb_sh,(1-g4));
    }
    
    g_sh     = cbind(g_sh,as.matrix(g_t));
    b_sh     = cbind(b_sh,as.matrix(b_t));
    
    
    pieffct  = heal_pi[as.matrix(ai2),3]*pisim96 + heal_pi[as.matrix(ai2),4]*(pisim96^2); # Women 
    bb_t     = logit(healswb[as.matrix(ai2)] + pieffct);
    gg_t     = 1 - logit(healswg[as.matrix(ai2)] + pieffct);
    g_t = get1yr(gg_t,bb_t)$g;
    b_t = get1yr(gg_t,bb_t)$b;
    
    if(ii<simyrs){
      g1 = onestep(g_t, g_t, gg_t, b_t, bb_t)$g1;
      g2 = onestep(g_t, g_t, gg_t, b_t, bb_t)$g2;
      g3 = onestep(g_t, g_t, gg_t, b_t, bb_t)$g3;
      g4 = onestep(g_t, g_t, gg_t, b_t, bb_t)$g4;
      
      B_gg_sw = cbind(B_gg_sw,(1-g1));
      B_gb_sw = cbind(B_gb_sw,(1-g2));
      B_bg_sw = cbind(B_bg_sw,(1-g3));
      B_bb_sw = cbind(B_bb_sw,(1-g4));
    }
    
    b_sw   = cbind(b_sw,as.matrix(b_t));
    g_sw   = cbind(g_sw,as.matrix(g_t)); 
    
    ii=ii+1
  }
  
  datlist = c( "B_gg_sh", "B_gb_sh", "B_bg_sh", "B_bb_sh", "B_gg_sw", "B_gb_sw", "B_bg_sw", "B_bb_sw",
               "g_sh", "b_sh", "g_sw", "b_sw");
  
  for (i in 1:length(datlist)) {
    
    dat = datlist[i];
    
    write.table(get(dat), file = paste0(shkpath,dat,".csv"), sep=",",row.names=FALSE,col.names = FALSE);
  }
}



#INITMX:  Generate medical expenditure shocks.  Here, "medex shocks" are  
#uniformly-distributed variables applied to Markov chain probabilities
#output datasets:  ztacdfsim96      Initial draws of AR(1)
#                  epscdfsim        Simulated innovation on AR(1)
#                  xicdfsim         Simulated white noises shock

initmx <- function(ns){
  
  ztacdfsim96 = runif(ns, min = 0, max = 1);          # Initial draw of AR(1)         
  epscdfsim   = matrix(runif(nn*(simyrs+1), min = 0, max = 1),ns,(simyrs+1) );   # simulated innovation on AR(1)  
  xicdfsim    = matrix(runif(nn*(simyrs+1), min = 0, max = 1),ns,(simyrs+1) );   #  simulated white noise shock
  
  if(simtype==1){
    
    datlist = c("ztacdfsim96", "xicdfsim", "epscdfsim");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      
      write.table(get(dat), file = paste0(shkpath,dat,".csv"), sep=",",row.names=FALSE, col.names = FALSE);
    }
    
  }else if(simtype==2){
    
    ztacdfx = ztacdfsim96;  xicdfx = xicdfsim;  epscdfx = epscdfsim;
    
    datlist = c("ztacdfx", "xicdfx", "epscdfx");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      
      write.table(get(dat), file = paste0(shkpath,dat,".csv"),sep=",", row.names=FALSE, col.names = FALSE);
    }
  }
}




#INITSIM2:  Simulate sequences of health, health cost and demographic shocks
#Here data are sorted by year, with different ages alive each year
#(Get simulated age, health probabilities and observation status.)
#This version uses health and mortality shocks observed in the data,
#with missing (no-wave) values imputed using Baye's Rule.

initsim2 <- function(morts, mort_pi, healdats, heal_pi){
  
  #---Construct health status imputation probabilities using Baye's Rule----#
  
  hsimpute(healdats, heal_pi, nn);
  
  if(simtype==1){
    
    datlist = c("hsimh96", "hsimw96", "mxsim96", "mstat96", "pisim96", "agesim96",
                "msdat", "hhdat", "hwdat", "hmissdat");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
      assign(dat,data) ;
    }
  }else if(simtype==2){
    
    datlist = c("hsimhx", "hsimwx", "mstatx", "mxsim96x", "pisimx", "agesimx", "mstatdx", "hhdx", "hwdx", "hmissdx");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
      assign(dat,data) ;
    }
    
    hsimh96  = hsimhx;    hsimw96 = hsimwx;  mstat96  = mstatx[.,1];  
    mxsim96  = mxsim96x;  pisim96 = pisimx;  agesim96 = agesimx[.,1];  
    msdat    = mstatdx;   hhdat   = hhdx;    hwdat    = hwdx; 
    hmissdat = hmissdx;
  }
  
  mstatsim = mstat96;                    # 0=>defunct, 1=>husband, 2=>wife 
  agesim   = agesim96;
  
  healsimh = hsimh96;                    # 0=>bad health, 1=>good health 
  healsimw = hsimw96;
  
  ageindx  = agesim96-bornage;           # Offset used to get age-appropriate probs 
  
  jj = mmtyrs;
  while (jj >1) {
    msmiss = (msdat[,jj-1]==0)*(msdat[,jj]>0);
    msdat[,jj-1] = msdat[,jj-1]*(1-msmiss) + msmiss*msdat[,jj];
    jj=jj-1;
  }
  
  jj=1;
  while (jj<=mmtyrs) {
    
    #------------First, do simulations to find missing values-------------#
    
    if(jj < mmtyrs){
      ii_L = mmtcols[jj];
      ii_H = mmtcols[jj+1]-1;
    }else{       #Tack on an extra year beyond the sample period 
      ii_L = ncol(healsimh);
      ii_H = ii_L;
    }
    
    ii = ii_L;
    while (ii <= ii_H) { # Update with transition probabilities 
      
      ai2   = ageindx+ii;                   # for appropriate age 
      too_old  = ai2>TR;
      ai2      = ai2*(1-too_old)+too_old*TR;
      
      outcome = gethms(ai2, mstatsim[,ii], healsimh[,ii], healsimw[,ii], pisim96, heal_pi, mort_pi, morts, healdats,nn);
      hsimpoh = outcome$hsimpoh;
      hsimpow = outcome$hsimpow;
      mstpo = outcome$mstpo;
      
      healsimh = cbind(healsimh,hsimpoh);
      healsimw = cbind(healsimw,hsimpow);
      mstatsim = cbind(mstatsim,mstpo);
      
      ii=ii+1;
    }
    
    #-------Replace simulations with observed values, when possible--------#
    
    if(jj >= mmtyrs){
      jj=jj+1     # At this point, no more data is available 
      next;
    }
    
    ii_H     = ii_H+1;
    hmiss    = hmissdat[,jj+1];
    healsimh[,ii_H] = hhdat[,jj+1]*(1-hmiss) + healsimh[,ii_H]*hmiss;
    healsimw[,ii_H] = hwdat[,jj+1]*(1-hmiss) + healsimw[,ii_H]*hmiss;
    mstatsim[,ii_H] = msdat[,jj+1];
    
    #----------Now fill in non-wave years, using imputation probs----------#
    
    obs_ggh = (healsimh[,ii_L]>0.5)*(healsimh[,ii_H]>0.5);
    obs_gbh = (healsimh[,ii_L]>0.5)*(healsimh[,ii_H]<0.5);
    obs_bgh = (healsimh[,ii_L]<0.5)*(healsimh[,ii_H]>0.5);
    obs_bbh = (healsimh[,ii_L]<0.5)*(healsimh[,ii_H]<0.5);
    obs_ggw = (healsimw[,ii_L]>0.5)*(healsimw[,ii_H]>0.5);
    obs_gbw = (healsimw[,ii_L]>0.5)*(healsimw[,ii_H]<0.5);
    obs_bgw = (healsimw[,ii_L]<0.5)*(healsimw[,ii_H]>0.5);
    obs_bbw = (healsimw[,ii_L]<0.5)*(healsimw[,ii_H]<0.5);
    healshk =  matrix(runif(nn*2, min = 0, max = 1),nn,2);
    
    mstatsim[,(ii_L+1)] = mstatsim[,ii_H];
    
    datlist = c("B_gg_sh", "B_gg_sw", "B_gb_sh", "B_gb_sw","B_bg_sh", "B_bg_sw", "B_bb_sh", "B_bb_sw");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
      assign(dat,data) ;
    }
    
    iprob  =  obs_ggh*B_gg_sh[,(ii_L+1)] + obs_gbh*B_gb_sh[,(ii_L+1)] +
      obs_bgh*B_bg_sh[,ii_L+1] + obs_bbh*B_bb_sh[,ii_L+1];
    healsimh[,(ii_L+1)] = (healshk[,1]>iprob);
    iprob  = obs_ggw*B_gg_sw[,(ii_L+1)] + obs_gbw*B_gb_sw[,(ii_L+1)] +
      obs_bgw*B_bg_sw[,(ii_L+1)] + obs_bbw*B_bb_sw[,(ii_L+1)];
    healsimw[,(ii_L+1)] = (healshk[,2]>iprob);
    
    jj = jj +1;
  }
  
  ii = 1;
  
  while (ii <= simyrs) {
    
    agesim   = cbind(agesim,(agesim[,ii]+(mstatsim[,ii+1]>0)));
    too_old  = (agesim[,ii+1]>dieage);
    agesim[,ii+1] = agesim[,ii+1]*(1-too_old) + dieage*too_old;
    mstatsim[,ii+1] = mstatsim[,ii+1]*(1-too_old);
    
    ii = ii +1;
  }
  
  agesim = agesim*(mstatsim>0);
  healsimh = healsimh*1;
  healsimw = healsimw*1;
  
  if(simtype==1){
    
    datlist = c("agesim", "healsimh", "healsimw", "mstatsim");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      
      write.table(get(dat), file = paste0(shkpath,dat,".csv"),sep=",", row.names=FALSE, col.names = FALSE);
    }
    
  }else if(simtype==2){
    
    agesimx = agesim;    mstatx  = mstatsim;  
    hlsimhx = healsimh;  hlsimwx = healsimw;
    
    datlist = c("agesimx", "hlsimhx", "hlsimwx", "mstatx");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      
      write.table(get(dat), file = paste0(shkpath,dat,".csv"),sep=",", row.names=FALSE, col.names = FALSE);
    }
  }
  
  initmx(nn);                         # Generate and save medex shocks
}




#sumshks: Summarize the health status of the simulated health datasets. 
sumshks <- function(){
  
  yearseq2 = c(yearseq,(max(yearseq)+1));
  
  
  datlist = c("healsimh", "healsimw", "mstatsim", "agesim");
  
  for (i in 1:length(datlist)) {
    
    dat = datlist[i];
    data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
    assign(dat,data) ;
  }
  
  frdead = (mstatsim==0)*1;
  frh    = (mstatsim==1)*1;
  frw    = (mstatsim==2)*1;
  frc    = (mstatsim==3)*1;
  
  print( "Fraction of Households that are husbands, wives, couples or dead");
  
  print(cbind(yearseq2,colMeans(frh),colMeans(frw), colMeans(frc), colMeans(frdead)));
  
  hsm = (mstatsim==1)+(mstatsim==3);  #ones(_nn,1); 
  hsw = (mstatsim==2)+(mstatsim==3);  #ones(_nn,1); 
  
  print("Average health status (1=>good) for men and women that are still alive");
  print(cbind(yearseq2, colMeans(healsimh*hsm)/colMeans(hsm), colMeans(healsimw*hsw)/colMeans(hsw), 
              as.vector(sapply(agesim, max))))
  
  cn  = ncol(mstatsim);
  cp1 = healsimh[,1:(cn-1)];
  cp2 = healsimh[,2:cn];
  cph = as.vector(colMeans(cp2*cp1)/colMeans(cp1));
  cp1 = 1-cp1;
  cph =  cbind(cph, as.vector(sapply(cp2*cp1, mean)/sapply(cp1, mean)));
  
  cp1 = healsimw[,1:(cn-1)];
  cp2 = healsimw[,2:cn];
  cpw = as.vector(sapply(cp2*cp1, mean)/sapply(cp1, mean));
  cp1 = 1-cp1;
  cpw = cbind(cpw, as.vector(sapply(cp2*cp1, mean)/sapply(cp1, mean)));
  
  print("Time-t Probability of Good Health at time-t+1")
  health = cbind(yearseq2[1:(cn-1)],cph,cpw);
  colnames(health) = c("year","hgood",  "hbad" , "wgood", "wbad");
  
  print(health);
}



#PIquant2: Generate tables of numbers & shares of the obs less than each quant. 
# e.g.     0.5       1288
#           1        1400
#Specify the quant level for each observation in the data.
#Output: qnts, cndmnum, qnttype

PIquant2 <- function(data,quants,chktie){
  
  qnum     = length(quants);
  rn       = length(data);
  
  if(quants[qnum]==0){
    
    qnts    = c(1,sapply(as.data.frame(data), max))
    cndmnum = c(1,rn);
    qntype  = rep(1,rn);
  } else {
    
    qnts  = matrix(rep(0,qnum*2),qnum,2); 
    qnts[,1]= quants;
    cndmnum  = matrix(rep(0,(qnum+1)*2), (qnum+1),2); 
    cndmnum[,1] = c(quants,1);
    
    #-----------------------Remove Missing Observations-----------------------#
    data[is.na(data)] <- mvcode;
    data  = cbind(data,seq(1,rn,by=1));
    srtddata = data[order(data[,1]), ];
    k        = (srtddata[,1]==mvcode);
    k        = t(k)%*%k; 
    
    if(k > 0){
      srtd2   = srtddata[1:k,];
    }else if(k == 0){
      srtd2   =  srtddata[FALSE,]
    }
    
    srtddata = srtddata[(k+1):rn,];
    
    rn1      = rn-k;                  # Number of non-missing observations 
    cdf      = srtddata[,1];
    qno      = 0;
    qntype   = matrix(1,rn1,1);
    
    i = 1;
    while (i <= qnum) {
      
      qn   = cdf<= quants[i]; 
      qn   = t(qn)%*%qn;
      qnt  = srtddata[qn,1];
      qnts[i,2] = qnt;
      
      #--------Adjust for ties:  shouldn't happen w/ continuous dist--------#
      
      if(chktie == 1){
        j  = (srtddata[qn:rn1,1]==qnt);
        j  = t(j)%*%j -1; 
        qn = qn+j;
      }
      
      if(qn==qno){
        i=i+1;
        next;
      }
      
      rn2 = qn-qno;                     # Counts for interval (q_i-1,q_i] 
      cndmnum[i,2] = rn2;
      qntype[(qno+1):qn] = i*rep(1,rn2);    # observations in (q_i-1,q_i] 
      
      qno = qn;
      
      i = i+1
    }
    #----------------Count observations above highest quantile----------------#
    
    if(qn < rn1){
      qntype[(qn+1):rn1] = (qnum+1)*rep(1,rn1-qn);
      cndmnum[(qnum+1),2] = rn1 - qn; 
    }
    
    if(k>0){
      qntype = c(rep(0,k),qntype);
    }
    
    data   = rbind(srtd2,srtddata);
    data   = cbind(data,qntype);
    data   = data[order(data[,2]), ];
    qntype = data[,3];
  }
  
  output = list("qnts"=qnts,"cndmnum" = cndmnum,"qntype" = qntype);
}

#Getmean:  Find weighted mean, after removing missing values
#output: mns,  cndmnum.

getmean <- function(data,wgts){
  mns     = rep(0,2);
  cndmnum = rep(0,2);
  data    = na.omit(cbind(data,wgts));              # Remove Missing Observations 
  if(anyNA(data)){
    mns[2] = -1e6;
  }else{
    wgts   = data[,2];
    wgts   = wgts/sum(wgts);
    data   = data[,1]*wgts;                     # Use weighted data S
    mns[2] = sum(data);
    cndmnum[2] = length(data);
  }
  
  output = list("mns" = mns,"cndmnum" = cndmnum);
}

#Getquant:  Finds (weighted) quantiles, after removing missing values 
#Can handle weighted data
#output: qnts, cndmnum

getquant <- function(data,wgts,quants,chktie){
  
  qnum     = length(quants);
  rn       = length(data);
  qnts     = matrix(0,qnum,2);
  qnts[,1] = quants;
  cndmnum  = matrix(0,qnum+1,2);
  cndmnum[,1] = c(quants,1);
  
  data     = na.omit(cbind(data,wgts));              # Remove Missing Observations 
  srtddata = data[order(data[,1]), ];
  rn1      = length(srtddata);         # Number of non-missing observations 
  cdf      = srtddata[,2];             # Base distribution on Weights 
  cdf      = cumsum(cdf)/sum(cdf);
  qno      = 0;
  
  i=1;
  while (i<=qnum) {
    qn   = cdf< quants[i]; 
    qn   = t(qn)%*%qn+1;
    qnt  = srtddata[qn,1];
    qnts[i,2] = qnt;
    
    #--------Adjust for ties:  shouldn't happen w/ continuous dist--------#
    if(chktie == 1){
      j  = (srtddata[qn:rn1,1]==qnt);
      j  = t(j)%*%j -1; 
      qn = qn+j;
    }
    
    if(qn==qno){
      i = i+1
      next
    }
    
    rn2 = qn-qno;                    # Counts for interval (q_i-1,q_i] 
    cndmnum[i,2] = rn2;
    qno = qn;
    
    i= i+1
    
  }
  
  if(qn < rn1){
    cndmnum[qnum+1,2] = rn1 - qn;
  }
  
  output = list("qnts" = qnts,"cndmnum"= cndmnum)
  
}

#Getqunt2:  Finds quantiles, after removing missing values
#output: qnts, cndmnum

getqunt2 <- function(data,quants){
  
  qnum      = length(quants);
  qnts      = matrix(0,qnum,2);
  qnts[,1]  = quants;
  cndmnum   = matrix(0,(qnum+1),2);
  cndmnum[,1] = c(quants,1);
  
  data      = na.omit(data);                  #Remove Missing Observations 
  rn        = length(data);
  if(rn>1){
    qnts[,2] = quantile(data,quants,type = 4);
  }else{
    qnts[,2] = -1e6;
  }
  
  qno      = 0;
  i = 1;
  while(i<=qnum){
    qn   = data <=qnts[i,2]; 
    qn   = t(qn)%*%qn;
    if(qn==qno){
      i = i+1;
      next
    }
    rn2 = qn-qno;                     # Counts for interval (q_i-1,q_i] 
    cndmnum[i,2] = rn2;
    qno = qn;
    i=i+1;
  }
  #----------------Count observations above highest quantile----------------#
  if(qn < rn){
    cndmnum[(qnum+1),2] = rn - qn;
  }
  output = list("qnts" = qnts,"cndmnum" = cndmnum);
}




#simqunts: Calculate quantiles & number of observations in each quantile of the data for given year given
#cohorts given PI quant. Save these quantiles & number of observations in two lists: dataprfs, datacnts.

simqunts <- function(pisim96, agesim96, datasim, mssim2, alivesim, pistate_j, 
                     cohorts_j, quants_j, mmtyrs_j, comptype, wgts){
  
  pinum_j   = length(pistate_j) + (identical(pistate_j,rep(0,length(pistate_j))) == FALSE);
  chrtnum_j = length(cohorts_j)-1;
  qnum_j    = length(quants_j); 
  
  outcome = PIquant2(pisim96,pistate_j,chktie);
  PIqnts = outcome$qnts;
  PIcnts = outcome$cndmnum;
  PItype = outcome$qntype;
  
  outcome = getchrt(agesim96,cohorts_j);
  chrtcnts = outcome$chrtcnts;
  chrttype = outcome$chrttype;
  
  cn = ncol(alivesim);
  
  if(comptype==1){
    print("Looking at Survivors Only!!");
    alivesim =as.data.frame(replicate(cn,as.vector(rowMins(alivesim[,c(1,cn)]))));
  }
  
  nummtx   = matrix(0,chrtnum_j,pinum_j);
  dataprfs = rep(list(rep(list(rep(list(matrix(mvcode,qnum_j,mmtyrs_j)),MSnum)),pinum_j)),chrtnum_j);
  datacnts = rep(list(rep(list(rep(list(matrix(0,qnum_j,mmtyrs_j)),MSnum)),pinum_j)),chrtnum_j);
  
  iChrt = 1;
  while(iChrt <= chrtnum_j){
    
    iPI=1;
    while (iPI <= pinum_j) {
      
      indicat0 = (chrttype==iChrt)*(PItype==iPI);
      indicat  = indicat0*alivesim[,1];
      nummtx[iChrt,iPI] = length(indicat)*mean(indicat);
      
      iMStat=1; 
      while (iMStat <= MSnum) {
        if(MSsplit==0){
          indicat = replicate(mmtyrs_j, as.vector(indicat0));
        } else if(MSsplit==1){
          
          indicat = replicate(mmtyrs_j, as.vector(indicat0))*(mssim2==iMStat);
        }
        
        indicat = indicat*alivesim;
        indicat[indicat==0] <- NA
        data    = datasim*indicat;
        
        iYear=1;
        while (iYear <= mmtyrs_j) {
          
          if(sum(indicat[,iYear],na.rm=TRUE)<quantmin){
            break;
          }
          if(quants_j==0){  # using means, rather than quantiles 
            outcome = getmean(data[,iYear],wgts);
            tempprf = outcome$mns;
            tempnum = outcome$cndmnum;
            tempprf = t(as.matrix(tempprf));
            tempnum = t(as.matrix(tempnum));
          }else{
            
            if(wgtddata==0){
              outcome = getqunt2(data[,iYear],quants_j);
              tempprf = outcome$qnts;
              tempnum = outcome$cndmnum;
            }else{
              outcome = getquant(data[.,iYear],wgts,quants_j,chktie);
              tempprf = outcome$qnts;
              tempnum = outcome$cndmnum;
            }
          }
          
          iQunt=1;
          while (iQunt <= qnum_j) {
            
            dataprfs[[iChrt]][[iPI]][[iMStat]][iQunt,iYear] = tempprf[iQunt,2];
            datacnts[[iChrt]][[iPI]][[iMStat]][iQunt,iYear] = tempnum[iQunt,2];
            
            iQunt = iQunt+1;
          }
          
          iYear = iYear+1;
        }
        
        iMStat = iMStat+1;
      }
      
      iPI = iPI+1;
    }
    
    iChrt = iChrt+1;
  }
  
  if(prnres>1){
    print( "Cohort-PI counts (from first wave in moment criterion)")
    nummtx
  }
  
  output = list("dataprfs" = dataprfs,"datacnts" = datacnts);
  
}

#simcrrl: Calculate means, standard deviations, number of obs and correlations between t and t+1, t and t+2. 
#Save these data in 7 lists: "meanprfs","meancnts","stdprfs","crrlprf1", "crrlcnt1","crrlprf2","crrlcnt2".

simcrrl <- function(pisim96, agesim96, datasim, mssim2, alivesim, pistate_j, cohorts_j, mmtyrs_j, comptype, wgts){
  
  pinum_j   = length(pistate_j) + (identical(pistate_j,rep(0,length(pistate_j))) == FALSE);
  chrtnum_j = length(cohorts_j)-1;
  
  outcome = PIquant2(pisim96,pistate_j,chktie);
  PIqnts = outcome$qnts;
  PIcnts = outcome$cndmnum;
  PItype = outcome$qntype;
  
  outcome = getchrt(agesim96,cohorts_j);
  chrtcnts = outcome$chrtcnts;
  chrttype = outcome$chrttype;
  
  cn = ncol(alivesim);
  
  if(comptype==1){
    print("Looking at Survivors Only!!")
    
    alivesim = as.data.frame(replicate(cn,as.vector(rowMins(alivesim[,c(1,cn)]))));
    
  }
  
  meanprfs = rep(list(rep(list(matrix(mvcode,MSnum,mmtyrs_j)),pinum_j)),chrtnum_j);
  meancnts = rep(list(rep(list(matrix(mvcode,MSnum,mmtyrs_j)),pinum_j)),chrtnum_j);
  stdprfs  = rep(list(rep(list(matrix(mvcode,MSnum,mmtyrs_j)),pinum_j)),chrtnum_j);
  crrlprf1 = rep(list(rep(list(matrix(mvcode,MSnum,mmtyrs_j-1)),pinum_j)),chrtnum_j);
  crrlcnt1 = rep(list(rep(list(matrix(0,MSnum,mmtyrs_j-1)),pinum_j)),chrtnum_j);
  crrlprf2 = rep(list(rep(list(matrix(mvcode,MSnum,mmtyrs_j-2)),pinum_j)),chrtnum_j);
  crrlcnt2 = rep(list(rep(list(matrix(0,MSnum,mmtyrs_j-2)),pinum_j)),chrtnum_j);
  
  iChrt=1;
  while(iChrt <= chrtnum_j){
    
    iPI = 1;
    while (iPI <= pinum_j) {
      
      indicat0 = (chrttype==iChrt)*(PItype==iPI);
      indicat  = indicat0*alivesim[,1];
      
      iMStat=1;
      while (iMStat <= MSnum) {
        
        if(MSsplit==0){
          indicat = replicate(mmtyrs_j, as.vector(indicat0));
        }else if(MSsplit==1){
          indicat = replicate(mmtyrs_j, as.vector(indicat0))*(mssim2==iMStat);
        }
        
        indicat = indicat*alivesim;
        
        #  "Cohort = ";;iCHrt;; "PI Quintile = ";;iPI;; meanc(indicat)';  
        
        indicat[indicat==0] <- NA
        
        data0   = datasim*indicat;
        dY0     = datasim[,1]*0;
        dY1     = dY0;
        dY2     = dY0;
        
        iYear=1;
        while (iYear <= mmtyrs_j) {
          
          if(sum(indicat[,iYear],na.rm = TRUE)<quantmin){
            break;
          }
          
          dY0 = data0[,iYear];
          
          outcome = getmean(dY0,wgts);
          tempprf = outcome$mns;
          tempnum = outcome$cndmnum;
          tempprf = t(as.matrix(tempprf));
          tempnum = t(as.matrix(tempnum))
          
          meanprfs[[iChrt]][[iPI]][iMStat,iYear]= tempprf[1,2];
          meancnts[[iChrt]][[iPI]][iMStat,iYear]= tempnum[1,2];
          
          dY0 = dY0 - tempprf[1,2];     # Convert into residual 
          
          outcome = getmean(dY0^2,wgts);
          tempprf = outcome$mns;
          tempnum = outcome$cndmnum;
          tempprf = t(as.matrix(tempprf))
          tempnum = t(as.matrix(tempnum));
          
          stdprfs[[iChrt]][[iPI]][iMStat,iYear] = sqrt(tempprf[1,2]);
          dY0 = dY0/sqrt(tempprf[1,2]);    # Normalize residual 
          
          if(iYear>1){
            outcome = getmean(dY0*dY1,wgts);
            tempprf = outcome$mns;
            tempnum = outcome$cndmnum;
            tempprf = t(as.matrix(tempprf));
            tempnum = t(as.matrix(tempnum));
            
            
            crrlprf1[[iChrt]][[iPI]][iMStat,iYear-1] = tempprf[1,2];
            crrlcnt1[[iChrt]][[iPI]][iMStat,iYear-1] = tempnum[1,2];
          }
          
          if(iYear>2){
            outcome = getmean(dY0*dY2,wgts);
            tempprf = outcome$mns;
            tempnum = outcome$cndmnum;
            tempprf = t(as.matrix(tempprf));
            tempnum = t(as.matrix(tempnum));
            
            crrlprf2[[iChrt]][[iPI]][iMStat,iYear-2] = tempprf[1,2];
            crrlcnt2[[iChrt]][[iPI]][iMStat,iYear-2] = tempnum[1,2];
          }
          dY2 = dY1;
          dY1 = dY0;       
          
          iYear=iYear+1;
        }
        iMStat=iMStat+1;
      }
      iPI=iPI+1;
    }
    iChrt=iChrt+1;
  }
  
  output = list("meanprfs" = meanprfs,"meancnts"= meancnts,"stdprfs"= stdprfs,"crrlprf1"= crrlprf1, 
                "crrlcnt1" = crrlcnt1,"crrlprf2"= crrlprf2,"crrlcnt2"= crrlcnt2)
}

# procedure to calculate a univariate kernel density estimate and its derivative
#usage:      {f, d, h} = ukernel(x, z, h, w, &kf);
#input:      x:      T vector where density is to be estimated
#z:      n vector with observed data points
#h:      scalar bandwidth, if h<=0 bandwidth is determined by
#procedure bandw1(z)
#w:      n vector with weights
#kf:    pointer to weighting function
#output:     f:      T vector with estimated density
#d:      T vector with estimated derivative
#h:      scaler bandwidth

ukernel <- function(x,z,h,w, kf){
  
  if(ncol(as.matrix(x))>1){
    warning("ukernel.g: x has too many columns");
    return(list("f"=-1,"d"=-1,"h"=-1));
  };
  
  if(ncol(as.matrix(z))>1){
    warning("ukernel.g: z has too many columns");
    return(list("f"=-1,"d"=-1,"h"=-1));
  };
  
  # initialization 
  
  i = 1;
  n = length(x);
  k = ncol(x);
  
  # determine bandwidth 
  
  if(h<=0){
    h=bandw1(z);
  };
  
  f=rep(0,n);
  d=rep(0,n);
  
  while (i<=n) {
    
    arg = (x[i]-z)/h;
    
    kff = get(kf)(arg)$kff;
    kfd = get(kf)(arg)$kfd;
    
    f[i] = mean(kff*w)/h;
    d[i] = mean(kfd*w)/(h^2);
    
    i=i+1;
  }
  
  return(list("f"=f,"d"=d,"h"=h));
  
}


# bandw1
#procedure to calculate the optimal bandwidth in kernel estimation of a density.
#The optimal bandwidth is calculated according to eq. 3.31 of Silverman (1986)

#usage:      h=bandw1(y);
#input:      y:      n-vector whose density will be estimated;
#output:     h:      scalar, optimal bandwidth choice;


bandw1 <- function(y){
  
  if(ncol(as.matrix(y))>1){
    warning("input error in bandw1.g: too many columns");
    return(-1);
  };
  
  s=sqrt(var(y));
  n=length(y);
  ys=sort(y);
  qi1=round(0.25*n);
  qi3=round(0.75*n);
  iqr=ys[qi3]-ys[qi1];
  
  return(0.9*min(c(s,(iqr/1.34)))/n^0.2)
}


#This part contains some kernel functions. All functions take an
#n x k matrix u as their argument and return an n x k matrix with the
#function evaluated in each point of x and an n x k matrix d with the
#derivative of the function in each point of x.
#k_bw:     biweight kernel function
#k_epan:   Epanechnikov kernel
#k_gauss:  Gaussian kernel
#k_triang: triangular kernel
#k_rect:   rectangular kernel

k_bw <- function(u){
  
  select = abs(u)<=1;
  
  return(list("kff" = 15/16*((1-u^2)^2)*select,"kfd" = -15/4*u*(1-u^2)*select));
}

k_gauss <- function(u){
  
  return(list("kff" = dnorm(u),"kfd" = -u*dnorm(u)));
}

k_epan <-function(u){
  
  s = abs(u)<sqrt(5);
  c = 0.75/sqrt(5);
  return(list("kff" = c*(1-0.2*u^2)*s, "kfd" = -0.4*c*u*s));
}

k_rect <- function(u){
  
  return(list("kff" = 0.5*(abs(u)<1), "kfd" = 0*u ));
}

k_trian <- function(u){
  
  a = abs(u);
  s = a<1;
  
  return(list("kff" = (1-a)*s, "kfd" = (-(u>=0) + (u<=0))*s ));
}

# getpdf:  Finds kernel densities, using non-missing data.
# Kernel density estimator written by Ruud Koenig
# return "pdfs

getpdf <- function(data,indicat,xvals){
  
  indicat = indicat>0;
  indicat[indicat==0] <- NA
  data    = data*indicat;
  pdfs    = vector();
  rn      = nrow(data);
  
  iWave=1;
  while (iWave <= ncol(data)) {
    
    kdns  = 0*xvals;
    datac = data[,iWave];
    datac = na.omit(datac)       # Remove Missing Observations 
    
    if(anyNA(datac)==FALSE){
      
      kdns = ukernel(xvals[iWave],datac,0,1,"k_gauss")$f;
      dkdns= ukernel(xvals[iWave],datac,0,1,"k_gauss")$d;
      bw   = ukernel(xvals[iWave],datac,0,1,"k_gauss")$h;
      
    }
    
    pdfs  = c(pdfs,kdns);
    
    iWave=iWave+1; 
  }
  
  return(pdfs)
}


#makemmts_q:Put qunatiles saved in a list in a vector "qntvec". Put a minus sign in front of "qntvec" to obtain "pdfvec".
#Generate a dataset indicating the realtionship between the data and the quantile, do this for the dataset of each given year
#given cohort, given PI quant. Then column combine them as "mmtmtx". 
#Column combine the indicate data for each cohort, PI quant and year as "obsmtx".
#"mmtskeep" save the column number of moments to save (not with too few observations).

makemmts_q <- function(PIdat, agedat, data, MStatdat, obsdat, wgts, dataprfs,
                       pistate_j, cohorts_j, quants_j, mmtyrs_j, keepinit, addpdfs){
  
  nobs      = length(PIdat);
  pinum_j   = length(pistate_j) + (identical(pistate_j,rep(0,length(pistate_j))) == FALSE);
  chrtnum_j = length(cohorts_j)-1;
  qnum_j    = length(quants_j);
  
  if(quants_j==0){
    addpdfs=0;  # Looking at means
  }
  
  outcome = PIquant2(PIdat,pistate_j,chktie);
  PIqnts = outcome$qnts;
  PIcnts = outcome$cndmnum;
  PItype = outcome$qntype;
  
  outcome = getchrt(agedat,cohorts_j);
  chrtcnts = outcome$chrtcnts;
  chrttype = outcome$chrttype;
  
  mmtmtx = vector();
  obsmtx = vector();
  qntvec = vector();
  pdfvec = vector();
  mmtskip= vector();
  mmtskeep=vector();
  iM       = 0;
  
  iChrt=1;
  while (iChrt <= chrtnum_j) {
    
    iPI=1;
    while (iPI<=pinum_j) {
      
      indicat0 = (chrttype ==iChrt)*(PItype==iPI);
      
      iMStat=1;
      while(iMStat <= MSnum){
        
        if(MSsplit==0){
          indicat = replicate(mmtyrs_j, as.vector(indicat0));
        }else if(MSsplit==1){
          indicat = replicate(mmtyrs_j, as.vector(indicat0))*(MStatdat==iMStat);
        };
        
        indicat = indicat*obsdat;
        indicat = indicat*wgts;
        
        iQunt=1;
        while (iQunt<=qnum_j) {
          
          tempprf = dataprfs[[iChrt]][[iPI]][[iMStat]][iQunt,]
          qntvec = c(qntvec,tempprf)
          
          if(quants_j==0){   #looking at means 
            mmtdata = sweep(data, 2, tempprf);
          }else{
            mmtdata =  sweep(data,2,tempprf,"<=") - quants_j[iQunt];
          }
          
          mmtdata = mmtdata*indicat;
          mmtmtx  = cbind(mmtmtx,as.matrix(mmtdata));
          obsmtx  = cbind(obsmtx,as.matrix(indicat));
          temppdf = -tempprf;
          
          if(addpdfs==1){
            temppdf = getpdf(data,indicat,tempprf);
          };
          
          pdfvec  = c(pdfvec,temppdf);
          cnum    = colSums(indicat);
          
          # Identify moments to drop:  These include moments using  
          # the initial asset distribution, or moments with too few observations 
          
          iYear=1;
          while (iYear<=mmtyrs_j) {
            
            iM = iM+1;
            mmtskeep = c(mmtskeep,iM);
            
            if(cnum[iYear] < cellmin){
              mmtskip = c(mmtskip,iM);
              cnum[iYear:mmtyrs_j] = rep(0,mmtyrs_j-iYear+1); # Skip all future years 
            }else if(iYear==1&keepinit==0){
              mmtskip = c(mmtskip,iM);
            };
            iYear = iYear+1;
          }
          iQunt=iQunt+1;
        }
        iMStat=iMStat+1;
      }
      iPI=iPI+1;
    }
    iChrt=iChrt+1;
  }
  
  if(length(mmtskip)>0){
    mmtskeep[mmtskip] = rep(0,length(mmtskip));
  };
  
  mmtskeep = sort(mmtskeep);
  zn       = sum(mmtskeep < 0.99);
  mmtskeep = mmtskeep[(zn+1):iM];
  
  mmtmtx = mmtmtx[,mmtskeep];     #Delete moments derived from initial 
  obsmtx = obsmtx[,mmtskeep];     #distribution of assets SS
  qntvec = qntvec[mmtskeep];
  pdfvec = pdfvec[mmtskeep];
  
  mmtvec = colMeans(mmtmtx); 
  obsvec = colMeans(obsmtx);
  
  output = list("qntvec" = qntvec,"pdfvec" = pdfvec,"mmtskeep"= mmtskeep,"iM" = iM,"mmtmtx"= mmtmtx,"obsmtx"= obsmtx);
}


makemmts_c <- function(PIdat, agedat, data0, MStatdat, obsdat, wgts, meanprfs, 
                       stdprfs, crrlprf1, crrlprf2, pistate_j, cohorts_j, mmtyrs_j){
  
  nobs      = length(PIdat);
  pinum_j   = length(pistate_j) + (identical(pistate_j,rep(0,length(pistate_j))) == FALSE);
  chrtnum_j = length(cohorts_j)-1;
  
  outcome = PIquant2(PIdat,pistate_j,chktie);
  PIqnts = outcome$qnts;
  PIcnts = outcome$cndmnum;
  PItype = outcome$qntype;
  
  outcome = getchrt(agedat,cohorts_j)
  chrtcnts = outcome$chrtcnts;
  chrttype = outcome$chrttype;
  
  mmtmtx   = vector();                       # rows->obs, cols->moment function 
  obsmtx   = vector();
  crrlvec  = vector();
  mmtskip  = vector();
  mmtskeep = vector();
  iM       = 0;
  
  iChrt=1;
  while (iChrt <= chrtnum_j) {
    
    iPI=1;
    while (iPI <= pinum_j) {
      
      indicat0 = (chrttype==iChrt)*(PItype==iPI);
      
      iMStat=1;
      while (iMStat<=MSnum) {
        
        if(MSsplit==0){
          indicat = replicate(mmtyrs_j, as.vector(indicat0));
        } else if(MSsplit==1){
          indicat = replicate(mmtyrs_j, as.vector(indicat0))*(MStatdat==iMStat);
        };
        
        indicat = indicat*obsdat;
        indicat = indicat*sqrt(wgts);      # terms are multiplied below 
        
        # meanvec and stdvec should be from data as the model is matching level means, not log means or deviations.                        
        
        meanvec = meanprfs[[iChrt]][[iPI]][iMStat,];  # a row vector 
        data1   = sweep(data0, 2, meanvec);           # Convert into residuals            
        stdvec  = stdprfs[[iChrt]][[iPI]][iMStat,];
        data1   = sweep(data1,2,stdvec,"/" )          # Normalize residual
        dY0     = data1[,1]*0;
        dY1     = dY0;
        dY2     = dY0;
        icY0    = dY0;
        icY1    = dY0;
        icY2    = dY0;
        
        iYear=1;
        while (iYear <= mmtyrs_j) {
          
          dY0  = data1[,iYear];
          icY0 = indicat[,iYear];
          
          if(iYear>2){   #We lack 2-year averages for year 1
            
            tempprf = crrlprf1[[iChrt]][[iPI]][iMStat,iYear-1]
            crrlvec  = c(crrlvec,tempprf);
            mmtdata  = dY0*dY1 - tempprf;
            indicati = icY0*icY1;
            mmtdata  = mmtdata*indicati;
            mmtmtx   = cbind(mmtmtx,as.matrix(mmtdata));
            obsmtx   = cbind(obsmtx,as.matrix(indicati));
            cnum     = sum(indicati>0);
            iM       = iM+1;
            mmtskeep = c(mmtskeep,iM);
            
            if(cnum < cellmin){   #Drop moments with too few observations
              mmtskip = c(mmtskip,iM);
            }
          }
          
          if(iYear>3){
            
            tempprf = crrlprf2[[iChrt]][[iPI]][iMStat,iYear-2];
            crrlvec  = c(crrlvec,tempprf);
            mmtdata  = dY0*dY2 - tempprf;
            indicati = icY0*icY2;
            mmtdata  = mmtdata*indicati;
            mmtmtx   = cbind(mmtmtx,as.matrix(mmtdata));
            obsmtx   = cbind(obsmtx,as.matrix(indicati));
            cnum     = sum(indicati>0);
            iM       = iM+1;
            mmtskeep = c(mmtskeep,iM);
            
            if(cnum < cellmin){
              
              mmtskip = c(mmtskip,iM);
            }
          }
          
          dY2  = dY1;
          dY1  = dY0; 
          icY2 = icY1;
          icY1 = icY0;       
          
          iYear=iYear+1;
        }
        
        iMStat=iMStat+1;
      }
      
      iPI=iPI+1;
    }
    
    iChrt=iChrt+1;
  }
  
  if(length(mmtskip)>0){
    
    mmtskeep[mmtskip] = rep(0,rows(mmtskip));
  }
  
  mmtskeep = sort(mmtskeep);
  zn       = sum(mmtskeep< 0.99);
  mmtskeep = mmtskeep[zn+1:iM];
  
  mmtmtx  = mmtmtx[,mmtskeep];  
  obsmtx  = obsmtx[,mmtskeep];  
  crrlvec = crrlvec[mmtskeep];
  
  mmtvec  = colMeans(mmtmtx);
  obsvec  = colMeans(obsmtx);
  
  output = list("crrlvec" = crrlvec,"mmtskeep" = mmtskeep,"iM" = iM,"mmtmtx" = mmtmtx,"obsmtx" = obsmtx);
}

makemmts <- function(PIdat, agedat, asstdat, MStatdat, obsdat, wgts, asstqnts,
                     mxdat, mxobsdat, mxquants, mxmeans, lmxmnsdat, lmxstddat, 
                     mxcrrls1, mxcrrls2, addpdfs){
  #  First, evaluate quantiles for assets 
  
  outcome = makemmts_q(PIdat,as.data.frame(agedat)[,1],asstdat,MStatdat,obsdat,datawgts,
                       asstqnts,pistate_a,cohorts_a,quants_a,mmtyrs,0,addpdfs);
  qntvec_a = outcome$qntvec;
  pdfvec_a = outcome$pdfvec;
  mmtskeep_a = outcome$mmtskeep;
  iM_a = outcome$iM;
  mmtmtx_a= outcome$mmtmtx;
  obsmtx_a = outcome$obsmtx;
  
  mmtmtx     = mmtmtx_a;
  obsmtx     = obsmtx_a;
  qntvec     = qntvec_a;
  pdfvec     = pdfvec_a;
  mmtskeep   = mmtskeep_a;
  iM         = iM_a;
  mmttype    = rep(1,iM_a);
  
  
  if(medexmmts==1){
    #  Next, evaluate quantiles for Medex 
    
    outcome = makemmts_q(PIdat,data.frame(agedat)[,1],mxdat,MStatdat,mxobsdat,datawgts,
                         mxquants,pistate_m,cohorts_m,quants_m,mmtyrs,0,addpdfs);
    qntvec_m = outcome$qntvec;
    pdfvec_m = outcome$pdfvec;
    mmtskeep_m = outcome$mmtskeep;
    iM_m = outcome$iM;
    mmtmtx_m = outcome$mmtmtx;
    obsmtx_m = outcome$obsmtx;
    
    mmtmtx     = cbind(mmtmtx,mmtmtx_m);
    obsmtx     = cbind(obsmtx,obsmtx_m);
    qntvec     = c(qntvec,qntvec_m);
    pdfvec     = c(pdfvec,pdfvec_m);
    mmtskeep_m = mmtskeep_m + iM_a;
    mmtskeep   = c(mmtskeep,mmtskeep_m);
    iM         = iM+iM_m;
    mmttype    = c(mmttype,(2*rep(1,iM_m)));
    
    # Next, evaluate Medex means
    outcome  = makemmts_q(PIdat,data.frame(agedat)[,1],mxdat,MStatdat,mxobsdat,datawgts,
                          mxmeans,pistate_m,cohorts_m,0,mmtyrs,0,0);
    meanvec_m = outcome$qntvec;
    pdfvec_m = outcome$pdfvec;
    mmtskeep_m = outcome$mmtskeep;
    iM_m = outcome$iM;
    mmtmtx_m = outcome$mmtmtx;
    obsmtx_m = outcome$obsmtx;
    
    mmtmtx     = cbind(mmtmtx,mmtmtx_m);
    obsmtx     = cbind(obsmtx,obsmtx_m);
    qntvec     = c(qntvec,meanvec_m);
    mmtskeep_m = mmtskeep_m + iM;
    mmtskeep   = c(mmtskeep,mmtskeep_m);
    iM         = iM+iM_m;
    mmttype    = c(mmttype,(3*rep(1,iM_m)));
    
    # Next, look at log Medex autocorrelations
    
    outcome = makemmts_c(PIdat,data.frame(agedat)[,1],log(mxdat+(1-mxobsdat)),MStatdat,
                         mxobsdat,datawgts,lmxmnsdat,lmxstddat,mxcrrls1,mxcrrls2,pistate_m,cohorts_m,mmtyrs);
    crrlvec_m = outcome$crrlvec;
    mmtskeep_m = outcome$mmtskeep;
    iM_m = outcome$iM;
    mmtmtx_m = outcome$mmtmtx;
    obsmtx_m = outcome$obsmtx;
    
    mmtmtx     = cbind(mmtmtx,mmtmtx_m);
    obsmtx     = cbind(obsmtx,obsmtx_m);
    qntvec     = c(qntvec,crrlvec_m);
    mmtskeep_m = mmtskeep_m + iM;
    mmtskeep   = c(mmtskeep,mmtskeep_m);
    mmttype    = c(mmttype,(4*rep(1,iM_m)));
  }
  
  mmtvec     = colMeans(mmtmtx);
  obsvec     = colMeans(obsmtx);
  nmom       = length(mmtvec);
  mmttype    = mmttype[mmtskeep];
  MSMval     = (totobs*(mmtvec^2)*diag(W));
  
  if(savemmts==1){     #N.B. This is a global
    
    datlist = c( "mmtmtx", "obsmtx");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      
      write.table(get(dat), file = paste0(datapath,dat,".csv"), sep=",", row.names=FALSE, col.names = FALSE);
    }
  }
  
  output = list("mmtvec"= mmtvec,"obsvec" = obsvec,"qntvec" = qntvec,"pdfvec" = pdfvec,"mmttype" = mmttype);
}


#compvinv: generate covariance matrix, inverse of the covariance matrix (as the weighting matrix.)
compvinv<-function(optwgts){
  
  mmtmtx = read.csv(file = paste0(datapath,"mmtmtx",".csv"), header = FALSE);
  
  mnvals  = 0;
  mnvals  = colMeans(mmtmtx);      # This data should already be zero-mean 
  mmtmtx  = sweep(mmtmtx, 2, mnvals);
  
  totobs  = nrow(mmtmtx);
  mmtmtx = as.matrix(mmtmtx)
  vcv     = t(mmtmtx)%*%mmtmtx;
  rm(mmtmtx);
  vcv     = vcv/totobs;
  
  #-----Take the principal diagonal and form a diagonal weighting matrix-----
  
  vinv    = solve(vcv);  
  rn      = nrow(vcv);
  vdiag   = diag(vcv);
  
  write.table(vcv, file = paste0(datapath,"vcv",".csv"), sep=",", row.names=FALSE, col.names = FALSE);
  vdiag   = diag(vdiag);
  vdiag   = solve(vdiag);
  
  if(optwgts==0){
    W  = diag(rn);
  }else if(optwgts==1){
    W = vdiag;
  }else if(optwgts==2){
    W = vinv;
  }
  
  output = list("VCV" = vcv,"vinv" = vinv,"W" = W);
  
}


fillin <- function(x,ageseq2){
  
  rn = length(x);
  
  i=2;
  while (i<=rn) {
    
    if(is.na(x[i])==TRUE & is.na(x[i-1])==FALSE){
      
      j=i+1;
      while (j<=rn) {
        
        if(is.na(x[j])==FALSE){
          break
        }
        
        j = j+1;
      }
      
      if(j<=rn){
        frac = (ageseq2[i]-ageseq2[i-1])/(ageseq2[j]-ageseq2[i-1]);
        x[i] = x[i-1] + frac*(x[j]-x[i-1]);
      }
    }
    
    i=i+1;
  }
  
  return(x) 
}

grphmtx <- function(dataprfs,vartype,comptype,datatype,qnum_j,pinum_j,chrtnum_j){
  
  if(vartype==1){        #consumption vs. assets
    name1 = "a";
  }else if(vartype==2){
    name1 = "c";
  }else if(vartype==3){
    name1 = "m";
  }
  
  if(comptype==0){       #All observations vs. survivors
    name2 = "all";
  }else if(comptype==1){
    name2 = "srv";
  }
  
  if(datatype==0){       #Data vs. simulations
    name3 = "dt";  
  }else if(datatype==1){
    name3 = "sm"; 
  }
  
  if(qnum_j==0){         #looking at means
    iQunt=0;
  }else{
    iQunt=1;
  }
  
  tr2 = max(avgage96)+simyrs - bornage;
  tr2 = max(c(TR,tr2));
  ageseq2 = seq(bornage,bornage+tr2-1,by=1); 
  
  while (iQunt <= qnum_j) {
    
    name5 = toString(round(iQunt));
    
    iMStat=1;
    while (iMStat <= MSnum) {
      
      if(MSsplit==0){      #No Marital Status Distinctions
        MSIndex = 1;
      }else if(MSsplit==1){    #Single M vs. Single F vs. Couples
        MSIndex = iMStat;
      }
      
      gmat = matrix(, nrow = tr2, ncol = chrtnum_j*pinum_j );
      gmat = cbind(ageseq2,gmat);
      gotsome = rep(0,tr2);        # Records ages with observations 
      cn   = 1;
      
      iChrt=1;
      while(iChrt <= chrtnum_j){
        
        iPI=1;
        while (iPI <= pinum_j) {
          
          if(iQunt==0){   #Means
            tempprf = dataprfs[[iChrt]][[iPI]][[MSIndex]][1,]  #a row vector
          }else{
            tempprf = dataprfs[[iChrt]][[iPI]][[MSIndex]][iQunt,]  #a row vector
          }
          
          skipem  = (tempprf==mvcode);
          skipem  = t(skipem)%*%skipem;
          cn      = cn+1;
          rn      = mmtcols+avgage96[iChrt]-bornage;
          rn2     = mmtyrs - skipem;
          
          if(rn2>0){
            
            rn = rn[1:rn2];
            rn  = rn[which(rn<=tr2)];
            rn2 = length(rn);
            gmat[rn,cn] = tempprf[1:rn2];
            gotsome[rn] = rep(1,rn2); 
            
          }
          
          iPI=iPI+1; 
        }
        
        iChrt=iChrt+1; 
      }
      
      iYear= 1;
      while (max(iYear) <= max(rn-1)) {
        
        if(gotsome[iYear]==1){
          
          jYear=tr2;
          while (jYear >iYear) {
            
            if(gotsome[jYear]==1){
              
              gotsome[iYear:jYear]=rep(1,jYear-iYear+1); 
              jYear=iYear+1;
            }
            
            jYear=jYear-1;
          }
          
          iYear=rn-1;
        }
        
        iYear=iYear+1; 
      }
      
      gmat = cbind(gmat,gotsome);
      gmat = gmat[gmat[,ncol(gmat)] == 1,-ncol(gmat)];   #Drop ages with no observations 
      ageseq3 = cbind(ageseq2,gotsome)
      ageseq3 = as.vector(ageseq3[ageseq3[,ncol(ageseq3)] == 1,-ncol(ageseq3)]);
      
      cn = 2;
      while (cn <= ncol(gmat)) {
        
        gmat[,cn] = fillin(gmat[,cn],ageseq3); 
        cn=cn+1;
      }
      
      if(MSsplit==0){
        
        name4 = "";
      }else if(MSsplit==1){
        if(iMStat==1){
          name4 = "m";
        }else if(IMStat==2){
          name4 = "f";
        }
      }
      
      fnamestr = paste0(name1,name2,name3,name4,name5)
      
      write.table(gmat, file = paste0(grphpath,fnamestr,".csv"),sep=",", row.names = FALSE,  col.names = FALSE)
      
      if((datatype*basecase) == 1){
        
        fnamestr = paste0(name1,name2,"bn",name4,name5)
        
        write.table(gmat, file = paste0(grphpath,fnamestr,".csv"), sep=",",row.names = FALSE,col.names = FALSE)
      }
      iMStat=iMStat+1;
    }
    iQunt=iQunt+1;
  }
}


getWmtx <- function(agedat96,PIdat,asstdat,MStatdat,obsdat,mxdat,mxobsdat,datawgts,optwgts,xtrasst){
  
  outcome = simqunts(PIdat, agedat96, asstdat, MStatdat, obsdat,pistate_a, cohorts_a, quants_a, mmtyrs, 0, datawgts);
  aqntdat = outcome$dataprfs;
  aqntcnts = outcome$datacnts;
  
  outcome = simqunts(PIdat, agedat96, mxdat, MStatdat, mxobsdat,pistate_m, cohorts_m, quants_m, mmtyrs, 0, datawgts);
  mxqntdat = outcome$dataprfs;
  mxqntcnts = outcome$datacnts;
  
  outcome = simqunts(PIdat, agedat96, mxdat, MStatdat, mxobsdat, pistate_m, cohorts_m, 0, mmtyrs, 0, datawgts);
  mxmnsdat = outcome$dataprfs;
  mxmnscnts = outcome$datacnts;
  
  outcome = simcrrl(PIdat, agedat96, log(mxdat+(1-mxobsdat)), MStatdat, mxobsdat, pistate_m, cohorts_m, mmtyrs, 0, datawgts);
  lmxmnsdat = outcome$meanprfs;
  lmxmnscnts= outcome$meancnts;
  lmxstddat= outcome$stdprfs;
  mxcrldat1= outcome$crrlprf1;
  mxcrlcnt1= outcome$crrlcnt1;
  mxcrldat2= outcome$crrlprf2;
  mxcrlcnt2= outcome$crrlcnt2;
  
  savemmts <<- 1;                     # Modify Global 
  
  outcome = makemmts(PIdat,agedat[,1],asstdat,MStatdat,obsdat,datawgts,aqntdat,
                     mxdat,mxobsdat,mxqntdat,mxmnsdat,lmxmnsdat,lmxstddat,mxcrldat1, mxcrldat2,0);
  mmtvec = outcome$mmtvec;
  obsvec = outcome$obsvec;
  qntvec = outcome$qntvec;
  pdfvec = outcome$pdfvec;
  mmttype = outcome$mmttype;
  
  savemmts <<- 0;
  
  VCV = compvinv(optwgts)$VCV;
  vinv = compvinv(optwgts)$vinv;
  W = compvinv(optwgts)$W;
  
  if(optwgts<2){
    
    Wdiag = diag(W);
    Wdiag = Wdiag + (xtrasst-1)*(mmttype==1)*Wdiag;
    diag(W) <- Wdiag
  }
  
  outcome = simqunts(PIdat, agedat96, asstdat, MStatdat, obsdat, pistate_a, cohorts_a, quants_a, mmtyrs, 1, datawgts);
  aqntdat = outcome$dataprfs;
  aqntcnts =outcome$datacnts;
  
  outcome = simqunts(PIdat, agedat96, mxdat, MStatdat, mxobsdat,pistate_m, cohorts_m, quants_m, mmtyrs, 1, datawgts);
  mxqntdat= outcome$dataprfs;
  mxqntcnts=outcome$datacnts;
  
  outcome = simqunts(PIdat, agedat96, mxdat, MStatdat, mxobsdat,pistate_m, cohorts_m, 0, mmtyrs, 1, datawgts);
  mxmnsdat= outcome$dataprfs;
  mxmnscnts=outcome$datacnts;
  
  grphmtx(aqntdat,1,1,0,qnum_a,pinum_a,chrtnum_a);   # Adjusts for composition bias
  grphmtx(mxqntdat,3,1,0,qnum_m,pinum_m,chrtnum_m);
  grphmtx(mxmnsdat,3,1,0,0,pinum_m,chrtnum_m);
  
  output = list("VCV" = VCV,"vinv"= vinv,"W" = W,"lmxmnsdat"= lmxmnsdat,"lmxstddat" = lmxstddat,
                "mmtvec" = mmtvec,"obsvec" = obsvec,"qntvec" = qntvec,"pdfvec" = pdfvec,"mmttype" = mmttype);
  
}



# RNDER:  Randomizes bigp to create a meaningful parameter simplex
# Rescales proportionally or additively according to rndtype

rnder <- function(isimplx,stder,stder2,parscale,rndtype){
  
  stder    = stder*parscale;
  cn       = length(stder);
  stder    = rep(1,cn)%*%t(stder);
  controll = (stder^2)/2;
  rndmat   = rbind(rep(0,cn),(stder*matrix(rnorm(cn^2, mean = 0, sd = 1),cn,cn)-controll));
  nusimplx = isimplx*exp(rndmat);
  
  stder2   = stder2*parscale;
  stder2   = rep(1,cn)%*%t(stder2);
  rndmat   = rbind(rep(0,cn),(stder2*matrix(rnorm(cn^2, mean = 0, sd = 1),cn,cn)));
  rndmat   = rndmat*(rep(1,cn+1)%*%t(rndtype));
  nusimplx = nusimplx+rndmat;
  
  return(nusimplx)
}

prnswtch <- function(){
  
  print("Switch Settings:");
  print(paste0("    Uncertain Mortality:             ", swchMort));
  print(paste0("    Beta fixed to equal 1:           ", swchBeta));
  print(paste0("    Non-asset income included:       ", swchY));
  print(paste0("    Deterministic Health Costs:      ", swchmxst));
  print(paste0("    Income Taxes:                    ", swchTax));
  print(paste0("    Persistent Health Cost Shocks:   ", swchZeta));
  print(paste0("    Transitory Health Cost Shocks:   ", swchXi));
  print(paste0("    Rate of Return Shocks:           ", swchROR));
  print(paste0("    Bequest Motive Activated:        ", swchBeq));
  print(paste0("    Gender-specific input profiles:  ", swchGdif));
  print(paste0("    Moments split by Marital Status: ", MSsplit));
  print(paste0("    Health Cost Variance Scaling:    ", mxvarscl));
  print(paste0("    Disaggregated Mortality Profiles:", mortdif));
}


#PUNSCALE_M:  Converts a vector of transformed parameters into the parameters used ARMA model of medex shocks
punscale_m <- function(shkparms){
  
  rhomx   = shkparms[1];
  fracar1 = shkparms[2];
  
  if(pscaled_m == 1){
    rhomx   = logit(rhomx);
    fracar1 = logit(fracar1);
  }
  output = list("rhomx" = rhomx, "fracar1" = fracar1);
}

#Print out the medex coefficients 
prnmedex <- function(rhomx,fracar1,mxcoef){
  
  print("Medex Coefficients");
  print(paste0("rho_mx          ",rhomx));
  print(paste0("fracar1         ",fracar1));
  
  print("                       Mean      Variance");
  j=1;
  while (j<=nrow(mxcoef)) {
    print(paste(mxlabel[j], "      ",paste(mxcoef[j,],collapse="           ")));
    j = j+1;
  }
}

#Calculate mean, standard deviation and 99.5% percentile of the medical expense.
mxprofs <- function(mnlnmxs,stdmxs,mnmx_pi,varmx_pi,pctle){
  
  pieffcts = replicate(2,(mnmx_pi[,1]*pctle+mnmx_pi[,2]*(pctle^2)));
  pieffcts = cbind(pieffcts,replicate(2,(mnmx_pi[,3]*pctle+mnmx_pi[,4]*(pctle^2))));
  mnmxs2   = mnlnmxs + pieffcts;
  pieffcts = replicate(2,(varmx_pi[,1]*pctle+varmx_pi[,2]*(pctle^2)));
  pieffcts = cbind(pieffcts,replicate(2,(varmx_pi[,3]*pctle+varmx_pi[,4]*(pctle^2))));
  stdmxs2  = sqrt(stdmxs^2 + pieffcts);
  
  mxmean   = exp(mnmxs2 + (stdmxs2^2)/2);
  mxstd    = sqrt( exp(2*mnmxs2 + stdmxs2^2)*(exp(stdmxs2^2)-1) );
  pct995   = 2.5758293;
  mx995    = exp(mnmxs2 + pct995*stdmxs2);
  
  output = list("mxmean"=mxmean,"mxstd"=mxstd,"mx995"=mx995);
}


#read in the medical expense data, clean and prepare the datasets. Then call function mxprofs to calculate the
#mean, standard deviation and 99.5% percentile of the medical expense, return these results. 
getmxtab <- function(mxcoef,pctle,smplyrs){
  
  if(smplyrs != 0){
    smplyrs = smplyrs-bornage+1;
  };
  
  #First for the means
  
  if(findmedex==0){
    
    datstr = paste0(datapath,"medexprofX.out");
    data <- read.delim(datstr, header=FALSE);
    outcome = getages(bornage,dieage,data);
    k = outcome$k;
    j = outcome$j;
    ib = outcome$ib;
    data = fixobs(data,bornage,dieage,k,j,0,-1e10);
    
    m_const = data[,2];
    m_badh  = data[,3];
    m_male  = data[,4];
    m_PI    = data[,5];
    m_PI2   = data[,6];
    v_const = data[,7];
    v_badh  = data[,8];
    v_male  = data[,9];
    v_PI    = data[,10];
    v_PI2   = data[,11];
    
  }else{
    
    m_const = mxcoef[1,1] + mxcoef[2,1]*ageseq + mxcoef[3,1]*(ageseq^2/100)+ mxcoef[4,1]*(ageseq^3/10000);
    m_badh  = mxcoef[5,1] + mxcoef[6,1]*ageseq;
    m_male  = mxcoef[7,1] + mxcoef[8,1]*ageseq;
    m_PI    = mxcoef[9,1] + mxcoef[10,1]*ageseq;
    m_PI2   = rep(mxcoef[11,1],TR);
    v_const = mxcoef[1,2] + mxcoef[2,2]*ageseq + mxcoef[3,2]*(ageseq^2/100)
    + mxcoef[4,2]*(ageseq^3/10000);
    v_badh  = mxcoef[5,2] + mxcoef[6,2]*ageseq;
    v_male  = mxcoef[7,2] + mxcoef[8,2]*ageseq;
    v_PI    = mxcoef[9,2] + mxcoef[10,2]*ageseq;
    v_PI2   = rep(mxcoef[11,2],TR);
  }
  
  mnmxs    = cbind((m_const[1:TR]+m_badh[1:TR]),m_const[1:TR]);
  mnmxs    = cbind((m_const+m_badh+m_male),(m_const+m_male),mnmxs);
  mnmx_pi  = cbind(m_PI[1:TR],m_PI2[1:TR]);
  mnmx_pi  = cbind(m_PI,m_PI2,mnmx_pi);
  stdmxs   = cbind((v_const[1:TR]+v_badh[1:TR]),v_const[1:TR]);
  stdmxs   = cbind((v_const+v_badh+v_male),(v_const+v_male),stdmxs);
  varmx_pi = cbind(v_PI[1:TR],v_PI2[1:TR]);
  varmx_pi = cbind(v_PI,v_PI2,varmx_pi);
  
  if(min(stdmxs)<0){     #penalize negative variances in GMM criterion
    negpen = 1;
  }else{
    negpen = 0;
  }
  
  mnmxs    = mnmxs + (mxvarscl==0)*discrtzn;   # Ad Hoc discretization adjustment 
  stdmxs   = abs(stdmxs);
  mnmxs    = mnmxs + (1-mxvarscl)*(stdmxs)/2;  # Mean-preserving adjustment 
  stdmxs   = sqrt(mxvarscl*stdmxs);
  mnmx_pi  = mnmx_pi + (1-mxvarscl)*(varmx_pi)/2;
  mnmx_pi  = mnmx_pi*exp((mxvarscl==0)*discrtzn);    # Ad Hoc discretization adjustment 
  varmx_pi = mxvarscl*varmx_pi;
  
  
  if(pctle != -1){
    
    outcome = mxprofs(mnmxs,stdmxs,mnmx_pi,varmx_pi,pctle);
    mxmean = outcome$mxmean;
    mxstd = outcome$mxstd;
    mx995 = outcome$mx995;
    
    # print(paste0(" Analytical Mean Health Care Costs at the ", pctle, " PI Percentile"));
    # print("  age     hbad      hgood      wbad      wgood");
    # if(smplyrs==0){
    #   print(cbind(ageseq,mxmean));
    # }else{
    #   print(cbind(ageseq[smplyrs],mxmean[smplyrs,]));
    # }
    # 
    # print(paste0(" Analytical Std. deviation Health Care Costs at the ", pctle, " PI Percentile"));
    # print("  age     hbad      hgood      wbad      wgood");
    # if(smplyrs==0){
    #   print(cbind(ageseq,mxstd));
    # }else{
    #   print(cbind(ageseq[smplyrs],mxstd[smplyrs,]));
    # }
    # 
    # print(paste0(" 99.5-th Percentile Annual Health Care Costs at the ", pctle, " PI Percentile"));
    # print("  age     hbad      hgood      wbad      wgood");
    # if(smplyrs==0){
    #   print(cbind(ageseq,mx995));
    # }else{
    #   print(cbind(ageseq[smplyrs],mx995[smplyrs,]));
    # }
  }
  
  output = list("mnmxs"=mnmxs,"stdmxs"=stdmxs,"mnmx_pi"=mnmx_pi,"varmx_pi"=varmx_pi);
}

#PUNSCALE_P:  Converts a vector of transformed parameters into the parameters used in the life-cycle model.

punscale_p <- function(parmvec){
  
  delta = parmvec[1];
  beta  = parmvec[2];
  nu     = parmvec[3];
  cfloor = exp(parmvec[4]);
  phi0   = parmvec[5];
  phi0   = max(c(min(c(phi0,1)),0.000001));
  K0     = parmvec[6]*1000;
  bigR   = 1+mu_r;
  phi0   = ((bigR*(1-phi0)/phi0)^nu)/(bigR*beta); 
  
  output= list("delta"= delta, "beta"= beta, "nu" = nu, "cfloor"= cfloor, "phi0" = phi0, "K0" = K0);
}




#bound the cash on hand at the consumption floor
fixcoh <- function(cohsim96){
  
  toopoor  = cohsim96<cfloor;               # Bound at consumption floor 
  cohsim96 = cohsim96*(1-toopoor) + toopoor*cfloor;
  
  print("Mean and std dev of cash-on-hand, year 1996 = ")
  print(cbind(colMeans(as.matrix(cohsim96)),colSds(as.matrix(cohsim96))));
  
  return(cohsim96);
}


#Prepare (stack rows/columns of matrixs to form vectors) and save relevent datasets in the iofiles folder,
#These files are going to be used in the C program. 
savevecs <- function(){
  
  agevec   = c(bornage,dieage,TR);
  simvec   = c(nn,simyrs,momyr1,momyr2);
  swchvec  = c(swchMort,swchBeta,swchY,swchmxst,swchTax,swchZeta,swchXi);
  swchvec  = c(swchvec,swchROR,swchBeq,swchGdif);
  medexvec = c(rhomx,fracar1,fracar1i,fracwn,minmedex);
  prefvec  = c(delta,beta,nu,phi0,K0);
  asstvec  = c(cfloor,tauBeq,exBeq,mu_r,sigma_r);
  
  
  mortprfs =  as.vector(t(as.matrix(mortprfs)));
  mort_pi =  as.vector(t(as.matrix(mort_pi)));
  hsprobs =  as.vector(t(as.matrix(hsprobs)));
  heal_pi =  as.vector(t(as.matrix(heal_pi)));
  yprof =  as.vector(t(as.matrix(yprof)));
  y_pi =  as.vector(t(as.matrix(y_pi)));
  mnlnmxs =  as.vector(t(as.matrix(mnlnmxs)));
  stdlnmxs =  as.vector(t(as.matrix(stdlnmxs)));
  mxpicoef =  c(as.vector(t(as.matrix(mnmx_pi))), as.vector(t(as.matrix(varmx_pi))));
  
  if(simtype==1){
    
    datlist = c("cohsim96", "ztacdfsim96", "pisim96", "agesim96", "healsimh", "healsimw","mstatsim","xicdfsim","epscdfsim");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
      assign(dat,data) ;
    }
  } else if(simtype==2){
    
    datlist = c("cohsimx", "ztacdfx", "pisimx", "agesimx", "hlsimhx", "hlsimwx","mstatx","xicdfx","epscdfx");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
      assign(dat,data) ;
    }
    
    cohsim96 = cohsimx[,1];  ztacdfsim96 = ztacdfx;  pisim96   = pisimx;
    agesim96 = agesimx[,1];  healsimh    = hlsimhx;  healsimw  = hlsimwx;
    mstatsim = mstatx;        xicdfsim    = xicdfx;   epscdfsim = epscdfx;
  }
  
  if(allalive==1){
    
    mstatsim = replicate(ncol(mstatsim),as.vector(mstatsim[,1]));
  }
  
  healsimh  = as.vector(as.matrix(healsimh));         #  simulated husband's health status 
  healsimw  = as.vector(as.matrix(healsimw));         #  simulated wife's health status   
  mstatsim  = as.vector(as.matrix(mstatsim));         #  simulated marital status         
  xicdfsim  = as.vector(as.matrix(xicdfsim));         #  simulated innovation on AR(1)    
  epscdfsim = as.vector(as.matrix(epscdfsim));        #  simulated white noise shock     
  cohsim96  = fixcoh(cohsim96);            #  Bound below by consumption floor  
  
  
  
  datlist = c("agevec", "simvec", "swchvec", "medexvec", "prefvec", "asstvec", "rorshk", "mortprfs", "mort_pi", "hsprobs",
              "heal_pi", "yprof", "y_pi", "mnlnmxs", "stdlnmxs","mxpicoef", "cohsim96", "ztacdfsim96","pisim96",
              "agesim96","healsimh","healsimw","mstatsim","xicdfsim","epscdfsim");
  
  for (i in 1:length(datlist)) {
    
    dat = datlist[i];
    
    write.table(get(dat), file = paste0(iopath,dat,".csv"),sep=",", row.names=FALSE, col.names =FALSE);
  }
}


getparms <- function(parmvec){
  
  if(pscaled_p==1){
    
    outcome = punscale_p(parmvec);
    delta <<- outcome$delta;
    beta <<- outcome$beta;
    nu <<- outcome$nu;
    cfloor <<- outcome$cfloor;
    phi0 <<-outcome$phi0;
    K0 <<- outcome$K0;
    #print('before swchBeq')
    #print('phi0')
    #print(phi0)
    #print('K0')
    #print(K0)
    
    
  }else if(pscaled_p==0){
    delta <<- parmvec[1];
    beta <<- parmvec[2];
    nu <<- parmvec[3];
    cfloor <<- parmvec[4];
    phi0 <<- parmvec[5];
    K0 <<- parmvec[6];
  };
  
  if(swchBeq==0){
    phi0 <<- 0;  
    K0 <<- 1; 
  };
  print('after swchBeq')
  print('phi0')
  print(phi0)
  print('K0')
  print(K0)
  
  if(prnres>0){
    uparams = c(delta,beta,nu,cfloor,phi0,K0);
    print("Transformed parameters, and parameters actually used:");
    rn=1;
    while (rn<=length(parmvec)) {
      print(paste(plabel[rn]," ",parmvec[rn], "      ", uparams[rn]));
      rn = rn+1;
    }
  }
  
}

getcrit <- function(parmvec){
  
  getparms(parmvec);
  savevecs();                 # Save input vector for C program 
  
}


oneloop <- function(allparms){
  
  allparms = allparms*zerovec + fixvals*(1-zerovec);
  
  outcome = punscale_m(allparms[1:2]);
  rhomx <<- outcome$rhomx;
  fracar1 <<- outcome$fracar1;
  
  fracar1i <<- fracar1*(1-rhomx^2);    # sigma of AR(1) innovations: global  
  fracwn   <<- 1-fracar1;              # sigma of white noise: global 
  mxmcoef  = allparms[3:13];
  mxvcoef  = allparms[14:24];
  mxcoef   = cbind(mxmcoef,mxvcoef);
  
  #  prnmedex(rhomx,fracar1,mxcoef);
  
  outcome = getmxtab(mxcoef,0.4,smplyrs);
  
  mnlnmxs <<- outcome$mnmxs;
  stdlnmxs <<- outcome$stdmxs;
  mnmx_pi <<- outcome$mnmx_pi;
  varmx_pi <<- outcome$varmx_pi;
  
  parmvec  = allparms[25:30]; 
  criter   = -100;
  
  if(job != 4){
    
    criter  = getcrit(parmvec);
  }else{
    print("WRONG JOB!!!");
  }
  
  return(criter);
}


loadsim <- function(){
  
  datlist = c("cohsim", "ztasim", "ztaindexsim", "xisim", "xiindexsim", "Medicaidsim", "medexsim", "conssim", 
              "beqsim", "mssim2", "asstsim","netIncomesim", "transfersim");
  
  for (i in 1:length(datlist)) {
    
    dat = datlist[i];
    data = read.csv(file = paste0(iopath,dat,".csv"), header = FALSE);
    assign(dat,data) ;
  }
  
  cohsim = matrix(as.matrix(cohsim),nn,simyrs+1);
  asstsim = matrix(as.matrix(asstsim),nn,simyrs+1);
  ztasim = matrix(as.matrix(ztasim),nn,simyrs+1);
  xisim = matrix(as.matrix(xisim),nn,simyrs+1);
  medexsim = matrix(as.matrix(medexsim),nn,simyrs+1);
  ztaindexsim = matrix(as.matrix(ztaindexsim),nn,simyrs+1);
  xiindexsim = matrix(as.matrix(xiindexsim),nn,simyrs+1);
  Medicaidsim = matrix(as.matrix(Medicaidsim),nn,simyrs+1);
  transfersim = matrix(as.matrix(transfersim),nn,simyrs+1);
  medex1yrsim = medexsim;
  conssim = matrix(as.matrix(conssim),nn,simyrs+1);
  beqsim = matrix(as.matrix(beqsim),nn,simyrs+1);
  netIncomesim = matrix(as.matrix(netIncomesim),nn,simyrs+1);
  mssim2 = matrix(as.matrix(mssim2),nn,simyrs+1);
  alivesim = mssim2>0;
  aliveavg = colMeans(alivesim);
  
  if(simtype==1){
    
    datlist = c("asim96", "mxsim96", "incsim96", "cohsim96");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
      assign(dat,data) ;
    }
  }else if(simtype==2){
    
    datlist = c("asimx", "incsim96x", "mxsim96x", "cohsimx");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
      assign(dat,data) ;
    }
    
    asim96 = asimx[,1];  
    mxsim96 = mxsim96x[,1];  
    incsim96 = incsim96x[,1];
    cohsim96 = cohsimx[,1];
  }
  
  asstsim[,1]  = as.matrix(asim96);
  medexsim[,1] = as.matrix(mxsim96);
  cohsim[,1]   = as.matrix(cohsim96);
  netIncomesim[,1] = as.matrix(incsim96);
  medexsim = (medexsim[,1:simyrs]+medexsim[,2:(simyrs+1)])/2; # 2-year averages 
  medexsim = cbind(mxsim96,medexsim);
  toosmall = medexsim<minmedex;   # bottom coding 
  medexsim = medexsim*(1-toosmall) + toosmall*minmedex;
  
  if(prnres>1){
    
    simavg = cbind(aliveavg,(colMeans(cohsim*alivesim)/aliveavg),(colMeans((ztasim+xisim)*alivesim)/aliveavg));
    simavg = cbind(simavg,(colMeans(Medicaidsim*alivesim)/aliveavg));
    simavg = cbind(simavg,(colMeans(medexsim*alivesim)/aliveavg),(colMeans(conssim*alivesim)/aliveavg));
    simavg = cbind(simavg,(colMeans(beqsim*alivesim)/aliveavg),(colMeans(asstsim*alivesim)/aliveavg));
    simavg = cbind(simavg,(colMeans(netIncomesim*alivesim)/aliveavg),(colMeans((Medicaidsim+medex1yrsim)*alivesim)/aliveavg));
    
    print("Year Survival cash-o-h zeta+xi Medicaid  h. costs  consumption bequests assets  net income tot 1yr mx");
    
    print(cbind(seq(momyr1,momyr1+simyrs, by=1),simavg))
    
  }
  
  cohsim   = cohsim[,1:simyrs];
  asstsim  = asstsim[,1:simyrs];
  ztasim   = ztasim[,1:simyrs];
  medexsim = medexsim[,1:simyrs];
  conssim  = conssim[,1:simyrs];
  mssim2   = mssim2[,1:simyrs];
  beqsim   = beqsim[,1:simyrs];
  alivesim = alivesim[,1:simyrs];
  netIncomesim = netIncomesim[,1:simyrs];
  medex1yrsim  = medex1yrsim[,1:simyrs];
  Medicaidsim  = Medicaidsim[,1:simyrs];
  transfersim  = transfersim[,1:simyrs];
  alivesim = alivesim*1
  
  datlist = c("cohsim", "asstsim", "ztasim", "ztaindexsim", "xisim", "xiindexsim", 
              "Medicaidsim", "medexsim", "medex1yrsim", "conssim", "beqsim", 
              "mssim2", "alivesim", "netIncomesim", "transfersim");
  
  for (i in 1:length(datlist)) {
    
    dat = datlist[i];
    
    write.table(get(dat), file = paste0(shkpath,dat,".csv"), sep=",", row.names=FALSE, col.names = FALSE);
  }
  
  if(simtype==1){
    
    datlist = c("cohsim", "asstsim", "ztasim", "ztaindexsim", "xisim", "xiindexsim", 
                "Medicaidsim", "medexsim", "medex1yrsim", "conssim", "beqsim", 
                "mssim2", "alivesim", "netIncomesim");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      
      write.table(get(dat), file = paste0(shkpath,dat,".csv"), sep=",", row.names=FALSE, col.names = FALSE);
    }
    
  }else if(simtype==2){
    
    cohsimx   = cohsim;       asstsimx = asstsim;   ztasimx = ztasim;
    ztaindx   = ztaindexsim;  xisimx   = xisim;     xiindx  = xiindexsim;
    Medicaidx = Medicaidsim;  mxsimx   = medexsim;  mx1simx = medex1yrsim;
    transferx = transfersim;  conssimx = conssim;   beqsimx = beqsim;    
    mssim2x   = mssim2;       alivex   = alivesim;  incsimx = netIncomesim;
    
    datlist = c("cohsimx", "asstsimx", "ztasimx", "ztaindx", "xisimx", "xiindx",
                "Medicaidx", "transferx", "mxsimx", "mx1simx", "conssimx", "beqsimx", 
                "mssim2x", "alivex", "incsimx");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      
      write.table(get(dat), file = paste0(shkpath,dat,".csv"),sep=",",row.names=FALSE, col.names = FALSE);
    }
    
  }
  
}

loadsim2 <- function(mean_asset,stdev_asset){
  
  datlist = c("cohsim", "ztasim", "ztaindexsim", "xisim", "xiindexsim", "Medicaidsim", "medexsim", "conssim", 
              "beqsim", "mssim2", "asstsim","netIncomesim", "transfersim");
  
  for (i in 1:length(datlist)) {
    
    dat = datlist[i];
    data = read.csv(file = paste0(iopath,dat,".csv"), header = FALSE);
    assign(dat,data) ;
  }
  
  cohsim = matrix(as.matrix(cohsim),nn,simyrs+1);
  asstsim = matrix(as.matrix(asstsim),nn,simyrs+1);
  ztasim = matrix(as.matrix(ztasim),nn,simyrs+1);
  xisim = matrix(as.matrix(xisim),nn,simyrs+1);
  medexsim = matrix(as.matrix(medexsim),nn,simyrs+1);
  ztaindexsim = matrix(as.matrix(ztaindexsim),nn,simyrs+1);
  xiindexsim = matrix(as.matrix(xiindexsim),nn,simyrs+1);
  Medicaidsim = matrix(as.matrix(Medicaidsim),nn,simyrs+1);
  transfersim = matrix(as.matrix(transfersim),nn,simyrs+1);
  medex1yrsim = medexsim;
  conssim = matrix(as.matrix(conssim),nn,simyrs+1);
  beqsim = matrix(as.matrix(beqsim),nn,simyrs+1);
  netIncomesim = matrix(as.matrix(netIncomesim),nn,simyrs+1);
  mssim2 = matrix(as.matrix(mssim2),nn,simyrs+1);
  alivesim = mssim2>0;
  aliveavg = colMeans(alivesim);
  
  if(simtype==1){
    
    datlist = c("asim96", "mxsim96", "incsim96", "cohsim96");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
      assign(dat,data) ;
    }
  }else if(simtype==2){
    
    datlist = c("asimx", "incsim96x", "mxsim96x", "cohsimx");
    
    for (i in 1:length(datlist)) {
      
      dat = datlist[i];
      data = read.csv(file = paste0(shkpath,dat,".csv"), header = FALSE);
      assign(dat,data) ;
    }
    
    asim96 = asimx[,1];  
    mxsim96 = mxsim96x[,1];  
    incsim96 = incsim96x[,1];
    cohsim96 = cohsimx[,1];
  }
  
  asstsim[,1]  = as.matrix(asim96);
  medexsim[,1] = as.matrix(mxsim96);
  cohsim[,1]   = as.matrix(cohsim96);
  netIncomesim[,1] = as.matrix(incsim96);
  medexsim = (medexsim[,1:simyrs]+medexsim[,2:(simyrs+1)])/2; # 2-year averages 
  medexsim = cbind(mxsim96,medexsim);
  toosmall = medexsim<minmedex;   # bottom coding 
  medexsim = medexsim*(1-toosmall) + toosmall*minmedex;
  
  if(prnres>1){
    
    simavg = cbind(aliveavg,(colMeans(cohsim*alivesim)/aliveavg),(colMeans((ztasim+xisim)*alivesim)/aliveavg));
    simavg = cbind(simavg,(colMeans(Medicaidsim*alivesim)/aliveavg));
    simavg = cbind(simavg,(colMeans(medexsim*alivesim)/aliveavg),(colMeans(conssim*alivesim)/aliveavg));
    simavg = cbind(simavg,(colMeans(beqsim*alivesim)/aliveavg),(colMeans(asstsim*alivesim)/aliveavg));
    simavg = cbind(simavg,(colMeans(netIncomesim*alivesim)/aliveavg),(colMeans((Medicaidsim+medex1yrsim)*alivesim)/aliveavg));
    
    #print("Year Survival cash-o-h zeta+xi Medicaid  h. costs  consumption bequests assets  net income tot 1yr mx");
    
    #print(cbind(seq(momyr1,momyr1+simyrs, by=1),simavg))
    
  }
  
  cohsim   = cohsim[,1:simyrs];
  asstsim  = asstsim[,1:simyrs];
  ztasim   = ztasim[,1:simyrs];
  medexsim = medexsim[,1:simyrs];
  conssim  = conssim[,1:simyrs];
  mssim2   = mssim2[,1:simyrs];
  beqsim   = beqsim[,1:simyrs];
  alivesim = alivesim[,1:simyrs];
  netIncomesim = netIncomesim[,1:simyrs];
  medex1yrsim  = medex1yrsim[,1:simyrs];
  Medicaidsim  = Medicaidsim[,1:simyrs];
  transfersim  = transfersim[,1:simyrs];
  alivesim = alivesim*1
  
# We set a data table with the simulated data
  ID <- seq(1,nn, by=1)
  DATASIM <- data.table(ID)
  
  # Variables by assetsim
  DATASIM[, asset_1996 := asstsim[,1]]
  DATASIM[, asset_1998 := asstsim[,3]]
  DATASIM[, asset_2000 := asstsim[,5]]
  DATASIM[, asset_2002 := asstsim[,7]]
  DATASIM[, asset_2004 := asstsim[,9]]
  DATASIM[, asset_2006 := asstsim[,11]]
  
  # Variables by mssim
  DATASIM[, alive_1996 := (mssim2[,1]>0)*1]
  DATASIM[, alive_1998 := (mssim2[,3]>0)*1]
  DATASIM[, alive_2000 := (mssim2[,5]>0)*1]
  DATASIM[, alive_2002 := (mssim2[,7]>0)*1]
  DATASIM[, alive_2004 := (mssim2[,9]>0)*1]
  DATASIM[, alive_2006 := (mssim2[,11]>0)*1]
  
  # Variables by assetsim
  DATASIM[, asset_1996 := ((asset_1996-mean_asset)/stdev_asset)*alive_1996]
  DATASIM[, asset_1998 := ((asset_1998-mean_asset)/stdev_asset)*alive_1998]
  DATASIM[, asset_2000 := ((asset_2000-mean_asset)/stdev_asset)*alive_2000]
  DATASIM[, asset_2002 := ((asset_2002-mean_asset)/stdev_asset)*alive_2002]
  DATASIM[, asset_2004 := ((asset_2004-mean_asset)/stdev_asset)*alive_2004]
  DATASIM[, asset_2006 := ((asset_2006-mean_asset)/stdev_asset)*alive_2006]
  
  
  # DATASIM[is.nan(asset_1996)==1, asset_1996 := 0]
  # DATASIM[is.nan(asset_1998)==1, asset_1998 := 0]
  # DATASIM[is.nan(asset_2000)==1, asset_2000 := 0]
  # DATASIM[is.nan(asset_2002)==1, asset_2002 := 0]
  # DATASIM[is.nan(asset_2004)==1, asset_2004 := 0]
  # DATASIM[is.nan(asset_2006)==1, asset_2006 := 0]
  
  
  # gender (female is 1)
  DATASIM[, gender:=(mssim2[,1]==2)*1]
  
  # we need to add ht, I, age
  agesim = read.csv(file = paste0(shkpath,"agesim",".csv"), header = FALSE);
  
  DATASIM[, age_1996 := log(agesim[,1])]
  DATASIM[, age_1998 := agesim[,3]]
  DATASIM[, age_2000 := agesim[,5]]
  DATASIM[, age_2002 := agesim[,7]]
  DATASIM[, age_2004 := agesim[,9]]
  DATASIM[, age_2006 := agesim[,11]]
  
  healthsimh = read.csv(file = paste0(shkpath,"healsimh",".csv"), header = FALSE);
  healthsimw = read.csv(file = paste0(shkpath,"healsimw",".csv"), header = FALSE);
  
  
  healthsimh2 <- cbind(healthsimh[,1],healthsimh[,3],healthsimh[,5],healthsimh[,7],healthsimh[,9],healthsimh[,11])
  healthsimw2 <- cbind(healthsimw[,1],healthsimw[,3],healthsimw[,5],healthsimw[,7],healthsimw[,9],healthsimw[,11])
  
  Isim<- dim(healthsimh2)[1]
  TT <-dim(healthsimh2)[2]
  TT2 <- TT-1
  
  # Fix up the observations for when the individuals die
  # We fix the health status by making that two periods after death the h's shoudl be 0 (one period after death seems to be populated in the real data)
  
  alive_mat <- cbind(DATASIM$alive_1996,DATASIM$alive_1998,DATASIM$alive_2000,DATASIM$alive_2002,DATASIM$alive_2004,DATASIM$alive_2006)
  obs_health_mat <- alive_mat
  # for(ii in 1:Isim){
  #   alive_aux = alive_mat[ii,1]
  #   tt = 1
  #   while((alive_aux == 1)*(tt<=TT2))
  #     {
  #       tt = (tt+1)
  #       alive_aux = alive_mat[ii,tt]
  #     }
  #   obs_health_mat[ii,tt] = 1
  # }
  
  healthsimw2 = healthsimw2*obs_health_mat
  healthsimh2 = healthsimh2*obs_health_mat
  
  
  
  DATASIM[, health_1996 := (healthsimw2[,1]*gender+healthsimh2[,1]*(1-gender))]
  DATASIM[, health_1998 := (healthsimw2[,2]*gender+healthsimh2[,2]*(1-gender))]
  DATASIM[, health_2000 := (healthsimw2[,3]*gender+healthsimh2[,3]*(1-gender))]
  DATASIM[, health_2002 := (healthsimw2[,4]*gender+healthsimh2[,4]*(1-gender))]
  DATASIM[, health_2004 := (healthsimw2[,5]*gender+healthsimh2[,5]*(1-gender))]
  DATASIM[, health_2006 := (healthsimw2[,6]*gender+healthsimh2[,6]*(1-gender))]
  
  DATASIM[, health_bar:= (health_1996 + health_1998 + health_2000 + health_2002 + health_2004 + health_2006)/(alive_1996 + alive_1998 + alive_2000 + alive_2002 + alive_2004 +alive_2006)]
  
  pisim = read.csv(file = paste0(shkpath,"pisim96",".csv"), header = FALSE);
  DATASIM[, PI  := pisim]
  
  
  DATASIM[,PI2:= 0]
  DATASIM[,PI3:= 0]
  DATASIM[,PI4:= 0]
  DATASIM[,PI5:= 0]
  
  DATASIM[PI > .2 & PI <= .4, PI2:= 1]
  DATASIM[PI > .4 & PI <= .6, PI3:= 1]
  DATASIM[PI > .6 & PI <= .8, PI4:= 1]
  DATASIM[PI > .8, PI5:= 1]
  
  DATASIM[PI < .2, PIq:= 1]
  DATASIM[PI > .2 & PI <= .4, PIq:= 2]
  DATASIM[PI > .4 & PI <= .6, PIq:= 3]
  DATASIM[PI > .6 & PI <= .8, PIq:= 4]
  DATASIM[PI > .8, PIq:= 5]
  
  
  return(DATASIM)
  
}
