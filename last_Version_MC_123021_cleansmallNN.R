library(tensorflow)
library(matrixStats)
library(parallel)
library(MASS)
require(Matrix)
require(boot)
require(data.table)
require(caTools)
require(magrittr)
require(keras)
require(pryr)
library(data.table)
library(gtools)
#library(tictoc)
library(data.table)

library(reticulate)
conda_list()
use_condaenv("python_empirical_cpu2")


#set.seed(2)
set.seed(3)

## This is a Monte Carlo exercise to assess the gains of the GAN approach
S = 100
#num_s = 1

num_s <- commandArgs(trailingOnly = TRUE)
num_s <- strtoi(num_s)

print("MC replication is")
print(num_s)


rep_num = 5



#for(num_s in 4:S){

mat_seed = sample.int(10000, S)
set.seed(mat_seed[num_s])

best_mat = matrix(0,rep_num, 15)

############################# Optimization parameters ##############################################

# compute derivative
#default: 0.1, zoom = 0.01
hstep = 0.01
# learning rate: default = 1, zoom = 0.05
step0 = .01
# momentum
gamma_g0 = 0
gamma_g = gamma_g0

# type of optimization
# stochastic gradient descent
stoch_grad = 1
# momentum
mom = 1
# accelerated momentum
acc_mom = 1


epoch_param_loss = 100
################   MAIN ##################

#change the working directory
useMPI = 1                 # 1 => mac; 0 => windows


#---------------------Set Optimization/Simulation Globals--------------------#







#---------------------Set Optimization/Simulation Globals--------------------#


W = 0;

job        =     5;        # 1 => estimating and getting se's, all parms;
# 2 => two-stage estimation and se's, all parms;
# 3 => getting se's only;                       
# 4 => estimating, medex parms only;            
# 5 => get graphs;  6 => experiments;
# 7 => single cohort, long time period.

fixparams  = c(1,1,1,1);   # fixparam[j]==0 => item[j] fixed in estimation; 
# item 1 <-> _delta;       item 2 <-> _beta;    
# item 3 <-> _nu (1/IES);  item 4 <-> _cloor.

basecase   =     1;        # Output saved as baseline for future graphs 
pscaled_p  =     1;        # 1=> parameters scaled to help optimization
pscaled_m  =     1;
rshktype   =     0;        # 0 => none, 1 => EP (2007) diffs , 2 => Blau
# (JOLE, 2008) diffs, 3 => EP levels, 4 => Blau levels
#nn        =  (2688*2);        # Number of simulations#
nn = 2688
newmat     =     1;        # 0 => use previous set of health and wage shocks;
# 1 => create new matrices. Set to 1, except when doing experiments
simtype    =     1;        # 1=> all ages; 2=> specific age
hsimtype   =     2;        # 1 => simulate entire health status/mortality trajectory from initial conditions
# 2=> use data trajectories w/ imputations
allalive   =     0;        # 1=> nobody dies in simulations (and only sims)
prnres     =     2;        # 1 => print simulation means; 0 => don't print
prngrph    =     1;        # 1 => update graphs; 0 => don't update
prnout     =     1;        # 1 => call output on/off
prnout = ifelse(useMPI >0,0,prnout)

optwgts    =     2;        # 0 => use identity matrix as the weighting matrix;
# 1 => use diagonal, non-identity matrix; 2 => GMM optimal matrix 
xtrasst    =     1;        # Weight on asset-related moments
savemmts   =     0;        # 1=> save matrix of moment variables 
wgtddata   =     0;        # 1=> use AHEAD data weights 
findmedex  =     0;        # 0=> load medex parms from data; 1=> estimate them 
medexmmts  =     1;        # 1=> include medex moments in criterion function
medexmmts  = medexmmts*findmedex;

onesim     = rep(1,nn);

pinum_a    =     5;

if(pinum_a>1) {            # PI quantiles used in analyses:  assets
  pistate_a = seq(1/pinum_a,1-1/pinum_a,1/pinum_a)
}

qnum_a     =     2;
quants_a   = seq(1/qnum_a,1-1/qnum_a,1/qnum_a);  # Asset quantiles within each PI quantile
quants_a   = 0.5;                                 
qnum_a     = length(quants_a);

pinum_m    =     2;

if(pinum_m>1) {            # PI quantiles used in analyses:  medex
  pistate_m = seq(1/pinum_m,1-1/pinum_m,1/pinum_m)
}

qnum_m     =     1;
quants_m   =   0.9;            # medex quantiles within Each PI quantile 
qnum_m     = length(quants_m);

chktie     =     0;        # 0=> quantiles might split obs with same values 
mvcode     =  -1e10;       # Missing value code used for quantile sorts 
MSsplit    =     0;        # 1=> moments split by male vs. female 
MSnum = ifelse(MSsplit==0, 1, 2)

mortdif    =     1;        # 0=> aggregated mort profiles; 1=> disaggregated 
coherr     =     0;        # Variance of measurement error in coh vs. assets 
assterr    =     0;        # Variance of measurement error in assets 
minmedex   =   250;        # Bottom-coding for medical expenditures 
discrtzn   = log(0.942);    # Discretization Adjustment for zero-var experiments 
smplyrs    =     0;        # Rows displayed in medex table; 0 => all rows 
mean_asset = 134818.9;
stdev_asset = 263023.6;

#--------------------------Age/Time Parameters---------------------------#
bornage    =    70;        # First year 
dieage     =   102;
TR        = dieage-bornage+1;  # Number of periods in a life 
momyr1     =  1996;        # Calendar year at which we start matching moments 
momyr2     =  2006;        # Calendar year at which we stop matching moments 
simyrs     = momyr2-momyr1+1;   # Number of periods for matching moments 
ageshft    = momyr1-1994;             # Data start at calendar year 1994 

ageseq     = seq(bornage,dieage,1); 
yearseq    = seq(momyr1,momyr2,1);

mmtcols    = c(1996,1998,2000,2002,2004,2006);  # Years where we have data 

if(ageshft==0){
  mmtcols = c(1994, mmtcols)
}
mmtcols    = mmtcols - momyr1 +1;             # convert to column indices 
mmtyrs     = length(mmtcols);

cohorts_a  = c(c(69,74,79,84,89),dieage);           # Cohorts defined in terms of 1994 ages 
chrtnum_a  = length(cohorts_a)-1;
cohorts_m  = c(c(69,74,79,84,89),dieage);
chrtnum_m  = length(cohorts_m)-1;
cohorts    = intersect(cohorts_a,cohorts_m);
chrtnum    = length(cohorts)-1;

exprage  = c(70,74);              # Initial 1994 range for policy experiments 
trexpr  = dieage-max(exprage+ageshft)+1;
simyrsx  = trexpr;
ageseqx  = seq(mean(exprage+ageshft),mean(exprage+ageshft)+trexpr-1,1);
exprcols = ageseqx-mean(exprage+ageshft)+1;
cohortsx = exprage - c(1,0);
chrtnumx = length(cohortsx)-1;

cellmin    = 10;          # minimum # of obs to be included in criterion 
quantmin   = 10;          # minimum # of obs needed for quantile calcs 

#-----------------------Switches for C executables-----------------------#
swchMort   = 1;                    # 0=> no shocks; 1=> mortality shocks 
swchBeta   = 0;                       # 0=> use GAUSS value; 1=> _beta=1 
swchY      = 1;           # 0=> noIncome; 1=> use loaded income profiles 
swchmxst   = 1;      # 0=> no deterministic health cost; 1=> Health cost 
swchTax    = 1;                   # 0=> noIncome taxes; 1=> income taxes 
swchZeta   = 1;        # 0=> no shocks; 1=> health cost shocks from zeta 
swchXi     = 1;          # 0=> no shocks; 1=> health cost shocks from zi 
swchROR    = 0;                # 0=> no shocks; 1=> interest rate shocks 
swchBeq    = 1;                             # 0=> no; 1=> bequest motive 
swchGdif   = 1;                                # 0=> no; 1=> gender diff 

#-----------------------Asset/Financial Variables------------------------#
#cfloor=8.406262
#cfloor     = 8.261902;

#cfloor = 7.88
#cfloor = 7.768956
asstmax    = 2000000;
asstmin    = 0;
tauBeq     = 0.15;                                     # Estate tax rate 
exBeq      = 6e5;                           # Estate tax exemption level 
mu_r       = 1/0.98-1;                                # Mean net ROR 
bigR   = 1+mu_r;

if(rshktype==3){
  mu_r = 0.0433                # Predicted ROR 1997-2005 (per EP, 2007) 
} else if(rshktype==4){
  mu_r = 0.0093                # Mean ROR 1976-2008 (per Blau, 2008) 
}

sigma_r    = 0.0295;          # Std deviation net ROR 1989-2008 (per Blau, 2008) 

if(rshktype%%2==1){
  sigma_r = 0.0229            # Std deviation net ROR 1960-1995 (per EP, 2007) 
}

#------------------------------------------------------------------------#
#              Discretize the state/choice variables                     #
#------------------------------------------------------------------------#

#-----------------------------Health Status------------------------------#

hsnum      = 2;                  # 0 => bad health, 1   => good health 



#---------------------------------Paths----------------------------------#
datapath =  "data/";
shkpath  =  "shk/";
iopath   =  "iofiles/";
#grphpath =  "C:/wealth_singles - Copy/graphs/";
#rulecall =  "C:/wealth_singles - Copy/ccode/wealth20_singles";


source("external1_fix_copy.R")
source("external2_gendersmallNN.R")
source("external_nacho.R")
source("f_ea_trainingsamples.R")


if(.Platform$OS.type == "windows"){
	path_exe = "C:/wealth_singles - MC/ccode/wealth21_singles/test1/x64/Release/test1.exe"
}else{
	path_exe = "./ccode/wealth21_singles/test1.exe" 
}


###################   MAIN ############################################
#------------------Load up the income, etc. profiles-------------------#

outcome = getprofs();
mortprfs = outcome$morts;
hsprobs = outcome$healdats;
yprof = outcome$yprof;
mort_pi = outcome$mort_pi;
heal_pi = outcome$heal_pi;
y_pi = outcome$y_pi;
rorshk = outcome$rorshk;

# mortprf*:   mortality rates for singles (s)                         
# hsprob*:    health status transition rates for singles (s)          
# yprob*:     income profiles                                          
# picoefs:    coefficients on Permanent Income and PI squared         
# rorshk:     rate of return shocks, by calendar year      


#getchrt: Specify the corhorts for each observation in the data, 
#and record in a table the number of observations in each cohort
#e.g. 1     70-74     730
#     2     75-79     764
#     3     80-84     658
#     4     85-89     372
#     5     90-102    159 



#------------------Load up the datasets-------------------#
outcome = getdata();
agedat = outcome$agedat;
PIdat = outcome$PIdat;
asstdat = outcome$asstdat;
beqdat = outcome$beqdat;
MStatdat = outcome$MStatdat;
obsdat = outcome$obsdat;
mxdat = outcome$mxdat;
mxobsdat = outcome$mxobsdat;
hsdat = outcome$hsdat;
hsobsdat = outcome$hsobsdat;
totobs = outcome$totobs;
datawgts = outcome$datawgts;
agedat96 = outcome$agedat96;
avgage96 = outcome$avgage96;
HHIDdat = outcome$HHIDdat;



# real parameters of MC
cfloor0  = 8.411833;
delta = 0;
beta  = 0.97106
nu = 3.8
phi0 = 0.3
K0 =  10

thres = 0 

param_set = list(cfloor0,delta,beta,nu,phi0,K0,thres)




############################# Generating all the uncertainty ############################
#Generate the health and mortality shocks
if(newmat==1){            # 1 => load in new matrices 
  initdist(0);
  if(hsimtype==1){
    initsim(mortprfs, mort_pi, hsprobs, heal_pi)
  }else if(hsimtype==2){
    initsim2(mortprfs, mort_pi, hsprobs, heal_pi);
  }
}

# print a summary of the shocks
sumshks();
grid_solution = matrix(0,9,9)
random_num = matrix(rnorm(10),10,1)



#-------------------Initialize preference parameters-------------------#
pscaled_p = 1;    # 1 => Parameters transformed to aid optimization 
pscaled_m = 1;

plabel = c("_delta (health shift)    ", "_beta (discount rate)    ", "nu (1/IES)               ", 
           "cfloor (cons. floor)     ", "phi (bequest scale)      ", "K  (bequest curvature)   ");

#delta = 0.009;       # Utility shift due to good health
#delta = 0;       # Utility shift due to good health

# beta0  = 0.97106;  # Discount factor
#beta0 = 0.99;
# beta = beta0;
# nu    = 5.587;  # 1/IES
#nu = 2.9301721 
#nu =7 
#nu    = 3.81
# 
# # Bequest utility parameters                                           
# # Phi_0(b_net) = phi_0*[(b_net+K_0)^(1-nu)]/[1-nu]                    
# 
#phi0   =    0.2999829;  # MPC out of terminal period wealth:  1=> NO bequest motive
#K0     = 1.2084281
#phi0 = 1;
#K0 = 100;
#phi0   =    0.15;  # MPC out of terminal period wealth:  1=> NO bequest motive
#K0 = 123;

#-----------------------Health shock parameters------------------------#
#The process for logged, detrended health costs is the sum of an AR(1)-#
#-and a white noise process.  The processes are normalized to a total--#
#--variance of one.  Note that we have to integrate over each process--#
#----------and that the AR(1) component is a state variable.-----------#

mxlabel = c("Intercept       ", "Age             ", "A^2/100         ","A^3/10000       ", 
            "Bad Health      ", "Bad Health x Age", "Male            ", "Male x Age      ", 
            "PI (percentile) ", "PI x Age        ", "PI^2            ");

mxvec=c(
  2.41663,
  -0.590504,
  -5.6576,
  0.371098,
  -0.423898,
  0.182864,
  -0.480791,
  0.00585711,
  0,
  0,
  0.491667,
  0.0236076,
  -1.42724,
  28.0464,
  -0.624084,
  0.355875,
  0,
  -2.113,
  0.0331839,
  0,
  0,
  -5.87487,
  0.101302,
  -1.62123);

rn_m = length(mxvec);

if(findmedex==0){
  rhomx    = logitrv(0.9215);   # From the `standard' AR(1), transformed for logit 
  fracar1  = logitrv(0.5236/(0.5236+1.039));    
  mxvarscl = (1.039+0.5236)/(1.0476^2);
  
  mxmcoef  = rep(0,11);
  mxvcoef  = rep(0,11);
  
  mxvec    = c(rhomx,fracar1,mxmcoef,mxvcoef);
  zerovec  = rep(0,rn_m);
}else{
  mxvarscl = 1;
  zerovec  = rep(1,rn_m);
  zerovec[9:10]  = rep(0,2);  #male coefficients, mean 
  zerovec[20:21] = rep(0,2);  #male coefficients, variance 
  zerovec[17] = 0;
}

parmvec  = c(delta,beta,nu,cfloor0,phi0,K0);
parmvec_grd = parmvec
#rn_p     = length(parmvec);

allparms = c(mxvec,parmvec);
fixvals  = allparms;
zerovec  = c(zerovec,fixparams);

if(swchBeq == 0){
  
  zerovec = c(zerovec,rep(0,2));
}else if(swchBeq == 1){
  
  zerovec = c(zerovec,rep(1,2)); 
}

zerovec0 = zerovec
#prnswtch();


oneloop(allparms)

system2(path_exe, invisible = FALSE)


data_real <- loadsim2(mean_asset,stdev_asset)

# we add the labeling
data_real[,d:=1]
N <- length(data_real$d)

N0 <- nn
data_real[,d:=1]

# How many simulated samples we need:
nn = nn*40


# We augment the real data to obtain balanced labels: we do a simple bootstrap approach: we resample with replacement
augm <-floor(nn/N0)
sample_real <- sample(N0,replace=TRUE)
data_aug = data_real[sample_real,]
#data_aug <- data_real
if(augm > 1){
  for(ll in (1:(augm-1))){
    print(ll)
    sample_real <- sample(N0,replace=TRUE)
    data_aug_aux = data_real[sample_real,]
    data_aug <- as.data.table(rbind(data_aug,data_aug_aux))
  }
}

data_real <- data_aug

## we select the population of individuals with at least thres periods
NN <- data_real[,.N]
data_real[,ID:=seq(1,NN,1)]
data_real[, sum_alive:= alive_1996+alive_1998 + alive_2000 + alive_2002 + alive_2004 + alive_2006, by=ID]
data_real <- data_real[sum_alive>thres]
data_real[,ID:=NULL]


#data_real <- data_real[PIq == 5]
data_real[,d:=1]
N <- length(data_real$d)




############################# Generating all the uncertainty ############################
#Generate the health and mortality shocks
if(newmat==1){            # 1 => load in new matrices 
  initdist(0);
  if(hsimtype==1){
    initsim(mortprfs, mort_pi, hsprobs, heal_pi)
  }else if(hsimtype==2){
    initsim2(mortprfs, mort_pi, hsprobs, heal_pi);
  }
}

# parameters for estimation search
#cfloor0 = 8.411833
#nu0_grid = c(3,4)
#phi0_grid = c(0.25,0.35)
#k0_grid = c(5,10)

#matrix_grid = matrix(0,8,3)

 	

#for(ll3 in 1:4){
#  for(ll1 in 1:4){
#    for(ll2 in 1:4){
#      matrix_grid[(ll1-1)*4+ll2 + 16*(ll3-1),3] = k0_grid[ll3]
#      matrix_grid[(ll1-1)*4+ll2 + 16*(ll3-1),1] = nu0_grid[ll1]
#      matrix_grid[(ll1-1)*4+ll2 + 16*(ll3-1),2] = phi0_grid[ll2]
#    }
#  }
#}


# print a summary of the shocks
sumshks();
#grid_solution = matrix(0,64,9)
#random_num = matrix(rnorm(15),5,3)

#print(random_num)


# ----------------- We generate the split to partition training and validation -------------#

fraction = 0.2
unif_vec <- runif(2*N)
split_vec <- 1*(unif_vec>fraction)


#print(LineNumber)

#grid_num = as.integer(LineNumber)


loss_f_mat = matrix(0,64,1)

data_real[,d:=1]




#grid_best = order(loss_f_mat)[1:rep_num]



##########################################################################################



max_kkk = 100

for(repp in 1:rep_num){


if(repp == 1){
cfloor = cfloor0
delta=0;
beta  = 0.97106
nu = 3.8
phi0 = 0.3
K0 = 10
phi00 = phi0
thres = 0
}else{
# best estimation for GAN long lived
cfloor = cfloor0
delta = 0;
beta  = 0.97106
nu = 3.8 + 0.4*rnorm(1)
phi0 = 0.3 + 0.04*rnorm(1)
K0 =  10 + 2*rnorm(1)
phi00 = phi0
thres = 0
}


#-------------------Initialize preference parameters-------------------#
pscaled_p = 1;    # 1 => Parameters transformed to aid optimization 
pscaled_m = 1;

plabel = c("_delta (health shift)    ", "_beta (discount rate)    ", "nu (1/IES)               ", 
           "cfloor (cons. floor)     ", "phi (bequest scale)      ", "K  (bequest curvature)   ");

#delta = 0.009;       # Utility shift due to good health
#delta = 0;       # Utility shift due to good health

# beta0  = 0.97106;  # Discount factor
#beta0 = 0.99;
# beta = beta0;
# nu    = 5.587;  # 1/IES
#nu = 2.9301721 
#nu =7 
#nu    = 3.81
# 
# # Bequest utility parameters                                           
# # Phi_0(b_net) = phi_0*[(b_net+K_0)^(1-nu)]/[1-nu]                    
# 
#phi0   =    0.2999829;  # MPC out of terminal period wealth:  1=> NO bequest motive
#K0     = 1.2084281
#phi0 = 1;
#K0 = 100;
#phi0   =    0.15;  # MPC out of terminal period wealth:  1=> NO bequest motive
#K0 = 123;

#-----------------------Health shock parameters------------------------#
#The process for logged, detrended health costs is the sum of an AR(1)-#
#-and a white noise process.  The processes are normalized to a total--#
#--variance of one.  Note that we have to integrate over each process--#
#----------and that the AR(1) component is a state variable.-----------#

mxlabel = c("Intercept       ", "Age             ", "A^2/100         ","A^3/10000       ", 
            "Bad Health      ", "Bad Health x Age", "Male            ", "Male x Age      ", 
            "PI (percentile) ", "PI x Age        ", "PI^2            ");

mxvec=c(
  2.41663,
  -0.590504,
  -5.6576,
  0.371098,
  -0.423898,
  0.182864,
  -0.480791,
  0.00585711,
  0,
  0,
  0.491667,
  0.0236076,
  -1.42724,
  28.0464,
  -0.624084,
  0.355875,
  0,
  -2.113,
  0.0331839,
  0,
  0,
  -5.87487,
  0.101302,
  -1.62123);

rn_m = length(mxvec);

if(findmedex==0){
  rhomx    = logitrv(0.9215);   # From the `standard' AR(1), transformed for logit 
  fracar1  = logitrv(0.5236/(0.5236+1.039));    
  mxvarscl = (1.039+0.5236)/(1.0476^2);
  
  mxmcoef  = rep(0,11);
  mxvcoef  = rep(0,11);
  
  mxvec    = c(rhomx,fracar1,mxmcoef,mxvcoef);
  zerovec  = rep(0,rn_m);
}else{
  mxvarscl = 1;
  zerovec  = rep(1,rn_m);
  zerovec[9:10]  = rep(0,2);  #male coefficients, mean 
  zerovec[20:21] = rep(0,2);  #male coefficients, variance 
  zerovec[17] = 0;
}

parmvec  = c(delta,beta,nu,cfloor,phi0,K0);
parmvec_grd = parmvec
#rn_p     = length(parmvec);

allparms = c(mxvec,parmvec);
fixvals  = allparms;
zerovec  = c(zerovec,fixparams);

if(swchBeq == 0){
  
  zerovec = c(zerovec,rep(0,2));
}else if(swchBeq == 1){
  
  zerovec = c(zerovec,rep(1,2)); 
}

zerovec0 = zerovec
prnswtch();

data_real[,d:=1]


grad_num = matrix(1,1,6)
v_vec <- matrix(0,1,6)
v_vec_old <-matrix(0,1,6)
#mini_batch_size = floor(N/2) # we assign a size to the mini-batch
#mini_batch_size = floor(N/20)
mini_batch_size = floor(N)
rounds <- floor(N/mini_batch_size)
print(rounds)


#Objective function
loss_f <- matrix(0,500*rounds)
chloss <- matrix(0,500*rounds)
loss_fsmooth <- matrix(0,500*rounds)



parmvec_mat<-matrix(0,500*rounds,6)
grad_mat <- matrix(0,500*rounds,6)

parmvec_mat[1,] <- parmvec

step_g = step0
small_ch = 0
small_ch2 = 0
#num_start <-0
#loc <- 0

kkk <-1 
#num_start2 <-1

#epoch_param = 30
epoch_param = epoch_param_loss
full_derivative = 0

first = 0
count_all <- 0
count_in <- 0
new_param <- 0


# -------------------------------------------------------------------------


# compute derivative
#default: 0.1, zoom = 0.01
hstep = 0.01
# learning rate: default = 1, zoom = 0.05
step0 = 0.1
#step_g = step0
# momentum
gamma_g0 = 0.5
gamma_g = gamma_g0

# type of optimization
# stochastic gradient descent
stoch_grad = 0
# momentum
mom = 1
# accelerated momentum
acc_mom = 1


# optimization as a function of certain parameters
num_param <- 6
# delta, beta, nu, cfloor, phi0, k0
optim_param <- c(0,0,1,0,1,1)
h_mat <- hstep*diag(c(1/10,1/10,1,1,1/10,100))
h_mat2 <- hstep*c(1/10,1/10,1,1,1/10,100)


step_g_mat0 <- step0*c(1/10,1/10,1/10,1/10,1/10,100)
step_g_mat <- step_g_mat0

mat_weight = matrix(0,max_kkk,4)

kkk_old = max_kkk + 1

kkk = 1

while(((max(abs(grad_num))>0.001) + (kkk<10))&(kkk < max_kkk)){
  
  step_g <- step0
  #set.seed(kkk*5)
  # we shuffle the data for each epoch
  sgd_v <- seq(1,N,1)
  
  # we start with the minibatch
  for(ll in (1:rounds)){
    
    if(kkk >2){
      cur = ((kkk-1)*rounds+ll-1)
      cur_1 = ((kkk-1)*rounds+ll-2)
      
      for(grad_ii in 1:num_param){
        if(optim_param[grad_ii]==1){
          if(grad_mat[cur,grad_ii]*grad_mat[cur_1,grad_ii]<0){
            step_g_mat[grad_ii] <- max(step_g_mat[grad_ii]*0.5,0.001*step_g_mat0[grad_ii])
            #step_g_mat[grad_ii] <- max(step_g_mat[grad_ii]*1,0.001*step_g_mat0[grad_ii])
          }else{
            step_g_mat[grad_ii] <- min(step_g_mat[grad_ii]*1.5,100*step_g_mat0[grad_ii])
            #step_g_mat[grad_ii] <- min(step_g_mat[grad_ii]*1,100*step_g_mat0[grad_ii])
          }
        }
      }
    }
    
    
    parmvec_mat[(kkk-1)*rounds+ll,] <- parmvec
    grad_mat[(kkk-1)*rounds+ll,] <- grad_num
    
    # counts going up
    #num_start <- num_start + 1
    # counts overall iterations up to when we need to look "locally".
    #num_start2 <-num_start2 + 1
    
    vec_sgd <- sgd_v[seq((ll-1)*mini_batch_size+1,ll*mini_batch_size,1)]
    
    
    if (mod(kkk,20)==0){
      epoch_param=min(floor(epoch_param*1.05),epoch_param_loss)
    }
    
    ### WE ONLY COMPUTE THE LOSS FUNCTION EVERY 10 ITERATIONS
    #if(mod(count_all,1)==0){
      #count_in = count_in + 1
      oneloop(allparms)
      
      system2(path_exe, invisible = FALSE)
               
      
      syn_data <- loadsim2(mean_asset,stdev_asset)
      NN_syn <- syn_data[,.N]
      #syn_data <- syn_data[PIq ==5]
      syn_data[,d:=0]
      # we select individuals with at least 3 periods
      syn_data[,ID:=seq(1,NN_syn,1)]
      syn_data[, sum_alive:= alive_1996+alive_1998 + alive_2000 + alive_2002 + alive_2004 + alive_2006, by=ID]
      syn_data <- syn_data[sum_alive>thres]
      syn_data[ ,ID:=NULL]
      
      N_syn <- syn_data[,.N]
      
      data_D <- rbind(data_real,syn_data)
      data_D_param <- data_D
      
      
      ##################################3 Train the disciminator ###########################
      
      
      ###############  COMPUTE LOSS FUNCTION #################################
      
      
      # we generate training and validation samples
      out <- training_samples(data_D_param,split_vec)
      
      x_train <- out$x_train
      y_train <- out$y_train
      
      x_test <- out$x_val
      y_test <- out$y_val
      
      # this is the file where we will store the best model
      
      aux_num = (kkk-1)*rounds+ll
      file_name = paste("model_h5_loss/checkpoint",aux_num,".h5",sep="")
      
      # we fix the seed and clear previous models
      #use_session_with_seed(142)
 #     rm(model_keras)
 #     gc()
      
      # we define the model
      model_keras <- define_model()
      
      accur = train_model2(x_train,y_train,x_test, y_test,epoch_param_loss,model_keras,file_name)
      
      # load the best model:
      best_model <- load_model_hdf5(file_name)
      
      
      #aux_weights <- get_weights(best_model)
      
      # we store them
      #mat_weight[kkk,1] = sum(abs(aux_weights[[1]]))
      #mat_weight[kkk,2] = sum(abs(aux_weights[[2]]))
      #mat_weight[kkk,3] = sum(abs(aux_weights[[3]]))
      #mat_weight[kkk,4] = sum(abs(aux_weights[[4]]))
      
      #print("sum weights")
      #print(mat_weight[kkk,])
      
      
      
      predict_aux <- predict_proba(object = best_model, x= as.matrix(rbind(x_train,x_test)))
      table_aux_loss <- as.data.table(cbind(as.matrix(rbind(y_train,y_test)),predict_aux))
      

      colnames(table_aux_loss)[2] <- "V1"
      
      predict_theta <- table_aux_loss[d ==0]$V1
      predict_theta_true <- table_aux_loss[d==1]$V1
      
   
      
      loss_f[kkk]<-Loss(table_aux_loss)
      print(loss_f[1:kkk])
    
    
    ################ Neural Network training for gradient (early stopping) #######################
    ## We compute the gradient in (theta)
    
      
      #use_session_with_seed(142)
      #aux_num = (kkk-1)*rounds+ll
      #file_name_t = paste("model_h5/checkpoint",aux_num,".h5",sep="")
      
      # we define the model
      #rm(model_keras)
      #gc()
      #model_keras <- define_model()
      
      #accur = train_model2(x_train,y_train,x_test, y_test,epoch_param_loss,model_keras,file_name_t)
      
      # load the best model:
      #best_model <- load_model_hdf5(file_name_t)
      
      #predict_aux <- predict_proba(object = best_model, x= as.matrix(rbind(x_train,x_test)))
      #table_aux_loss <- as.data.table(cbind(as.matrix(rbind(y_train,y_test)),predict_aux))
      #predict_theta <- table_aux_loss[d ==0]$V2
    
    
    ###############  COMPUTE NUMERICAL DERIVATIVE #######################
    grad_num_plus <- matrix(0,1,6)
    grad_num_minus <- matrix(0,1,6)
    grad_num_center <- matrix(0,1,6)
    
    #### delta forward
    #if (mom == 1){
    #  parmvec0 <- parmvec_grd
    #}else{
    #}
    
    parmvec0 <- parmvec

    for(grad_ii in (1:num_param)){
      #if(stoch_grad == 1){
      
      if(optim_param[grad_ii]==1){
        parmvec_aux <- parmvec0 + h_mat[grad_ii,]
        
        allparms_aux <-c(mxvec,parmvec_aux);
        fixvals_aux = allparms_aux;
        zerovec  = zerovec0
        
        # save input for C
        oneloop(allparms_aux)
        
        # Generate data
        system2(path_exe, invisible = FALSE)
        
        syn_data_aux <- loadsim2(mean_asset,stdev_asset)
        NN_syn <- syn_data_aux[,.N]
        #syn_data <- syn_data[PIq ==5]
        syn_data_aux[,d:=0]
        # we select individuals with at least 3 periods
        syn_data_aux[,ID:=seq(1,NN_syn,1)]
        syn_data_aux[, sum_alive:= alive_1996 + alive_1998 + alive_2000 + alive_2002 + alive_2004 + alive_2006, by=ID]
        syn_data_aux <- syn_data_aux[sum_alive>thres]
        syn_data_aux[,ID:=NULL]
        
        N_syn <- syn_data_aux[,.N]
        
        syn_data_aux <- syn_data_aux[,.(gender,PI,age_1996,asset_1996,asset_1998,asset_2000,asset_2002,asset_2004,asset_2006,alive_1998,alive_2000,alive_2002,alive_2004,alive_2006,health_1996,health_1998,health_2000,health_2002,health_2004,health_2006)]
        
        
        vec_sgd[vec_sgd>N_syn] = N_syn
        pr_plus <- as.double(predict_proba(object = best_model,x=as.matrix(syn_data_aux[vec_sgd])))
        aux <- data.table(D_p=pr_plus,D_theta=as.double(predict_theta[vec_sgd]))
        aux[,id:=seq(1,mini_batch_size,by=1)]
        aux[,D_plus_star:= max(0.01,min(D_p,0.99)), by=id]
        aux[,D_theta_star:= max(0.01,min(D_theta,0.99)), by=id]
        pr_plus_c <- aux$D_plus_star
        
        #browser()
        grad_num_plus[1,grad_ii] = (mean(log(1-aux$D_plus_star)) - mean(log(1-aux$D_theta_star)))/h_mat2[grad_ii]
        
        
        
      }
      
      print(grad_num_plus)
      # print(grad_num_minus)
      # print(grad_num_center)
      #}
      
    }
    
    
    ####################### update parameters ###############################
    ## We first limit the size of the gradient and then update 
    v_vec_old <- v_vec
    
    
    grad_num = matrix(0,1,6)
    for (ii in 1:6){
      #grad_num[1,ii] = grad_num_plus[1,ii]*(grad_num_plus[1,ii]*grad_num_minus[1,ii]<0)
      grad_num[1,ii] = min(abs(grad_num_plus[1,ii]),.1)*sign(grad_num_plus[1,ii])
    }
    
   
    
    # compute the update
    #if(mom==1){
    #  v_vec = gamma_g*v_vec_old + step_g_mat*grad_num;
    #}else{
      v_vec = 
    #}
    
    parmvec_old <- parmvec
    
    
    # Update parameters
    delta <- parmvec[1] - step_g_mat[1]*grad_num[1,1]
    beta <- parmvec[2] - step_g_mat[2]*grad_num[1,2]
    nu <- parmvec[3] - step_g_mat[3]*grad_num[1,3]
    cfloor <- parmvec[4] - step_g_mat[4]*grad_num[1,4]
    phi0 <- parmvec[5] - step_g_mat[5]*grad_num[1,5]
    K0 <- parmvec[6] - step_g_mat[6]*grad_num[1,6]
    
    parmvec <- c(delta,beta,nu,cfloor,phi0,K0)
    
    #print(c(loss_f[kkk],mean(predict_theta),sd(predict_theta)))
    print(grad_num)
    print(parmvec)
    
    
    allparms = c(mxvec,parmvec);
    fixvals = allparms;
    zerovec  = zerovec0
    
    ## We update the parameter for the next gradient evaluation, if using accelerated momentum
    #if(mom == 1){
    #delta_grd <- parmvec[1] - gamma_g*v_vec[1,1]
    #beta_grd <- parmvec[2] - gamma_g*v_vec[1,2]
    #nu_grd <- parmvec[3] - gamma_g*v_vec[1,3]
    #cfloor_grd <- parmvec[4] - gamma_g*v_vec[1,4]
    #phi0_grd <- parmvec[5] - gamma_g*v_vec[1,5]
    #K0_grd <- parmvec[6] - gamma_g*v_vec[1,6]
    
    #parmvec_grd <- c(delta_grd,beta_grd,nu_grd,cfloor_grd,phi0_grd,K0_grd)
    
    
    kkk <- kkk + 1
    
#    keras::k_clear_session()
    rm(model_keras)
    gc()
    
    #print('v_vec')
    #print(v_vec)
    #print('v_vec_old')
    #print(v_vec_old)
    
    
    #if(kkk>15){
    #  aux = (max(abs(v_vec))<0.0005)==1
    #  if(aux == TRUE){
    #    kkk_old = kkk
    #    kkk = max_kkk
    #  }
    #}
    
    
    
    print(kkk)
    print(repp)
  }
  
  
  
  
}

kkk_old = min(kkk_old,kkk)

uu = which.min(loss_f)
best_mat[repp,1] = loss_f[uu]
best_mat[repp,2] = loss_f[1]
best_mat[repp,3:8] = parmvec_mat[uu,1:6]
best_mat[repp,9] = kkk_old
best_mat[repp,10:15] = parmvec_mat[1,1:6]


#best_weight[rep,] = mat_weight[uu,]


#output = list(best_mat,best_weight)

save(best_mat,file=paste("MC_",num_s,"_",repp,".Rdata",sep=""))

}



#}















