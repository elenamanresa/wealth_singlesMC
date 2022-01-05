/*  WEALTH45.C

  - Finds value fn and policy fns for model in De Nardi, French and Jones (2006) 
  - Written by Fang Yang, with Cristina De Nardi 10/06 
    (cosmetic edits by John Jones 1/07) 
  - Doss: version with singles only (2008)
  - GAUSS-C I/O procedures written by K. Housinger (2003)

  - Endogenizing medical expenditures, De Nardi and Doss 6/1/08

    ~ Note to John: we have changed the GetProfiles code, which interacts with GAUSS.
    ~ We will need both prefShifter_shocksCoeffs and prefShifter_mean to be passed in, 
      holding coefficients for computing mu and the precomputed mean values for mu, 
     the preference shifter.
    ~ We removed references to 'sigma' from the old code. 
    ~ Preference shifter is composed of two parts, deterministic, D, and stochastic, S
    ~ The actual shifter is:   exp(D+S)
    ~ GetProfiles produces the matrix for D, stored in lnMuM_det, provides standard 
      deviations and autocorrelation coefficient for the shocks, and the mus for the 
      stochastic part, S.  Then main discretizes the stochastic variables.  Then GetMu computes
      S and puts D and S together to compute muM.

  - Modifed by JBJ in Sept 08 to simulate discrete-valued transitory shocks.
    Lots of new output matrices created and saved
    New naming convention:
   ~ *cdfsim => simulations of U[0,1] draws
   ~ *indexim => simulations of index numbers
   ~ *sim => actual values

  - I/O instructions
    ~ Put the (GAUSS) *.fmt files in a subdirectory called \iofiles\
    ~ create a folder called \data
    ~ create a folder called \output

  - There are many subroutines.
    ~ GetRulesSingle does the maximization for the single cases. 
    ~ GetContinuationUtility calculates value for each consumption choice and is used by GetRulesSingle. 
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <limits.h>

/*  MPI-specific Code:  Uncomment if used in UNIX */
//#include <mpi.h>  
//#define max(x,y) ((x) > (y) ? (x) : (y)) // i am scared of pre-compilatoin macros (ex: max(x++, y) ).

// Disable warning messages 4996 
#pragma warning( once : 4996 )

#define switchJohn 1     /* 0: Cristina in Chicago, 1: John in Albany; 2 Charles Doss in Chicago */ 
#define useMPI 0         /* 0 => PC, 1 => Unix, 2 => cluster  */
#define PI 3.14159265
#define rounder 0.0001
#define ADDRESS_LEN 300

/* These constants are used to size various matrices */
#define TDIMS 33       /* dim. of age for singles  Programmer must ensure compatibility */
#define IDIM 5         /* dim. of permanent income, from lowest to highest */
#define HSDIM 2        /* dim. of health status: bad=0, good=1 */
#define ZETADIM 8      /* dim. of persistent component of medical shock */
#define XIDIM 8        /* dim. of temporary component of medical shock */

//Note: Assets also define consumption! So Amin fixes the maximum that can be consumed/expended for another fixed asset level.
#define AMIN 0
#define AMAX  2000000        
#define AMAX1 200000         /* maximum assets for segment 1, old grid algorithm */ 
#define AMAX2 1000000        /* maximum assets for segment 2 */
#define AMAX3 2000000        /* maximum assets for segment 3 */

// need Extra parameter for new grid; because the consumption grid DNE anymore, but is defined by the asset grid's differences from itself
#define ADIM  200          /* dim. of assets */
#define LAMBDA 2
#define ADIM1 50 // 250   /* dim. of assets for segment 1, old grid algorithm  */
#define ADIM2 10 // 100   /* dim. of assets for segment 2; in old grid, ADIM > ADIM1+ADIM2*/

#define taxDim 7       /* dim. of vector of marginal tax rates */

#define BASIC_HEADER_LEN 128 /* Length of a scalar header */
#define BYTE_POS 5           /* Int offset in header of byte order */
#define BIT_POS 6            /* Int offset in header of bit order */
#define TYPE_POS 15          /* Int offset in header of matrix type */
#define HEADER_LEN_POS 18    /* Int offset in header of header length */
#define M_POS 32             /* Int offset in header of value for m */
#define N_POS 33             /* Int offset in header of value for n */
#define SCALAR 0             /* Value of type that indicates a scalar */
#define MATRIX 2             /* Value of type that indicates a matrix */
#define BIT_SYSTEM 0         /* Value of bit order (0=backwards/i386) */
#define BYTE_SYSTEM 0        /* Value of byte order (0=backwards/i386) */

// To make this not a global would either have to pass it to all routines or have a 
// 'universe' class of which this is a member
FILE *errorOutput; //initialized at the beginning of main and closed' at the end of main.  Use for all output.
//Log file; for knowing what you've done.
FILE *logOutput;
int printOn; //'bool' for whether to print or not.  (will probably depend on job (whether estimating or not), MPI, etc.
int errorPrintOn; //want separate tests for printing errors and printing other things.
/*
**This is my GMatrix opaque datatype.  The gauss fmt format has implied 
**dimensionality, but is just a list of doubles.  I could rewrite the
**functions to use pointers-to-pointers to improve the interface with
**Eric's functions, but that would add to the complexity of my code.  I
**suggest writing an intermediate function that turns *data into **data
**to cope with 2-dimensional arrays.
*/ 
//
typedef struct {   /* this defines a matrix with given # of rows and */
  unsigned int m;  /* columns and gives the address of the first element */
  unsigned int n;
  double *data;} GMatrix;

///*--------------------------------------------------------------------------------*/
///*----------------------Global PARAMETERS read in from GAUSS----------------------*/
//
///* Preference Parameters */
   double delta, beta, nu, omega, IESfrac, phi0, K0; /* **** dropped eta, phi1 and phi2, k1 and 2 */

/* Simulation Parameters */
   int nsims, TDimSims, TSTART;  

///* Switches */
   int switchMor;     /* 0: no shocks, 1: shocks of mortality risk */
   int switchBeta;    /* 0: from GAUSS, 1: beta=1 */
   int switchY;       /* 0: noIncome 1: use income loaded */
   int switchPrefShifter;  /* 0: no deterministic health cost, 1:  Health cost */
   int switchTax;     /* 0: noIncome, 1: use tax  loaded  */
   int switchZeta;    /* 0: no shocks, 1: shocks of health costs zeta,*/
   int switchXi;      /* 0: no shocks, 1: shocks of health costs zi */
   int switchBeq;     /* 0: no, 1: bequest motive  */
   double switchCMin; /* 0: no, 1: force everyone to consume xMin */ 
   int switchFloor;   /* 0: expenditure floor; 1: utility floor */   
  
   int job; //read from GAUSS ( 1 => estimating and getting se's;2 => getting se's only; 3 => get graphs;  4 => experiments. )
   int rank;   /* Index number of this node  */
   int size;   /*  Number of nodes in cluster */

///* Miscellaneous */
   double eMin, uMin, uMinAbs, tauBeq, exBeq, mu_r, medex_bottomcode;  
                                                       //Also note: 'sigma_r' for interest rate eliminated
///* Non-asset income:  y_(t)  =  y(f,IDIM,TDIMS+1), a deterministic fn.
//   INDEXES:  
//      f: family structure
//      f=0 single male
//      f=1 single female
//      IDIM: permanent income quantile, ordered from the lowest to highest 
//      TDIMS: age
//*/
   double yM[2][IDIM][TDIMS];

///* Health status uncertainty:
//   Health status at time t, m_(t)(i), is a Markov process taking 
//   two values, good and bad. The transition probabilities for health status 
//   depend on current health status and age. The elements of the health status 
//   transition matrix are pi_(kjt)(i) = Pr(m_(t+1)(i) =  j|m_(t)(i)  =  k),   
//   k, j~(good, bad).
//   Singles first:
//   INDEXES: 
//      SEX: 0=husband, 1=wife
//      HSDIM: health state today 0=bad, 1=good
//      HSDIM: health state tomorrow 0=bad, 1=good 
//*/
   double hsProbM[2][TDIMS][IDIM][HSDIM][HSDIM];  /* Singles */


///* Survival uncertainty:
//   s_(m, I, t)(i) denotes the probability of individual i being alive at age t
//   conditional on current health status, permanent income, and being alive at age t-1. 
//   0=male, 1=female 
//*/
   double survivalProbM[2][TDIMS][HSDIM][IDIM]; 
//
/* Health Costs
   ln hc_(t)  =  hc(f_(t), m_(I,t)(h), m_(I,t)(w), t, I)+
                 sigma(f_(t), m_(I,t)(h), m_(I,t)(w), t, I) * psi_(t).
   psi_(t)    =  zeta_(t)+xi_(t),   xi_(t) ~ N(0, sigma_(xi)^2),  
   zeta_(t)   =  rho_(hc)zeta_(t-1) + epsilon_(t),   epsilon_(t) ~ N(0, sigma_(epsilon)^2) 
*/
   double rhoHc, sigma_xi, sigma_epsilon; 

/* health cost functions */
//   double hcSingleLogMean[2][HSDIM][TDIMS+1][IDIM]; 
//   double hcSingleSigma[2][HSDIM][TDIMS+1][IDIM]; 

/* Income tax structure, from French's code
   - taxBrk gives tax brackets
   - taxMar gives marginal tax rates  
*/
   double taxBrk[taxDim-1] = {6250, 40200, 68400, 93950, 148250, 284700}; 
   double taxMar[taxDim] = {0.0765, 0.2616, 0.4119, 0.3499, 0.3834, 0.4360, 0.4761}; 

/* Strings of directory names */
   char rootdir[ADDRESS_LEN]; 
   char outputdir[ADDRESS_LEN]; 
   char datadir[ADDRESS_LEN]; 
   
/*----------------Global VARIABLES determined inside this program-----------------*/

double aA[ADIM];       /* asset grid */
double hsA[HSDIM];     /* health status grid */
double IncomeA[IDIM];  /* Permanent income grid, income as percentage ranking. */
int tA[TDIMS];         /* age grid */
int gA[2];             /* gender grid, 0==male, 1==female */
                       /* (See subroutine Grid for detail.) */

//* Used to store markov chains */
double xiA[XIDIM], xiProbA[XIDIM];  /* temporary medical shock  */
double xiProbAcdf[XIDIM+1];         
double tempInv[XIDIM];              /* Placeholder              */
double zetaA[ZETADIM];              /* persistent medical shock */
double zetaProbM[ZETADIM][ZETADIM];
double zetaInv[ZETADIM];            /* stationary probabilities */
double zetaInvcdf[ZETADIM+1];
int zetaIndN; 
double zetaProbMcdf[ZETADIM][ZETADIM+1]; 

//* income after tax at brackets  */
double incomeBrk[taxDim-1];   

//* expected utility from bequest */
double bequestUM[ADIM];  /* single */
//precompute medical preference shock and store it in muM
double muM[2][TDIMS][IDIM][HSDIM][ZETADIM][XIDIM];
double qM[2][TDIMS][HSDIM];

//Note: the "TDIMS" index refers to _age_ not _time_ (or year). (TDIMS is the number of allowed ages.)
//* value function matrix */
double valueFunM[2][TDIMS+1][IDIM][HSDIM][ZETADIM][ADIM][XIDIM];  
//* consumption policy function matrix */
double consumFunM[2][TDIMS][IDIM][HSDIM][ZETADIM][ADIM][XIDIM]; 
//* total medical expenditures policy function matrix */
double medFunM[2][TDIMS][IDIM][HSDIM][ZETADIM][ADIM][XIDIM]; 
double outofpocketmedFunM[2][TDIMS][IDIM][HSDIM][ZETADIM][ADIM][XIDIM]; 
//* bequest policy function matrix*/
double savingsFunM[2][TDIMS][IDIM][HSDIM][ZETADIM][ADIM][XIDIM];  /* singles */

struct result
{ int Ind1;
double weight;
};
//
///*--------------------------------------------------------------------------------*/
///*--------------------------------------------------------------------------------*/
///* functions used in main*/
//
///*  Converts vectors to matrices and vice-versa*/
GMatrix zeromat(unsigned int recsize);
void switchem(double **dataPtr, double **tempPtr);
void copyit(double *permvec, double *tempvec, int recsize);
double **SetUpSim(double *dataVec, int extrayears);

///*  GAUSS-C++ I/O programs written by K. Housinger */
unsigned char * gread(unsigned char *inbuf, int bytes, int byte_reverse, int bit_reverse);
GMatrix gau5read(char *fmt);
void gau5write(char *fmt, GMatrix mat); /* reads vector from hard drive*/

//*   Update global parameters with values passed in from GAUSS*/
int globderef(double *prefvecPtr, double *asstvecPtr, double *medexvecPtr, 
              double *simvecPtr, double *agevecPtr, double *switchvecPtr);

//Gets data from GAUSS.
//Fills in lnMuM_det, ie "D," see the comments at the top of this file about how the shifter is computed.
int GetProfiles(double *agevecPtr, double *prefShifter_mean,
                double *prefShifter_shockCoeffs,
                double *muPIPtr, double *qMPtr, double *mortrateSPtr,  
                double *mortratePIPtr, double *hstranSPtr, 
                double *hstranPIPtr, double *yprofPtr, double *yprofPIPtr,
                double lnMuM_det[2][TDIMS][IDIM][HSDIM],
                double lnMu_shockCoeffs[2][TDIMS][IDIM][HSDIM]);

double LogitSQRT(double x);
double Logit(double x);
void TwoYearToOneYear(double gg, double bb, double *_g, double *_b); 

// create grid for cash on hand and consumption; new version of grid allows for a general parameter lambda
void Grid(double lambda, double aA[ADIM], double hsA[HSDIM], double IncomeA[IDIM]);
void Grid_Old(double aA[], double hsA[HSDIM], double IncomeA[IDIM]);

void Discretization(int n,double rho, double mu, double sigma, double zArray[], 
                    double *piMatrixP, double *piInvarV); 

// Define Utility from consumption matrix 
double getUtility(double cons, double totalMedex, double muRealization); 

// Find expenditures that satisfy utility floor
double getXMin(double qRealization, double muRealization);  

// q, the out-of-pocket share of medical costs
double qShare(int gInd, int tInd, int iInd, int hsInd);

// Not discounted
double GetContinuationUtility(int aTomorrowInd, int tInd, int gInd, int hsInd, 
                              int iInd, int zetaInd);

// Define Utility from bequest matrix 
void GetUtilityBeq(double bequestUM[ADIM]);

//precompute mu (medical preference shock)
//  GetMu leaves muM in its final form.
//  lnMuM_det and lnMuM_shockCoeffs are filled by getProfiles.  The former contains the entire deterministic
//  portion of lnMu.  The latter contains coefficients to multiply zeta+xi by, including interaction terms
//  (i.e. terms like coeff*t*g).
void GetMu(double lnMuM_det[2][TDIMS][IDIM][HSDIM], double lnMuM_shockCoeffs[2][TDIMS][IDIM][HSDIM], 
           double xiA[XIDIM], double zetaA[ZETADIM], double muM[2][TDIMS][IDIM][HSDIM][ZETADIM][XIDIM]);

//* after-tax income at bracket points  */
void IncomeAtBrk(double taxBrk[], double taxMar[], double incomeBrk[]);

//* Interpolation and extrapolation */
//Interpolation(): fp is the array on which we interpolate; xP is the grid giving the values associated with
// each element in fP; that is, the appropriate weight is calculated using xP and x (and DIM).

double Interpolation(double *fP,  double *xP,  double x, int DIM);
int Locate(double *Xarray, double x, int DIM);
struct result GetLocation(double *xP, double x, int DIM);

void WriteNonSimData();
void WriteSimData(double **assetsimMat, double **netIncomesimMat, double **consumptionsimMat, 
                  double **totalmedexsimMat, double **oopmedexsimMat,  double **zetaindexsimMat, 
                  double **marstatsimMat, GMatrix agesim96Ptr, GMatrix PIsim96Ptr, 
                  double **healthsimhMat, double **healthsimwMat);

double findConsumption(double allofrac, double totalfunds);

void GetRulesSingle(int iAssetsmin,int iAssetsmax);

//Global tdimSims is used for number of periods the sim will run.  The state and decision grids
//are allocated to be the same dimension.  The last decision variable is extra;
//(i.e. start with a[0]; choose c[0]; get a[1]).  We make tdimsims+1 choices but the last is irrelevant.
// Thus we store tdimsims+1 state variables (and implicitly tdimsims+2 exist).
// The timing and use of Epsilon is analagous to a decision variable.

void simulation(GMatrix zetacdfsim96Ptr,  GMatrix PIsim96Ptr,  GMatrix agesim96Ptr, GMatrix assetsim96Ptr,
                double **assetsimMat, //storing levels not indices; indices are done only internally
                double **cohsimMat, double **netIncomesimMat, double **consumptionsimMat, //empty on entry; for output to GAUSS
                double **healthindexsimhMat, double **healthindexsimwMat, double **marstatsimMat, //store indices
                double **xicdfsimMat, double **epsiloncdfsimMat,  //epsilon is analagous to decision variable; has one less draw than xi; we would like it to be numbered from 1 to tsimdims; C does start-at-0 indexing, so we call it "T+1" since the index is one behind where we want it.
                double **totalmedexsimMat, double **oopmedexsimMat, double **MedicaidsimMat, 
                double **grossbeqsimMat, double **zetaindexsimMat, double **zetasimMat, 
                double **xiindexsimMat, double **xisimMat, double **musimMat, double **qsimMat,
                double **marstatsim2Mat, double **transfersimMat, int iSimmin, int iSimmax);
            
void GetCdf(int nRows, int nCols, double *piMatrixP, double *piMatrixCDF); 

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* main program, note that everything before this is global */
int main(int argc, char *argv[])
{
  double *tmpTempV, *tmpPtr, *tempvec;

   GMatrix prefvecPtr, asstvecPtr, agevecPtr, simvecPtr, switchvecPtr, medexvecPtr, 
           prefShifter_meanPtr, prefShifter_shockCoeffsPtr, qMPtr, muPIPtr, mortrateSPtr, 
           mortratePIPtr, hstranSPtr, hstranPIPtr, yprofPtr, yprofPIPtr, //rorsimPtr, 
           assetsim96Ptr, zetacdfsim96Ptr, PIsim96Ptr, agesim96Ptr, healthindexsimhPtr, 
           healthindexsimwPtr, marstatsimPtr, xicdfsimPtr, epsiloncdfsimPtr, 
           assetsimPtr, zetasimPtr, zetaindexsimPtr, xisimPtr, xiindexsimPtr, musimPtr, 
           totalmedexsimPtr, oopmedexsimPtr, MedicaidsimPtr, consumptionsimPtr, beqsimPtr, 
           netIncomesimPtr, qsimPtr, cohsimPtr, GaussJobPtr, marstatsim2Ptr, transfersimPtr;

   double **assetsimMat, **netIncomesimMat, **consumptionsimMat, **healthindexsimhMat, 
          **healthindexsimwMat, **marstatsimMat, //0==dead, 1==male, 2==female, 3==couple
          **xicdfsimMat, **epsiloncdfsimMat, **zetasimMat, **zetaindexsimMat, **xisimMat, 
          **xiindexsimMat, **qsimMat, **totalmedexsimMat, **oopmedexsimMat, **MedicaidsimMat,
          **grossbeqsimMat, **cohsimMat, **musimMat, **marstatsim2Mat, **transfersimMat;

   double lnMuM_det[2][TDIMS][IDIM][HSDIM], lnMu_shockCoeffs[2][TDIMS][IDIM][HSDIM];

   int gotderef, gotprofiles, recsize, mpierr;  //yearInd, 
   int zetaInd, xiInd;  /*loop indices */

   char fullpath[ADDRESS_LEN];

   //double *rorsim; /* define ror sequence in main. remember to change it.  */

   clock_t start, end;  /* recode time */

   int simspermachine, iAssetsmin, iAssetsmax, iSimmin, iSimmax, statespermachine;
   double spm2;

   FILE *parameterP, /*point to file containing parameters */
      *timeStamp, //an ID for each run and has length of run.
      *fileP, *fileP2; //some file output.

   if (useMPI==0)  /*  PC */
   {
      if (switchJohn==1)  /* John is using in Albany */
      {
         strcpy(rootdir,"c:\\wealth_singles\\iofiles\\");
         strcpy(outputdir, "c:\\wealth_singles\\output\\"); 
         strcpy(datadir, "c:\\wealth_singles\\data\\"); 
      }
      else if (switchJohn==0)
      {
         strcpy(rootdir,"c:\\cristina\\eric\\Charles_July08\\iofiles\\");
         strcpy(outputdir, "c:\\cristina\\eric\\Charles_July08\\output\\");
         strcpy(datadir, "c:\\cristina\\eric\\Charles_July08\\data\\");
      }
     else if (switchJohn==2) //charles
     {
         strcpy(rootdir,"c:\\Documents and Settings\\g1crd01\\My Documents\\wealth\\fullcode\\iofiles\\");
         strcpy(outputdir, "c:\\Documents and Settings\\g1crd01\\My Documents\\wealth\\fullcode\\output\\");
         strcpy(datadir, "c:\\Documents and Settings\\g1crd01\\My Documents\\wealth\\fullcode\\data\\");
     }
   }
   else if (useMPI==1){ //unix
    if (switchJohn == 2) //charles
      {
         strcpy(rootdir,"/home/cdoss/Charles_July08/iofiles_nompi/"); 
         strcpy(outputdir,"/home/cdoss/Charles_July08/output_nompi/"); 
         strcpy(datadir,"/home/cdoss/Charles_July08/data/"); 
      }
   }
   else if (useMPI==2) /*   Cluster*/
   {
      if (switchJohn==1)
      {
         strcpy(rootdir,"/home/jjones/wealth_singles/job_021009/iofiles/"); 
         strcpy(outputdir,"/home/jjones/wealth_singles/job_021009/output/"); 
         strcpy(datadir,"/home/jjones/wealth_singles/data/"); 
      }
      else if (switchJohn == 2)
      {
         strcpy(rootdir,"/home/cdoss/Charles_July08/iofiles/"); 
         strcpy(outputdir,"/home/cdoss/Charles_July08/output/"); 
         strcpy(datadir,"/home/cdoss/Charles_July08/data/"); 
      }
   }
   GaussJobPtr = gau5read(strcat(strcpy(fullpath,rootdir),"job.fmt"));
   job   = (int) floor(rounder+GaussJobPtr.data[0]); //job is global
   
   start = clock();
   rank  = 1; /* Global */
   size  = 1;

// MPI-specific Code:
/*
   if (useMPI==2)
   {
      MPI_Init(&argc, &argv);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      printf("size=%3d rank=%3d\n", size, rank);
      MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
   }
*/
   if( (useMPI<2) || (rank==(size-1)) ) { 
      printOn=1;
      logOutput         =   fopen(strcat(strcpy(fullpath,outputdir), "c_logoutput.txt"), "w");
      timeStamp = fopen(strcat(strcpy(fullpath,outputdir), "timeStamp.txt"), "w");
   }
   else {//else don't need to open the output files because they should never be written to.
      printOn=0;
   }
   if ((useMPI<2)||(rank==(size-1))){
      errorPrintOn=1;
      errorOutput=fopen(strcat(strcpy(fullpath,outputdir), "erroroutput.txt"), "w");
   }
   else{
      errorPrintOn=0;
   }

   if (printOn==1)
   {
   fprintf(timeStamp, "ID: %f\n", start); //ID==start is the time at which the program was run.
   }

// Read in parameter vectors from the *.fmt files (read into memory as arrays)
   prefvecPtr    = gau5read(strcat(strcpy(fullpath,rootdir),"prefvec.fmt"));  
   asstvecPtr    = gau5read(strcat(strcpy(fullpath,rootdir),"asstvec.fmt"));
   agevecPtr     = gau5read(strcat(strcpy(fullpath,rootdir),"agevec.fmt")); 
   switchvecPtr  = gau5read(strcat(strcpy(fullpath,rootdir),"swchvec.fmt")); 
   simvecPtr     = gau5read(strcat(strcpy(fullpath,rootdir),"simvec.fmt")); 
   medexvecPtr   = gau5read(strcat(strcpy(fullpath,rootdir),"medexvec.fmt"));   
   prefShifter_meanPtr        = gau5read(strcat(strcpy(fullpath,rootdir),"mu_mean.fmt"));
   prefShifter_shockCoeffsPtr = gau5read(strcat(strcpy(fullpath,rootdir),"mu_stdev.fmt")); 
   muPIPtr       = gau5read(strcat(strcpy(fullpath,rootdir),"mupicoef.fmt"));
   mortrateSPtr  = gau5read(strcat(strcpy(fullpath,rootdir),"mortprfs.fmt"));   
   mortratePIPtr = gau5read(strcat(strcpy(fullpath,rootdir),"mort_pi.fmt"));
   hstranSPtr    = gau5read(strcat(strcpy(fullpath,rootdir),"hsprobs.fmt")); 
   hstranPIPtr   = gau5read(strcat(strcpy(fullpath,rootdir),"heal_pi.fmt"));
   yprofPtr      = gau5read(strcat(strcpy(fullpath,rootdir),"yprof.fmt")); 
   yprofPIPtr    = gau5read(strcat(strcpy(fullpath,rootdir),"y_pi.fmt"));
   qMPtr         = gau5read(strcat(strcpy(fullpath,rootdir),"qcoefs.fmt"));

// Initialize parameters by assigning array elements to appropriate parameters
   gotderef     = globderef(prefvecPtr.data, asstvecPtr.data, medexvecPtr.data,
                            simvecPtr.data, agevecPtr.data, switchvecPtr.data);

// Create grids for cash on hand, consumption, health status and income 
   Grid(LAMBDA, aA, hsA, IncomeA);

// Initialize profiles 
   gotprofiles = GetProfiles(agevecPtr.data, prefShifter_meanPtr.data,
                             prefShifter_shockCoeffsPtr.data, muPIPtr.data, 
                             qMPtr.data, mortrateSPtr.data, mortratePIPtr.data, 
                             hstranSPtr.data, hstranPIPtr.data, yprofPtr.data, 
                             yprofPIPtr.data, lnMuM_det, lnMu_shockCoeffs);

//* Save dimensions to ASCII files so that matlab can transform output files into matrices */
   if (printOn==1)
   {
      parameterP=fopen(strcat(strcpy(fullpath,outputdir),"parameter.txt"),"w" );
      fprintf(parameterP, "%5d\n %5d\n %5d\n %5d\n %5d\n %5d\n  %10.8f\n %10.8f\n %10.8f\n %10.8f\n %10.8f\n %10.8f\n %10.8f\n %10.8f\n",
              IDIM,  TDIMS, HSDIM, ZETADIM, XIDIM, ADIM, eMin, delta, beta, nu, omega, phi0, K0, mu_r);

     //CCC == modifying the COH args to be asset instead, to fprintf above 
     /* **** note I dropped some arguments above */
    /* Items are
      - dimension of permanent income
      - dimension of age  
      - dimension of health status: good = 1,  bad = 0
      - dimension of persistent component of medical shock 
      - dimension of temporary component of medical shock 
      - dimension of assets
      - dimension of bequest to the child  
      - eMin, expenditure floor
      - delta, health dependence utility param 
      - beta, discount factor 
      - nu, coefficient of relative risk aversion
      - phi0, beq parameter 
      - K0, beq parameter       
   */
      fclose(parameterP);   
   }

   if (switchZeta==1)
   {
   /* discretize persistent AR(1) process zeta into Markov chain */
      Discretization(ZETADIM, rhoHc, 0.0, sigma_epsilon, zetaA, &zetaProbM[0][0], &zetaInv[0]);  
   }
   else
   {
      for (zetaInd=0; zetaInd<ZETADIM; zetaInd++)
      { 
         zetaA[zetaInd]=0.0;
         zetaInv[zetaInd]=1.0/ZETADIM; 
         for (zetaIndN=0; zetaIndN<ZETADIM; zetaIndN++)
         {
            zetaProbM[zetaInd][zetaIndN]=1.0/ZETADIM;
         }
      }
   }

    GetCdf(ZETADIM, ZETADIM, &zetaProbM[0][0], &zetaProbMcdf[0][0]);  /* get conditional CDF of zeta */
    GetCdf(1, ZETADIM, &zetaInv[0], &zetaInvcdf[0]);  /* get invariant CDF of zeta */

/* Transitory component of medex preference shock */
   if (switchXi==1)
   {
      Discretization(XIDIM,0.0,0.0,sigma_xi, xiA, &xiProbA[0], &tempInv[0]); 
   }
   else /* no shocks */
   {   
      for (xiInd=0; xiInd<XIDIM; xiInd++)
      {
         xiA[xiInd]=0.0;
         xiProbA[xiInd]=1.0/XIDIM;
      }
   }

   GetCdf(1, XIDIM, &xiProbA[0], &xiProbAcdf[0]);  /* get CDF of xi */


   if (printOn==1)
   {
      fileP=fopen(strcat(strcpy(fullpath,outputdir),"xiA.txt"),"w");
      for (xiInd=0; xiInd<XIDIM; xiInd++)
      { 
         fprintf(fileP,"%15.5f\n", xiA[xiInd]);
      }
      fclose(fileP);
     
      fileP=fopen(strcat(strcpy(fullpath,outputdir),"xiProbA.txt"),"w");
      for (xiInd=0; xiInd<XIDIM; xiInd++)
      { 
         fprintf(fileP,"%15.5f\n", xiProbA[xiInd]);
      }
      fclose(fileP);

      fileP2=fopen(strcat(strcpy(fullpath,outputdir),"zetaA.txt"),"w");
      fileP=fopen(strcat(strcpy(fullpath,outputdir),"zetaProbM.txt"),"w");

      for (zetaInd=0; zetaInd<ZETADIM; zetaInd++)
      { 
        fprintf(fileP2,"%15.5f\n", zetaA[zetaInd]);
        for (zetaIndN=0; zetaIndN<ZETADIM; zetaIndN++){
         fprintf(fileP,"%15.5f ", zetaProbM[zetaInd][zetaIndN]);
        }
        fprintf(fileP,"\n");
      }
      fclose(fileP);
     fclose(fileP2); 
   }

// Initialize matrices
   GetUtilityBeq(bequestUM); /* Find utility from bequest matrix */
   GetMu(lnMuM_det, lnMu_shockCoeffs,  xiA, zetaA, muM); //precompute mu (medical preference shock); requires discretization
   IncomeAtBrk(taxBrk, taxMar, incomeBrk); /* after-tax income at bracket points  */
   
   if (useMPI==2)
   {
      spm2 = ((double) ADIM)/ ((double) size);
      statespermachine = (int) ceil(spm2);
      iAssetsmin = rank*statespermachine;
      iAssetsmax = (rank+1)*statespermachine;
      if (iAssetsmax > ADIM) iAssetsmax = ADIM;

      spm2 = ((double) nsims)/ ((double) size);
      statespermachine = (int) ceil(spm2);
      iSimmin = rank*statespermachine;
      iSimmax = (rank+1)*statespermachine;
      if (iSimmax > nsims) iSimmax = nsims;
   }
   else
   {
      iAssetsmin=0;
      iAssetsmax=ADIM;

      iSimmin=0;
      iSimmax=nsims;
   }

   GetRulesSingle(iAssetsmin,iAssetsmax);
   printf("Finished with decision rules for singles (rank==%d).\n", rank);

// Load inputs simulated in GAUSS

// Initial year
   PIsim96Ptr      = gau5read(strcat(strcpy(fullpath,rootdir),"pisim96.fmt"));  /* simulated permanent income in 1996 */
   agesim96Ptr     = gau5read(strcat(strcpy(fullpath,rootdir),"agesim96.fmt")); /* simulated age in 1996 */
   assetsim96Ptr   = gau5read(strcat(strcpy(fullpath,rootdir),"asim96.fmt"));   /* simulated assets in 1996 */
   zetacdfsim96Ptr = gau5read(strcat(strcpy(fullpath,rootdir),"ztacdfsim96.fmt")); /* simulated U[0,1] shocks for zeta in 1996. */

// Read in whole history */  
   healthindexsimhPtr = gau5read(strcat(strcpy(fullpath,rootdir),"healsimh.fmt")); /* simulation of health status for males */
   healthindexsimwPtr = gau5read(strcat(strcpy(fullpath,rootdir),"healsimw.fmt")); /* simulation of health status for females */
   marstatsimPtr      = gau5read(strcat(strcpy(fullpath,rootdir),"mstatsim.fmt")); /* simulation of marital status */
   xicdfsimPtr        = gau5read(strcat(strcpy(fullpath,rootdir),"xicdfsim.fmt")); /* simulated U[0,1] shocks for transitory variable xi. */
   epsiloncdfsimPtr   = gau5read(strcat(strcpy(fullpath,rootdir),"epscdfsim.fmt"));/* simulated U[0,1] shocks for innovation epsilon. */

// Reshape array to matrix
// If there is a '0' passed, matrix has TDSIMS rows; '1' => TDSIMS+1
   healthindexsimhMat = SetUpSim(healthindexsimhPtr.data,1); /* fill in entire history */
   healthindexsimwMat = SetUpSim(healthindexsimwPtr.data,1); /* fill in entire history */
   marstatsimMat      = SetUpSim(marstatsimPtr.data,1);      /* fill in entire history */
   xicdfsimMat        = SetUpSim(xicdfsimPtr.data,1);        /* Fill in entire history. */
   epsiloncdfsimMat   = SetUpSim(epsiloncdfsimPtr.data,1);   /* Fill in entire history. */

// Set up matrices for storing simulation results
// Need zero matrix (rather than undefined values) so that MPI_SUM can be used.
   recsize           = (TDimSims+1)*nsims;
   zetaindexsimPtr   = zeromat(recsize);
   zetaindexsimPtr   = zeromat(recsize);
   assetsimPtr       = zeromat(recsize);
   marstatsim2Ptr    = zeromat(recsize);  
   consumptionsimPtr = zeromat(recsize);
   beqsimPtr         = zeromat(recsize);
   netIncomesimPtr   = zeromat(recsize);
   totalmedexsimPtr  = zeromat(recsize);
   oopmedexsimPtr    = zeromat(recsize);
   cohsimPtr         = zeromat(recsize);
   musimPtr          = zeromat(recsize);
   qsimPtr           = zeromat(recsize);
   xiindexsimPtr     = zeromat(recsize); 
   zetasimPtr        = zeromat(recsize); 
   xisimPtr          = zeromat(recsize); 
   MedicaidsimPtr    = zeromat(recsize); 
   transfersimPtr    = zeromat(recsize); 

   tempvec = (double *)calloc(recsize,sizeof(double));

   zetaindexsimMat   = SetUpSim(zetaindexsimPtr.data,1);  // zeta, marstat and assets get computed for T+1;
   zetasimMat        = SetUpSim(zetasimPtr.data,1);       
   assetsimMat       = SetUpSim(assetsimPtr.data,1);
   marstatsim2Mat    = SetUpSim(marstatsim2Ptr.data,1);    
   grossbeqsimMat    = SetUpSim(beqsimPtr.data,1);        // everything else gets computed up to T
   consumptionsimMat = SetUpSim(consumptionsimPtr.data,1);// empty column is a placeholders to enable switchem
   netIncomesimMat   = SetUpSim(netIncomesimPtr.data,1);
   totalmedexsimMat  = SetUpSim(totalmedexsimPtr.data,1);
   oopmedexsimMat    = SetUpSim(oopmedexsimPtr.data,1);
   MedicaidsimMat    = SetUpSim(MedicaidsimPtr.data,1);  
   transfersimMat    = SetUpSim(transfersimPtr.data,1);  
   xiindexsimMat     = SetUpSim(xiindexsimPtr.data,1);
   xisimMat          = SetUpSim(xisimPtr.data,1);        
   qsimMat           = SetUpSim(qsimPtr.data,1);
   musimMat          = SetUpSim(musimPtr.data,1);
   cohsimMat         = SetUpSim(cohsimPtr.data,1);

   printf("simulation starts for rank=%d\n", rank);

   simulation(zetacdfsim96Ptr, PIsim96Ptr, agesim96Ptr, assetsim96Ptr, assetsimMat,
              cohsimMat, netIncomesimMat, consumptionsimMat, healthindexsimhMat, 
              healthindexsimwMat, marstatsimMat, xicdfsimMat, epsiloncdfsimMat, 
              totalmedexsimMat, oopmedexsimMat, MedicaidsimMat, grossbeqsimMat, 
              zetaindexsimMat, zetasimMat, xiindexsimMat, xisimMat, musimMat, 
              qsimMat, marstatsim2Mat, transfersimMat, iSimmin, iSimmax);  

   printf("Simulation ends for rank=%d\n", rank);

// MPI-specific Code:
/*
   if (useMPI==2)
   {
    MPI_Reduce((void *) cohsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //CARE: if you "switchem()" which switches the ADDRESSES then your *Mat
    //will not point to the correct location!!!
    //     switchem(&cohsimPtr.data, &tempvec); //tempvec won't be 0s, but it will be ignored.
    copyit(cohsimPtr.data, tempvec, recsize);

    mpierr = MPI_Reduce((void *) assetsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                        MPI_SUM, size - 1, MPI_COMM_WORLD);
    //    switchem(&assetsimPtr.data, &tempvec);
    copyit(assetsimPtr.data, tempvec, recsize);


    MPI_Reduce((void *) netIncomesimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&netIncomesimPtr.data, &tempvec);
    copyit(netIncomesimPtr.data, tempvec, recsize);

    MPI_Reduce((void *) zetaindexsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&zetaindexsimPtr.data, &tempvec);
    copyit(zetaindexsimPtr.data, tempvec, recsize);

    MPI_Reduce((void *) consumptionsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&consumptionsimPtr.data, &tempvec);
    copyit(consumptionsimPtr.data, tempvec, recsize);

    MPI_Reduce((void *) beqsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&beqsimPtr.data, &tempvec);
    copyit(beqsimPtr.data, tempvec, recsize);

    MPI_Reduce((void *) marstatsim2Ptr.data, (void *) tempvec, recsize, MPI_DOUBLE,
               MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&marstatsim2Ptr.data, &tempvec);
    copyit(marstatsim2Ptr.data, tempvec, recsize);

 // new mpi addtions -- doss -- yet to be checked
   
    MPI_Reduce((void *) oopmedexsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&oopmedexsimPtr.data, &tempvec);
    copyit(oopmedexsimPtr.data, tempvec, recsize);

    MPI_Reduce((void *) totalmedexsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&totalmedexsimPtr.data, &tempvec);
    copyit(totalmedexsimPtr.data, tempvec, recsize);

    MPI_Reduce((void *) xiindexsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&xiindexsimPtr.data, &tempvec);
    copyit(xiindexsimPtr.data, tempvec, recsize);

 // mu and q don't depend on anything determined within the model.  
 // We calculate them in the simulation out of laziness
	
	MPI_Reduce((void *) musimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&musimPtr.data, &tempvec);
    copyit(musimPtr.data, tempvec, recsize);

    MPI_Reduce((void *) qsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&qsimPtr.data, &tempvec);
    copyit(qsimPtr.data, tempvec, recsize);

    MPI_Reduce((void *) MedicaidsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&MedicaidsimPtr.data, &tempvec);
    copyit(MedicaidsimPtr.data, tempvec, recsize);

    MPI_Reduce((void *) zetasimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&zetasimPtr.data, &tempvec);
    copyit(zetasimPtr.data, tempvec, recsize);

    MPI_Reduce((void *) xisimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&xisimPtr.data, &tempvec);
    copyit(xisimPtr.data, tempvec, recsize);

    MPI_Reduce((void *) transfersimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                MPI_SUM, size - 1, MPI_COMM_WORLD);
    //     switchem(&transfersimPtr.data, &tempvec);
    copyit(transfersimPtr.data, tempvec, recsize);
   }
*/

   if (printOn==1){
    printf("Postreduce.\n");
    printf("marstatsimptr.data[0] is %f\n", marstatsimPtr.data[0]);
    printf("marstatsimMat[0][0] is %f\n", marstatsimMat[0][0]);
   }

   //   free(tempvec); //tempveec points to qsimptr now ...
   if (printOn==1)
   {
    WriteNonSimData();
    WriteSimData(assetsimMat, netIncomesimMat, consumptionsimMat, totalmedexsimMat, oopmedexsimMat, zetaindexsimMat,
                 marstatsimMat, agesim96Ptr,PIsim96Ptr, healthindexsimhMat, healthindexsimwMat);
    
    gau5write(strcat(strcpy(fullpath,rootdir),"asstsim.fmt"), assetsimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"cohsim.fmt"), cohsimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"netIncomesim.fmt"), netIncomesimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"ztaindxsim.fmt"), zetaindexsimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"xiindxsim.fmt"), xiindexsimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"ztasim.fmt"), zetasimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"xisim.fmt"), xisimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"totmedexsim.fmt"), totalmedexsimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"oopmedexsim.fmt"), oopmedexsimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"Medicaidsim.fmt"), MedicaidsimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"conssim.fmt"), consumptionsimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"beqsim.fmt"), beqsimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"mssim2.fmt"), marstatsim2Ptr);  
    gau5write(strcat(strcpy(fullpath,rootdir),"_qsim.fmt"), qsimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"_musim.fmt"), musimPtr);
    gau5write(strcat(strcpy(fullpath,rootdir),"transfersim.fmt"), transfersimPtr);
    
    end = clock();
    printf("This program ends in %5d minutes %5d seconds \n ",(end-start)/CLOCKS_PER_SEC/60,
         ((end-start)/CLOCKS_PER_SEC)%60);
    fprintf(timeStamp, "This program ends in %5d minutes %5d seconds \n ",(end-start)/CLOCKS_PER_SEC/60,
          ((end-start)/CLOCKS_PER_SEC)%60);
   }
   
   if (printOn==1)
   {
    printf("Postprinting.\n");
    printf("marstatsimptr.data[0] is %f\n", marstatsimPtr.data[0]);
    printf("marstatsimMat[0][0] is %f\n", marstatsimMat[0][0]);
    fclose(timeStamp); //didn't open unless printOn==1.
    fclose(logOutput);
   }

   if (errorPrintOn==1) fclose(errorOutput); //didnt open unless errorprintOn==1

// MPI-specific Code:
/*
   if (useMPI==2) MPI_Finalize();
*/   
   return 0;
}   /* End of main*/

/*--------------------------------------------------------------------------------*/
/*----------------------------------SUBROUTINES-----------------------------------*/

void Grid(double lambda, double aA[ADIM], double hsA[HSDIM], double IncomeA[IDIM])
{

   FILE *fileP;
   char fullpath[ADDRESS_LEN];
   int IInd, hsInd;
   int aInd;

   double lowValArith, highValArith, currValArith, incrArith;
   lowValArith = pow(AMIN,1/lambda);
   highValArith = pow(AMAX, 1/lambda);
   incrArith = (highValArith-lowValArith)/(ADIM-1);

   currValArith=lowValArith;
   for (aInd=0; aInd<ADIM; aInd++){ 
      aA[aInd] = pow(currValArith, lambda);
      currValArith += incrArith;
   }

   for (IInd = 0; IInd<IDIM; IInd++)
   {
      IncomeA[IInd]= ((double) IInd)/(IDIM-1); /* Income expressed as percentile ranking */
   }
      
   for (hsInd = 0; hsInd<HSDIM; hsInd++)
   {
      hsA[hsInd]=hsInd;
   }

   if (printOn==1)
   {
      fileP=fopen(strcat(strcpy(fullpath,outputdir),"aA.txt"),"w");
      for (aInd=0; aInd<ADIM; aInd++)
      { 
         fprintf(fileP,"%10.3f\n", aA[aInd]);
      }
      fclose(fileP);

      fileP=fopen(strcat(strcpy(fullpath,outputdir),"incomeA.txt"),"w");
      for (IInd=0; IInd<IDIM; IInd++)
      { 
         fprintf(fileP,"%10.3f\n", IncomeA[IInd]);
      }
      fclose(fileP);

      fileP=fopen(strcat(strcpy(fullpath,outputdir),"hsA.txt"),"w");
      for (hsInd=0; hsInd<HSDIM; hsInd++)
      { 
         fprintf(fileP,"%10.3f\n", hsA[hsInd]);
      }
      fclose(fileP);
   }
 
 
   if (printOn==1){
       fprintf(logOutput, "AMIN: %d\nAMAX: %d\nADIM: %d\nLAMBDA: %d\nIDIM: %d\nHSDIM: %d\nTDIMS: %d\nZETADIM: %d\nXIDIM: %d\n", 
               AMIN, AMAX, ADIM, LAMBDA, IDIM, HSDIM, TDIMS, ZETADIM, XIDIM);
   }
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* create grid for cash on hand   */ 

void Grid_Old(double aA[], double hsA[HSDIM], double IncomeA[IDIM])
{
   int aInd, IInd, hsInd;
   FILE *cashP; 
   char fullpath[ADDRESS_LEN];
   double aMin2, aMin3;

//* create asset grid, aMin<=a<=CASHMAX */

   for (aInd=0; aInd<ADIM1; aInd++)
   { 
      aA[aInd]=pow(sqrt(AMIN)+(sqrt(AMAX1)-sqrt(AMIN))*aInd/(ADIM1-1),2);
   }

   aMin2=aA[ADIM1-1];

   for (aInd=0; aInd<ADIM2; aInd++)
   { 
      aA[ADIM1+aInd]=pow(sqrt(aMin2)+(sqrt(AMAX2)-sqrt(aMin2))*(aInd+1)/(ADIM2),2);
   }

   aMin3=aA[ADIM1+ADIM2-1];

   for (aInd=1; aInd<ADIM-ADIM1-ADIM2+1; aInd++)
   {
      aA[aInd+ADIM1+ADIM2-1]
         = pow(sqrt(aMin3)+(sqrt(AMAX3)-sqrt(aMin3))*aInd/(ADIM-ADIM1-ADIM2),2);
   }

   for (IInd = 0; IInd<IDIM; IInd++)
   {
      IncomeA[IInd]= ((double) IInd)/(IDIM-1); /* Income expressed as percentile ranking */
   }
      
   for ( hsInd = 0; hsInd<HSDIM; hsInd++)
   {
      hsA[hsInd]=hsInd;
   }

   if (printOn==1)
   {
      cashP=fopen(strcat(strcpy(fullpath,outputdir),"xA.txt"),"w");

      for (aInd=0; aInd<ADIM1; aInd++)
      { 
         fprintf(cashP,"%10.3f\n", aA[aInd]);
      }

      aMin2=aA[ADIM1-1];

      for (aInd=0; aInd<ADIM2; aInd++)
      { 
         fprintf(cashP,"%10.3f\n", aA[ADIM1+aInd]);
      }

      aMin3=aA[ADIM1+ADIM2-1];

      for (aInd=1; aInd<ADIM-ADIM1-ADIM2+1; aInd++)
      {
         fprintf(cashP,"%10.3f\n", aA[aInd+ADIM1+ADIM2-1]);
      }
  
      fclose(cashP);
     
   }
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
 //Records nodes and weights for Gauss-Hermite Quadrature
 //     n: number of Gauss points to be used
 //     xArray[j]: j-th point
 //     wArray[j]: weight of the j-th point
 //     Nodes and weights are from Judd(1999), page 266 

void Quadrature(int n, double xArray[], double wArray[])
{
   int j;

   if (n==1)
   {
      xArray[0] = 0.0;
      wArray[0] = 1.0;
   }
   
   if (n==2)
   {
      xArray[1] = 0.7071067811;
      wArray[1] = 0.8862269254;
   }

   if (n==3)
   {
      xArray[1] = 0;
      wArray[1] = 1.181635900;
      xArray[2] = 1.224744871;
      wArray[2] = 0.2954089751;
   }

   if (n==4)
   {
      xArray[3] = 1.650680123;
      wArray[3] = 0.08131293544;
      xArray[2] = 0.5246476232;
      wArray[2] = 0.8049140900;}

   if (n==5)
   {
      xArray[4] = 2.020182870;
      wArray[4] = 0.01995324205;
      xArray[3] = 0.9585724646;
      wArray[3] = 0.3936193231;
      xArray[2] = 0;
      wArray[2] = 0.9453087204; 
   }

   if (n==6)
   {
      xArray[5]=0.2350604973e1;
      wArray[5]=0.4530009905e-2;
      xArray[4]=0.1335849074e1;
      wArray[4]=0.1570673203;
      xArray[3]=0.4360774119;
      wArray[3]=0.7246295952;
   }

   if (n==7)
   {
      xArray[3] = 0;
      wArray[3] = 0.8102646175;
      xArray[4] = 0.8162878828;
      wArray[4] = 0.4256072526;
      xArray[5] = 1.673551628;
      wArray[5] = 0.05451558281;
      xArray[6] = 2.651961356;
      wArray[6] = 0.0009717812450;
   }

   if (n==8)
   {
      xArray[7]=0.2930637420e1;
      wArray[7]=0.1996040722e-3;
      xArray[6]=0.1981656756e1;
      wArray[6]=0.1707798300e-1;
      xArray[5]=0.1157193712e1;
      wArray[5]=0.2078023258;
      xArray[4]=0.3811869902;
      wArray[4]=0.6611470125;
   }

   if (n==10)
   {
      xArray[5] = 0.3429013272;
      wArray[5] = 0.6108626337;
      xArray[6] = 1.036610829;
      wArray[6] = 0.2401386110;
      xArray[7] = 1.756683649;
      wArray[7] = 0.03387439445;
      xArray[8] = 2.532731674;
      wArray[8] = 0.001343645746;
      xArray[9] = 3.436159118;
      wArray[9] = 0.000007640432855;
   }

   for (j = 0;j<(n/2);j++)
   {
      xArray[j] = -xArray[n-1-j];
      wArray[j] = wArray[n-1-j];
   }
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* Discretization using Gauss Quadrature  
   This program discretizes a AR(1) process with autoregressive coefficient 
   rho and standard normal innoviation with standard deviation sigma into a 
   n-point first-order markov chain following Tauchen and Hussey(1991).
   z(t) = z(t-1)*rho+epsilon    
   epsilon ~ N(mu,sigma^2) 
   Also returns the invariant probability vector.
*/

void Discretization(int n, double rho, double mu, double sigma, double zArray[], 
                    double *piMatrixP, double *piInvarV)
{
   int h,i,j,k;  /* indices */

   double xArray[10];  /* xArray[j]: j-th point,used in Quadrature */
   double wArray[10];  /* wArray[j]: weight of the j-th point,used in Quadrature */
   double sum[10];
   double **ProductMat, **subtot;

   int m;  /* number of rows. If rho = 0 then AR(1) process is a iid series, and 
              we define m = 1, since every row is the same for transition matrix  */
    
   if ((sigma==0)||(rho==1)) /* not random */
   {
      for (i = 0;i<n;i++)
      {    
         zArray[i] =mu;

         for (j = 0; j< n; j++)
         {     
            *(piMatrixP+i*n+j) =0;
         }
         *(piMatrixP+i*n+i) =1;
      }
   }

/* detect errors */
   else if ((!((n<9)||(n==10)))||((rho>=1)||(rho<0))||(sigma<0))
   {
      printf("Error! The requirement for the AR(1) process is that the number of Guass points n=1,2,3,4,5,7,10, the autoregession coefficiency 0<=rho<1 and the standard deviation sigma is positive.\n");     
   }
   else
   {
   /* find discrete points and weights associated with each point */
      Quadrature(n,xArray,wArray);

   /* transformation 
      x = (z+mu/(rho-1))/(sqrt(2)*sigma) 
      thus z = sqrt(2.0)*sigma*x-mu/(rho-1) */ 
   
      for (i = 0;i<n;i++)
      {
         zArray[i] = sqrt(2.0)*sigma*xArray[i]-mu/(rho-1);
      }
      
   /* calculate  and print transition matrix */
      if ( rho ==0 ) m = 1;          /* iid process  */
      else m = n;     /* AR(1) process  */
            
      for (i = 0; i< m; i++)
      {     
         sum[i] = 0.0;
         for (j = 0; j<n; j++)
         {
         /* sum is used to normalize transition matrix */
            sum[i] = sum[i] +wArray[j]*exp(pow(xArray[j],2)-pow(xArray[j]-rho*xArray[i],2))/sqrt(PI);
         }         
         for (j = 0; j< n; j++)
         {     
            *(piMatrixP+i*n+j) = wArray[j]*exp(pow(xArray[j],2)-pow(xArray[j]-rho*xArray[i],2))/sqrt(PI)/sum[i];    
         }
      }
   } /* end of else */ 

/* Find stationary distribution using brute force, by repeatedly multiplying Pi by itself */

   ProductMat = (double **)malloc(n*sizeof(double *));
   subtot     = (double **)malloc(n*sizeof(double *));
   for(i=0; i<n; i++)
   {
      ProductMat[i] = (double *)malloc(n*sizeof(double));
      subtot[i]     = (double *)calloc(n,sizeof(double));
     for(j=0; j<n; j++)
     {
         ProductMat[i][j] = *(piMatrixP+i*n+j);
     }
   }

   for(h=0; h<500; h++)
   {
      for(i=0; i<n; i++)
      {
         for(j=0; j<n; j++)
          {
             subtot[i][j]=0;
            for(k=0; k<n; k++)          
                subtot[i][j] += ProductMat[i][k]*(*(piMatrixP+k*n+j));
        }
      }

      for(i=0; i<n; i++)
      {
         for(j=0; j<n; j++)
          {
            ProductMat[i][j] = subtot[i][j];
         }
      }
   } 

   for(i=0; i<n; i++)
   {
      piInvarV[i] = ProductMat[0][i];
   }
} 

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
// Allows IES_C and IES_M to differ

double getUtility(double cons, double totalMedex, double muRealization)
{
   double utility_C, utility_M, utility;
   if (nu==0) {printf("getUtility(): passed nu==0, exiting.\n"); exit(1);} // divide by zero error coming up
   if (omega==0) {printf("getUtility(): passed omega==0, exiting.\n"); exit(1);} // divide by zero error coming up

   if (nu==1.0)     /* log utility */
      utility_C = log(cons);
   else
      utility_C = pow(cons,1-nu)/(1-nu);

   if (omega==1.0)     /* log utility */
      utility_M = log(totalMedex);
   else
      utility_M = pow(totalMedex,1-omega)/(1-omega);

   utility = utility_C + muRealization*utility_M;
   
   return utility;
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* getXmin:  Uses bisection to find value of consumption that allows individual to
             enjoy guaranteed level of utility.  Utilizes optimal INTRA-temporal  
             allocation of funds between consumption and medex.
 */
double getXMin(double qRealization, double muRealization)
{
   double medexfrac, totalMedex, c_low, c_hi, c_mid, diff, diff2, thisUtility, xMin;
   if (nu==0) {printf("getUtility(): passed nu==0, exiting.\n"); exit(1);} // divide by zero error coming up
   if (omega==0) {printf("getUtility(): passed omega==0, exiting.\n"); exit(1);} // divide by zero error coming up

   if (switchFloor==0) // expenditure floor
   {
      xMin = eMin;
   }

   else
   {
      medexfrac = pow(muRealization/qRealization,1/omega); 	   
      c_low     = eMin/1000000;
      c_hi      = AMAX;
      diff      = 1;

      while (diff>(1e-9))
      {
         c_mid        = (c_hi+c_low)/2;
         totalMedex   = medexfrac*pow(c_mid,IESfrac);
         thisUtility  = getUtility(c_mid, totalMedex, muRealization);
         diff         = (uMin-thisUtility)/uMinAbs;
         diff2        = (c_hi-c_low)/c_mid;

         if (diff>0)
            c_low = c_mid;
         else
            c_hi  = c_mid;
 
         diff  = sqrt(diff*diff);
         if (diff > diff2) diff=diff2; // don't get stuck at corners
      }

      xMin = c_mid + qRealization*totalMedex;
   }

   return xMin;
}


/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/

void GetMu(double lnMuM_det[2][TDIMS][IDIM][HSDIM], double lnMuM_shockCoeffs[2][TDIMS][IDIM][HSDIM], 
           double xiA[XIDIM], double zetaA[ZETADIM], double muM[2][TDIMS][IDIM][HSDIM][ZETADIM][XIDIM])
{
   int tInd, gInd, hsInd, iInd, zetaInd, xiInd;

   if (switchPrefShifter==1)
   {   
      for (tInd=0; tInd<TDIMS; tInd++)
      {
         for (gInd=0; gInd<2; gInd++)
         {
            for (hsInd=0; hsInd<HSDIM; hsInd++)
            {
               for (iInd=0; iInd<IDIM; iInd++)
               {
                  for (zetaInd=0; zetaInd<ZETADIM; zetaInd++)
                  {
                     for (xiInd=0; xiInd<XIDIM; xiInd++)
                     {   
                        muM[gInd][tInd][iInd][hsInd][zetaInd][xiInd] = exp(lnMuM_det[gInd][tInd][iInd][hsInd]
                         + lnMuM_shockCoeffs[gInd][tInd][iInd][hsInd]*(zetaA[zetaInd]+xiA[xiInd]));
                      }
                  }
               }
            }
         } 
      }
   }
   else
   {
      for (tInd=0; tInd<TDIMS; tInd++)
      {
         for (gInd=0; gInd<2; gInd++)
         {
            for (hsInd=0; hsInd<HSDIM; hsInd++)
            {
               for (iInd=0; iInd<IDIM; iInd++)
               {
                  for (zetaInd=0; zetaInd<ZETADIM; zetaInd++)
                  {
                     for (xiInd=0; xiInd<XIDIM; xiInd++)
                     {
                        muM[gInd][tInd][iInd][hsInd][zetaInd][xiInd] = 1e-20; // Need non-zero value to prevent /0 error.
                     }
                  }
               }
            }
         }
      }
   }
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/

double GetContinuationUtility(int aTomorrowInd, int tInd, int gInd, int hsInd, 
                              int iInd, int zetaInd)
{
   int hsIndN, zetaIndN, xiIndN;
   double sum; /* sum is the expected value next period */   
   sum = 0.0; 
   
/* Loop through possible next period states */   
   for (hsIndN = 0; hsIndN<HSDIM; hsIndN++)
   {
       for (zetaIndN = 0; zetaIndN<ZETADIM; zetaIndN++)
       {
          for (xiIndN = 0; xiIndN<XIDIM; xiIndN++)   
          {
             sum+= hsProbM[gInd][tInd][iInd][hsInd][hsIndN]*zetaProbM[zetaInd][zetaIndN]*xiProbA[xiIndN]
             *valueFunM[gInd][tInd+1][iInd][hsIndN][zetaIndN][aTomorrowInd][xiIndN];
          }
       }
   } /* end looping through possible next period states */

   return sum; 
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  net bequest function  */

double NetBequest(double amountBequestedP)
{
   double netBequest = 0; 
   if (amountBequestedP>exBeq)
      netBequest = exBeq+(1-tauBeq)*(amountBequestedP-exBeq); 
   else netBequest = amountBequestedP; 
   return netBequest; 
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* utility from a net bequest for a single person
   phi_j(b_net) = phi_j*( (b_net+K_j)^(1-nu) )/(1-nu) 
   not discounted */

double UBeqSingle(double netBequestP)
{
   double utils;    
   if (nu  ==  1)   /* log utility */
      utils  =  phi0*log(netBequestP+K0); 
   else
      utils  =  phi0*pow(netBequestP+K0, 1-nu)/(1-nu);         
   return utils; 
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* Expected utility from leaving bequest matrix for a single, not discounted */

void GetUtilityBeq(double bequestUM[ADIM])
{
   int aInd;
   for (aInd = 0; aInd<ADIM; aInd++)
   {
      bequestUM[aInd] = UBeqSingle(NetBequest(aA[aInd]));
   }
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* After-tax income at bracket points, used to calculate income after tax*/

void IncomeAtBrk(double taxBrk[], double taxMar[], double incomeBrk[])
{
   int j; 

   incomeBrk[0] = (1-taxMar[0])*taxBrk[0]; /* The leftmost interval  */
   for (j = 1; j<(taxDim-1); j++)
   {
      incomeBrk[j] = incomeBrk[j-1]+(1-taxMar[j])*(taxBrk[j]-taxBrk[j-1]); 
   }
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* Calculate income after tax 
   netIncome(j) = netIncome(j-1)+(1-taxMar(j))*(taxBrk(j)-taxBrk(j-1))  */

double AfterTaxIncome(double y)
{
   if (y<0)
   {
   /* printf("Error! negtive gross income!\n "); */
      return -1;  /* this case will be ruled out in maximization  */
   }
   else if (y<taxBrk[0])
      return ((1-taxMar[0])*(y)); 
   else if (y<taxBrk[1])
      return (incomeBrk[0]+(1-taxMar[1])*(y-taxBrk[0])); 
   else if (y<taxBrk[2])
      return (incomeBrk[1]+(1-taxMar[2])*(y-taxBrk[1])); 
   else if (y<taxBrk[3])
      return (incomeBrk[2]+(1-taxMar[3])*(y-taxBrk[2])); 
   else if (y<taxBrk[4])
      return (incomeBrk[3]+(1-taxMar[4])*(y-taxBrk[3])); 
   else if (y<taxBrk[5])
      return (incomeBrk[4]+(1-taxMar[5])*(y-taxBrk[4])); 
   else 
      return (incomeBrk[5]+(1-taxMar[6])*(y-taxBrk[5])); 
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  Locates nearest point _below_ x in a sorted array
    From Numerical Recipes in C, p. 117
*/

int Locate(double *Xarray, double x, int DIM)
{
   int j_L, j_U, j_M, ascend, dif, j_Star;

   ascend = 1;
   if (Xarray[DIM-1]<Xarray[0]) ascend = 0;

   j_L = 0;
   j_U = DIM-1;
   dif = j_U-j_L;

   if (ascend==1)
   {
      while (dif>1)
      {
         j_M = (int) (j_U+j_L)/2;
         if (x>Xarray[j_M]) 
            j_L = j_M;
         else 
            j_U = j_M;

         dif = j_U-j_L;
      }
      j_Star = j_L;
   }

   else
   {
      while (dif>1);
      {
         j_M = (int) (j_U+j_L)/2;
         if (x<Xarray[j_M]) 
            j_L = j_M;
         else 
            j_U = j_M;
         dif = j_U-j_L;
      }
      j_Star = j_L;
   }

   return j_Star;
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* Returns linearly interpolated or extrapolated value at x */

double Interpolation(double *fP,  double *xP,  double x, int DIM)
{
   int j;

   j = Locate(xP,x,DIM) + 1;
   //return (x-*(xP+j))*(*(fP+j-1)-*(fP+j))/(*(xP+j-1)-*(xP+j))+(*(fP+j));  /*x[j-1]<x<x[j] */
   return (x-*(xP+j))/(*(xP+j-1)-*(xP+j)) * (*(fP+j-1)-*(fP+j)) + (*(fP+j));  /*x[j-1]<x<x[j] */
   //return ((x-xP[j])*(fP[j-1]) - fP[j]) / (xP[j-1]-xP[j]) + fP[j];
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* Returns grid point and linear interpolation weight */
struct result GetLocation(double *xP,  double x, int DIM)
{
   int j;
   struct result result1;

   j = Locate(xP,x,DIM)+1;
   result1.Ind1=j-1;
   result1.weight=(*(xP+j)-x)/(*(xP+j)-*(xP+j-1));
   return result1;
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/

void GetCdf(int nRows, int nCols, double *piMatrixP, double *piMatrixCDF) 
{
   int iRow, iCol;
   double sum;
   
   for (iRow=0; iRow<nRows; iRow++)
   {
      sum = 0;   
      for (iCol=0; iCol<nCols+1; iCol++)
      {
         *(piMatrixCDF + iRow*(nCols+1) + iCol) = sum;
         sum += *(piMatrixP + iRow*nCols + iCol);
      }
   }
return;
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/

//Write data as column vectors to ascii files;
//Looping over indices in same order they are written so they can be read in using
// a row-major order index function.
void WriteNonSimData(){
   FILE *valueFP;  /*point to files related to singles case */
   FILE *consumFP;
   FILE *medFP;
   FILE *saveFP;  
   FILE *oopMedFP;
   FILE *fp;
   int i, sexInd, tInd, IInd, hsInd, hsIndN, zetaInd, aInd, xiInd; //i and sexInd are used for same thing...
   char fullpath[ADDRESS_LEN];


   fp = fopen(strcat(strcpy(fullpath,outputdir),"imcprof1.txt"), "w");
   for (i = 0; i<2; i++) /* for (i = 0; i<3; i++) ***** */
   {
      for (IInd = 0; IInd<IDIM; IInd++)
      {
         for (tInd = 0; tInd<TDIMS; tInd++)
         {   
            fprintf(fp,"%20.8f\n",yM[i][IInd][tInd]);
         }
      }
   }    
   fclose(fp);

//*-------------------------------------------------------------------------------*/

   fp = fopen(strcat(strcpy(fullpath,outputdir),"survivalPro1.txt"), "w"); /* open file*/

   for (i = 0; i<2; i++)
   {
      for (tInd = 0; tInd<TDIMS; tInd++)
      {
         for ( hsInd = 0; hsInd<HSDIM; hsInd++)
         {
            for (IInd = 0; IInd<IDIM; IInd++)
            {
               fprintf(fp, "%20.8f\n",survivalProbM[i][tInd][hsInd][IInd]);
            }
         }
      }
   } 
   fclose(fp);

//*--------------------------------------------------------------------------------*/
   fp = fopen(strcat(strcpy(fullpath,outputdir),"profheal1.txt"), "w"); 
   for (i = 0; i<2; i++)
   {
      for (tInd = 0; tInd<TDIMS; tInd++)
      {
         for (IInd = 0; IInd<IDIM; IInd++)
         {
            for ( hsInd = 0; hsInd<HSDIM; hsInd++)
            {
               for ( hsIndN = 0; hsIndN<HSDIM; hsIndN++)
               {
                  fprintf(fp,"%20.8f\n",hsProbM[i][tInd][IInd][hsInd][hsIndN]);
               }
            }
         }
      } 
   }
   fclose(fp);
   //*--------------------------------------------------------------------------------*/
//* mu and q, the pref shifter and the shareof med expenses */   
   fp = fopen(strcat(strcpy(fullpath,outputdir),"muM.txt"), "w");
   for (sexInd = 0; sexInd<2; sexInd++)
   {
     for (tInd = 0; tInd<TDIMS; tInd++)
      {
         for (IInd = 0; IInd<IDIM; IInd++)
         {
            for (hsInd = 0; hsInd<HSDIM; hsInd++)
            {
            for (zetaInd=0; zetaInd<ZETADIM; zetaInd++)
            {
               for (xiInd=0; xiInd<XIDIM; xiInd++)
               {
                  fprintf(fp, "%16.4f \n", muM[sexInd][tInd][IInd][hsInd][zetaInd][xiInd]);
               }
            }
            }
         }
      }
   }
   fclose(fp);

   fp = fopen(strcat(strcpy(fullpath,outputdir),"qM.txt"), "w");
   for (sexInd = 0; sexInd<2; sexInd++)
   {
      for ( hsInd = 0; hsInd<HSDIM; hsInd++)
      {
        for (tInd = 0; tInd<TDIMS; tInd++)   /*JJJJJJJ  time reordered to be column dimension */
        {
           fprintf(fp, "%20.8f\n", qM[sexInd][tInd][hsInd]);
        }
      }
   }
   fclose(fp);

   
//* value and policy functions */
//*--------------------------------------------------------------------------------*/

   valueFP  = fopen(strcat(strcpy(fullpath,outputdir),"valueF.txt"),"w");  /* value function */
   consumFP = fopen(strcat(strcpy(fullpath,outputdir),"consumptionF.txt"),"w");  /* consumption policy function */
   medFP    = fopen(strcat(strcpy(fullpath,outputdir),"medicalF.txt"),"w");  /* medical expenditures policy function */
   oopMedFP = fopen(strcat(strcpy(fullpath,outputdir),"outofpocketmedF.txt"),"w");  /* outofpocket medical expenditures function */
   saveFP   = fopen(strcat(strcpy(fullpath,outputdir),"savingsF.txt"),"w");  /* bequest policy function */
   
//* Print value function for t = 1, ..., T  where value function for T+1 is defined as 0 matrix */
//* single case */
   for (i = 0; i<2; i++)
   {
      for (tInd=0; tInd<TDIMS; tInd++) //for (tInd = (TDIMS-1); tInd>= 0; tInd--) //6/9/08 time was written in reverse order and read (in MATLAB) in reverse order; we switched both
      {
         for (IInd = 0; IInd<IDIM; IInd++)
         {
            for (hsInd = 0; hsInd<HSDIM; hsInd++)
            {
               for (zetaInd = 0; zetaInd<ZETADIM; zetaInd++)
               {
                  for (aInd = 0; aInd<ADIM; aInd++)
              {
                 for (xiInd=0; xiInd<XIDIM; xiInd++)
                 {
                 /* write value function and policy function to files */
                   fprintf(valueFP, "%22.20f \n", valueFunM[i][tInd][IInd][hsInd][zetaInd][aInd][xiInd]); //note: we dont write valuefun_t+1.
                   fprintf(consumFP, "%22.20f \n", consumFunM[i][tInd][IInd][hsInd][zetaInd][aInd][xiInd]);
                   fprintf(medFP, "%22.20f \n", medFunM[i][tInd][IInd][hsInd][zetaInd][aInd][xiInd]);
                   fprintf(saveFP, "%22.20f \n", savingsFunM[i][tInd][IInd][hsInd][zetaInd][aInd][xiInd]);
                   fprintf(oopMedFP, "%22.20f \n", outofpocketmedFunM[i][tInd][IInd][hsInd][zetaInd][aInd][xiInd]);

                 }//xi
                 //fprintf(valueFP, "\n");
                 //fprintf(consumFP, "\n");
                 //fprintf(medFP, "\n");
                 //fprintf(saveFP, "\n");
                 //fprintf(oopMedFP, "\n");
              }//assets
               //fprintf(valueFP, "\n");
               //fprintf(consumFP, "\n");
               //fprintf(medFP, "\n");
               //fprintf(saveFP, "\n");   
               //fprintf(oopMedFP, "\n");
               }
            }/* end loop through health status */
         } /* end loop through permanent income */   
      } /* end loop through age */
   } /* end loop through sex */

//* End writing to files  */
   fclose(valueFP);
   fclose(consumFP);
   fclose(saveFP);
   fclose(medFP);
   fclose(oopMedFP);
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
// Need to check all time indices are done correctly for new code.


void WriteSimData(double **assetsimMat, double **netIncomesimMat, double **consumptionsimMat, 
                  double **totalmedexsimMat, double **oopmedexsimMat,  double **zetaindexsimMat, 
                  double **marstatsimMat, GMatrix agesim96Ptr, GMatrix PIsim96Ptr, 
                  double **healthsimhMat, double **healthsimwMat)

{ 
   FILE *fp;
   char fullpath[ADDRESS_LEN];
   int yearInd, personInd;

//* simulation results */
   fp = fopen(strcat(strcpy(fullpath,outputdir),"assetsimMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
    for (yearInd=0; yearInd<TDimSims+1; yearInd++)  /*calendar year */
      {
       fprintf(fp, "%16.4f ", assetsimMat[yearInd][personInd]);
      }
    fprintf(fp, "\n");
   }
   fclose(fp);


   fp = fopen(strcat(strcpy(fullpath,outputdir),"zetaindsimMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
   // for (yearInd=0; yearInd<TDimSims; yearInd++)   // calendar year
      for (yearInd=0; yearInd<TDimSims+1; yearInd++) // calendar year/
      {
         fprintf(fp, "%16.4f ", zetaindexsimMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);

   fp = fopen(strcat(strcpy(fullpath,outputdir),"consumptionsimMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {

      for (yearInd=0; yearInd<TDimSims; yearInd++)  /*calendar year */
      {
         fprintf(fp, "%16.4f ", consumptionsimMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);

   fp = fopen(strcat(strcpy(fullpath,outputdir),"totalmedexsimMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims; yearInd++)  /*calendar year */
      {
         fprintf(fp, "%16.4f ", totalmedexsimMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);

   fp = fopen(strcat(strcpy(fullpath,outputdir),"oopmedexsimMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims; yearInd++)  /*calendar year */
      {
         fprintf(fp, "%16.4f ", oopmedexsimMat[yearInd][personInd]); //extra space is important.
      }
      fprintf(fp, "\n");
   }
   fclose(fp);

   fp = fopen(strcat(strcpy(fullpath,outputdir),"marstatsimMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims; yearInd++)
      {
         fprintf(fp, "%16.4f ", marstatsimMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);


   fp = fopen(strcat(strcpy(fullpath,outputdir),"PI96.txt"), "w");

   for(personInd=0; personInd<nsims; personInd++)
   {
      fprintf(fp, "%16.4f ", PIsim96Ptr.data[personInd] );
   }
   fclose(fp);

   fp = fopen(strcat(strcpy(fullpath,outputdir),"age96.txt"), "w");

   for(personInd=0; personInd<nsims; personInd++)
   {
      fprintf(fp, "%16.4f ", agesim96Ptr.data[personInd] );
   }
   fclose(fp);


   fp = fopen(strcat(strcpy(fullpath,outputdir),"healthsimhMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims; yearInd++)  /*calendar year */
      {
         fprintf(fp, "%16.4f ", healthsimhMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);


   fp = fopen(strcat(strcpy(fullpath,outputdir),"healthsimwMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims; yearInd++)  /*calendar year */
      {
         fprintf(fp, "%16.4f ", healthsimwMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);


   fp = fopen(strcat(strcpy(fullpath,outputdir),"PIMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims; yearInd++)  /*calendar year */
      {
         fprintf(fp, "%16.4f ", netIncomesimMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);

}


/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/

int GetProfiles(double *agevecPtr, double *prefShifter_mean,
                double *prefShifter_shockCoeffs,
                double *muPIPtr, double *qMPtr, double *mortrateSPtr,  
                double *mortratePIPtr, double *hstranSPtr, 
                double *hstranPIPtr, double *yprofPtr, double *yprofPIPtr,
                double lnMuM_det[2][TDIMS][IDIM][HSDIM],
                double lnMu_shockCoeffs[2][TDIMS][IDIM][HSDIM]) 				           
{
   int i, j, k, varshift, MSInd, sexInd, tInd, IInd, hsInd, tdimGAUSS; 

   double income, survprobS, PI1, PI2,   
          PIC1, PIC2, PIC11, PIC12, PIC21, PIC22, PIC31, PIC32, PIC41, PIC42,
          gg, bb, _g, _b, bhprob;
           
   tdimGAUSS = (int) floor(rounder+agevecPtr[2]);

   if (tdimGAUSS != TDIMS)
   {
      fprintf(errorOutput, "Incompatible timespans!!! \n");      
   }

//* Income */
   j=0;
   i=0;
   for (tInd = 0; tInd<TDIMS; tInd++)
   {
      for (MSInd = 0; MSInd<2; MSInd++) /* [MSind == sexInd?] Single male, single female */  
      {
         income = yprofPtr[i];
         PIC1   = yprofPIPtr[j];
         PIC2   = yprofPIPtr[j+1];
         i++;
         j+=2;
         for (IInd = 0; IInd<IDIM; IInd++)
         {
            PI1 = IncomeA[IInd];
            PI2 = pow(PI1,2);
            yM[MSInd][IInd][tInd] = exp(income + PIC1*PI1 + PIC2*PI2);
         }
      }
   }
  
//* Health Status transition probabilities for singles */
//* Health expense parameters for singles */
//* Survival Probabilities for singles */

   i=0;
   j=0;
   k=0;
   varshift  = tdimGAUSS*4; /* coefficients on PI and PI^2 for male and female, by age*/

   for (tInd = 0; tInd<TDIMS; tInd++)
   {
      for (sexInd = 0; sexInd<2; sexInd++)
      {
         PIC11 = mortratePIPtr[j];
         PIC12 = mortratePIPtr[j+1];

         PIC21 = hstranPIPtr[j];
         PIC22 = hstranPIPtr[j+1];

         PIC31 = muPIPtr[k];//PI coefficients affecting stochastic part of pref-shifter (mu); note they include interaction term [e.g. coeff*t*g or some such]
         PIC32 = muPIPtr[k+1];

         PIC41 = muPIPtr[k+varshift]; //PI coefficients affecting stochastic part of pref-shifter (mu); note they include interaction term [e.g. coeff*t*g or some such]
         PIC42 = muPIPtr[k+1+varshift];

         j+=2;
         k+=2;

         for (IInd = 0; IInd<IDIM; IInd++)
         {
            PI1 = IncomeA[IInd];
            PI2 = pow(PI1,2);

            bb  = Logit(hstranSPtr[i] + PIC21*PI1 + PIC22*PI2);   /* Pr(hs_t+2=bad|hs_t=bad) */
            gg  = Logit(hstranSPtr[i+1] + PIC21*PI1 + PIC22*PI2); /* Pr(hs_t+2=bad|hs_t=good) */
            gg  = 1-gg;  /* Pr(hs_t+2=good|hs_t=good) */

            TwoYearToOneYear(gg, bb, &_g, &_b);
      
            for ( hsInd = 0; hsInd<HSDIM; hsInd++)  /* Bad health, then good*/
            {
               survprobS = mortrateSPtr[i+hsInd];
               survivalProbM[sexInd][tInd][hsInd][IInd] = LogitSQRT(survprobS + PIC11*PI1 + PIC12*PI2);
               
               if (hsInd==0) bhprob = _b;
               else bhprob = 1-_g;

               hsProbM[sexInd][tInd][IInd][hsInd][0] = bhprob;   /* Bad health at t+1 */
               hsProbM[sexInd][tInd][IInd][hsInd][1] = 1-bhprob; /* Good health at t+1 */

         // We calculate the deterministic portion of mu (pref shifter) here;
         // PIC31 and PIC32 are indexed by t and g now, which is how interactions are handled.
            lnMuM_det[sexInd][tInd][IInd][hsInd] = prefShifter_mean[i+hsInd] + PIC31*PI1 + PIC32*PI2;
         // calculate coefficients for the stochastic portion of mu here; these multiply against (Zeta+Xi)
             lnMu_shockCoeffs[sexInd][tInd][IInd][hsInd] = sqrt(pow(prefShifter_shockCoeffs[i+hsInd],2) + PIC41*PI1+PIC42*PI2); 
            }
         }
         i+=2;
      }
   }

   /* kill men and women by 100 */
   for ( hsInd = 0; hsInd<HSDIM; hsInd++)  /* Bad health, then good*/
    {
      for (IInd = 0; IInd<IDIM; IInd++)
      {
         survivalProbM[1][TDIMS-1][hsInd][IInd] = 0; // women 
         survivalProbM[0][TDIMS-1][hsInd][IInd] = 0; // men
      }
    }   
    
   i=0;
   for (tInd = 0; tInd<TDIMS; tInd++)
   {
      for (sexInd = 0; sexInd<2; sexInd++)
     {
        for (hsInd=0; hsInd<HSDIM; hsInd++)
        {
           qM[sexInd][tInd][hsInd] = qMPtr[i]; // qMPtr is identical format to prefShifter_mean
           i++;
        }
     }
   }


   if (switchTax==0) 
   {
      memset(taxMar,0, taxDim*sizeof(double));  /* no tax */
   }

   if (switchY==0) 
   { /*   Annual Income */
      for (tInd = 0; tInd<TDIMS; tInd++)
      {
         for (MSInd = 0; MSInd<2; MSInd++) /* Single male, single female */ 
         {
            for (IInd = 0; IInd<IDIM; IInd++)
            {
               yM[MSInd][IInd][tInd] = 0.0;            
            }
         }
      }
   }

   if (switchMor==0) 
   { 
   /* Survival Probabilities */

      for (tInd = 0; tInd<TDIMS-1; tInd++) 
      {
         for (sexInd = 0; sexInd<2; sexInd++)
         {
            for ( hsInd = 0; hsInd<HSDIM; hsInd++)  /* Bad health, then good */
            {
               for (IInd = 0; IInd<IDIM; IInd++)
               {
                  survivalProbM[sexInd][tInd][hsInd][IInd] = 1.0;                  
               }
            }
         }
      }  
   }   
   return 1;
}//end GetProfiles()

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  LOGIT  */

double Logit(double x)
{
   return exp(x)/(1+exp(x));
}
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  LOGITSQRT  */

double LogitSQRT(double x)
{
   return sqrt(exp(x)/(1+exp(x)));
}
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  TWOYEARTOONEYEAR:  Converts 2-year transition probes to 1-year probs
                       Using formula by O. Nartova
                       gg = pr(h_t+2=good|h_t=good); bb = pr(h_t+2=bad|h_t=bad)
                       _g = pr(h_t+1=good|h_t=good); _b = pr(h_t+1=bad|h_t=bad)
*/

void TwoYearToOneYear(double gg, double bb, double *_g, double *_b)
{
   double big_A, big_B, big_C;

   big_A = 2-gg-bb;
   big_B = 2*(gg-1);
   big_C = 1-bb-gg+(bb*bb);

   *_b   = -big_B + sqrt(pow(big_B,2) - 4*big_A*big_C);
   *_b   = (*_b)/(2*big_A);
   *_g   = ( (1-bb) - (*_b)*(1-(*_b)) )/(1-(*_b));

   return;
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  GLOBDEREF:  Dereference a bunch of pointers to C++ globals  */

int globderef(double *prefvecPtr, double *asstvecPtr, double *medexvecPtr, 
           double *simvecPtr, double *agevecPtr, double *switchvecPtr)
{
//* Flow utility parameters:
//   beta = discount factor
//   d(m(i)) = 1+d*1{m(i) = good}
//   u(c, m(i)) = d(m(i))*(c^(1-nu))/(1-nu)
//*/
   delta   = prefvecPtr[0];
   beta    = prefvecPtr[1];
   nu      = prefvecPtr[2];
   omega   = prefvecPtr[3];
   IESfrac = nu/omega;

   eMin    = asstvecPtr[0]; /* Minimum expenditure provided by the government */
   uMin    = getUtility(eMin, eMin, 1); /* Utility implied by minimum expenditure */  
   uMinAbs = sqrt(uMin*uMin);
   tauBeq  = asstvecPtr[1]; /* Estate tax rate */
   exBeq   = asstvecPtr[2]; /* Estate tax exemption level */

// interest rate r is an i.i.d. random variable, with mean mu_r and variance (sigma_r)^2.
   mu_r          = asstvecPtr[3];
// sigma_r       = asstvecPtr[4]; 

   rhoHc         = medexvecPtr[0];
   sigma_epsilon = sqrt(medexvecPtr[2]);
   sigma_xi      = sqrt(medexvecPtr[3]);
   medex_bottomcode = medexvecPtr[4]; 


   nsims         = (int) floor(rounder+simvecPtr[0]);
   TDimSims      = (int) floor(rounder+simvecPtr[1]);
// if (TDimSims >= TDIMS) {exit(1);}

   TSTART        = (int) floor(rounder+agevecPtr[0]);

   switchMor     = (int) floor(rounder+switchvecPtr[0]);
   switchBeta    = (int) floor(rounder+switchvecPtr[1]);
   switchY       = (int) floor(rounder+switchvecPtr[2]);
   switchPrefShifter   = (int) floor(rounder+switchvecPtr[3]);
   switchTax     = (int) floor(rounder+switchvecPtr[4]);
   switchZeta    = (int) floor(rounder+switchvecPtr[5]);
   switchXi      = (int) floor(rounder+switchvecPtr[6]);
// switchR       = (int) floor(rounder+switchvecPtr[7]);
   switchBeq     = (int) floor(rounder+switchvecPtr[8]);
   switchCMin    = switchvecPtr[9]; 
   switchFloor   = (int) floor(rounder+switchvecPtr[10]); 

   if (switchBeta==1) beta=1;

// Bequest utility parameters for singles:
// Phi_j(b_net) = phi_j*[(b_net+K_j)^(1-nu)]/[1-nu] */

   phi0 = prefvecPtr[4]*((double) switchBeq);
   K0   = prefvecPtr[5];
   
   if (printOn==1)
   {
      fprintf(logOutput, "delta: %f\nbeta: %f\nnu: %f\nomega: %f\neMin: %f\ntauBeq: %f\nexBeq: %f\nmu_r: %f\nrhoHC: %f\nsigma_epsilon: %f\nsigma_xi: %f\nphi0: %f\nK0: %f\n", 
              delta, beta, nu, omega, eMin, tauBeq, exBeq, mu_r, rhoHc, sigma_epsilon, sigma_xi, phi0, K0);
      fprintf(logOutput, "switchMor: %d\nswitchBeta: %d\nswitchY: %d\nswitchPrefShifter: %d\nswitchTax: %d\n", 
              switchMor, switchBeta, switchY, switchPrefShifter, switchTax);
      fprintf(logOutput, "switchZeta: %d\nswitchXi: %d\nswitchBeq: %d\n", 
              switchZeta, switchXi, switchBeq);
      fprintf(logOutput, "switchCMin: %f\nswitchCMin: %f\n", switchCMin); 
      fprintf(logOutput, "switchFloor: %d\nswitchCMin: %d\n", switchFloor); 
   }

   return 1;
} 

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  SETUPSIM:  Initialize matrices that hold simulation results  
               Each row of dataMat is a subsection of dataVec    
               This economizes on memory and reduces copying  
*/
double **SetUpSim(double *dataVec, int extrayears)

{
   double **dataMat;
   int iYear;

   dataMat = (double **)malloc((TDimSims+extrayears)*sizeof(double *));
   for(iYear=0; iYear<(TDimSims+extrayears); iYear++)
   {
      dataMat[iYear] = &dataVec[iYear*nsims];
   }

   return dataMat;
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  ZEROMAT:  Initializes matrix structure and fills it with zeros  */

GMatrix zeromat(unsigned int recsize)
{
   GMatrix mat = {1, 1};
   mat.m = recsize;
   mat.n = 1;
   mat.data = (double *)calloc(recsize,sizeof(double));
   return mat;
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  SWITCHEM:  Switches two pointer addresses  */

void switchem(double **dataPtr, double **tempPtr)
{
   double *tempPtr2;
   tempPtr2 = *dataPtr;  /* dataPtr, tempPtr are pointers to pointers.  */
   *dataPtr = *tempPtr;  /* Dereferencing them yields memory addresses, */
   *tempPtr = tempPtr2;  /* rather than values. */
}
//*--------------------------------------------------------------------------------*/
//*--------------------------------------------------------------------------------*/
void copyit(double *permvec, double *tempvec, int recsize)
{
   int i;
   for(i=0;i<recsize;i++) permvec[i]=tempvec[i];
   return;
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  findConsumption:  Finds value of consumption that satisfies INTRA-temporal  
                      allocation of funds between consumption and medex using 
                      bisection
*/
double findConsumption(double allocfrac, double totalfunds)
{
   double c_low, c_hi, c_mid, diff, impliedtotal;

   c_low = 0;
   c_hi  = totalfunds;
   diff  = 1;

   while (diff>(1e-8))
   {
      c_mid = (c_hi+c_low)/2;
	  impliedtotal = c_mid + allocfrac*pow(c_mid,IESfrac);
	  diff  = (totalfunds-impliedtotal)/totalfunds;

      if (diff>0)
         c_low = c_mid;
      else
         c_hi  = c_mid;

      diff = sqrt(diff*diff);
   }

   return c_mid;
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* Solve value function for t = 1, ..., T   
   Singles case 
*/
void GetRulesSingle(int iAssetsmin,int iAssetsmax)
{
   int gInd, tInd, iInd, hsInd,  zetaInd, xiInd, aInd, savingsInd,
       jMax; //jMax is max index for feasible assets tomorrow
   
   double diff, todaysUtility, cons, maxCons, value, totalMedex, maxMedex, OOPMedex,
          maxOOP, continuationUtilityAlive;
   double muRealization, qRealization, allocfrac, xMin;
   double maximum = 0.0, x;     /* used in maximization to record maximum, x is cash-on-hand*/
   int maxSavingsInd;            /* used in maximization to record index of maximizer */

   double *tempvec; //this may need to be changed? lookup MPI code
   int recsize = IDIM*HSDIM*ZETADIM*ADIM*XIDIM; // all dimensions except gender and time
   tempvec = (double *)calloc(recsize,sizeof(double));

   for (gInd = 0; gInd<2; gInd++)
   {
      if  (gInd == 0)
      {
         if (printOn==1)
            printf("Solving policy functions. Case of a male (rank==%d).\n", rank); 
      }
      else 
      {
         if (printOn==1)
            printf("Solving policy functions. Case of a female (rank==%d).\n", rank); 
      }

      for (tInd = (TDIMS-1); tInd>= 0; tInd--)
      {
         for (iInd = 0; iInd<IDIM; iInd++)
         {
            for (hsInd = 0; hsInd<HSDIM; hsInd++)
            {
               qRealization = qM[gInd][tInd][hsInd]; // out-of-pocket share for given state (doesn't rely on stochastic elements)
               for (aInd = iAssetsmin; aInd<iAssetsmax; aInd++) 
               {   
               // Compute cash-on-hand
                  x = aA[aInd] + AfterTaxIncome(mu_r*aA[aInd] + yM[gInd][iInd][tInd]);
               
                  for (zetaInd=0; zetaInd<ZETADIM; zetaInd++)
                  {
                     for (xiInd=0; xiInd<XIDIM; xiInd++)
                     {//last state variable for today
                        maximum       = -10000000000;
                        maxSavingsInd = -1; //if you end with -1, get an appropriate error.
                        maxCons       = -1;
                        maxMedex      = -1;
                        maxOOP        = -1;
                        muRealization = muM[gInd][tInd][iInd][hsInd][zetaInd][xiInd]; //current medical-needs shock
                        allocfrac     = qRealization*pow(muRealization/qRealization,1/omega); 
                     
                        xMin          = getXMin(qRealization, muRealization);
                  
                     // People receiving transfers have to consume xMin; however, other people can be flexible
                        if (x<xMin)
                        {
                           cons          = findConsumption(allocfrac, x); 
                           OOPMedex      = allocfrac*pow(cons,IESfrac); 
                           cons          = findConsumption(allocfrac, xMin);
                           totalMedex    = allocfrac*pow(cons,IESfrac)/qRealization; 
                           savingsInd    = 0;
                           todaysUtility = getUtility(cons, totalMedex, muRealization);

                           if (tInd==(TDIMS-1))
                              continuationUtilityAlive = 0;
                           else
                              continuationUtilityAlive = GetContinuationUtility(savingsInd, tInd, gInd, hsInd, iInd, zetaInd);

                           value = todaysUtility + beta*survivalProbM[gInd][tInd][hsInd][iInd]*continuationUtilityAlive 
                                                 + beta*(1-survivalProbM[gInd][tInd][hsInd][iInd])*bequestUM[savingsInd];
                           maximum  = value;
                           maxSavingsInd = savingsInd;
                           maxCons  = cons;
                           maxMedex = totalMedex;
                           maxOOP   = OOPMedex;  // Omit transfer-funded spending. 
                        }

                        else
                        {	 
                        // savings for tomorrow (a_t+1) cannot exceed COH today
                        // we find the index for the highest feasible savings
                           for (savingsInd=0; savingsInd<ADIM; savingsInd++)
                           {// aA is increasing in savingsInd
                              diff = x - aA[savingsInd] - xMin*switchCMin; 
                              if (diff<=0) break; //sign(diff) = sign(consumption); and if diff==0 then cons==0 which we dont want.
                           }
                           jMax = savingsInd-1; 

                           for (savingsInd=0; savingsInd<=jMax; savingsInd++)
                           {
                           // loop over decision variable, all possible savings levels tomorrow 
                           // compute consumption, utility(consumption)
                              cons          = findConsumption(allocfrac, x-aA[savingsInd]);
                              OOPMedex      = allocfrac*pow(cons,IESfrac); 
                              totalMedex    = OOPMedex/qRealization; 
                              todaysUtility = getUtility(cons, totalMedex, muRealization); 
                              if (tInd==(TDIMS-1))
                                 continuationUtilityAlive = 0;
                              else
                                 continuationUtilityAlive = GetContinuationUtility(savingsInd, tInd, gInd, hsInd, iInd, zetaInd);

                              value = todaysUtility + beta*survivalProbM[gInd][tInd][hsInd][iInd]*continuationUtilityAlive 
                                                    + beta*(1-survivalProbM[gInd][tInd][hsInd][iInd])*bequestUM[savingsInd];
                              if (value>maximum)
                              {
                                 maximum  = value;
                                 maxSavingsInd = savingsInd;
                                 maxCons  = cons;
                                 maxMedex = totalMedex;
                                 maxOOP   = OOPMedex;
                              }
                           } // end loop over decision of savings for tomorrow
                        } // end else branch

                     // Record value function for this period and policy function
                        valueFunM[gInd][tInd][iInd][hsInd][zetaInd][aInd][xiInd]   = maximum; 
                        consumFunM[gInd][tInd][iInd][hsInd][zetaInd][aInd][xiInd]  = maxCons; 
                        savingsFunM[gInd][tInd][iInd][hsInd][zetaInd][aInd][xiInd] = aA[maxSavingsInd];
                        outofpocketmedFunM[gInd][tInd][iInd][hsInd][zetaInd][aInd][xiInd] = maxOOP;
                        medFunM[gInd][tInd][iInd][hsInd][zetaInd][aInd][xiInd]     = maxMedex;                     
                     }//end loop over xi, transitory medical preference shock
                  }//end loop over zeta, persistent medical preference shock
               }//end loop over todays assets
            }//end loop over todays health status
         }//end loop over todays permanent income

      // MPI-specific Code:  Merge output from different processors 

         if (useMPI==2) //Why merge bequest and conumtion now? Each processor only needs the value function to continue.
         {
         // MPI_Allreduce and copyit take pointers for their first two args, so that's why valueFunM is missing an index.
/*
            MPI_Allreduce( (void*)&valueFunM[gInd][tInd][0][0][0][0][0], (void *)tempvec, recsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            copyit(valueFunM[gInd][tInd][0][0][0][0], tempvec,recsize);
            MPI_Allreduce( (void *)&consumFunM[gInd][tInd][0][0][0][0][0], (void *)tempvec, recsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            copyit( consumFunM[gInd][tInd][0][0][0][0], tempvec,recsize);
            MPI_Allreduce( (void *)&savingsFunM[gInd][tInd][0][0][0][0][0], (void *)tempvec, recsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            copyit( savingsFunM[gInd][tInd][0][0][0][0], tempvec,recsize);
            MPI_Allreduce( (void *)&medFunM[gInd][tInd][0][0][0][0][0], (void *)tempvec, recsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            copyit( medFunM[gInd][tInd][0][0][0][0], tempvec,recsize);
            MPI_Allreduce( (void *)&outofpocketmedFunM[gInd][tInd][0][0][0][0][0], (void *)tempvec, recsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            copyit( outofpocketmedFunM[gInd][tInd][0][0][0][0], tempvec,recsize);
*/
         }
      }//end loop over age
   }//end loop over gender

// printf("valueFunM_0s = %5.2lf\n", valueFunM[0][0][0][0][0][0][0]);

   free(tempvec);
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* Simulation of assets sequence 
   Subroutines needed:
      1 dimensional interpolation
      Discretization into markov chain
*/

//Note: The index matrices are still of type double; [for GAUSS interfacing]
void simulation(GMatrix zetacdfsim96Ptr,  GMatrix PIsim96Ptr,  GMatrix agesim96Ptr, GMatrix assetsim96Ptr,
                double **assetsimMat, // storing levels not indices; indices are done only internally
                double **cohsimMat, double **netIncomesimMat, double **consumptionsimMat, //empty on entry; for output to GAUSS
                double **healthindexsimhMat, double **healthindexsimwMat, double **marstatsimMat, //store indices
                double **xicdfsimMat, double **epsiloncdfsimMat,  //epsilon is analogous to decision variable; has one less draw than xi; we would like it to be numbered from 1 to tsimdims; C does start-at-0 indexing, so we call it "T+1" since the index is one behind where we want it.
                double **totalmedexsimMat, double **oopmedexsimMat, double **MedicaidsimMat, 
                double **grossbeqsimMat, double **zetaindexsimMat, double **zetasimMat, 
                double **xiindexsimMat, double **xisimMat, double **musimMat, double **qsimMat,
                double **marstatsim2Mat, double **transfersimMat, int iSimmin, int iSimmax)
{
   int personInd, yearInd, tInd, hsInd, zetaInd, xiInd, PIInd1,PIInd2, aIndLow, aIndHigh;
   double PIweightLow, AweightLow;    /* Contains weight used to interpolate on PI dimension. */
   struct result result1;             /* return of GetLocation, for interpolation*/

   double coh, coh_pre_transfer, zeta, assets, xi, epsilon, muRealization, qRealization, 
          allocfrac, xMin, assetsNextPeriod, cons, OOPMedex, trueoopmedex;
   double grossIncomePIInd1, grossIncomePIInd2,netIncomePIInd1,netIncomePIInd2;
            
   int sex, age96;
   double *tempvec;
   int recsize = nsims;
   
   tempvec = (double *)calloc(recsize,sizeof(double));

//* timing: at the begining of each period, know assets, shocks, marital status, age, health statues, PI*/
//      printf("switchMor: %f\n", switchMor);

   if (switchMor==0) /* no mortality risk*/
   {
      for(personInd=iSimmin; personInd<iSimmax; personInd++)
      {
         age96=(int)floor(rounder+agesim96Ptr.data[personInd]);
         for (yearInd=1; yearInd<TDimSims+1; yearInd++) 
         { 
            if ((age96+yearInd)<(TSTART+TDIMS)) /* Not too old for model */
            {
               marstatsimMat[yearInd][personInd]=marstatsimMat[yearInd-1][personInd]; 
            }/* Alive, for sure */
            else 
            {
               marstatsimMat[yearInd][personInd]=0; 
            }/* Too old for model, treat as dead */

         }
      }
   }
   printf("Rank=%d, iSimmin=%d iSimmax=%d \n", rank, iSimmin, iSimmax);

   for(personInd=iSimmin; personInd<iSimmax; personInd++)  
   {   /* if sample household is couple, skip */
      if ( (marstatsimMat[0][personInd]==3) || (marstatsimMat[0][personInd]==0)) 
      {
         for (yearInd=0; yearInd<TDimSims; yearInd++)
         {
            marstatsimMat[yearInd][personInd]     =  0;
            totalmedexsimMat[yearInd][personInd]  = -1e5;
            oopmedexsimMat[yearInd][personInd]    = -1e5;
            MedicaidsimMat[yearInd][personInd]    = -1e5; 
            transfersimMat[yearInd][personInd]    = -1e5;  
            zetaindexsimMat[yearInd][personInd]   = -1;   // Index number
            zetasimMat[yearInd][personInd]        = -1e5; // Value, from grid. 
            xiindexsimMat[yearInd][personInd]     = -1;   // Index number
            xisimMat[yearInd][personInd]          = -1e5; // Value, from grid. 
            musimMat[yearInd][personInd]          = -1;   // today's pref shock
            qsimMat[yearInd][personInd]           = -1;   // today's medexshare
            consumptionsimMat[yearInd][personInd] = -1e5;
            grossbeqsimMat[yearInd][personInd]    = -1e5;
            cohsimMat[yearInd][personInd]         = -1e5; // compute coh for time between today and tomorrow
            assetsimMat[yearInd][personInd]       = -1e5;
            netIncomesimMat[yearInd][personInd]   = -1e5;
         }
         continue;
      }

   /* Initialize zetaindexsimMat[0][personInd].  GAUSS provides a draw     
      from a U[0,1] distribution, which is combined with the invariant cdf 
      to produce an index on the discretized chain.  
      Also simulate initial draws of xiindexsim and xi.
    */
      zeta    = zetacdfsim96Ptr.data[personInd];
      zetaInd = GetLocation(&zetaInvcdf[0], zeta, ZETADIM+1).Ind1;
      zeta    = zetaA[zetaInd];  /* actual value of zeta */
      xi      = xicdfsimMat[0][personInd];
      xiInd   = GetLocation(&xiProbAcdf[0], xi, XIDIM+1).Ind1;
      xi      = xiA[xiInd];  /* actual value of xi */
      zetaindexsimMat[0][personInd] = zetaInd;
      zetasimMat[0][personInd] = zeta;
      xiindexsimMat[0][personInd] = xiInd;
      xisimMat[0][personInd] = xi;

      assetsimMat[0][personInd]=assetsim96Ptr.data[personInd];

   /* To interpolate along PI dimension, need to find indexes that this point 
      lies between and distance to each index. */
      result1     = GetLocation(IncomeA, PIsim96Ptr.data[personInd], IDIM); 
      PIInd1      = result1.Ind1;
      PIweightLow = result1.weight;
      PIInd2      = PIInd1+1;
      age96       = (int)floor(rounder+agesim96Ptr.data[personInd])-TSTART;

      for (yearInd=0; yearInd<TDimSims; yearInd++)  /*calendar year*/ //TdimSims loops yields TDimSims+1 asset variables [calculate for yearInd+1]
      {
 /* **************************************************************************************** */
      /* Decide whether it is a couple or single case by looking at marital status */
      /* in beq10.gau, 1-> male, 2->wife, (3->used to be couple, we should not get here if there are any) */
         if (marstatsimMat[yearInd][personInd]==0) //never check marstat[0][]; no one dies at time zero, so its always >0.
         {
        // household members all dead 
            grossbeqsimMat[yearInd][personInd]    =  assetsimMat[yearInd][personInd]; //0s; then a[t]; then -1e5s;
            marstatsimMat[yearInd][personInd]     =  0;
            totalmedexsimMat[yearInd][personInd]  = -1e5;
            oopmedexsimMat[yearInd][personInd]    = -1e5;
            MedicaidsimMat[yearInd][personInd]    = -1e5; 
            transfersimMat[yearInd][personInd]    = -1e5; 
            zetaindexsimMat[yearInd+1][personInd] = -1;   // Index number
            zetasimMat[yearInd+1][personInd]      = -1e5; // Value, from grid.
            xiindexsimMat[yearInd][personInd]     = -1;   // Index number
            xisimMat[yearInd][personInd]          = -1e5; // Value, from grid. 
            musimMat[yearInd][personInd]          = -1;   // todays pref shock
            qsimMat[yearInd][personInd]           = -1;   // todays medexshare
            consumptionsimMat[yearInd][personInd] = -1e5;
            cohsimMat[yearInd][personInd]         = -1e5; // compute coh for time between today and tomorrow
            assetsimMat[yearInd+1][personInd]     = -1e5;
            netIncomesimMat[yearInd][personInd]   = -1e5;
         }     
         else /* Single */
         {
            tInd=age96+yearInd; /* age*/
            if (marstatsimMat[yearInd][personInd]==1)
            {
               sex=0; /*husband */
               hsInd=(int)floor(rounder+healthindexsimhMat[yearInd][personInd]);                  
            } 

            if (marstatsimMat[yearInd][personInd]==2)
            { 
               sex=1; /*wife */
               hsInd=(int)floor(rounder+healthindexsimwMat[yearInd][personInd]);                  
            } 

         /* Zeta and xi are discretized into Markov chains.  xicdfsimMat 
              contains U[0,1] random variables that are combined with xi's  
              cdf to produce discrete-valued draws. 
            Write xi's index values into xiindexsimMat
            Write values into zetasimMat and xisimMat
         */
            zeta    = zetaA[zetaInd];  /* actual value of zeta */
            xi      = xicdfsimMat[yearInd][personInd];
            xiInd   = GetLocation(&xiProbAcdf[0], xi, XIDIM+1).Ind1;
            xi      = xiA[xiInd];  /* actual value of xi */
            zetasimMat[yearInd][personInd] = zeta;
            xiindexsimMat[yearInd][personInd] = xiInd;
            xisimMat[yearInd][personInd] = xi;

            assets     = assetsimMat[yearInd][personInd]; //grab asset values for beginning of period
            result1    = GetLocation(aA, assets, ADIM);
            aIndLow    = result1.Ind1;
            aIndHigh   = aIndLow+1;
            AweightLow = result1.weight;

         // Calculate net after tax income for realized interest rate
         // Need to interpolate over PI to get after-tax income
            grossIncomePIInd1 = mu_r*assetsimMat[yearInd][personInd] + yM[sex][PIInd1][tInd];
            grossIncomePIInd2 = mu_r*assetsimMat[yearInd][personInd] + yM[sex][PIInd2][tInd];
            netIncomePIInd1   = AfterTaxIncome(grossIncomePIInd1);
            netIncomePIInd2   = AfterTaxIncome(grossIncomePIInd2);
            netIncomesimMat[yearInd][personInd] = PIweightLow*netIncomePIInd1+(1-PIweightLow)*netIncomePIInd2;
            coh_pre_transfer  = assets+netIncomesimMat[yearInd][personInd];

         // We now have a shock-specific expenditure floor.
            muRealization     = muM[sex][tInd][PIInd1][hsInd][zetaInd][xiInd]*PIweightLow 
                                + muM[sex][tInd][PIInd2][hsInd][zetaInd][xiInd]*(1-PIweightLow);
            qRealization      = qM[sex][tInd][hsInd];        
            allocfrac         = qRealization*pow(muRealization/qRealization,1/omega);
            xMin              = getXMin(qRealization, muRealization);

            if (coh_pre_transfer < xMin)
               coh = xMin;
            else
               coh = coh_pre_transfer;

         // Calculate tomorrows decision for savings:
         // To get saving decision, interpolate over PI and assets 
            assetsNextPeriod = savingsFunM[sex][tInd][PIInd1][hsInd][zetaInd][aIndLow][xiInd]*PIweightLow*AweightLow
                               + savingsFunM[sex][tInd][PIInd2][hsInd][zetaInd][aIndLow][xiInd]*(1-PIweightLow)*AweightLow
                               + savingsFunM[sex][tInd][PIInd1][hsInd][zetaInd][aIndHigh][xiInd]*PIweightLow*(1-AweightLow)
                               + savingsFunM[sex][tInd][PIInd2][hsInd][zetaInd][aIndHigh][xiInd]*(1-PIweightLow)*(1-AweightLow);

         // Rather than interpolate highly non-linear decision rules for consumption, etc., redo calcs
         // Check asset accumulation equation to ensure internal consistency. 
            if (assetsNextPeriod > (coh-switchCMin*xMin)) assetsNextPeriod = (coh-switchCMin*xMin);  
            if (assetsNextPeriod < 0) assetsNextPeriod = 0;
            if (coh_pre_transfer < xMin) assetsNextPeriod = 0; // consumer forced to consume at least xMin. 

            cons     = findConsumption(allocfrac, coh-assetsNextPeriod); 
            OOPMedex = allocfrac*pow(cons,IESfrac);

            musimMat[yearInd][personInd]          = muRealization;
            qsimMat[yearInd][personInd]           = qRealization;
            cohsimMat[yearInd][personInd]         = coh;
            consumptionsimMat[yearInd][personInd] = cons;    
            totalmedexsimMat[yearInd][personInd]  = OOPMedex/qRealization; 
            assetsimMat[yearInd+1][personInd]     = assetsNextPeriod;
            transfersimMat[yearInd][personInd]    = coh - coh_pre_transfer;      

            if (coh_pre_transfer < xMin)
            {
            // find the OOP medex the consumer would make in the absence of any transfers
               trueoopmedex = 0;
               if (coh_pre_transfer>0)
               {
                  trueoopmedex = findConsumption(allocfrac, coh_pre_transfer);
                  trueoopmedex = allocfrac*pow(trueoopmedex,IESfrac);
               }
               MedicaidsimMat[yearInd][personInd] = OOPMedex - trueoopmedex;
               oopmedexsimMat[yearInd][personInd] = trueoopmedex;
            }
            else
            {
               MedicaidsimMat[yearInd][personInd] = 0;
               oopmedexsimMat[yearInd][personInd] = OOPMedex;
            }

         /* Calculate tomorrow's zeta (in parallel with calculating assets:
            Start with value for time t, calculate value for t+1
            zeta_t=rhoHc*zeta_t-1+epsilon_t;  
            Zeta is discretized into a Markov chain.  Epsiloncdfsim contains
              U[0,1] random variables that are combined with the zeta's 
              conditional cdf to produce discrete-valued draws. 
            Write zeta's index values into zetaindexsimMat 
         */
            epsilon = epsiloncdfsimMat[yearInd+1][personInd];
            zetaInd = GetLocation(&zetaProbMcdf[zetaInd][0], epsilon, ZETADIM+1).Ind1;
            zetaindexsimMat[yearInd+1][personInd] = zetaInd;

         }/* end else */

        marstatsim2Mat[yearInd][personInd] = marstatsimMat[yearInd][personInd]; 
      } /* end loop of year */

      marstatsim2Mat[TDimSims][personInd] = marstatsimMat[TDimSims][personInd]; 
   }  /* end loop of household */
   
   free(tempvec);
} //end simulation

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  Below is C/GAUSS I/O code written by Ken Housinger                       */
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  This function takes a char* and either reverses bits within bytes or bytes */
/*  themselves or both.  It returns a pointer to the original buffer which is  */
/*  altered. */

unsigned char * gread(unsigned char *inbuf, int bytes, int byte_rev, int bit_rev) 
   {
    unsigned char *tempbuf;
    unsigned char tempbyte, tempbit;
    int i, j;

    tempbuf = (unsigned char *) malloc(bytes);
    for (i = 0; i < bytes; i++) 
       {
        if (byte_rev) *(tempbuf + i) = *(inbuf + bytes - i - 1);
        else *(tempbuf + i) = *(inbuf + i);

        if (bit_rev) 
           {
            tempbyte = 0;
            for (j = 0; j < CHAR_BIT; j++) 
               {
                tempbit = *(tempbuf + i) >> (CHAR_BIT - j - 1);
                tempbit = tempbit << (CHAR_BIT - 1);
                tempbit = tempbit >> (CHAR_BIT - j - 1);
                tempbyte = tempbyte | tempbit;
               }
            *(tempbuf + i) = tempbyte;
           }
       }

    for (i = 0; i < bytes; i++)
        *(inbuf + i) = *(tempbuf + i);
    free(tempbuf);

    return(inbuf);
   }

//*--------------------------------------------------------------------------------*/
//*--------------------------------------------------------------------------------*/
//*  This function reads a Gauss v5.0 fmt file into a Matrix structure and  */
//*  returns the structure.  */

GMatrix gau5read(char *fmt) 
   {
//*  Initialize the matrix to be 1x1 and the byte/bit order to 0.  */
    GMatrix mat = {1, 1}; 
    unsigned int i;
    int type, byte, bit;
    unsigned char *cread;
    int bit_rev = 0, byte_rev = 0;
    FILE *fd;

    if (sizeof(int) != 4 || sizeof(double) != 8 || CHAR_BIT != 8) 
       {
        printf("Incompatable machine architecture.\n");
        return (mat);
       }

//*  Allocate enough space to store the header.  */
    cread = (unsigned char *) malloc(BASIC_HEADER_LEN); 
//*  Open *fmt for reading only.  */
    fd = fopen(fmt, "rb"); 
  
//*  Read the basic header (128 bytes) all at once.  */

    fread((void *) cread, 1, BASIC_HEADER_LEN, fd);
    byte = (int) *(cread + (BYTE_POS * sizeof(int)));  /* (0=Backward) */
    bit = (int) *(cread + (BIT_POS * sizeof(int)));    /* (0=Backward) */

//*  To get some system independence, we detect whether we have to reverse */
//*  the bytes or bits or both.  If x and x_SYSTEM match, no reverse is  */
//*  necessary. */

    if ((bit || BIT_SYSTEM) && !(bit && BIT_SYSTEM)) bit_rev=1;
    if ((byte || BYTE_SYSTEM) && !(byte && BYTE_SYSTEM)) byte_rev=1;

    type = *( (int *) gread((cread + (TYPE_POS * sizeof(int))), sizeof(int), 
                             byte_rev, bit_rev) );

//*  If the fmt file type is not a scalar, there are another two */
//*  ints of header giving the values of m and n.  If a matrix, also reset n. */

    if (type > SCALAR) 
       { 
        fread((void *) cread, 1, sizeof(int), fd);
        mat.m = *((unsigned int *) gread(cread, sizeof(int), byte_rev, bit_rev));
        fread((void *) cread, 1, sizeof(int), fd);
        if (type == MATRIX)
          mat.n = *((unsigned int *) gread(cread, sizeof(int), byte_rev, bit_rev));
      } 

//*  Allocate memory for the matrix.  The amount needed is m * n * sizeof(double). */
//*  Next, read in the data all at once.  Then use gread to reverse */
//*  bytes/bits if necessary. */

    free(cread);

    mat.data = (double *) malloc(mat.m * mat.n * sizeof(double));
    fread((void *) mat.data, sizeof(double), mat.m * mat.n, fd);
    if (byte_rev || bit_rev)
      for(i = 0; i < mat.m * mat.n; i++)
        gread((unsigned char *) mat.data + (i * sizeof(double)), sizeof(double), 
               byte_rev, bit_rev);

    fclose(fd);

    return (mat);
   }

//*--------------------------------------------------------------------------------*/
//*--------------------------------------------------------------------------------*/
//*  This function writes a Gauss v5.0 fmt file from a Matrix structure. */

void gau5write(char *fmt, GMatrix mat) 
{
//*  This ugly mess is the basic header. */

    unsigned int header[(BASIC_HEADER_LEN / sizeof(int)) + 2] = 
        {0xffffffff, 0, 0xffffffff, 0, 0xffffffff, 0, 0,
         0xabcdef01,1, 0, 1, 1008, sizeof(double), 0, 1,
         SCALAR, 1, 0, BASIC_HEADER_LEN};

    FILE *fd;

    if (sizeof(int) != 4 || sizeof(double) != 8 || CHAR_BIT != 8) 
       {
        printf("Incompatible machine architecture.\n");
        return;
       }

//*  If forward byte, make 6th int 0xffffffff. */
//*  If forward bit, make 7th int 0xffffffff. */

    if (BYTE_SYSTEM) header[BYTE_POS] = 0xffffffff;
    if (BIT_SYSTEM) header[BIT_POS] = 0xffffffff;

//*  If not a scalar, increase the 16th int by 1 and the 19th int (header */ 
//*  length) by 8 (2 * sizeof(int)).  Also, set m in int 33. */

    if (!(mat.m * mat.n == 1)) 
       {
        header[TYPE_POS] += 1;
        header[HEADER_LEN_POS] += (2 * sizeof(int));
        header[M_POS] = mat.m;

    /*  If not a vector (and not a scalar), increase the 16th int by 1 again */
    /*  and set m in int 34. */
        if (!(mat.n == 1)) 
           {
            header[TYPE_POS] += 1;
            header[N_POS] = mat.n;
           }
       }
  /*
  **Open fmt for writing and create if it does not exist.  If you create it,
  **make it a regular file with permissions 0640.  See comment in gau5read
  **for detail on how read (and similarly, write) work.
  **
  **Order: Write the basic header
  **       If not a scalar, write the other 2 ints of header 
  **       Write the m * n elements of data
  */

    fd = fopen(fmt, "wb"); 
    if ((mat.m * mat.n == 1))
      fwrite((void *) header, 1, BASIC_HEADER_LEN, fd);
    else
      fwrite((void *) header, 1, BASIC_HEADER_LEN + (2 * sizeof(int)), fd);
    fwrite((void *) mat.data, sizeof(double), mat.m * mat.n, fd);
    fclose(fd);
   }  
//*--------------------------------------------------------------------------------*/
//*--------------------------------------------------------------------------------*/
