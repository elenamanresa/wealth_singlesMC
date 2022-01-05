/*  WEALTH21.C
  - Finds value fn and policy fns for model in De Nardi, French and Jones (2006) 
  - Written by Fang Yang, with Cristina De Nardi 10/06 
    (cosmetic edits by John Jones 1/07) 
  - Doss: version with singles only 
  - GAUSS-C I/O procedures written by K. Housinger (2003)
  
  - I/O instructions
    ~ Put the (GAUSS) *.fmt files in a subdirectory called \iofiles\
    ~ create a folder called \data
    ~ create a folder called \output

  - We need to check if the difference in indexes exists also for singles (as it did 
    for couples)
   Note the different indexes for males and females: for a male age t=1, 
    means age=70,  for a female t=1 means age=67. Both 
    live up to age 100 (t=31 for males and 34 for females).

  - This version redefines the grid for cash on hand (cashA) and 
    consumption (consumA), putting more points near minimum level. 
    Also consumption (consumA) does not depend on cash on hand. 
    This could speed up searching.

  - There are two subroutines.
    ~ GetRulesSingle does the maximization for the single cases. 
    ~ GetValueSingle calculates value for each consumption choice and is used by GetRulesSingle. 

  - Medical expense shocks are discretized, both in the model and in the simulations.

  - Modifed by JBJ in July 08 to simulate discrete-valued transitory shocks.
    Lots of new output matrices created and saved
    New naming convention:
    ~ *cdfsim => simulations of U[0,1] draws
    ~ *indexim => simulations of index numbers
    ~ *sim => actual values

  - Future modifications: 
    ~ Solve policy function at T separately. This could save some time.
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <limits.h>

// Disable warning messages 4996 
#pragma warning( once : 4996 )

#define NUM_THREADS 52  // PC has 56 threads, but save a thread for other work
// #define NUM_THREADS 1  // PC has 56 threads, but save a thread for other work

#define useMPI 1         /* 0 => PC, 1 => OpenMP, 2 => cluster  */
// #define useMPI 0         /* 0 => PC, 1 => OpenMP, 2 => cluster  */


// MPI_COMMENT: 0 => single thread; 1 => openMP; 2 => MPI
// per Rory McGee in couples project

#define MPI_COMMENT 1

#if(MPI_COMMENT==0)
 #define useMPI 0 
 #define NUM_THREADS 1
#elif(MPI_COMMENT==1) 
 #include <omp.h>   
#elif(MPI_COMMENT==2)
 #include <mpi.h>
 #define useMPI 2    
 #define max(x,y) ((x) > (y) ? (x) : (y))
#endif(MPI_COMMENT==0)

#define switchJohn 1     /* 0: Cristina in Chicago, 1: John in Albany; 2 Charles Doss in Chicago */ 
#define PI 3.14159265
#define rounder 0.0001
#define ADDRESS_LEN 90
#define PRINTOUT 1

/* These constants are used to size various matrices */
#define TDIMS 33       /* dim. of age for singles  Programmer must ensure compatibility */
#define IDIM 5         /* dim. of permanent income, from lowest to highest */
#define HSDIM 2        /* dim. of health status: bad=0, good=1 */
#define RDIM 2         /* dim. of interest rate */
#define ZETADIM 8      /* dim. of persistent component of medical shock */
#define XIDIM 8        /* dim. of temporary component of medical shock */

/*  CASHDIM and CONSUMDIM should be chosen so that the grid of coh is a 
    subset of the grid for consumption.  Otherwise accidental  
    bequests might arise. */

#define CASHMAX1 200000         /* maximum cash on hand for segment 1 */ 
#define CASHMAX2 1000000        /* maximum cash on hand for segment 2 */
#define CASHMAX3 1600000        /* maximum cash on hand for segment 2 */
#define CONSUMPTIONMAX1 60000   /* maximum consumption for segment 1  */
#define CONSUMPTIONMAX2 200000  /* maximum consumption for segment 2  */
#define CONSUMPTIONMAX3 1600000 /* maximum consumption for segment 2  */

#define CASHDIM1 55   /* dim. of cash on hand for segment 1 :  up to 13/*00 if bequestdim=2 */
#define CASHDIM2 50   /* dim. of cash on hand:  up to 1300 if bequestdim=2 */
#define CASHDIM  112   /* dim. of cash on hand:  up to 1300 if bequestdim=2 */
#define CONSUMDIM1 55 /* dim. of consumption for segment 1:  up to 1300 if bequestdim=2 */
#define CONSUMDIM2 50 /* dim. of consumption for segment 1:  up to 1300 if bequestdim=2 */
#define CONSUMDIM  112 /* dim. of consumption:  up to 1300 if bequestdim=2 */

/*#define CASHDIM1 110   /* dim. of cash on hand for segment 1 :  up to 1300 if bequestdim=2 */
/*#define CASHDIM2 100   /* dim. of cash on hand:  up to 1300 if bequestdim=2 */
/*#define CASHDIM  235   /* dim. of cash on hand:  up to 1300 if bequestdim=2 */
/*#define CONSUMDIM1 110 /* dim. of consumption for segment 1:  up to 1300 if bequestdim=2 */
/*#define CONSUMDIM2 100 /* dim. of consumption for segment 1:  up to 1300 if bequestdim=2 */
/*#define CONSUMDIM  235 /* dim. of consumption:  up to 1300 if bequestdim=2 */

#define BEQUESTDIM 2   /* dim. of bequest to heirs */
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

/*
**This is my GMatrix opaque datatype.  The gauss fmt format has implied 
**dimensionality, but is just a list of doubles.  I could rewrite the
**functions to use pointers-to-pointers to improve the interface with
**Eric's functions, but that would add to the complexity of my code.  I
**suggest writing an intermediate function that turns *data into **data
**to cope with 2-dimensional arrays.
*/ 

typedef struct {   /* this defines a matrix with given # of rows and */
  unsigned int m;  /* columns and gives the address of the first element */
  unsigned int n;
  double *data;} GMatrix;

/*--------------------------------------------------------------------------------*/
/*----------------------Global PARAMETERS read in from GAUSS----------------------*/

/* Preference Parameters */
   double delta, beta, nu, phi0, K0; /* **** dropped eta, phi1 and phi2, k1 and 2 */

/* Simulation Parameters */
   int nsims, TDimSims, TSTART;  

/* Switches */
   int switchMor;    /* 0: no shocks, 1: shocks of mortality risk */
   int switchBeta;   /* 0: from GAUSS, 1: beta=1 */
   int switchY;      /* 0: noIncome 1: use income loaded */
   int switchHCost;  /* 0: no deterministic health cost, 1:  Health cost */
   int switchTax;    /* 0: noIncome, 1: use tax  loaded  */
   int switchZeta;   /* 0: no shocks, 1: shocks of health costs zeta,*/
   int switchXi;     /* 0: no shocks, 1: shocks of health costs zi */
   int switchR;      /* 0: no shocks, 1: shocks of interest rate,   */
   int switchBeq;    /* 0: no, 1: bequest motive  */
   int switchGender; /* 0: no, 1: gender diff */
   
   int rank;   /* Index number of this node  */
   int size;   /*  Number of nodes in cluster */

/* Miscellaneous */
   double cMin, tauBeq, exBeq, mu_r, sigma_r, medex_bottomcode; 

/* Non-asset income:  y_(t)  =  y(f,IDIM,TDIMS+1), a deterministic fn.
   INDEXES:  
      f: family structure
      f=0 single male
      f=1 single female
      IDIM: permanent income quantile, ordered from the lowest to highest 
      TDIMS: age, where +1 is just a place holder 
*/
   double yM[2][IDIM][TDIMS+1]; /* **** double yM[3][IDIM][TDIMS+1]; */

/* Health status uncertainty:
   Health status at time t, m_(t)(i), is a Markov process taking 
   two values, good and bad. The transition probabilities for health status 
   depend on current health status and age. The elements of the health status 
   transition matrix are pi_(kjt)(i) = Pr(m_(t+1)(i) =  j|m_(t)(i)  =  k),   
   k, j~(good, bad).
   Singles first:
   INDEXES: 
      SEX: 0=husband, 1=wife
      HSDIM: health state today 0=bad, 1=good
      HSDIM: health state tomorrow 0=bad, 1=good 
*/
      double hsProbM[2][TDIMS][IDIM][HSDIM][HSDIM];  /* Singles */


/* Survival uncertainty:
   s_(m, I, t)(i) denotes the probability of individual i being alive at age t
   conditional on current health status, permanent income, and being alive at age t-1. 
   0=male, 1=female 
*/
   double survivalProbM[2][TDIMS][HSDIM][IDIM]; 

/* Health Costs
   ln hc_(t)  =  hc(f_(t), m_(I,t)(h), m_(I,t)(w), t, I)+
                 sigma(f_(t), m_(I,t)(h), m_(I,t)(w), t, I) * psi_(t).
   psi_(t)    =  zeta_(t)+xi_(t),   xi_(t) ~ N(0, sigma_(xi)^2),  
   zeta_(t)   =  rho_(hc)zeta_(t-1) + epsilon_(t),   epsilon_(t) ~ N(0, sigma_(epsilon)^2) 
*/
   double rhoHc, sigma_xi, sigma_epsilon; 

/* health cost functions */
   double hcSingleLogMean[2][HSDIM][TDIMS+1][IDIM]; 
   double hcSingleSigma[2][HSDIM][TDIMS+1][IDIM]; 

/* Income tax structure, from French's code
      taxBrk gives tax brackets
      taxMar gives marginal tax rates  
*/
   double taxBrk[taxDim-1] = {6250, 40200, 68400, 93950, 148250, 284700}; 
   double taxMar[taxDim] = {0.0765, 0.2616, 0.4119, 0.3499, 0.3834, 0.4360, 0.4761}; 

/* Strings of directory names */
   char rootdir[ADDRESS_LEN]; 
   char outputdir[ADDRESS_LEN]; 
   char datadir[ADDRESS_LEN]; 
   
/*----------------Global VARIABLES determined inside this program-----------------*/

double consumA[CONSUMDIM];   /* consumption grid */
double cashA[CASHDIM];       /* cash on hand grid */
double hsA[HSDIM];           /* health status grid */
double IncomeA[IDIM];        /* Permanent income grid, income as percentage ranking. */
                             /* (See subroutine Grid for detail.) */

/* Used to store markov chains */
double xiA[XIDIM], xiProbA[XIDIM];  /* temporary medical shock */
double xiProbAcdf[XIDIM+1];         
double tempInv[XIDIM+RDIM];         /* Placeholder */
double zetaA[ZETADIM];              /* persistent medical shock */
double zetaProbM[ZETADIM][ZETADIM];
double zetaInv[ZETADIM];            /* stationary dist.  */
double zetaInvcdf[ZETADIM+1]; 
double rA[RDIM], rProbA[RDIM];      /* interest rate shock */
int zetaIndN; 
double zetaProbMcdf[ZETADIM][ZETADIM+1]; 

/* income after tax at brackets  */
double incomeBrk[taxDim-1];   

/* utility matrix  */
double utilityM[HSDIM][CONSUMDIM];  /* single */

/* expected utility from bequest */
double bequestUM[CASHDIM][CONSUMDIM];  /* single */

/* value function matrix */
double valueFunM[2][TDIMS+1][IDIM][HSDIM][ZETADIM][CASHDIM];  

/* consumption policy function matrix */
double consumFunM[2][TDIMS][IDIM][HSDIM][ZETADIM][CASHDIM]; 

/* bequest policy function matrix*/
double bequestFunSM[2][TDIMS][IDIM][HSDIM][ZETADIM][CASHDIM];  /* singles */

struct result
{ int Ind1;
double weight;
};

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* functions used in main*/

/*  Converts vectors to matrices and vice-versa*/
GMatrix zeromat(unsigned int recsize);
void switchem(double **dataPtr, double **tempPtr);
double **SetUpSim(double *dataVec, int extrayears);

/*  GAUSS-C++ I/O programs written by K. Housinger */
unsigned char * gread(unsigned char *inbuf, int bytes, int byte_reverse, int bit_reverse);
GMatrix gau5read(char *fmt);
void gau5write(char *fmt, GMatrix mat); /* reads vector from hard drive*/

/*  GAUSS-C++ I/O programs added by Elena Manresa */
GMatrix getcsvdat(char *csv, int nRows, int nCols);
void writecsvdat(char *csv, GMatrix mat);





/*   Update global parameters with values passed in from GAUSS*/
int globderef(double *prefvecPtr, double *asstvecPtr, double *hcostvecPtr, 
              double *simvecPtr, double *agevecPtr, double *switchvecPtr);
int GetProfiles(double *agevecPtr, double *mnlnhcSPtr, double *stdlnhcSPtr,  
                double *hcostPIPtr, double *mortrateSPtr,  
                double *mortratePIPtr, double *hstranSPtr, 
                double *hstranPIPtr, double *yprofPtr, double *yprofPIPtr);
double LogitSQRT(double x);
double Logit(double x);
void TwoYearToOneYear(double gg, double bb, double *_g, double *_b); 

/* create grid for cash on hand and consumption   */ 
void Grid(double cashA[], double consumA[CONSUMDIM], double hsA[HSDIM], double IncomeA[IDIM]);

/*discretize AR(1)and iid processes into Markov chain */
void Discretization(int n,double rho, double mu, double sigma, double zArray[], 
                    double *piMatrixP, double *piInvarV); 

/* define Utility from consumption matrix   */
void GetUtility(double utilityM[HSDIM][CONSUMDIM]);

/* define Utility from bequest matrix   */
void GetUtilityBeq(double bequestUM[][CONSUMDIM]);

/* after-tax income at bracket points  */
void IncomeAtBrk(double taxBrk[], double taxMar[], double incomeBrk[]);

/* Interpolation and extrapolation */
double Interpolation(double *fP,  double *xP,  double x, int DIM);
int Locate(double *Xarray, double x, int DIM);
struct result GetLocation(double *xP,  double x, int DIM);

void WriteData(double **cohsimMat, double **netIncomesimMat, double **consumptionsimMat, 
               double **healthcostsimMat, double **zetaindexsimMat, double **marstatsimMat, 
               GMatrix agesim96Ptr, GMatrix PIsim96Ptr, double **healthsimhMat, 
               double **healthsimwMat);
int *getAssignmentVec(int numPoints, int numNodes);

void GetRulesSingle(int iAssetsmin, int iAssetsmax, int *assetAssignments);
double GetValueSingle(int i,int tInd, int IInd, int hsInd, int zetaInd, 
                 int cashInd, int consumInd);

void simulation(GMatrix zetacdfsim96Ptr, GMatrix PIsim96Ptr, GMatrix agesim96Ptr, 
                double **cohsimMat, double **assetsimMat,  double **netIncomesimMat, 
                double **consumptionsimMat, double **healthsimhMat, double **healthsimwMat, 
                double **marstatsimMat, double **xicdfsimMat, double **epsiloncdfsimMat, 
                double **healthcostsimMat, double **zetaindexsimMat, double **zetasimMat, 
                double **xiindexsimMat, double **xisimMat, double **MedicaidsimMat, 
                double *rorsim, double **beqsimMat, double **marstatsim2Mat, 
                double **transfersimMat, int iSimsmin, int iSimsmax); 
            
void GetCdf(int nRows, int nCols, double *piMatrixP, double *piMatrixCDF); 

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* main program, note that everything before this is global */
void main(int argc, char *argv[])
{
   GMatrix prefvecPtr, asstvecPtr, agevecPtr, simvecPtr, switchvecPtr, hcostvecPtr, 
           mnlnhcSPtr, stdlnhcSPtr, hcostPIPtr, mortrateSPtr, mortratePIPtr, 
           hstranSPtr, hstranPIPtr, yprofPtr, yprofPIPtr, rorsimPtr, cohsim96Ptr, 
           zetacdfsim96Ptr, PIsim96Ptr, agesim96Ptr, healthsimhPtr, healthsimwPtr, 
           marstatsimPtr, xicdfsimPtr, epsiloncdfsimPtr, cohsimPtr, assetsimPtr, 
           healthcostsimPtr, zetasimPtr, zetaindexsimPtr, xisimPtr, xiindexsimPtr, 
           MedicaidsimPtr, consumptionsimPtr, beqsimPtr, netIncomesimPtr, 
           marstatsim2Ptr, transfersimPtr;

   double **cohsimMat, **assetsimMat, **netIncomesimMat, **consumptionsimMat, 
          **healthsimhMat, **healthsimwMat, **marstatsimMat, **xicdfsimMat, 
          **epsiloncdfsimMat, **healthcostsimMat, **zetasimMat, **zetaindexsimMat, 
          **xisimMat, **xiindexsimMat, **MedicaidsimMat, **beqsimMat, 
          **marstatsim2Mat, **transfersimMat;

   double *tempvec;

   int gotderef, gotprofiles, yearInd, recsize;
   char fullpath[ADDRESS_LEN];

   double *rorsim; /* define ror sequence in main. remember to change it.  */

   clock_t start, end;  /* recode time */

   int personInd, rInd, zetaInd, xiInd;  /*loop indices */
   int *assetAssignments;
   int simspermachine, iAssetsmin, iAssetsmax, iSimsmin, iSimsmax, statespermachine;
   double spm2;

   FILE *parameterP; /*point to file containing parameters */

   if (useMPI<2)  /*  PC */
   {
      if (switchJohn==1)  /* John is using in Albany */
      {
         /*strcpy(rootdir,"e:\\Users\\John\\wealth_singles\\iofiles\\");
         strcpy(outputdir, "e:\\Users\\John\\wealth_singles\\output\\"); 
         strcpy(datadir, "e:\\Users\\John\\wealth_singles\\data\\"); */
	  
		 strcpy(rootdir, "C:\\wealth_singles - MC\\iofiles\\");
		 strcpy(outputdir, "C:\\wealth_singles - MC\\output\\");
		 strcpy(datadir, "C:\\wealth_singles - MC\\data\\");
      }
      else
      {
         strcpy(rootdir,"c:\\cristina\\eric\\wealth_singles\\iofiles\\");
         strcpy(outputdir, "c:\\cristina\\eric\\wealth_singles\\output\\");
         strcpy(datadir, "c:\\cristina\\eric\\wealth_singles\\data\\");
      }
   }
   else  /*  Unix or Cluster*/
   {
      if (switchJohn==1)
      {
         strcpy(rootdir,"/home/jjones/wealth_singles/job_110308/iofiles/"); 
         strcpy(outputdir,"/home/jjones/wealth_singles/job_110308/output/"); 
         strcpy(datadir,"/home/jjones/wealth_singles/data/"); 
      }
      else if (switchJohn == 2)
      {
         strcpy(rootdir,"/home/cdoss/wealth_singles/iofiles/"); 
         strcpy(outputdir,"/home/cdoss/wealth_singles/output/"); 
         strcpy(datadir,"/home/cdoss/wealth_singles/data/"); 
      }
   }

   start = clock();
   rank  = 1; /* Global */
   size  = 1;

#if(MPI_COMMENT==1)
 omp_set_num_threads(NUM_THREADS); // Get number of processors. 
 #pragma omp parallel 
 {   size = omp_get_num_threads(); }
 printf("Number of threads = %d\n", size);
#elif(MPI_COMMENT==2)
 MPI_Init(&argc, &argv);
 MPI_Comm_size(MPI_COMM_WORLD, &size);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 printf("size=%3d rank=%3d\n", size, rank);
#endif(MPI_COMMENT==1)

/* Read in parameter vectors from the *.fmt files (read into memory as arrays) */
  // prefvecPtr    = gau5read(strcat(strcpy(fullpath,rootdir),"prefvec.fmt"));  
  // prefvecPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "prefvec.csv"),6,1);
  // asstvecPtr    = gau5read(strcat(strcpy(fullpath,rootdir),"asstvec.fmt"));
  // agevecPtr     = gau5read(strcat(strcpy(fullpath,rootdir),"agevec.fmt")); 
  // switchvecPtr  = gau5read(strcat(strcpy(fullpath,rootdir),"swchvec.fmt")); 
  // simvecPtr     = gau5read(strcat(strcpy(fullpath,rootdir),"simvec.fmt")); 
  // hcostvecPtr   = gau5read(strcat(strcpy(fullpath,rootdir),"medexvec.fmt"));   
  // mnlnhcSPtr    = gau5read(strcat(strcpy(fullpath,rootdir),"mnlnmxs.fmt"));   
  // stdlnhcSPtr   = gau5read(strcat(strcpy(fullpath,rootdir),"stdlnmxs.fmt")); 
  // hcostPIPtr    = gau5read(strcat(strcpy(fullpath,rootdir),"mxpicoef.fmt"));
  // mortrateSPtr  = gau5read(strcat(strcpy(fullpath,rootdir),"mortprfs.fmt"));   
  // mortratePIPtr = gau5read(strcat(strcpy(fullpath,rootdir),"mort_pi.fmt"));
  // hstranSPtr    = gau5read(strcat(strcpy(fullpath,rootdir),"hsprobs.fmt")); 
  // hstranPIPtr   = gau5read(strcat(strcpy(fullpath,rootdir),"heal_pi.fmt"));
  // yprofPtr      = gau5read(strcat(strcpy(fullpath,rootdir),"yprof.fmt")); 
  // yprofPIPtr    = gau5read(strcat(strcpy(fullpath,rootdir),"y_pi.fmt"));

   prefvecPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "prefvec.csv"), 5, 1);
   asstvecPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "asstvec.csv"),5,1);
   agevecPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "agevec.csv"),3,1);
   switchvecPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "swchvec.csv"),10,1);
   simvecPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "simvec.csv"),4,1);
   hcostvecPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "medexvec.csv"),5,1);
   mnlnhcSPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "mnlnmxs.csv"),132,1);
   stdlnhcSPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "stdlnmxs.csv"),132,1);
   hcostPIPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "mxpicoef.csv"),264,1);
   mortrateSPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "mortprfs.csv"),132,1);
   mortratePIPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "mort_pi.csv"),132,1);
   hstranSPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "hsprobs.csv"),132,1);
   hstranPIPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "heal_pi.csv"),132,1);
   yprofPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "yprof.csv"),66,1);
   yprofPIPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "y_pi.csv"),132,1);



/* Initialize parameters by assigning array elements to appropriate parameters */
   gotderef     = globderef(prefvecPtr.data,asstvecPtr.data,hcostvecPtr.data,
                            simvecPtr.data, agevecPtr.data, switchvecPtr.data);

/* Create grids for cash on hand, consumption, health status and income   */
   Grid(cashA, consumA, hsA, IncomeA);

/* Initialize profiles */
   gotprofiles = GetProfiles(agevecPtr.data, mnlnhcSPtr.data,  
                             stdlnhcSPtr.data, hcostPIPtr.data, 
                             mortrateSPtr.data, mortratePIPtr.data, 
                             hstranSPtr.data, hstranPIPtr.data, 
                             yprofPtr.data, yprofPIPtr.data);

/* Save dimensions to ASCII files so that matlab can transform output files into matrices */
   if ( (useMPI<2)||(rank==(size-1)) )
   {
      parameterP=fopen(strcat(strcpy(fullpath,outputdir),"parameter.txt"),"w" );
      fprintf(parameterP,"%5d\n %5d\n %5d\n %5d\n %5d\n %5d\n %5d\n %5d\n %5d\n %5d\n %5d\n %10d\n %10d\n %5.2lf\n %5.2lf\n %5.2f\n %5.2f\n %5.2f\n %5.2f\n",
              IDIM, TDIMS, HSDIM, RDIM, ZETADIM, XIDIM, CASHDIM1, CASHDIM ,CONSUMDIM1, CONSUMDIM ,BEQUESTDIM,CASHMAX1,CASHMAX2, cMin, delta, beta, nu, phi0, K0);
     /* **** note I dropped some arguments above */
     /* Items are
      - dimension of permanent income
      - dimension of age  
      - dimension of health status: good = 1,  bad = 0
      - dimension of interest rate  
      - dimension of persistent component of medical shock 
      - dimension of temporary component of medical shock 
      - dimension of cash on hand 
      - dimension of consumption 
      - maximum cash on hand 
      - dimension of bequest to the child  
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
   /* printf(" discretized persistent component zeta\n"); */
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

/* On transitory component of health costs */
   if (switchXi==1)
   {
   /* printf(" discretized transitory component xi\n"); */
      Discretization(XIDIM,0.0,0.0,sigma_xi, xiA, &xiProbA[0], &tempInv[0]); 
   }
   else /* no shocks */
   {   
      for (xiInd=0; xiInd<XIDIM; xiInd++)
      {
         xiA[xiInd]=0.0;
      /* xiA[xiInd]=1.0;  changed 08/31/05*/
         xiProbA[xiInd]=1.0/XIDIM;
      }
   }

   GetCdf(1, XIDIM, &xiProbA[0], &xiProbAcdf[0]);  /* get CDF of xi */


   if (switchR==1) /* on interest rate  */
   {
   /* printf(" discretized interest rate r\n"); */
      Discretization(RDIM,0.0,mu_r,sigma_r, rA, &rProbA[0], &tempInv[XIDIM]);  
   }   
   else
   {
      for (rInd=0; rInd<RDIM; rInd++)
      {
         rA[rInd]=mu_r;
         rProbA[rInd]=1.0/RDIM;
      }
   }


/* initialize matrices */
   GetUtility(utilityM);     /* Find utility from consumption matrix */
   GetUtilityBeq(bequestUM); /* Find utility from bequest matrix */
   IncomeAtBrk(taxBrk, taxMar, incomeBrk); /* after-tax income at bracket points  */
   
   iAssetsmin = 0;
   iAssetsmax = CASHDIM;
   iSimsmin = 0;
   iSimsmax = nsims;
   assetAssignments = getAssignmentVec(CASHDIM, size);   // Improves efficiency when parallelizing.  JBJ:  5/30/14 *?#?* Added RM 08/17

   if (useMPI==1)
   {
      spm2 = ((double)CASHDIM) / ((double)size);
      statespermachine = (int)ceil(spm2);
      spm2 = ((double)nsims) / ((double)size);
      simspermachine = (int)ceil(spm2);
   }
   else if (useMPI==2)
   {
      spm2 = ((double) CASHDIM)/ ((double) size);
      statespermachine = (int) ceil(spm2);
      iAssetsmin = rank*statespermachine;
      iAssetsmax = (rank+1)*statespermachine;
      if (iAssetsmax > CASHDIM) iAssetsmax = CASHDIM;

      spm2 = ((double) nsims)/ ((double) size);
      simspermachine = (int) ceil(spm2);
      iSimsmin = rank*simspermachine;
      iSimsmax = (rank+1)*simspermachine;
      if (iSimsmax > nsims) iSimsmax = nsims;
   }

#if(MPI_COMMENT==1)
 #pragma omp parallel private(iAssetsmin,iAssetsmax,rank)
 {
   rank = omp_get_thread_num();
   iAssetsmin = rank * statespermachine;
   iAssetsmax = (rank + 1)*statespermachine;
   if (iAssetsmax > CASHDIM) iAssetsmax = CASHDIM;

// *?#? removed precomputation of single-->death bequest utility and instead compute at each t in getrulessingle
   GetRulesSingle(iAssetsmin, iAssetsmax, assetAssignments); //*?#?* scrambling of asset grid in parallel RM 08/17
   }
#else
 GetRulesSingle(iAssetsmin,iAssetsmax,assetAssignments);
#endif
 
   if ((useMPI<2) || (rank == (size - 1)))printf("Finished with decision rules for singles. rank=%d \n ", rank);
   fflush(stdout);

/* load simulated inputs */
/* initial year */
//   PIsim96Ptr    = gau5read(strcat(strcpy(fullpath,rootdir),"pisim96.fmt"));  /* simulated permanent income in 1996 */
 //  agesim96Ptr   = gau5read(strcat(strcpy(fullpath,rootdir),"agesim96.fmt")); /* simulated age in 1996 */
 //  cohsim96Ptr   = gau5read(strcat(strcpy(fullpath,rootdir),"cohsim96.fmt")); /* simulated cash-on-hand in 1996 */  
/* read in whole history */  
 //  healthsimhPtr = gau5read(strcat(strcpy(fullpath,rootdir),"healsimh.fmt")); /* simulation of health status for males */
 //  healthsimwPtr = gau5read(strcat(strcpy(fullpath,rootdir),"healsimw.fmt")); /* simulation of health status for females */
 //  marstatsimPtr = gau5read(strcat(strcpy(fullpath,rootdir),"mstatsim.fmt")); /* simulation of marital status */
 //  zetacdfsim96Ptr  = gau5read(strcat(strcpy(fullpath,rootdir),"ztacdfsim96.fmt")); /* simulated zeta in 1996 */
 //  xicdfsimPtr      = gau5read(strcat(strcpy(fullpath,rootdir),"xicdfsim.fmt"));    /* simulation of transitory part of health cost */
 //  epsiloncdfsimPtr = gau5read(strcat(strcpy(fullpath,rootdir),"epscdfsim.fmt"));   /* simulation of epsilon from uniform distribution */


   /* initial year */
   PIsim96Ptr = getcsvdat(strcat(strcpy(fullpath, rootdir), "pisim96.csv"),150000,1);  /* simulated permanent income in 1996 */
   agesim96Ptr = getcsvdat(strcat(strcpy(fullpath, rootdir), "agesim96.csv"),150000,1); /* simulated age in 1996 */
   cohsim96Ptr = getcsvdat(strcat(strcpy(fullpath, rootdir), "cohsim96.csv"),150000,1); /* simulated cash-on-hand in 1996 */
/* read in whole history */
   healthsimhPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "healsimh.csv"), 1800000,1); /* simulation of health status for males */
   healthsimwPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "healsimw.csv"), 1800000,1); /* simulation of health status for females */
   marstatsimPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "mstatsim.csv"), 1800000,1); /* simulation of marital status */
   zetacdfsim96Ptr = getcsvdat(strcat(strcpy(fullpath, rootdir), "ztacdfsim96.csv"),150000,1); /* simulated zeta in 1996 */
   xicdfsimPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "xicdfsim.csv"), 1800000,1);    /* simulation of transitory part of health cost */
   epsiloncdfsimPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "epscdfsim.csv"), 1800000,1);   /* simulation of epsilon from uniform distribution */



/* Reshape array to matrix */
   healthsimhMat = SetUpSim(healthsimhPtr.data,1); /* fill in all history */
   healthsimwMat = SetUpSim(healthsimwPtr.data,1); /* fill in all history */
   marstatsimMat = SetUpSim(marstatsimPtr.data,1); /* fill in all history */
   xicdfsimMat   = SetUpSim(xicdfsimPtr.data,1);      /* transitory shock fill in entire history.  0->1 */
   epsiloncdfsimMat = SetUpSim(epsiloncdfsimPtr.data,1); /* shocks fed to the markov process, for entire history. 0->1 */

/* Set up matrices for storing simulation results */
   recsize           = (TDimSims+1)*nsims;
   beqsimPtr         = zeromat(recsize);
   cohsimPtr         = zeromat(recsize);
   healthcostsimPtr  = zeromat(recsize);
   zetaindexsimPtr   = zeromat(recsize); 
   zetasimPtr        = zeromat(recsize); 
   xiindexsimPtr     = zeromat(recsize); 
   xisimPtr          = zeromat(recsize); 
   MedicaidsimPtr    = zeromat(recsize); 
   consumptionsimPtr = zeromat(recsize);
   assetsimPtr       = zeromat(recsize);
   netIncomesimPtr   = zeromat(recsize);
   marstatsim2Ptr    = zeromat(recsize); 
   transfersimPtr    = zeromat(recsize); 

   tempvec = (double *)calloc(recsize,sizeof(double));

/* Reshape array to matrix */
   beqsimMat         = SetUpSim(beqsimPtr.data,1);  
   healthcostsimMat  = SetUpSim(healthcostsimPtr.data,1);
   zetaindexsimMat   = SetUpSim(zetaindexsimPtr.data,1);
   zetasimMat        = SetUpSim(zetasimPtr.data,1);     
   xiindexsimMat     = SetUpSim(xiindexsimPtr.data,1);  
   xisimMat          = SetUpSim(xisimPtr.data,1);       
   MedicaidsimMat    = SetUpSim(MedicaidsimPtr.data,1); 
   consumptionsimMat = SetUpSim(consumptionsimPtr.data,1);
   cohsimMat         = SetUpSim(cohsimPtr.data,1);
   assetsimMat       = SetUpSim(assetsimPtr.data,1);
   netIncomesimMat   = SetUpSim(netIncomesimPtr.data,1);
   marstatsim2Mat    = SetUpSim(marstatsim2Ptr.data,1);  
   transfersimMat    = SetUpSim(transfersimPtr.data,1);  


/* initialize first column of cohsimMat for initial condition */
   for(personInd=0; personInd<nsims; personInd++) 
      cohsimMat[0][personInd]=cohsim96Ptr.data[personInd];

   // rorsimPtr = gau5read(strcat(strcpy(fullpath,rootdir),"rorshk.fmt"));
   rorsimPtr = getcsvdat(strcat(strcpy(fullpath, rootdir), "rorshk.csv"),12,1);

   rorsim    = rorsimPtr.data;
   for (yearInd=0; yearInd<TDimSims+1; yearInd++)
   {
      if (switchR==0) rorsim[yearInd]=mu_r;
      else rorsim[yearInd] += mu_r;
   }

   printf("simulation starts\n");

#if(MPI_COMMENT==1) //Use openMP
 #pragma omp parallel private(iSimsmin,iSimsmax,rank)
 {
    rank = omp_get_thread_num();
    iSimsmin = rank * simspermachine;
    iSimsmax = (rank + 1)*simspermachine;
    if (iSimsmax > nsims) iSimsmax = nsims;
    /*printf("simulation starts for rank=%d\n", rank);*/

    simulation(zetacdfsim96Ptr, PIsim96Ptr, agesim96Ptr, cohsimMat, assetsimMat,
               netIncomesimMat, consumptionsimMat, healthsimhMat, healthsimwMat,
               marstatsimMat, xicdfsimMat, epsiloncdfsimMat, healthcostsimMat,
               zetaindexsimMat, zetasimMat, xiindexsimMat, xisimMat, MedicaidsimMat,
               rorsim, beqsimMat, marstatsim2Mat, transfersimMat, iSimsmin, iSimsmax);

    /*printf("simulation ends for rank=%d\n", rank);*/
    #pragma omp barrier //wait until every processor has hit this point
 } //On openMP pragma
#else //not using openmp
   simulation(zetacdfsim96Ptr, PIsim96Ptr, agesim96Ptr, cohsimMat, assetsimMat,
              netIncomesimMat, consumptionsimMat, healthsimhMat, healthsimwMat, 
              marstatsimMat, xicdfsimMat, epsiloncdfsimMat, healthcostsimMat, 
              zetaindexsimMat, zetasimMat, xiindexsimMat, xisimMat, MedicaidsimMat, 
              rorsim, beqsimMat, marstatsim2Mat, transfersimMat, iSimsmin, iSimsmax); 
#endif

#if(MPI_COMMENT==2)
   if (useMPI==2)
   {
      MPI_Reduce((void *) cohsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&cohsimPtr.data, &tempvec);
   
      MPI_Reduce((void *) assetsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&assetsimPtr.data, &tempvec);

      MPI_Reduce((void *) netIncomesimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&netIncomesimPtr.data, &tempvec);

      MPI_Reduce((void *) zetaindexsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&zetaindexsimPtr.data, &tempvec);

      MPI_Reduce((void *) zetasimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&zetasimPtr.data, &tempvec);

      MPI_Reduce((void *) xiindexsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&xiindexsimPtr.data, &tempvec);

      MPI_Reduce((void *) xisimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&xisimPtr.data, &tempvec);

      MPI_Reduce((void *) MedicaidsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&MedicaidsimPtr.data, &tempvec);

      MPI_Reduce((void *) healthcostsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&healthcostsimPtr.data, &tempvec);

      MPI_Reduce((void *) consumptionsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&consumptionsimPtr.data, &tempvec);

      MPI_Reduce((void *) beqsimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&beqsimPtr.data, &tempvec);

      MPI_Reduce((void *) marstatsim2Ptr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&marstatsim2Ptr.data, &tempvec);  

      MPI_Reduce((void *) transfersimPtr.data, (void *) tempvec, recsize, MPI_DOUBLE,
                 MPI_SUM, size - 1, MPI_COMM_WORLD);
      switchem(&transfersimPtr.data, &tempvec);  
   }
#endif

   free(tempvec);

   if ( (useMPI<2)||(rank==(size-1)) )
   {
      if (PRINTOUT==1)
	  {
         WriteData(cohsimMat, netIncomesimMat, consumptionsimMat, healthcostsimMat, zetaindexsimMat, 
                   marstatsimMat, agesim96Ptr, PIsim96Ptr, healthsimhMat, healthsimwMat );
	  }
	  
    //  gau5write(strcat(strcpy(fullpath,rootdir),"cohsim.fmt"), cohsimPtr); 
	//  gau5write(strcat(strcpy(fullpath, rootdir), "ztasim.fmt"), zetasimPtr);
	//  gau5write(strcat(strcpy(fullpath,rootdir),"asstsim.fmt"), assetsimPtr);
	//  gau5write(strcat(strcpy(fullpath,rootdir),"netIncomesim.fmt"), netIncomesimPtr); 
	//  gau5write(strcat(strcpy(fullpath,rootdir),"ztasim.fmt"), zetasimPtr)//;
	//  gau5write(strcat(strcpy(fullpath,rootdir),"ztaindexsim.fmt"), zetaindexsimPtr);
	//  gau5write(strcat(strcpy(fullpath,rootdir),"xisim.fmt"), xisimPtr);
	//  gau5write(strcat(strcpy(fullpath,rootdir),"xiindexsim.fmt"), xiindexsimPtr);
	//  gau5write(strcat(strcpy(fullpath,rootdir),"Medicaidsim.fmt"), MedicaidsimPtr); 
	//  gau5write(strcat(strcpy(fullpath,rootdir),"medexsim.fmt"), healthcostsimPtr); 
	//  gau5write(strcat(strcpy(fullpath,rootdir),"conssim.fmt"), consumptionsimPtr); 
	//  gau5write(strcat(strcpy(fullpath,rootdir),"beqsim.fmt"), beqsimPtr); 
	//  gau5write(strcat(strcpy(fullpath,rootdir),"mssim2.fmt"), marstatsim2Ptr); 
	//  gau5write(strcat(strcpy(fullpath,rootdir),"transfersim.fmt"), transfersimPtr); 

	  writecsvdat(strcat(strcpy(fullpath, rootdir), "cohsim.csv"), cohsimPtr);
	  writecsvdat(strcat(strcpy(fullpath, rootdir), "asstsim.csv"), assetsimPtr);
	  writecsvdat(strcat(strcpy(fullpath, rootdir), "netIncomesim.csv"), netIncomesimPtr);
	  writecsvdat(strcat(strcpy(fullpath, rootdir), "ztasim.csv"), zetasimPtr);
	  writecsvdat(strcat(strcpy(fullpath, rootdir), "ztaindexsim.csv"), zetaindexsimPtr);
	  writecsvdat(strcat(strcpy(fullpath, rootdir), "xisim.csv"), xisimPtr);
	  writecsvdat(strcat(strcpy(fullpath, rootdir), "xiindexsim.csv"), xiindexsimPtr);
	  writecsvdat(strcat(strcpy(fullpath, rootdir), "Medicaidsim.csv"), MedicaidsimPtr);
	  writecsvdat(strcat(strcpy(fullpath, rootdir), "medexsim.csv"), healthcostsimPtr);
	  writecsvdat(strcat(strcpy(fullpath, rootdir), "conssim.csv"), consumptionsimPtr);
	  writecsvdat(strcat(strcpy(fullpath, rootdir), "beqsim.csv"), beqsimPtr);
	  writecsvdat(strcat(strcpy(fullpath, rootdir), "mssim2.csv"), marstatsim2Ptr);
	  writecsvdat(strcat(strcpy(fullpath, rootdir), "transfersim.csv"), transfersimPtr);

      end = clock();
      printf("This program ends in %5d minutes %5d seconds \n ",(end-start)/CLOCKS_PER_SEC/60,
             ((end-start)/CLOCKS_PER_SEC)%60);
   }

#if(MPI_COMMENT==2)
   MPI_Finalize();
#endif  
}   /* End of main*/

/*--------------------------------------------------------------------------------*/
/*----------------------------------SUBROUTINES-----------------------------------*/

/* create grid for cash on hand and consumption   */ 

void Grid(double cashA[], double consumA[CONSUMDIM], double hsA[HSDIM], double IncomeA[IDIM])
{
   int cashInd, consumInd,IInd, hsInd;
   FILE *cashP; 
   char fullpath[ADDRESS_LEN];
   double cMin2, cMin3;

/* create cash grid, cMin<=x<=CASHMAX */

   for (cashInd=0; cashInd<CASHDIM1; cashInd++)
   { 
      cashA[cashInd]=pow(sqrt(cMin)+(sqrt(CASHMAX1)-sqrt(cMin))*cashInd/(CASHDIM1-1),2);
   }

   cMin2=cashA[CASHDIM1-1];

   for (cashInd=0; cashInd<CASHDIM2; cashInd++)
   { 
      cashA[CASHDIM1+cashInd]=pow(sqrt(cMin2)+(sqrt(CASHMAX2)-sqrt(cMin2))*(cashInd+1)/(CASHDIM2),2);
   }

   cMin3=cashA[CASHDIM1+CASHDIM2-1];

   for (cashInd=1; cashInd<CASHDIM-CASHDIM1-CASHDIM2+1; cashInd++)
   {
      cashA[cashInd+CASHDIM1+CASHDIM2-1]
         = pow(sqrt(cMin3)+(sqrt(CASHMAX3)-sqrt(cMin3))*cashInd/(CASHDIM-CASHDIM1-CASHDIM2),2);
   }

   for (consumInd=0; consumInd<CONSUMDIM1; consumInd++)
   { 
      consumA[consumInd]=pow(sqrt(cMin)+(sqrt(CONSUMPTIONMAX1)-sqrt(cMin))*consumInd/(CONSUMDIM1-1),2);
   }

   cMin2=consumA[CONSUMDIM1-1];

   for (consumInd=0; consumInd<CONSUMDIM2; consumInd++)
   { 
      consumA[CONSUMDIM1+consumInd]
         = pow(sqrt(cMin2)+(sqrt(CONSUMPTIONMAX2)-sqrt(cMin2))*(consumInd+1)/(CONSUMDIM2),2);
   }

   cMin3=consumA[CONSUMDIM1+CONSUMDIM2-1];

   for (consumInd=1; consumInd<CONSUMDIM-CONSUMDIM1-CONSUMDIM2+1; consumInd++)
   {
      consumA[consumInd+CONSUMDIM1+CONSUMDIM2-1]
         = pow(sqrt(cMin3)+(sqrt(CONSUMPTIONMAX3)-sqrt(cMin3))*consumInd/(CONSUMDIM-CONSUMDIM1-CONSUMDIM2),2);
   }

   for (IInd = 0; IInd<IDIM; IInd++)
   {
      IncomeA[IInd]= ((double) IInd)/(IDIM-1); /* Income expressed as percentile ranking */
   }
      
   for ( hsInd = 0; hsInd<HSDIM; hsInd++)
   {
      hsA[hsInd]=hsInd;
   }

   if ( (useMPI<2)||(rank==(size-1)) )
   {
      cashP=fopen(strcat(strcpy(fullpath,outputdir),"xA.txt"),"w");

      for (cashInd=0; cashInd<CASHDIM1; cashInd++)
      { 
         fprintf(cashP,"%10.3f\n", cashA[cashInd]);
      }

      cMin2=cashA[CASHDIM1-1];

      for (cashInd=0; cashInd<CASHDIM2; cashInd++)
      { 
         fprintf(cashP,"%10.3f\n", cashA[CASHDIM1+cashInd]);
      }

      cMin3=cashA[CASHDIM1+CASHDIM2-1];

      for (cashInd=1; cashInd<CASHDIM-CASHDIM1-CASHDIM2+1; cashInd++)
      {
         fprintf(cashP,"%10.3f\n", cashA[cashInd+CASHDIM1+CASHDIM2-1]);
      }
  
      fclose(cashP);
      cashP=fopen(strcat(strcpy(fullpath,outputdir),"consumA.txt"),"w");

      for (consumInd=0; consumInd<CONSUMDIM1; consumInd++)
      { 
         fprintf(cashP,"%10.3f\n", consumA[consumInd]);
      }

      for (consumInd=0; consumInd<CONSUMDIM2; consumInd++)
      { 
         fprintf(cashP,"%10.3f\n", consumA[CONSUMDIM1+consumInd]);
      }

      for (consumInd=1; consumInd<CONSUMDIM-CONSUMDIM1-CONSUMDIM2+1; consumInd++)
      {
         fprintf(cashP,"%10.3f\n", consumA[consumInd+CONSUMDIM1+CONSUMDIM2-1]);
      }

      fclose(cashP);
   }
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* Records nodes and weights for Gauss-Hermite Quadrature
      n: number of Gauss points to be used
      xArray[j]: j-th point
      wArray[j]: weight of the j-th point
      Nodes and weights are from Judd(1999), page 266 
*/
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
    
   if ((sigma==0)||(rho==1)) /* not random: added 09/15.04 */
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
   else if ((!((n==1)||(n==2)||(n==3)||(n==4)||(n==5)||(n==6)||(n==7)||(n==8)||(n==10)))||((rho>=1)||(rho<0))||(sigma<0))
   {
      printf("Error! The requirement for the AR(1) process is that the number of Guass points n=1, 2,3,4,5,7,10, the autoregession coefficiency 0<=rho<1 and the standard deviation sigma is positive.\n");     
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
/* Defines utility matrix for singles */

void GetUtility(double utilityM[HSDIM][CONSUMDIM])
{
   int hsInd=0;
   int consumInd=0; /*  index */

   if (nu  ==  1.0)     /* log utility */
   {
      for (hsInd = 0; hsInd<HSDIM; hsInd++)
      {
         for (consumInd = 0; consumInd<CONSUMDIM; consumInd++)
         {
            utilityM[hsInd][consumInd] = (1+delta*hsA[hsInd])*log(consumA[consumInd]); 
         }
      }
   }
   else 
   {
      for (hsInd = 0; hsInd<HSDIM; hsInd++)
      {
         for (consumInd = 0; consumInd<CONSUMDIM; consumInd++)
         {
            utilityM[hsInd][consumInd] = (1+delta*hsA[hsInd])*pow(consumA[consumInd], 1-nu)/(1-nu);  
			
         }
      }
   }
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/

/*  net bequest function  */
double NetBequest(double *amountBequestedP)
{
   double netBequest = 0; 
   if (*amountBequestedP>exBeq)
      netBequest = exBeq+(1-tauBeq)*(*amountBequestedP-exBeq); 
   else netBequest = *amountBequestedP; 
   return netBequest; 
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* utility from bequest for a single person
   phi_j(b_net) = phi_j*( (b_net+K_j)^(1-nu) )/(1-nu) */

double UBeqSingle(double *netBequestP)
{
   double utils;    
   if (nu  ==  1)   /* log utility */
      utils  =  phi0*log(*netBequestP+K0); 
   else
      utils  =  phi0*pow(*netBequestP+K0, 1-nu)/(1-nu);         
   return utils; 
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* Expected utility from leaving bequest matrix for a single, not discounted */

void GetUtilityBeq(double bequestUM[][CONSUMDIM])
{
   int  cashInd, consumInd; /*  index */
   double temp, temp2;      
   double *tempP = &temp;      
   double *tempP2 = &temp2;  

   for (cashInd = 0; cashInd<CASHDIM; cashInd++)
   {
      for (consumInd = 0; consumInd<CONSUMDIM; consumInd++)
      {  
         temp = (cashA[cashInd]-consumA[consumInd]); /*gross bequest */
         if (temp>=0)
         {
            temp2 = NetBequest(tempP);  /*net bequest */
            bequestUM[cashInd][consumInd] = UBeqSingle(tempP2);                         
         }
         else break;      
      }
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
/*  Locates nearest point in an array
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
      while (dif>1)
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
   return (x-*(xP+j))*(*(fP+j-1)-*(fP+j))/(*(xP+j-1)-*(xP+j))+(*(fP+j));  /*x[j-1]<x<x[j] */
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

void WriteData(double **cohsimMat,double **netIncomesimMat, double **consumptionsimMat, 
               double **healthcostsimMat, double **zetaindexsimMat, double **marstatsimMat, 
               GMatrix agesim96Ptr, GMatrix PIsim96Ptr, double **healthsimhMat,
               double **healthsimwMat)
{
   int i, tInd, IInd, hsInd, hsIndN; 
   
   FILE *fp;
   char fullpath[ADDRESS_LEN];
   int yearInd, personInd;

   FILE *valueFP;  /*point to files related to singles case */
   FILE *consumFP;
   FILE *beqFP;  
   int zetaInd,  cashInd;

   fp = fopen(strcat(strcpy(fullpath,outputdir),"imcprof1.txt"), "w");
  
   for (i = 0; i<2; i++) /* for (i = 0; i<3; i++) ***** */
   {
      for (IInd = 0; IInd<IDIM; IInd++)
      {
         for (tInd = 0; tInd<TDIMS+1; tInd++)
         {   
            fprintf(fp,"%20.8lf",yM[i][IInd][tInd]);               
         }
         fprintf(fp,"\n");
      }
   }    

   fclose(fp);
 
/*-------------------------------------------------------------------------------*/

   fp = fopen(strcat(strcpy(fullpath,outputdir),"survivalPro1.txt"), "w"); /* open file*/

   for (i = 0; i<2; i++)
   {
      for ( hsInd = 0; hsInd<HSDIM; hsInd++)
      {
         for (IInd = 0; IInd<IDIM; IInd++)
         {
            for (tInd = 0; tInd<TDIMS; tInd++)
            {
               fprintf(fp, "%20.8lf",survivalProbM[i][tInd][hsInd][IInd]);
            }
            fprintf(fp,"\n");
         }
      }
   } 
   fclose(fp);

/*--------------------------------------------------------------------------------*/
/* hsProbM[2][TDIMS][IDIM][HSDIM][HSDIM]; */

   fp = fopen(strcat(strcpy(fullpath,outputdir),"profheal1.txt"), "w"); 

   for (i = 0; i<2; i++)
   {
      for ( hsInd = 0; hsInd<HSDIM; hsInd++)
      {
         for ( hsIndN = 0; hsIndN<HSDIM; hsIndN++)
         {
            for (IInd = 0; IInd<IDIM; IInd++)
            {
               for (tInd = 0; tInd<TDIMS; tInd++)
               {
                  fprintf(fp,"%20.8lf",hsProbM[i][tInd][IInd][hsInd][hsIndN]);
               }
            fprintf(fp,"\n");
            }
         }
      } 
   }

   fclose(fp);

/*--------------------------------------------------------------------------------*/
/* health status transition matrix singles */   

   fp = fopen(strcat(strcpy(fullpath,outputdir),"hcSingleLogMean1.txt"), "w");
   for (i = 0; i<2; i++)
   {
      for ( hsInd = 0; hsInd<HSDIM; hsInd++)
      {
         for (IInd = 0; IInd<IDIM; IInd++)
         {
            for (tInd = 0; tInd<TDIMS+1; tInd++)
            {
               fprintf(fp, "%20.8lf", hcSingleLogMean[i][hsInd][tInd][IInd]);
            }
            fprintf(fp,"\n");
         }
      }
   }
   fclose(fp);

   fp = fopen(strcat(strcpy(fullpath,outputdir),"hcSingleSigma1.txt"), "w");
   for (i = 0; i<2; i++)
   {
      for ( hsInd = 0; hsInd<HSDIM; hsInd++)
      {
         for (IInd = 0; IInd<IDIM; IInd++)
         {            
            for (tInd = 0; tInd<TDIMS+1; tInd++) 
            {
               fprintf(fp, "%20.8lf", hcSingleSigma[i][hsInd][tInd][IInd]);
            }
            fprintf(fp,"\n");
         }
      }
   }
   fclose(fp);

/* value and policy functions */
/*--------------------------------------------------------------------------------*/

   valueFP  = fopen(strcat(strcpy(fullpath,outputdir),"valueF.txt"),"w");  /* value function */
   consumFP = fopen(strcat(strcpy(fullpath,outputdir),"consumptionF.txt"),"w");  /* consumption policy function */
   beqFP    = fopen(strcat(strcpy(fullpath,outputdir),"bequestSF.txt"),"w");  /* bequest policy function */

/* Print value function for t = 1, ..., T  where value function for T+1 is defined as 0 matrix */
/* single case */
   for (i = 0; i<2; i++)
   {
      for (tInd = (TDIMS-1); tInd>= 0; tInd--)
      {
         for (IInd = 0; IInd<IDIM; IInd++)
         {
            for ( hsInd = 0; hsInd<HSDIM; hsInd++)
            {
               for (zetaInd = 0; zetaInd<ZETADIM; zetaInd++)
               {
                  for (cashInd = 0; cashInd<CASHDIM; cashInd++)   
                  {   
                  /* write value function and policy function to files */
                     fprintf(valueFP, "%12.4f  ", valueFunM[i][tInd][IInd][hsInd][zetaInd][cashInd]);
                     fprintf(consumFP, "%12.4f  ", consumFunM[i][tInd][IInd][hsInd][zetaInd][cashInd]);
                     fprintf(beqFP, "%12.4f  ", bequestFunSM[i][tInd][IInd][hsInd][zetaInd][cashInd]);
                  }
                  
                  fprintf(valueFP, "\n");
                  fprintf(consumFP, "\n");
                  fprintf(beqFP, "\n");
               }

               fprintf(valueFP, "\n");
               fprintf(consumFP, "\n");
               fprintf(beqFP, "\n");
            }/* end loop through health status */
         } /* end loop through permanent income */   
      } /* end loop through age */
   } /* end loop through sex */

/* End writing to files  */
   fclose(valueFP);
   fclose(consumFP);
   fclose(beqFP);

/*--------------------------------------------------------------------------------*/
/* simulation results */

   fp = fopen(strcat(strcpy(fullpath,outputdir),"zetaindsimMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims; yearInd++)  /*calendar year */
      {
         fprintf(fp, "%10.1lf", zetaindexsimMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);


   fp = fopen(strcat(strcpy(fullpath,outputdir),"cohsimMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims+1; yearInd++)  /*calendar year*/
      {
         fprintf(fp, "%20.3lf", cohsimMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);

   fp = fopen(strcat(strcpy(fullpath,outputdir),"consumptionsimMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {

      for (yearInd=0; yearInd<TDimSims+1; yearInd++)  /*calendar year */
      {
         fprintf(fp, "%20.3lf", consumptionsimMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);

   fp = fopen(strcat(strcpy(fullpath,outputdir),"healthcostsimMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims+1; yearInd++)  /*calendar year */
      {
         fprintf(fp, "%20.5lf", healthcostsimMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);


   fp = fopen(strcat(strcpy(fullpath,outputdir),"marstatsimMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims; yearInd++)  
      {
         fprintf(fp, "%20.5lf", marstatsimMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);


   fp = fopen(strcat(strcpy(fullpath,outputdir),"PI93.txt"), "w");

   for(personInd=0; personInd<nsims; personInd++)
   {
      fprintf(fp, "%10.5lf", PIsim96Ptr.data[personInd] );
   }
   fclose(fp);

   fp = fopen(strcat(strcpy(fullpath,outputdir),"age96.txt"), "w");

   for(personInd=0; personInd<nsims; personInd++)
   {
      fprintf(fp, "%10.5lf", agesim96Ptr.data[personInd] );
   }
   fclose(fp);


   fp = fopen(strcat(strcpy(fullpath,outputdir),"healthsimhMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims; yearInd++)  /*calendar year */
      {
         fprintf(fp, "%20.5lf", healthsimhMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);


   fp = fopen(strcat(strcpy(fullpath,outputdir),"healthsimwMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims; yearInd++)  /*calendar year */
      {
         fprintf(fp, "%20.5lf", healthsimwMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);


   fp = fopen(strcat(strcpy(fullpath,outputdir),"PIMat.txt"), "w");
   for(personInd=0; personInd<nsims; personInd++)
   {
      for (yearInd=0; yearInd<TDimSims+1; yearInd++)  /*calendar year */
      {
         fprintf(fp, "%20.5lf", netIncomesimMat[yearInd][personInd]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);

}


/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
int GetProfiles(double *agevecPtr, double *mnlnhcSPtr, double *stdlnhcSPtr,  
                double *hcostPIPtr, double *mortrateSPtr, double *mortratePIPtr,
                double *hstranSPtr, double *hstranPIPtr, double *yprofPtr, 
                double *yprofPIPtr)
{
   int i, j, k, varshift, MSInd, sexInd, tInd, IInd, hsInd, tdimGAUSS; 

   double income, survprobS, PI1, PI2,   
          PIC1, PIC2, PIC11, PIC12, PIC21, PIC22, PIC31, PIC32, PIC41, PIC42,
          gg, bb, _g, _b, bhprob, hcInt, hcSig;
           
   tdimGAUSS = (int) floor(rounder+agevecPtr[2]);

   if (tdimGAUSS != TDIMS)
   {
      printf("Incompatible timespans!!! \n");      
   }

/* Income */
   j=0;
   i=0;
   for (tInd = 0; tInd<TDIMS; tInd++)
   {
      for (MSInd = 0; MSInd<2; MSInd++) /* Single male, single female */  
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
  
/* Health Status transition probabilities for singles */
/* Health expense parameters for singles */
/* Survival Probabilities for singles */

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

         PIC31 = hcostPIPtr[k];
         PIC32 = hcostPIPtr[k+1];

         PIC41 = hcostPIPtr[k+varshift];
         PIC42 = hcostPIPtr[k+1+varshift];

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
               hcInt     = mnlnhcSPtr[i+hsInd];
               hcSig     = stdlnhcSPtr[i+hsInd];

               survivalProbM[sexInd][tInd][hsInd][IInd] = 
                  LogitSQRT(survprobS + PIC11*PI1 + PIC12*PI2);
               
               if (hsInd==0) bhprob = _b;
               else bhprob = 1-_g;

               hsProbM[sexInd][tInd][IInd][hsInd][0] = bhprob;   /* Bad health at t+1 */
               hsProbM[sexInd][tInd][IInd][hsInd][1] = 1-bhprob; /* Good health at t+1 */

               hcSingleLogMean[sexInd][hsInd][tInd][IInd] = hcInt + PIC31*PI1 + PIC32*PI2;
               hcSingleSigma[sexInd][hsInd][tInd][IInd]   = sqrt(pow(hcSig,2)+PIC41*PI1+PIC42*PI2); 
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
    
/* Fill in placeholder values for period T+1 */

   for (IInd = 0; IInd<IDIM; IInd++)
   {
      for (MSInd = 0; MSInd<2; MSInd++)   /* Income */ 
      {
         yM[MSInd][IInd][TDIMS] = 1;
      }

      for ( hsInd = 0; hsInd<HSDIM; hsInd++)   /* Health cost parameters for singles */
      {
         for (sexInd = 0; sexInd<2; sexInd++)
         {
            hcSingleLogMean[sexInd][hsInd][TDIMS][IInd] = 0;
            hcSingleSigma[sexInd][hsInd][TDIMS][IInd] = 0; 
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

   if (switchHCost==0) 
   {
   /* Health expense parameters */
   /* Survival Probabilities */

      for (tInd = 0; tInd<TDIMS; tInd++)
      {
         for (sexInd = 0; sexInd<2; sexInd++)
         {
            for ( hsInd = 0; hsInd<HSDIM; hsInd++)  /* Bad health, then good*/
            {
               for (IInd = 0; IInd<IDIM; IInd++)
               {
                  hcSingleLogMean[sexInd][hsInd][tInd][IInd] = 0.0;
                  hcSingleSigma[sexInd][hsInd][tInd][IInd] = 0.0; 
               }
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
}

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

int globderef(double *prefvecPtr, double *asstvecPtr, double *hcostvecPtr, 
           double *simvecPtr, double *agevecPtr, double *switchvecPtr)
{
/* Flow utility parameters:
   beta = discount factor
   d(m(i)) = 1+d*1{m(i) = good}
   u(c, m(i)) = d(m(i))*(c^(1-nu))/(1-nu)
*/
   delta  = prefvecPtr[0];

   beta   = prefvecPtr[1];
   nu     = prefvecPtr[2];

   cMin   = asstvecPtr[0]; /* Minimun consumption provided by the government */
   tauBeq = asstvecPtr[1]; /* Estate tax rate */
   exBeq  = asstvecPtr[2]; /* Estate tax exemption level */

/* interest rate r is an i.i.d. random variable, with mean mu_r and variance (sigma_r)^2. */
   mu_r          = asstvecPtr[3];
   sigma_r       = asstvecPtr[4]; 

   rhoHc         = hcostvecPtr[0];
   sigma_epsilon = sqrt(hcostvecPtr[2]);
   sigma_xi      = sqrt(hcostvecPtr[3]);
   medex_bottomcode = hcostvecPtr[4];

   nsims         = (int) floor(rounder+simvecPtr[0]);
   TDimSims      = (int) floor(rounder+simvecPtr[1]);

   TSTART        = (int) floor(rounder+agevecPtr[0]);

   switchMor     = (int) floor(rounder+switchvecPtr[0]);
   switchBeta    = (int) floor(rounder+switchvecPtr[1]);
   switchY       = (int) floor(rounder+switchvecPtr[2]);
   switchHCost   = (int) floor(rounder+switchvecPtr[3]);
   switchTax     = (int) floor(rounder+switchvecPtr[4]);
   switchZeta    = (int) floor(rounder+switchvecPtr[5]);
   switchXi      = (int) floor(rounder+switchvecPtr[6]);
   switchR       = (int) floor(rounder+switchvecPtr[7]);
   switchBeq     = (int) floor(rounder+switchvecPtr[8]);
   switchGender  = (int) floor(rounder+switchvecPtr[9]);
 
   if (switchBeta==1) beta=1;

/* Bequest utility parameters for singles
   Phi_j(b_net) = phi_j*[(b_net+K_j)^(1-nu)]/[1-nu] */

   phi0 = prefvecPtr[3]*((double) switchBeq);
   K0   = prefvecPtr[4];
   
   return 1;
} 
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  GETASSIGNMENTVEC:  Allocate points on a state vector across nodes             */
/*                     Uses Sangeeta Pratap's allocation trick                    */

int *getAssignmentVec(int numPoints, int numNodes)
{
   int numBlocks, iBlock, iNode, iPoint, thisPoint;
   double nb2;
   int *assignmentVec;

   assignmentVec = (int *)calloc(numPoints, sizeof(int));

   // for (iPoint=0; iPoint<numPoints; iPoint++)
   //    assignmentVec[iPoint] = iPoint;

   nb2 = ((double)numPoints) / ((double)numNodes);
   numBlocks = (int)ceil(nb2);
   iPoint = 0;

   // Now scramble, by pulling numbers from different "blocks" 
   for (iNode = 0; iNode<numNodes; iNode++)
   {
      for (iBlock = 0; iBlock<numBlocks; iBlock++)
      {
         thisPoint = iBlock * numNodes + iNode + 1;
         if (thisPoint>numPoints) { continue; }   // Recall that indexing starts at zero
         assignmentVec[iPoint] = thisPoint - 1;
         iPoint += 1;
      }
   }

   return assignmentVec;
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
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
void copyit(double *permvec, double *tempvec, int recsize)
{
   int i;
   for(i=0;i<recsize;i++) permvec[i]=tempvec[i];
   return;
}
/*--------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
/* Solve value function for t = 1, ..., T   
   Singles case 
*/
void GetRulesSingle(int iAssetsmin,int iAssetsmax, int *assetAssignments)
{
   int i, tInd, IInd, hsInd,  zetaInd,  cashInd2, cashInd, consumInd;
   
   double maximum = 0.0;     /* used in maximization to record maximum*/
   int index = 0;            /* used in maximization to record index of maximizer */
   double valueM[CONSUMDIM]; /* intermediate value function, used in maximization */

   double *tempvec;
   int recsize = IDIM*HSDIM*ZETADIM*CASHDIM;
   tempvec = (double *)calloc(recsize,sizeof(double));

   for (i = 0; i<2; i++)
   {
      if  (i == 0)
      {
         if ( (useMPI<2)||(rank==(size-1)) )
            /*printf("Solving policy functions. Case of a single male \n ")*/; 
      }
      else 
      {
         if ( (useMPI<2)||(rank==(size-1)) )
            /*printf("Solving policy functions. Case of a single female \n ")*/; 
      }

      for (tInd = (TDIMS-1); tInd>= 0; tInd--)
      {
         for (IInd = 0; IInd<IDIM; IInd++)
         {
            for (hsInd = 0; hsInd<HSDIM; hsInd++)
            {
               for (zetaInd = 0; zetaInd<ZETADIM; zetaInd++)
               {
                  for (cashInd2 = iAssetsmin; cashInd2<iAssetsmax; cashInd2++) 
                  {   
                  /* maximization on each node */
                     cashInd = assetAssignments[cashInd2];  // "Scramble" gridpoints across nodes/processors
                     maximum = -1000000.0;  /* maximum utility */
                     index   = 0;   /* Index of consumption giving maximum utility */
                     
                 /* loop through possible consumption choice */
                     for ( consumInd = 0; consumInd<CONSUMDIM; consumInd++)
                     { 
                        if (consumA[consumInd]<=cashA[cashInd])
                        {    
                           valueM[consumInd]=GetValueSingle(i, tInd, IInd, hsInd, zetaInd, cashInd, consumInd);

/*                         if ((tInd==28)&(hsInd==0)&(IInd==0))   
                           {
                              fprintf(sumP,"%10.5f", sum);
                              fprintf(sumP,"%10.5f %10.5f\n",utilityM[hsInd][consumInd],valueM[consumInd] );
                           }
*/
                        /* maximization over all possible consumption choices */
                           if (maximum<valueM[consumInd])
                           {
                              maximum = valueM[consumInd]; 
                              index = consumInd; 
                           } /*maximization */
                        }
                        else  /* consumption choice is bigger than available cash on hand */
                        {
                           break;
                        }

                     }/* end loop through consumption grid */
 
                  /* record value function for this period and policy function */
                     valueFunM[i][tInd][IInd][hsInd][zetaInd][cashInd]    = maximum; 
                     consumFunM[i][tInd][IInd][hsInd][zetaInd][cashInd]   = consumA[index]; 
                     bequestFunSM[i][tInd][IInd][hsInd][zetaInd][cashInd] = cashA[cashInd] -consumA[index]; 

                  }/* end loop through cash on hand  */      
               }/* end loop through persistent medical shock */   

            }/* end loop through health status */
         }/* end loop through permanent income */   

      /* Merge output from different processors */

#if(MPI_COMMENT==1)
 #pragma omp barrier //wait until every processor has hit this point
#elif(MPI_COMMENT==2)
         if (useMPI==2)  // Update value function across all nodes
         {
            MPI_Allreduce( valueFunM[i][tInd][0][0][0], tempvec, recsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            copyit(valueFunM[i][tInd][0][0][0], tempvec,recsize); 
            MPI_Allreduce(bequestFunSM[i][tInd][0][0][0], tempvec, recsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
             copyit(bequestFunSM[i][tInd][0][0][0] , tempvec,recsize);
            MPI_Allreduce( consumFunM[i][tInd][0][0][0], tempvec, recsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            copyit( consumFunM[i][tInd][0][0][0], tempvec,recsize);
         }
#endif

       } /* end loop through age */
   } /* end loop through sex */
   free(tempvec);
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* Calculate value function at each choice of consumption.
   Singles case 
*/
double GetValueSingle(int i,int tInd, int  IInd, int hsInd, int zetaInd, 
                      int cashInd, int consumInd)
{
   int  rInd, xiInd, zetaIndN, hsIndN; 
   double cashOnHand; /* cash on hand next period. used in maximization */

   double sum; /* sum is the expected value next period */   
   sum = 0.0; 
   
/* Loop through possible next period states */   
   
   for (rInd = 0; rInd<RDIM; rInd++)
   {
      for (hsIndN = 0; hsIndN<HSDIM; hsIndN++)
      {
         for (zetaIndN = 0; zetaIndN<ZETADIM; zetaIndN++)
         {
            for (xiInd = 0; xiInd<XIDIM; xiInd++)   
            {
               cashOnHand = max(cashA[cashInd]-consumA[consumInd]
                           +AfterTaxIncome(rA[rInd]*(cashA[cashInd]-consumA[consumInd])
                           +(yM[i][IInd][tInd+1]))-exp((hcSingleLogMean[i][hsIndN][tInd+1][IInd]
                           +hcSingleSigma[i][hsIndN][tInd+1][IInd]*(zetaA[zetaIndN]+xiA[xiInd]))), cMin); 
               sum+= hsProbM[i][tInd][IInd][hsInd][hsIndN]*zetaProbM[zetaInd][zetaIndN]*rProbA[rInd]*xiProbA[xiInd]*Interpolation(&valueFunM[i][tInd+1][IInd][hsIndN][zetaIndN][0], cashA, cashOnHand, CASHDIM); 
            }
         }
      }
   } /* end looping through possible next period states */

/* valueM stores the utility for each possible choice of consumption */   
   
   return utilityM[hsInd][consumInd]+beta*(1-survivalProbM[i][tInd][hsInd][IInd])*bequestUM[cashInd][consumInd]
         + beta*survivalProbM[i][tInd][hsInd][IInd]*sum; 
}


/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* Simulation of assets sequence 
   Subroutines needed:
      1 dimensional interpolation
      Discretization into markov chain
*/

void simulation(GMatrix zetacdfsim96Ptr, GMatrix PIsim96Ptr, GMatrix agesim96Ptr, 
                double **cohsimMat, double **assetsimMat,  double **netIncomesimMat, 
                double **consumptionsimMat, double **healthsimhMat, double **healthsimwMat, 
                double **marstatsimMat, double **xicdfsimMat, double **epsiloncdfsimMat, 
                double **healthcostsimMat, double **zetaindexsimMat, double **zetasimMat, 
                double **xiindexsimMat, double **xisimMat, double **MedicaidsimMat, 
                double *rorsim, double **beqsimMat, double **marstatsim2Mat, 
                double **transfersimMat, int iSimsmin, int iSimsmax) 

{
   int personInd, yearInd, tInd, zetaInd_before, hsInd, zetaInd, xiInd, 
       PIInd1, PIInd2;
   double weight;          /* Contains weight used to interpolate on PI dimension. */
   struct result result1;  /* return of GetLocation, for interpolation*/

   double coh, coh_pre_transfer, assets, consumption, grossIncome, netIncome, 
          xi, zeta, epsilon;
   double healthCost, Medicaid;      /* health costs, not their logs.  */
   int sex, age96;

   double *tempvec, *tempV;
   
   int recsize = nsims;

   tempvec = (double *)calloc(recsize,sizeof(double));
   tempV = (double *)calloc(recsize,sizeof(double));

/* timing: at the begining of each period, know coh, marital status, age, health statues, PI, zeta*/
   if (switchMor==0) /* no mortality risk*/
   {
      for(personInd=iSimsmin; personInd<iSimsmax; personInd++)
      {
         age96=(int)floor(rounder+agesim96Ptr.data[personInd]);
         for (yearInd=1; yearInd<TDimSims+1; yearInd++) 
         { 
            if ((age96+yearInd)<(TSTART+TDIMS)) /* Not too old for model */
            {
                marstatsimMat[yearInd][personInd]=marstatsimMat[yearInd-1][personInd]; 
            }  /* Alive, for sure */
            else 
            {
                marstatsimMat[yearInd][personInd]=0; 
            }  /* Too old for model, treat as dead */

         }
      }
   }
   /*printf("iSimsmin=%d iSimsmax=%d \n", iSimsmin, iSimsmax);*/
/* for(personInd=150016; personInd<150017; personInd++) */

   for(personInd=iSimsmin; personInd<iSimsmax; personInd++)  
   {

   /* if sample household is couple, skip */
      if ( (marstatsimMat[0][personInd]==3) || (marstatsimMat[0][personInd]==0)) 
      {
         for (yearInd=0; yearInd<TDimSims; yearInd++)
         {
            marstatsimMat[yearInd][personInd]     =  0;
            healthcostsimMat[yearInd][personInd]  = -1e5; 
            MedicaidsimMat[yearInd][personInd]    = -1e5;
            zetaindexsimMat[yearInd][personInd]   = -1; 
            zetasimMat[yearInd][personInd]        = -1e5;
            xiindexsimMat[yearInd][personInd]     = -1; 
            xisimMat[yearInd][personInd]          = -1e5;
            consumptionsimMat[yearInd][personInd] = -1e5;
            beqsimMat[yearInd][personInd]         = -1e5;
            cohsimMat[yearInd][personInd]         = -1e5;
            assetsimMat[yearInd][personInd]       = -1e5;
            netIncomesimMat[yearInd][personInd]   = -1e5;
            transfersimMat[yearInd][personInd]    = -1e5;
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

   /* To interpolate along PI dimension, need to find indexes that this point 
      lies between and distance to each index. */
      result1=GetLocation(IncomeA, PIsim96Ptr.data[personInd], IDIM); 
      PIInd1=result1.Ind1;
      weight=result1.weight;
      PIInd2=PIInd1+1;

      age96=(int)floor(rounder+agesim96Ptr.data[personInd])-TSTART;

   /* timing for year 0:
         marstatsimMat[0]->healthsimhMat[0], healthsimwMat[0], zetaindexsimMat[0],
         cohsimMat[0]->consumptionsimMat[0]
   */
      yearInd=0;

      tInd=age96+yearInd;  /* age*/

   /* Decide gender by looking at marital status ******* */
   /* in beq10.gau, 1->male, 2->female */
     if (marstatsimMat[yearInd][personInd]==1) /* male */
      {
        sex=0; /* male */        
        hsInd=(int)floor(rounder+healthsimhMat[yearInd][personInd]);
      } 
     if (marstatsimMat[yearInd][personInd]==2) /* female */
      { 
        sex=1;       /*female */
        hsInd=(int)floor(rounder+healthsimwMat[yearInd][personInd]);
      } 

      /* To get consumption function, interpolate over PI and COH*/
      /* First interpolate on the last dimension */

      coh=cohsimMat[yearInd][personInd];
      consumption = Interpolation(&consumFunM[sex][tInd][PIInd1][hsInd][zetaInd][0], 
                                  cashA, coh, CASHDIM)*weight
                  + Interpolation(&consumFunM[sex][tInd][PIInd2][hsInd][zetaInd][0], 
                                  cashA, coh, CASHDIM)*(1-weight);
      consumptionsimMat[yearInd][personInd]=consumption;

   /* Timing for year yearInd: 
         marstatsimMat[yearInd]-> healthsimhMat[yearInd], 
         healthsimwMat[yearInd] or beqsimMat[yearInd-1]-> xisimMat[yearInd-1],
         epsilonsimMat[yearInd-1], zetaindexsimMat[yearInd]-> healthcostsimMat[yearInd]-> 
         cohsimMat[yearInd]->consumptionsimMat[yearInd] 
   */
      for (yearInd=1; yearInd<TDimSims+1; yearInd++)  /*calendar year*/
      {
 /* **************************************************************************************** */
      /* Decide whether it is a couple or single case by looking at marital status */
      /* in beq10.gau, 1-> male, 2->wife, (3->used to couple, we should not get here is there are any) */
         if (marstatsimMat[yearInd][personInd]==0)
         { 
         /* household members all dead */
            healthcostsimMat[yearInd][personInd]  = -1e5;
            MedicaidsimMat[yearInd][personInd]    = -1e5;
            zetaindexsimMat[yearInd][personInd]   = -1; 
            zetasimMat[yearInd][personInd]        = -1e5;
            xiindexsimMat[yearInd][personInd]     = -1; 
            xisimMat[yearInd][personInd]          = -1e5;
            consumptionsimMat[yearInd][personInd] = -1e5;
            beqsimMat[yearInd][personInd]         = cohsimMat[yearInd-1][personInd]-consumptionsimMat[yearInd-1][personInd]; 
            cohsimMat[yearInd][personInd]         = -1e5;
            assetsimMat[yearInd][personInd]       = -1e5;
            netIncomesimMat[yearInd][personInd]   = -1e5;
            transfersimMat[yearInd][personInd]    = -1e5;

            /* Timing for year yearInd: 
                  marstatsimMat[yearInd]->healthsimhMat[yearInd],
                  healthsimwMat[yearInd] or beqsimMat[yearInd-1]-> xisimMat[yearInd-1],
                  epsilonsimMat[yearInd-1], zetaindexsimMat[yearInd]-> healthcostsimMat[yearInd]->
                  cohsimMat[yearInd]->consumptionsimMat[yearInd] */
        }     
        else /* Single */
        {
            zetaInd_before=(int)floor(rounder+zetaindexsimMat[yearInd-1][personInd]);
            tInd=age96+yearInd; /* age*/

            if (marstatsimMat[yearInd][personInd]==1)
            {
               sex=0; /*husband */
               hsInd=(int)floor(rounder+healthsimhMat[yearInd][personInd]);                  
            } 

            if (marstatsimMat[yearInd][personInd]==2)
            { 
               sex=1; /*wife */
               hsInd=(int)floor(rounder+healthsimwMat[yearInd][personInd]);                  
             } 

         /* zeta=rhoHc*zetaindexsimMat[yearInd][personInd]+epsilon;  
            Zeta and xi are discretized into Markov chains.  Epsilon and xi are
            U[0,1] random variables that are combined with the chains' cdfs to 
            produce discrete-valued draws. 
            Write zeta and xi's index values into zetaindexsimMat, xiindexsimMat 
         */
            epsilon = epsiloncdfsimMat[yearInd][personInd];
            zetaInd = GetLocation(&zetaProbMcdf[zetaInd_before][0], epsilon, ZETADIM+1).Ind1;
            zeta    = zetaA[zetaInd];  /* actual value of zeta */
            xi      = xicdfsimMat[yearInd][personInd];
            xiInd   = GetLocation(&xiProbAcdf[0], xi, XIDIM+1).Ind1;
            xi      = xiA[xiInd];  /* actual value of xi */
            zetaindexsimMat[yearInd][personInd] = zetaInd;
            zetasimMat[yearInd][personInd]      = zeta;
            xiindexsimMat[yearInd][personInd]   = xiInd;
            xisimMat[yearInd][personInd]        = xi;

         /* need to interpolate over PI to get health cost intercept and variance */
            healthCost=exp((hcSingleLogMean[sex][hsInd][tInd][PIInd1]*weight
                            +hcSingleLogMean[sex][hsInd][tInd][PIInd2]*(1-weight))
                            +(hcSingleSigma[sex][hsInd][tInd][PIInd1]*weight
                            +hcSingleSigma[sex][hsInd][tInd][PIInd2]*(1-weight))*(zeta+xi)); 

         /* Calculate net after tax income for realized interest rate shock */
         /* Need to interpolate over PI to get before tax income  */
            grossIncome=rorsim[yearInd]*(cohsimMat[yearInd-1][personInd]
                                         -consumptionsimMat[yearInd-1][personInd])
                        +yM[sex][PIInd1][tInd]*weight+yM[sex][PIInd2][tInd]*(1-weight);
            netIncome=AfterTaxIncome(grossIncome);

         /* update assets and COH */
            assets = cohsimMat[yearInd-1][personInd] - consumptionsimMat[yearInd-1][personInd];
            coh_pre_transfer = assets+netIncome; 
            Medicaid = 0;
            if (healthCost>coh_pre_transfer)     /* Code in Medicaid */
            {
               if (coh_pre_transfer>0)
               {
                  Medicaid   = healthCost - coh_pre_transfer;
                  healthCost = coh_pre_transfer;                  
               }
               else
               {
                  Medicaid   = healthCost;
                  healthCost = 0;                  
               }                 
            }              
            coh_pre_transfer = coh_pre_transfer - healthCost;
            healthcostsimMat[yearInd][personInd]=healthCost;
            MedicaidsimMat[yearInd][personInd]=Medicaid;

            coh=max(coh_pre_transfer,cMin);

         /* To get consumption function, interpolate over PI and COH */
            netIncomesimMat[yearInd][personInd]= netIncome;
            assetsimMat[yearInd][personInd]    = assets;
            cohsimMat[yearInd][personInd]      = coh;
            transfersimMat[yearInd][personInd] = coh - (coh_pre_transfer-Medicaid); 
            consumption = Interpolation(&consumFunM[sex][tInd][PIInd1][hsInd][zetaInd][0],
                                        cashA, coh, CASHDIM)*weight
                        + Interpolation(&consumFunM[sex][tInd][PIInd2][hsInd][zetaInd][0],
                                        cashA, coh, CASHDIM)*(1-weight);
            consumptionsimMat[yearInd][personInd]=consumption;

         }/* end else */

        marstatsim2Mat[yearInd][personInd] = marstatsimMat[yearInd][personInd];
      } /* end loop of year */

      marstatsim2Mat[0][personInd]        = marstatsimMat[0][personInd];
      marstatsim2Mat[TDimSims][personInd] = marstatsimMat[TDimSims][personInd];
   }  /* end loop of household */
   
   free(tempvec);
   free(tempV);
}

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

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  This function reads a Gauss v5.0 fmt file into a Matrix structure and  */
/*  returns the structure.  */


GMatrix getcsvdat(char *csv, int nRows, int nCols)
{
	/*  Initialize the matrix to be 1x1 and the byte/bit order to 0.  */
	GMatrix mat = { 1, 1 };
	int BSIZE = 2000;  // Meant to exceed the number of characters in a row of a standard dataset
	double *csvvec = malloc(nRows*nCols * sizeof(double));
	char *rowBuffer = malloc(BSIZE * sizeof(char));
	char *oneval;
	int iRow, iCol, count;
	FILE *csvfile;

	mat.m = nRows;
	mat.n = nCols;

	csvfile = fopen(csv, "rb");
	count = 0;

	for (iRow = 0; iRow < nRows; iRow++)
	{
		while (fgets(rowBuffer, BSIZE, csvfile))
		{
			oneval = strtok(rowBuffer, ",\n"); // first item
			csvvec[count] = atof(oneval);
			count++;

			for (iCol = 1; iCol < nCols; iCol++)
			{
				while (oneval = strtok(NULL, ",\n")) // NULL allows us to walk through rowBuffer
				{
					csvvec[count] = atof(oneval);
					count++;
				}
			}
		}
	}

	fclose(csvfile);
	mat.data = csvvec;
	return (mat);
}

void writecsvdat(char *csv, GMatrix mat)
{
FILE *file_ptr;
	int iCount;
	double aux;
	
	file_ptr = fopen(csv, "wb");
	for (iCount = 0; iCount< mat.m; iCount++) {
		aux = mat.data[iCount];
			fprintf(file_ptr, "%20.20lf \n", aux);
	}
	fclose(file_ptr);
}



GMatrix gau5read(char *fmt) 
   {
/*  Initialize the matrix to be 1x1 and the byte/bit order to 0.  */
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

/*  Allocate enough space to store the header.  */
    cread = (unsigned char *) malloc(BASIC_HEADER_LEN); 
/*  Open *fmt for reading only.  */
    fd = fopen(fmt, "rb"); 
  
/*  Read the basic header (128 bytes) all at once.  */

    fread((void *) cread, 1, BASIC_HEADER_LEN, fd);
    byte = (int) *(cread + (BYTE_POS * sizeof(int)));  /* (0=Backward) */
    bit = (int) *(cread + (BIT_POS * sizeof(int)));    /* (0=Backward) */

/*  To get some system independence, we detect whether we have to reverse */
/*  the bytes or bits or both.  If x and x_SYSTEM match, no reverse is  */
/*  necessary. */

    if ((bit || BIT_SYSTEM) && !(bit && BIT_SYSTEM)) bit_rev=1;
    if ((byte || BYTE_SYSTEM) && !(byte && BYTE_SYSTEM)) byte_rev=1;

    type = *( (int *) gread((cread + (TYPE_POS * sizeof(int))), sizeof(int), 
                             byte_rev, bit_rev) );

/*  If the fmt file type is not a scalar, there are another two */
/*  ints of header giving the values of m and n.  If a matrix, also reset n. */

    if (type > SCALAR) 
       { 
        fread((void *) cread, 1, sizeof(int), fd);
        mat.m = *((unsigned int *) gread(cread, sizeof(int), byte_rev, bit_rev));
        fread((void *) cread, 1, sizeof(int), fd);
        if (type == MATRIX)
          mat.n = *((unsigned int *) gread(cread, sizeof(int), byte_rev, bit_rev));
      } 

/*  Allocate memory for the matrix.  The amount needed is m * n * sizeof(double). */
/*  Next, read in the data all at once.  Then use gread to reverse */
/*  bytes/bits if necessary. */

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

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*  This function writes a Gauss v5.0 fmt file from a Matrix structure. */

void gau5write(char *fmt, GMatrix mat) 
   {
/*  This ugly mess is the basic header. */

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

/*  If forward byte, make 6th int 0xffffffff. */
/*  If forward bit, make 7th int 0xffffffff. */

    if (BYTE_SYSTEM) header[BYTE_POS] = 0xffffffff;
    if (BIT_SYSTEM) header[BIT_POS] = 0xffffffff;

/*  If not a scalar, increase the 16th int by 1 and the 19th int (header */ 
/*  length) by 8 (2 * sizeof(int)).  Also, set m in int 33. */

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
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* create grid for cash on hand and consumption   */ 

void Grid2(double cashA[], double consumA[CONSUMDIM], double hsA[HSDIM], double IncomeA[IDIM])
{
   int cashInd, consumInd,IInd, hsInd;
   FILE *cashP; 
   char fullpath[ADDRESS_LEN];
   double cMin2, cMin3;

/* create cash grid, cMin<=x<=CASHMAX */

   for (cashInd=0; cashInd<CASHDIM1; cashInd++)
   { 
      cashA[cashInd]=pow(sqrt(cMin)+(sqrt(CASHMAX1)-sqrt(cMin))*cashInd/(CASHDIM1-1),2);
   }

   cMin2=cashA[CASHDIM1-1];

   for (cashInd=0; cashInd<CASHDIM2; cashInd++)
   { 
      cashA[CASHDIM1+cashInd]=pow(sqrt(cMin2)+(sqrt(CASHMAX2)-sqrt(cMin2))*(cashInd+1)/(CASHDIM2),2);
   }

   cMin3=cashA[CASHDIM1+CASHDIM2-1];

   for (cashInd=1; cashInd<CASHDIM-CASHDIM1-CASHDIM2+1; cashInd++)
   {
      cashA[cashInd+CASHDIM1+CASHDIM2-1]
         = pow(sqrt(cMin3)+(sqrt(CASHMAX3)-sqrt(cMin3))*cashInd/(CASHDIM-CASHDIM1-CASHDIM2),2);
   }

   for (consumInd=0; consumInd<CONSUMDIM1; consumInd++)
   { 
      consumA[consumInd]=pow(sqrt(cMin)+(sqrt(CONSUMPTIONMAX1)-sqrt(cMin))*consumInd/(CONSUMDIM1-1),2);
   }

   cMin2=consumA[CONSUMDIM1-1];

   for (consumInd=0; consumInd<CONSUMDIM2; consumInd++)
   { 
      consumA[CONSUMDIM1+consumInd]
         = pow(sqrt(cMin2)+(sqrt(CONSUMPTIONMAX2)-sqrt(cMin2))*(consumInd+1)/(CONSUMDIM2),2);
   }

   cMin3=consumA[CONSUMDIM1+CONSUMDIM2-1];

   for (consumInd=1; consumInd<CONSUMDIM-CONSUMDIM1-CONSUMDIM2+1; consumInd++)
   {
      consumA[consumInd+CONSUMDIM1+CONSUMDIM2-1]
         = pow(sqrt(cMin3)+(sqrt(CONSUMPTIONMAX3)-sqrt(cMin3))*consumInd/(CONSUMDIM-CONSUMDIM1-CONSUMDIM2),2);
   }

   for (IInd = 0; IInd<IDIM; IInd++)
   {
      IncomeA[IInd]= ((double) IInd)/(IDIM-1); /* Income expressed as percentile ranking */
   }
      
   for ( hsInd = 0; hsInd<HSDIM; hsInd++)
   {
      hsA[hsInd]=hsInd;
   }

   if ( (useMPI<2)||(rank==(size-1)) )
   {
      cashP=fopen(strcat(strcpy(fullpath,outputdir),"xA.txt"),"w");

      for (cashInd=0; cashInd<CASHDIM1; cashInd++)
      { 
         fprintf(cashP,"%10.3f\n", cashA[cashInd]);
      }

      cMin2=cashA[CASHDIM1-1];

      for (cashInd=0; cashInd<CASHDIM2; cashInd++)
      { 
         fprintf(cashP,"%10.3f\n", cashA[CASHDIM1+cashInd]);
      }

      cMin3=cashA[CASHDIM1+CASHDIM2-1];

      for (cashInd=1; cashInd<CASHDIM-CASHDIM1-CASHDIM2+1; cashInd++)
      {
         fprintf(cashP,"%10.3f\n", cashA[cashInd+CASHDIM1+CASHDIM2-1]);
      }
  
      fclose(cashP);
      cashP=fopen(strcat(strcpy(fullpath,outputdir),"consumA.txt"),"w");

      for (consumInd=0; consumInd<CONSUMDIM1; consumInd++)
      { 
         fprintf(cashP,"%10.3f\n", consumA[consumInd]);
      }

      for (consumInd=0; consumInd<CONSUMDIM2; consumInd++)
      { 
         fprintf(cashP,"%10.3f\n", consumA[CONSUMDIM1+consumInd]);
      }

      for (consumInd=1; consumInd<CONSUMDIM-CONSUMDIM1-CONSUMDIM2+1; consumInd++)
      {
         fprintf(cashP,"%10.3f\n", consumA[consumInd+CONSUMDIM1+CONSUMDIM2-1]);
      }

      fclose(cashP);
   }
}
