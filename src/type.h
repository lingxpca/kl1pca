#define PATHLENGTH 500 /* length of path name describing location of data file*/
#define VERBOSITY 1 /* 1-write rotation file, print k, l , 2- print objective value  */
#define OBJ_INIT 100000000000000.0 
#define SCORES 0 /* create file for scores */
#define PROJPTS 0 /* create file for reconstructions */


struct entityinfo {
  int numentities_n;/*# of rows/points in your data*/
  int numattributes_m;/*# of columns*/
  double *points_XT;/*actual data points in column major format (i.e., transposed) - in a single array*/
  double **K;/*centered kernel matrix*/
  double **KK;//orginal kernel 
  double *a; 
  double *XBZX;
  double *BX;
  char *kernel;
  double rho; 
  char *norm;
  double vox;
};
typedef struct entityinfo ENTITYINFO, *ENTITYINFOptr;

/* input/output info */
struct ioinfo {
  int        status; 
  char       *problem;
  char       *location;
  char       *datafname;
  FILE       *datafile; 
  char       *rotationfname;
  FILE       *rotationfile;
  char       *projfname;
  FILE       *projfile;
  char       *scoresfname;
  FILE       *scoresfile;
  char       *sparse;
};

typedef struct ioinfo IOINFO, *IOINFOptr;

/* problem info */
struct probleminfo {
  int i;
  int j;
  int k;
  int l;
  int l1;
  int l2;
  double *c;
  double *ck_plus1; 
  double *PCA;
  double **cstar; 
  
  double *kc; 
  
  double *d;/*steepest descent direction at nbos*/
  
  int q;
  double minobjective;
  double objective;
  int    index;
  int    jstar;
  double *PCs;
  double *ratios;
  double **tosort;
  double *weights;
  double sumWeights;
  double medWeight;
  int    origIndex;
  double *v;
  double innerprod;
  double temp;
  double normv;

  double *lambdas;
  int num_lambdas;
  double lambda_max;
  int **num_lambdas_lj;
  double sum_below;
  double sum_above;
  double lambdaU;
  double lambdaL;
  double ***v_lj;
  double ***lambdas_lj;
  double ****lambdas_lj_sort;
  int **curr_lambda_lj;
  int t;
  double *lambdas_out;
  int num_distinct_lambdas;
  int max_memory;
  int **max_memory_lj;
  int max_memory_lambdas;
  char *sparse;
  double *vs;
  double *zs;
  double maxL;
  double minU;
  int lambda_curr;
  int lambda_next;
  int check_lambda;
  int v_index;
};
typedef struct probleminfo PROBLEMINFO, *PROBLEMINFOptr;
 
