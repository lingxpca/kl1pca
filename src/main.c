#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "type.h" 
#include "mkl.h"


int allocateMemory (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
int kernelMatrix (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
int pfkpca(IOINFOptr ioinfo,ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
int readdata  (IOINFOptr ioinfo, ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
int kernell1pca (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
int pca(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);

static void
   free_and_null (char **ptr);

int main (int argc, char *argv[])
{
  int status; 
  IOINFO      ioinfo;
  ENTITYINFO  entityinfo;
  PROBLEMINFO probleminfo;

  entityinfo.points_XT = NULL;
  entityinfo.K=NULL;
  ioinfo.location = NULL;
  ioinfo.problem  = NULL;
  ioinfo.datafname = NULL;
  ioinfo.datafile  = NULL;
  entityinfo.kernel = NULL; 
  if (argc < 4) {
    fprintf (stderr, "Usage \n krpca <filename> <q> <kernel> <norm> <para>\n");
    fprintf (stderr, "Example: \n krpca toy 2 rbf l1 2.1\n"); 
    return 1;
  }
  
  ioinfo.location = (char *) malloc ((PATHLENGTH) * sizeof (char));
  ioinfo.problem  = (char *) malloc ((PATHLENGTH) * sizeof (char));
  entityinfo.kernel = (char *) malloc ((PATHLENGTH)* sizeof(char));
  entityinfo.norm = (char *) malloc ((PATHLENGTH)* sizeof(char));

  strcpy (ioinfo.location, "./");
  strcpy (ioinfo.problem, argv[1]);
  probleminfo.q=atoi(argv[2]);
	strcpy(entityinfo.kernel,argv[3]); 
	strcpy(entityinfo.norm,argv[4]);
  probleminfo.var=atof(argv[5]);
  
 
  ioinfo.datafname = (char *) malloc ((PATHLENGTH) * sizeof (char));
  strcpy (ioinfo.datafname, ioinfo.location);
  strcat (ioinfo.datafname, ioinfo.problem);
  ioinfo.datafile = fopen (ioinfo.datafname, "r");
  if (ioinfo.datafile == NULL) {
    fprintf (stderr, "Unable to open data file.  Terminating...\n");
    return 1;
  }
  
  status = readdata(&ioinfo, &entityinfo, &probleminfo);  /* in readdata.c */
  if (status) {
    fprintf (stderr, "Unable to read data.  status %d.  Terminating...\n", status);
    goto TERMINATE;
  }
  
  status = allocateMemory(&entityinfo, &probleminfo); /* at the end of this file */
  if (status) {
    fprintf (stderr, "Unable to allocate memory\n");
    goto TERMINATE;
  } 

  status = kernelMatrix(&entityinfo, &probleminfo);
  if (status) {
    fprintf (stderr, "Unable to kernelMatrix\n");
    goto TERMINATE;
  }
 
  if (strcmp(entityinfo.norm,"l2")==0) {
    status = pca(&entityinfo, &probleminfo);
    if (status) {
      fprintf (stderr, "Unable to solve.  Terminating...; or done\n");
      goto TERMINATE;
    }
  }
  else {
    status = kernell1pca(&entityinfo, &probleminfo);// in krpca.c
    if (status) {
      fprintf (stderr, "Unable to solve.  Terminating...; or done\n");
      goto TERMINATE;
    }
  }
 
  status = pfkpca(&ioinfo,&entityinfo, &probleminfo);
  if (status) {
    fprintf (stderr, "Unable to pfkpca\n");
    goto TERMINATE;
  }

TERMINATE:

  free_and_null ((char **) &ioinfo.datafname);
  free_and_null ((char **) &entityinfo.K);
  return (status);
}
   															 /* END free_and_null */
static void
free_and_null (char **ptr)
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
}

																	//allocate memory
int allocateMemory (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo)  {   
	entityinfo->a= (double *)calloc((entityinfo->numentities_n*entityinfo->numentities_n),sizeof(double)); 
  
	entityinfo->KK=(double **)calloc(entityinfo->numentities_n,sizeof(double *)); 
	for (int i=0;i<entityinfo->numentities_n;++i){ 
		entityinfo->KK[i]=(double *)calloc(entityinfo->numentities_n,sizeof(double));
	}

  entityinfo->svd =(double *)calloc((entityinfo->numentities_n*entityinfo->numentities_n),sizeof(double)); 
  
	return 0;
}

 
int kernelMatrix (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo)
{ 
	int numattributes_m = entityinfo->numattributes_m;
	int numentities_n = entityinfo->numentities_n;   
	int i,j,l; 
	double *ki =(double *)calloc(numentities_n,sizeof(double)); 
  double s0,s1,s4,rho; 
  
	if (strcmp(entityinfo->kernel,"qua")==0) {   
		loop	(
      for (l=0,s1=probleminfo->var;l<numattributes_m;++l) s1 += entityinfo->points_XT[i*numattributes_m+l]*entityinfo->points_XT[j*numattributes_m+l];
      entityinfo->KK[i][j] = (s1)*(s1);
			s4 += entityinfo->KK[i][j];
		)
    
    loop (ki[i] += entityinfo->KK[i][j];)
    loop (entityinfo->svd[i*numentities_n+j]= entityinfo->KK[i][j]-ki[i]/numentities_n-ki[j]/numentities_n+s4/(numentities_n*numentities_n);)
    
 	}
  else if (strcmp(entityinfo->kernel,"rbf")==0) //Gaussian kernel
  { 
    entityinfo->var=0.0;
    for (j=0;j<numattributes_m;++j)  {  
      for (i = 0,s1=0.0;i<numentities_n;++i) s1 += entityinfo->points_XT[i*numattributes_m+j]/numentities_n;
      for (i = 0,s0=0.0;i<numentities_n;++i) s0 += pow(entityinfo->points_XT[i*numattributes_m+j]-s1,2)/(numentities_n-1); 
      entityinfo->var += s0;  
    }
     
    rho = probleminfo->var;//* entityinfo->var ;
    printf("%8.5g\n",entityinfo->var );
    loop (
          for (l=0;l<numattributes_m;++l) entityinfo->KK[i][j] += -pow(entityinfo->points_XT[i*numattributes_m + l]-entityinfo->points_XT[j*numattributes_m+l],2); 
          entityinfo->KK[i][j] = exp(entityinfo->KK[i][j]/rho);
          s4 += entityinfo->KK[i][j]; 
    )	

    loop (ki[i] += entityinfo->KK[i][j];)
    loop (entityinfo->svd[i*numentities_n+j] = entityinfo->KK[i][j]-ki[i]/numentities_n-ki[j]/numentities_n+s4/pow(numentities_n,2);)     
  }
  
  free(ki);
  return 0; 
}

int pca (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo){
 
	int numentities_n = entityinfo->numentities_n,i,j;  
  double *svd = (double *)calloc(numentities_n*numentities_n,sizeof(double)); 	
  for (i=0;i<numentities_n;++i){
    for (j=0;j<numentities_n;++j){
      svd[i*numentities_n+j] = (i>=j) ?entityinfo->svd[i*numentities_n+j] : 0.0;
      }
  } 

  MKL_INT n = numentities_n, lda = numentities_n;
  double w[numentities_n];
  LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'L', n, svd, lda,w);
   
  double *svd2 = (double *)calloc(probleminfo->q*numentities_n,sizeof(double)); 
  for (i =0;i<numentities_n;++i) {
    for (j = 0; j<probleminfo->q;++j) {
      svd2[i*probleminfo->q+j] = svd[i*numentities_n+numentities_n-1-j]/sqrt(w[numentities_n-j-1]);
    }
  } 

  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,numentities_n,numentities_n,probleminfo->q,1.0,svd2,probleminfo->q,svd2,probleminfo->q,0.0,entityinfo->a,numentities_n);
  
  free(svd2);
  free(svd);
	return 0;
}