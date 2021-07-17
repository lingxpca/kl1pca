#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "type.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>

#define pt(A) for(i=0.0;i<numentities_n;++i){for(j=0;j<numentities_n;++j) {(printf("%7.4g\t",A));}printf("\n");}
#define loop(A) for(i=0.0;i<numentities_n;++i){for(j=0;j<numentities_n;++j) {A}}

int allocateMemory (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
int kernelMatrix (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
int pfkpca(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
int readdata  (IOINFOptr ioinfo, ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
int kernell1pca (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);

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

  if (argc < 3) {
    fprintf (stderr, "Usage \n kernel_l1pca <filename> <q> <kernel> <parameter>\n");
    fprintf (stderr, "Example: \n kernel_l1pca toy 2 dot 1\n");
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
	entityinfo.rho=atoi(argv[4]);
	strcpy(entityinfo.norm,argv[5]);

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

 if (strcmp(entityinfo.norm,"l2")==0) {
   status = kernelMatrix(&entityinfo, &probleminfo);
  if (status) {
    fprintf (stderr, "Unable to kernelMatrix\n");
    goto TERMINATE;
  }
 }
 else{
    status = kernell1pca(&entityinfo, &probleminfo);
  if (status) {
    fprintf (stderr, "Unable to solve.  Terminating...; or done\n");
    goto TERMINATE;
  }
  }

  status = pfkpca(&entityinfo, &probleminfo);
  if (status) {
    fprintf (stderr, "Unable to pfkpca\n");
    goto TERMINATE;
  }

TERMINATE:

  free_and_null ((char **) &ioinfo.datafname);
  free_and_null ((char **) &entityinfo.K);
  free_and_null ((char **) &probleminfo.PCA);
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
 	entityinfo->BX = (double *)malloc( entityinfo->numattributes_m*sizeof(double));
	entityinfo->XBZX = (double *)malloc( entityinfo->numattributes_m*sizeof(double));
	entityinfo->a= (double *)calloc((entityinfo->numentities_n*entityinfo->numentities_n),sizeof(double)); 

	entityinfo->KK=(double **)calloc(entityinfo->numentities_n,sizeof(double *));
	entityinfo->K=(double **)calloc(entityinfo->numentities_n,sizeof(double *));
	for (int i=0;i<entityinfo->numentities_n;++i){
		entityinfo->K[i]=(double *)calloc(entityinfo->numentities_n,sizeof(double));
		entityinfo->KK[i]=(double *)calloc(entityinfo->numentities_n,sizeof(double));
	}
	return 0;
}

 
int kernelMatrix (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo)
{
	int numattributes_m = entityinfo->numattributes_m;
	int numentities_n = entityinfo->numentities_n;   
	int i,j,l,m; 
	double s0,s1,s4;
	double dd[numentities_n*numentities_n]; 


 															// polynomial kernel
	if (strcmp(entityinfo->kernel,"pol")==0) {  
		s4=1.0;
		// for(i=0;i<numentities_n;++i) {
		// 	for(j=0;j<numentities_n;++j) {
		// 		for (l=0;l<numattributes_m;++l) entityinfo->KK[i][j] += entityinfo->points_XT[i*numattributes_m+l]*entityinfo->points_XT[j*numattributes_m+l];
		// 		s4 += entityinfo->KK[i][j];
		// 	}
		// }
		s4=1.0;
		loop(for (l=0;l<numattributes_m;++l) entityinfo->KK[i][j] += entityinfo->points_XT[i*numattributes_m+l]*entityinfo->points_XT[j*numattributes_m+l];s4 += entityinfo->KK[i][j];)


   	

		for (i =0;i<numentities_n;++i) {
			for (j =0; j<numentities_n;++j)	{
				s0 = 0.0;
				s1 = 0.0;
				for(m=0;m<numentities_n;++m) {
					s0 += entityinfo->KK[i][m];
					s1 += entityinfo->KK[m][j];
				}
				entityinfo->K[i][j] = entityinfo->KK[i][j]-s0/numentities_n-s1/numentities_n+s4/pow(numentities_n,2); //centered POL kernel
			} 
		}
		pt(entityinfo->K[i][j]);
	}
		

             // Gaussian Kernel
	else if (strcmp(entityinfo->kernel,"rbf")==0) {

		entityinfo->vox=0.0;
		for (j=0;j<numattributes_m;++j)  {
			s0 = 0.0;
			s1 = 0.0;
			for (i = 0;i<numentities_n;++i) s1 += entityinfo->points_XT[i*numattributes_m+j]/numentities_n;
			for (i = 0;i<numentities_n;++i) s0 += pow(entityinfo->points_XT[i*numattributes_m+j]-s1,2)/(numentities_n-1);
			entityinfo->vox += s0;
		} 
		
		s4=0.0;
		for (i = 0; i<numentities_n;++i) {
			for (j = 0; j<numentities_n;++j) { 
				for (l=0;l<numattributes_m;++l) entityinfo->KK[i][j] += -pow(entityinfo->points_XT[i*numattributes_m + l]-entityinfo->points_XT[j*numattributes_m+l],2);
				entityinfo->KK[i][j] = exp(entityinfo->KK[i][j]/(entityinfo->vox*entityinfo->rho));
				s4 += entityinfo->KK[i][j];
			}
		}
											//centered kernel
		for (i =0;i<numentities_n;++i) 	{
			for (j =0; j<numentities_n;++j)	{
				s0 = 0.0;
				s1 = 0.0;
				for(m=0;m<numentities_n;++m) 	{
					s0 += entityinfo->KK[i][m];
					s1 += entityinfo->KK[m][j];
				}
				entityinfo->K[i][j] = entityinfo->KK[i][j]-s0/numentities_n-s1/numentities_n+s4/pow(numentities_n,2); //centered RBF kernel
			} 
		}	
	}

	for (i=0;i<numentities_n;++i)
 		for (j=0;j<numentities_n;++j)
			dd[i*numentities_n+j] = entityinfo->K[i][j];
						//the eigenvector has been centerlized
	gsl_matrix_view kernelmatrix = gsl_matrix_view_array (dd, numentities_n, numentities_n);
	gsl_vector *eval = gsl_vector_alloc (numentities_n);
	gsl_matrix *evec = gsl_matrix_alloc (numentities_n, numentities_n);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (numentities_n);
	gsl_eigen_symmv (&kernelmatrix.matrix, eval, evec, w);
	gsl_eigen_symmv_sort (eval, evec,GSL_EIGEN_SORT_ABS_DESC);

   for (i=0;i<numentities_n;++i)  
  		for(j=0;j<numentities_n;++j)   
	  		for(l=0;l<probleminfo->q;l++) entityinfo->a[i*numentities_n+j] += (gsl_matrix_get(evec,i,l)*pow(gsl_vector_get(eval,l),-0.5))*(gsl_matrix_get(evec,j,l)*pow(gsl_vector_get(eval,l),-0.5));

	gsl_vector_free (eval);
 	gsl_eigen_symmv_free (w);
 	gsl_matrix_free (evec);
	return 0;
}
