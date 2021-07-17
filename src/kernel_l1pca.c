#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "type.h"

int kernell1pca (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);  

int kernell1pca (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
	int i,j,k,l,m;
	int q = probleminfo->q;
 	int numentities_n = entityinfo->numentities_n; 
	int numattributes_m=entityinfo->numattributes_m;
 	double *points_XT = entityinfo->points_XT; 
	double rho = entityinfo->rho;
		
	double *c0;
	int jstar=0; 
 	double s0,s1,s2; 
	double *KK;
	double **K;  
	
	double **cstar; 
	double **cnorm; 
	double *kckc;
  																//allocate memory 
	kckc=  (double *)calloc(numentities_n,sizeof(double));   
	c0 = (double *)calloc(numentities_n,sizeof(double)); 
	cstar=(double **)calloc(numentities_n,sizeof(double *));
	KK = (double *)calloc(numentities_n*numentities_n,sizeof(double));  
	K=(double **)malloc(numentities_n*sizeof(double *));
	cnorm = (double **)calloc(numentities_n,sizeof(double)); 
	for (i=0;i<numentities_n;++i){ 
		K[i]=(double *)malloc(numentities_n*sizeof(double));
		cstar[i]=(double *)calloc(q,sizeof(double));
		cnorm[i]=(double *)calloc(q,sizeof(double));
	}
	
  
///////////////////////////////////////////////// polnomial kernel//////////////////////////////////////////////// 
	if (strcmp(entityinfo->kernel,"pol")==0) {	

		s2=1.0;
		for(i=0;i<numentities_n;++i) {
			for(j=0;j<numentities_n;++j) {
				for (l=0;l<numattributes_m;++l) KK[i*numentities_n+j]+= entityinfo->points_XT[i*numattributes_m+l]*entityinfo->points_XT[j*numattributes_m+l];
				s2 += KK[i*numentities_n+j];
			}
		}

		for (i =0;i<numentities_n;++i) {
			for (j =0; j<numentities_n;++j)	{
				s0 = s1 = 0.0; 
				for(m=0;m<numentities_n;++m) {
					s0 += KK[i*numentities_n+m];
					s1 += KK[m*numentities_n+j];
				}
				K[i][j] = KK[i*numentities_n+j]-s0/numentities_n-s1/numentities_n+s2/pow(numentities_n,2); //centered POL kernel
			} 
		}
	}	
	
	
	
//////////////////////////////////////////////// Gaussian Kernel//////////////////////////////////////////////// 
	else if (strcmp(entityinfo->kernel,"rbf")==0)      {

		s1=0.0;
		for (j=0;j<numattributes_m;++j)  {
			s0 = s2 = 0.0; 
			for (i = 0;i<numentities_n;++i) s2 += points_XT[i*numattributes_m+j]/numentities_n;
			for (i = 0;i<numentities_n;++i) s0 += pow(points_XT[i*numattributes_m+j]-s2,2)/(numentities_n-1);
			s1 += s0;
		}	
		entityinfo->vox = s1;
	
		s2=0.0;
		for (i = 0; i<numentities_n;++i)	{
			for (j = 0; j<numentities_n;++j) {
				for (l=0;l<numattributes_m;++l) KK[i*numentities_n+j]+= -pow(points_XT[i*numattributes_m + l]-points_XT[j*numattributes_m+l],2);
				KK[i*numentities_n+j] = exp(KK[i*numentities_n+j]/(s1*rho));
				entityinfo->KK[i][j]=KK[i*numentities_n+j]; //original kernel w/o centerlized
				s2 += KK[i*numentities_n+j]; 
			}
		}
					//centered kernel
 
		for (i =0;i<numentities_n;++i)		{
			for (j =0; j<numentities_n;++j)		{
				s0 = s1 = 0.0; 
				for(m=0;m<numentities_n;++m) 				{
					s0 += KK[i*numentities_n+m];
					s1 += KK[m*numentities_n+j];
				}
				K[i][j] = KK[i*numentities_n+j]-s0/numentities_n-s1/numentities_n+s2/pow(numentities_n,2);
			} 
		}
	}
  for (i =0;i<numentities_n;++i)		{
			{ for (j =0; j<numentities_n;++j)		
				printf("%10.5g\t",K[i][j]);			}printf("\n");}


	////////////////////////////////////////////////Kim L1-norm KPCA//////////////////////////////////////////////////////////////////////

   for (l=0;l<q;++l)	{ // find each loading 

 		s0=0.0;	
		for (j=0;j<numentities_n;++j)		{
			s1=0.0;
			for (i=0;i<numentities_n;++i) s1 += fabs(K[i][j]);
			if(s0 < s1/sqrt(K[j][j]))	{
				s0 = s1/sqrt(K[j][j]);
				jstar=j;
			}
		} 
	   for (i=0;i<numentities_n;++i) c0[i]= (K[i][jstar] > 0.0) ? 1.0 : -1.0;
		
	 	do	 {  
			for (i=0;i<numentities_n;++i)   { 
				s0 = 0.0;
			 	for (j=0; j<numentities_n;++j) s0 +=  K[i][j]*c0[j];  
				cstar[i][l] = (s0 > 0.0) ? 1.0 : -1.0;  
			}		 
			
			s1=0.0;
			for (i=0;i<numentities_n;++i)    {
				s0 = 0.0;
				for (j=0;j<numentities_n;++j) s0 += (c0[j]-cstar[j][l])*K[j][i];
				s1 += s0*(c0[i]-cstar[i][l]); 
			}
			 
			for (i=0;i<numentities_n;++i) c0[i] = cstar[i][l]; 
		} while(s1 != 0.0) ; 
		
		for (i=0;i<numentities_n;++i) printf("%g\n",cstar[i][l]);

 	 	s1 =0.0;
 	 	for (i=0;i<numentities_n;++i) 	{ 
			s0 =0.0; 
			for (j=0; j<numentities_n;++j) s0 += cstar[j][l]*K[j][i]; 
			s1 += s0*cstar[i][l];  
		}

		for (i = 0;i<numentities_n;++i) {
			cnorm[i][l] = cstar[i][l]/sqrt(s1);
		}

		/*update kernelmatrix for next loading*/	

		for (i = 0;i<numentities_n;++i) 	{	
			kckc[i] = 0.0; 
			for (k=0;k<numentities_n;++k) kckc[i] += K[i][k]*cstar[k][l]; 
	   }
 
	   for (i = 0;i<numentities_n;++i)   
	   	for (j=0;j<numentities_n;++j)  
		  	K[i][j]  =  K[i][j] - (kckc[i]*kckc[j])/s1;  

	}

	 for (i=0;i<numentities_n;++i)  
  		for(j=0;j<numentities_n;++j)   
	  		for(l=0;l<probleminfo->q;l++) entityinfo->a[i*numentities_n+j] +=  cnorm[i][l]*cnorm[j][l]; 

   
  	return 0;
}
  
 
    
 








































 








