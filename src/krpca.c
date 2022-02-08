#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "type.h"
#include <mkl.h>  

int kernell1pca (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);  

int kernell1pca (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
	int i,j,l; 
 	int numentities_n = entityinfo->numentities_n;  
		
	double *c0;
	int jstar=0; 
 	double s1;   
	double *findjstar;  
	double *cnorm;  
	double *ck;   
    
	c0 = (double *)calloc(numentities_n,sizeof(double));  
	
	cnorm = (double *)calloc(numentities_n*probleminfo->q,sizeof(double)); 
	ck = (double *)calloc(numentities_n,sizeof(double));  
	findjstar = (double *)calloc(numentities_n,sizeof(double));  
  
	double *A=(double *)calloc(numentities_n*numentities_n,sizeof(double)); 
	double *B=(double *)calloc(numentities_n*numentities_n,sizeof(double)); 
	double *svd =  (double *)calloc(numentities_n*numentities_n,sizeof(double)); 
	double *Aprime=(double *)calloc(numentities_n*numentities_n,sizeof(double)); 
	double *Y = (double *)calloc(numentities_n,sizeof(double)); 
	double *kc = (double *)calloc(numentities_n,sizeof(double)); 
	double *c = (double *)calloc(numentities_n,sizeof(double)); 
	double *AC = (double *)calloc(numentities_n,sizeof(double));  

 	for (i=0;i<numentities_n;++i) A[i*numentities_n+i] = 1.0; 
   	for (l=0;l<probleminfo->q;++l)	{
		 	   	
		for (j=0;j<numentities_n;++j) {
			findjstar[j]=0.0;   
			for (i=0;i<numentities_n;++i) findjstar[j] += fabs(entityinfo->svd[i*numentities_n+j])/sqrt(entityinfo->svd[j*numentities_n+j]);   
		}
		jstar = cblas_idamax(numentities_n,findjstar,1); 
		 
	   	for (i=0;i<numentities_n;++i) c0[i]= (entityinfo->svd[i*numentities_n+jstar] > 0.0) ? 1.0 : -1.0;  
		
		do	{
			cblas_dgemv(CblasRowMajor,CblasNoTrans,numentities_n,numentities_n,1.0,entityinfo->svd,numentities_n,c0,1,0.0,Y,1); 		 
			for (i=0;i<numentities_n;++i) c[i] = (Y[i] > 0.0) ? 1.0 : -1.0;
			for (i=0;i<numentities_n;++i) ck[i] = c0[i] - c[i];  
			s1 = cblas_ddot(numentities_n,ck,1,ck,1); 
			//cblas_dsymv(CblasRowMajor,CblasLower,numentities_n,1.0,entityinfo->svd,numentities_n,ck,1,0.0,Y,1);
			//s1 = cblas_ddot(numentities_n,Y,1,ck,1);	
			for (i=0;i<numentities_n;++i)   c0[i] =c[i];
		} while(fabs(s1) >10e-7);    
 
		 
		cblas_dsymv(CblasRowMajor,CblasLower,numentities_n,1.0,entityinfo->svd,numentities_n,c,1,0.0,kc,1);  //KC
		cblas_dgemv(CblasRowMajor,CblasNoTrans,numentities_n,numentities_n,1.0,A,numentities_n,c,1,0.0,AC,1); //Ac
		s1 = cblas_ddot(numentities_n,kc,1,c,1); 
		
		for (i=0;i<numentities_n;++i) cnorm[i*probleminfo->q+l] = AC[i]/sqrt(s1); 
		
		if (l == (probleminfo->q-1)) break;
		loop(Aprime[i*numentities_n+j] = (i==j) ? 1.0 - c[i]*kc[j]/s1 : -c[i]*kc[j]/s1;) 

		cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,numentities_n,numentities_n,numentities_n,1.0,entityinfo->svd,numentities_n,Aprime,numentities_n,0.0,svd,numentities_n);
		cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,numentities_n,numentities_n,numentities_n,1.0,A,numentities_n,Aprime,numentities_n,0.0,B,numentities_n);
		memcpy(A,B,numentities_n*numentities_n*sizeof(double));
	 	memcpy(entityinfo->svd,svd,numentities_n*numentities_n*sizeof(double));
		 
	} 


 	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,numentities_n,numentities_n,probleminfo->q,1.0,cnorm,probleminfo->q,cnorm,probleminfo->q,0.0,entityinfo->a,numentities_n);
	 

	free(c0); 
	free(cnorm);
	free(ck);
	free(findjstar); 
	free(A); 
	free(B);
	free(svd);
	free(Aprime);
  	return 0;
}
  
 
    
 








































 








