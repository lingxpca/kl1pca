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
	double **cstar; 
	double *cnorm;  
	double *ck;   
    
	c0 = (double *)calloc(numentities_n,sizeof(double)); 
	cstar=(double **)calloc(numentities_n,sizeof(double *));  
	
	cnorm = (double *)calloc(numentities_n*probleminfo->q,sizeof(double)); 
	ck = (double *)calloc(numentities_n,sizeof(double));  
	findjstar = (double *)calloc(numentities_n,sizeof(double));  

	for (i=0;i<numentities_n;++i){  
		cstar[i]=(double *)calloc(probleminfo->q,sizeof(double)); 
	} 

	double *svd = (double *)calloc(numentities_n*numentities_n,sizeof(double)); 	
	double *A=(double *)calloc(numentities_n*numentities_n,sizeof(double)); 
	double *Y = (double *)calloc(numentities_n,sizeof(double)); 
	double *kc = (double *)calloc(numentities_n,sizeof(double)); 
	double *c = (double *)calloc(numentities_n,sizeof(double)); 
	
	for (i=0;i<numentities_n;++i) A[i*numentities_n+i] = 1.0;
   for (l=0;l<probleminfo->q;++l)	{
	   
	   	for (i=0;i<numentities_n;++i){
			for (j=0;j<numentities_n;++j){
				svd[i*numentities_n+j] =  entityinfo->K_tilte[l][i][j];
			}
		} 
	  	
		for (j=0;j<numentities_n;++j) {
			findjstar[j]=0.0;   
			for (i=0;i<numentities_n;++i) findjstar[j] += fabs(svd[i*numentities_n+j])/sqrt(svd[j*numentities_n+j]);   
		}
		jstar = cblas_idamax(numentities_n,findjstar,1); 
	   	for (i=0;i<numentities_n;++i) c0[i]= (svd[i*numentities_n+jstar] > 0.0) ? 1.0 : -1.0;  
		

	 	do	{			
			memset(Y,0,numentities_n*sizeof(double));
			cblas_dgemv(CblasRowMajor,CblasNoTrans,numentities_n,numentities_n,1.0,svd,numentities_n,c0,1,0.0,Y,1); 
			 
			for (i=0;i<numentities_n;++i)   cstar[i][l] = (Y[i] > 0.0) ? 1.0 : -1.0;
			for (i=0;i<numentities_n;++i) ck[i] =c0[i] - cstar [i][l]; 
			memset(Y,0,numentities_n*sizeof(double));
			cblas_dsymv(CblasRowMajor,CblasLower,numentities_n,1.0,svd,numentities_n,ck,1,0.0,Y,1);
			s1 = cblas_ddot(numentities_n,Y,1,ck,1);			
		} while(fabs(s1) >10e-7);     
		
		
		for (i=0;i<numentities_n;++i) c[i] = cstar[i][l]; 	 
		memset(kc,0,numentities_n*sizeof(double));
		cblas_dsymv(CblasRowMajor,CblasLower,numentities_n,1.0,svd,numentities_n,c,1,0.0,kc,1);
		s1 = cblas_ddot(numentities_n,kc,1,c,1); 
		for (i=0;i<numentities_n;++i) cnorm[i*probleminfo->q+l] = c[i]/sqrt(s1);
		 
	   	if (l == (probleminfo->q-1)) break; 
		loop(entityinfo->K_tilte[l+1][i][j]  = entityinfo->K_tilte[l][i][j] - (kc[i]*kc[j])/s1;)  
	} 
 	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,numentities_n,numentities_n,probleminfo->q,1.0,cnorm,probleminfo->q,cnorm,probleminfo->q,0.0,entityinfo->a,numentities_n);
	 

	free(c0);
	free(cstar);
	free(cnorm);
	free(ck);
	free(findjstar);
	free(entityinfo->K_tilte);
  	return 0;
}
  
 
    
 
{
	if ( f > 10e-9) //calculate steepest descent 
			{  			 				
				memcpy(xbzx,jkxcolumsum,sizeof(double)*numattributes_m);			
				cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,numattributes_m,numentities_n,numentities_n,1.0,jt,numattributes_m,entityinfo->a,numentities_n,0.0,C,numentities_n);
				cblas_dgemv(CblasRowMajor,CblasNoTrans,numattributes_m,numentities_n,2.0,C,numentities_n,kx_tilde,1,2.0/numentities_n,xbzx,1); //xbzx
				
				sum4 = cblas_dnrm2(numattributes_m,xbzx,1);
				for (i = 0;i<numattributes_m;++i) xbzx[i] = xbzx[i]/sum4;    
								
				cblas_dgemv(CblasRowMajor,CblasNoTrans,numentities_n,numattributes_m,1.0,points_XT,numattributes_m,xbzx,1,0.0,maxA,1);
				A = maxA[0];
				for (i=0;i<numentities_n;++i) if (maxA[i]> A && i !=z) A = maxA[i];
				A = A-maxA[z];   	 
				 

				for (l=0;l<numattributes_m;++l) xhat[l] = points_XT[z*numattributes_m+l];
			 
				for (m=0;m<100;++m){
					f_exh = 1000000.0; 
					for (ci=0;ci<2*(numattributes_m);++ci) // 2*(numattributes_m-1)
					{					
						memcpy(xbzx,xhat,sizeof(double)*numattributes_m);  
						if (ci > numattributes_m) { //numattributes_m-1
						xbzx[ci-(numattributes_m)] = xhat[ci-(numattributes_m)] + A*h/400.0;} //numattributes_m-1
						else{ xbzx[ci]= xhat[ci] - A*h/400.0;}

						for (i=0;i<numentities_n;++i)  
							for (j = 0;j<numattributes_m;++j)  
								zmiusx[i*numattributes_m+j] = xbzx[j]-points_XT[i*numattributes_m+j];


						for (j =0,sum2=0.0;j<numentities_n;++j) {  
							kx[j]=1.0;
							for (l=0;l<numattributes_m;++l) kx[j] *= pow(p,zmiusx[j*numattributes_m+l]*zmiusx[j*numattributes_m+l]); //Kz
							sum2 += kx[j]; 
						}  
						for (i=0;i<numentities_n;++i) kx_tilde[i]= kx[i] - (kone[i]+sum2)/numentities_n+konesum/(numentities_n*numentities_n); //centralized kz  

						for (i=0;i<numentities_n;++i) Y[i] = 0.0;
						cblas_dgemv(CblasRowMajor,CblasNoTrans,numentities_n,numentities_n,1.0,entityinfo->a,numentities_n,kx_tilde,1,0.0,Y,1); 
						sum = cblas_ddot(numentities_n,kx_tilde,1,Y,1);	 
						
					 	if (1-(2*sum2/numentities_n)+konesum/(numentities_n*numentities_n)-sum < f_exh) {
							f_exh=1-(2*sum2/numentities_n)+konesum/(numentities_n*numentities_n)-sum; 
							k = (ci > numattributes_m-1)?(ci-(numattributes_m)):ci;  //k = (ci > numattributes_m-2)?(ci-(numattributes_m-1)):ci
							sig = (ci > numattributes_m-1)?1:-1; //sig = (ci > numattributes_m-2)?1:-1;
						}  

					} 
					xhat[k]  = xhat[k] + sig*A*h/400.0;
				} 
			}
}







































 








