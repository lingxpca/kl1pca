#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h> 
#include "type.h" 
#include <mkl.h> 
 
int pfkpca (IOINFOptr ioinfo,ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);  

int pfkpca (IOINFOptr ioinfo,ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo)
{
	int i,j,l,m,z;  
	int numattributes_m = entityinfo->numattributes_m;
	int numentities_n = entityinfo->numentities_n;
	double *points_XT = entityinfo->points_XT;  

	double *xbzx;  
	double sum,sum1,sum2,sum4, A;
	double* alpha;  
	double h = 1.2;
	double alphastar; 
	double xhat[numattributes_m];
	int jstar=0; 
	double *kone;
	double *kx;
	double *kx_tilde; 
	double **jkx;
	double **jkx_tilde;

	double f,p;
	double konesum =0.0;
	char output[40];
	char err[30]; 
	//int ci,k,sig;
	
	snprintf(output,sizeof(output),"./%s_%s%s_%d_%g",ioinfo->problem,entityinfo->kernel,entityinfo->norm,probleminfo->q,probleminfo->var);  
	FILE *myoutput=fopen(output,"w"); 

   	snprintf(err,sizeof(err),"./err_%s_%s%s_%d",ioinfo->problem,entityinfo->kernel,entityinfo->norm,probleminfo->q);  
	FILE *myerr=fopen(err,"w");   	
	
  	kx=(double *)malloc(numentities_n*sizeof(double));
	kone = (double *)calloc(numentities_n,sizeof(double));
  	kx_tilde = (double *)calloc(numentities_n,sizeof(double));
	xbzx =  (double *)calloc(numattributes_m,sizeof(double));
  
  	jkx=(double **)calloc(numentities_n,sizeof(double *)); 
	jkx_tilde=(double **)calloc(numentities_n,sizeof(double *));
	for (i=0;i<numentities_n;++i)  {
		jkx[i]=(double *)calloc(numattributes_m,sizeof(double));
		jkx_tilde[i]=(double *)calloc(numattributes_m,sizeof(double));
	}
	
	double *Y,*jt,*C,*zmiusx,*jkxcolumsum,*maxA,*zmiusxalpha,*gprime,f_exh,gp_exh;
	jt= (double *)calloc((entityinfo->numentities_n*entityinfo->numattributes_m),sizeof(double));
	C=  (double *)calloc((entityinfo->numentities_n*entityinfo->numattributes_m),sizeof(double));
	Y =  (double *)calloc(numentities_n,sizeof(double));  
	zmiusx = (double *)calloc((entityinfo->numentities_n*entityinfo->numattributes_m),sizeof(double));
	zmiusxalpha = (double *)calloc((entityinfo->numentities_n*entityinfo->numattributes_m),sizeof(double));
	alpha =(double *)malloc(600*sizeof(double)); 
	jkxcolumsum = (double *)calloc(numattributes_m,sizeof(double));  
	maxA =  (double *)calloc(numentities_n,sizeof(double));
	gprime = (double *)calloc(numattributes_m,sizeof(double));     
/*////////////////////////////////////////////////////////////////Gaussian Kernel////////////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////Gaussian Kernel/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////Gaussian Kernel////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////Gaussian Kernel////////////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////Gaussian Kernel////////////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////Gaussian Kernel/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////Gaussian Kernel////////////////////////////////////////////////////////////*/
		p = exp(-1.0/((probleminfo->var)));	 //*entityinfo->var
		konesum=0.0;
		for(i=0;i<numentities_n;++i){
			kone[i] =0.0;
			for (j=0;j<numentities_n;++j) kone[i] += entityinfo->KK[i][j];
			konesum += kone[i]; 
		}

		for (z=0;z<numentities_n;z++){   		
			
			for (i=0;i<numentities_n;++i)  
				for (j = 0;j<numattributes_m;++j)  
					zmiusx[i*numattributes_m+j] = points_XT[z*numattributes_m+j]-points_XT[i*numattributes_m+j];
			 

			for (j =0,sum2=0.0;j<numentities_n;++j) 	
			{  
				kx[j]=1.0;
 				for (l=0;l<numattributes_m;++l) kx[j] *= pow(p,zmiusx[j*numattributes_m+l]*zmiusx[j*numattributes_m+l]); //Kz
				sum2 += kx[j]; 
			} 
			
			for (i=0;i<numentities_n;++i)	kx_tilde[i]=kx[i] - (kone[i]+sum2)/numentities_n+konesum/(numentities_n*numentities_n); //centralized kz  
			
			for (i=0;i<numattributes_m;++i) 
			{  
				jkxcolumsum[i] =0.0;
				for (j = 0,sum1=0.0;j<numentities_n;++j)	 
				{   
 					jkx[j][i] = kx[j] *2*log(p)*(zmiusx[j*numattributes_m+i]); //Kz deriviative wrt z
					jkxcolumsum[i] += jkx[j][i]; 
				}
				for (j=0;j<numentities_n;++j) jt[j*numattributes_m+i]=jkx[j][i]-jkxcolumsum[i]/numentities_n;
			}
		 
			for (i=0;i<numentities_n;++i) Y[i] = 0.0;
			cblas_dgemv(CblasRowMajor,CblasNoTrans,numentities_n,numentities_n,1.0,entityinfo->a,numentities_n,kx_tilde,1,0.0,Y,1); 
         	sum = cblas_ddot(numentities_n,kx_tilde,1,Y,1);	//kx_tilde' * v * kx_tilde	
			
		 
			f = 1-(2*sum2/numentities_n)+konesum/(numentities_n*numentities_n)-sum;
			fprintf(myerr,"%g %g %g\n",points_XT[z*numattributes_m+0],points_XT[z*numattributes_m+1],f); 

			

			// if ( f > 10e-9) //calculate steepest descent 
			// {  			 				
			// 	memcpy(xbzx,jkxcolumsum,sizeof(double)*numattributes_m);			
			// 	cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,numattributes_m,numentities_n,numentities_n,1.0,jt,numattributes_m,entityinfo->a,numentities_n,0.0,C,numentities_n);
			// 	cblas_dgemv(CblasRowMajor,CblasNoTrans,numattributes_m,numentities_n,2.0,C,numentities_n,kx_tilde,1,2.0/numentities_n,xbzx,1); //xbzx
				
			// 	sum4 = cblas_dnrm2(numattributes_m,xbzx,1);
			// 	for (i = 0;i<numattributes_m;++i) xbzx[i] = xbzx[i]/sum4;    
								
			// 	cblas_dgemv(CblasRowMajor,CblasNoTrans,numentities_n,numattributes_m,1.0,points_XT,numattributes_m,xbzx,1,0.0,maxA,1);
			// 	A = maxA[0];
			// 	for (i=0;i<numentities_n;++i) if (maxA[i]> A && i !=z) A = maxA[i];
			// 	A = A-maxA[z];   	 
				 

			// 	for (l=0;l<numattributes_m;++l) xhat[l] = points_XT[z*numattributes_m+l];
			 
			// 	for (m=0;m<100;++m){
			// 		f_exh = 1000000.0; 
			// 		for (ci=0;ci<2*(numattributes_m);++ci) // 2*(numattributes_m-1)
			// 		{					
			// 			memcpy(xbzx,xhat,sizeof(double)*numattributes_m);  
			// 			if (ci > numattributes_m) { //numattributes_m-1
			// 			xbzx[ci-(numattributes_m)] = xhat[ci-(numattributes_m)] + A*h/100.0;} //numattributes_m-1
			// 			else{ xbzx[ci]= xhat[ci] - A*h/100.0;}

			// 			for (i=0;i<numentities_n;++i)  
			// 				for (j = 0;j<numattributes_m;++j)  
			// 					zmiusx[i*numattributes_m+j] = xbzx[j]-points_XT[i*numattributes_m+j];


			// 			for (j =0,sum2=0.0;j<numentities_n;++j) {  
			// 				kx[j]=1.0;
			// 				for (l=0;l<numattributes_m;++l) kx[j] *= pow(p,zmiusx[j*numattributes_m+l]*zmiusx[j*numattributes_m+l]); //Kz
			// 				sum2 += kx[j]; 
			// 			}  
			// 			for (i=0;i<numentities_n;++i) kx_tilde[i]= kx[i] - (kone[i]+sum2)/numentities_n+konesum/(numentities_n*numentities_n); //centralized kz  

			// 			for (i=0;i<numentities_n;++i) Y[i] = 0.0;
			// 			cblas_dgemv(CblasRowMajor,CblasNoTrans,numentities_n,numentities_n,1.0,entityinfo->a,numentities_n,kx_tilde,1,0.0,Y,1); 
			// 			sum = cblas_ddot(numentities_n,kx_tilde,1,Y,1);	 
						
			// 		 	if (1-(2*sum2/numentities_n)+konesum/(numentities_n*numentities_n)-sum < f_exh) {
			// 				f_exh=1-(2*sum2/numentities_n)+konesum/(numentities_n*numentities_n)-sum; 
			// 				k = (ci > numattributes_m-1)?(ci-(numattributes_m)):ci;  //k = (ci > numattributes_m-2)?(ci-(numattributes_m-1)):ci
			// 				sig = (ci > numattributes_m-1)?1:-1; //sig = (ci > numattributes_m-2)?1:-1;
			// 			}  

			// 		} 
			// 		xhat[k]  = xhat[k] + sig*A*h/100.0;
			// 	} 
			// } 
			
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
				 
				int r = 0; 
				alpha[r] =0.0;  
				while(alpha[r]< A ) { 		 
					alpha[r+1] = alpha[r] + A*h/499.0; 
					r = r + 1;
				};  
				alpha[r]=A;

     			alphastar=f;  
				for (jstar=1;jstar<r+1;++jstar) {
					
					for (l= 0;l<numentities_n;++l){
						for (j=0;j<numattributes_m;++j) 
							zmiusxalpha[l*numattributes_m+j] = zmiusx[l*numattributes_m+j]+alpha[jstar]*xbzx[j];}

					sum1 =0.0;
					for (j=0;j<numentities_n;++j){
						kx[j]=1.0;
						for (l=0;l<numattributes_m;++l) kx[j] *= pow(p,zmiusxalpha[j*numattributes_m+l]*zmiusxalpha[j*numattributes_m+l]);
						sum1 += kx[j];
					} 

					for (m=0;m<numentities_n;++m)	kx_tilde[m]=kx[m] - (kone[m]-sum1)/numentities_n+konesum/(numentities_n*numentities_n); //kx_tilde

					memset(Y,0,sizeof(double)*numentities_n);
		 			cblas_dgemv(CblasRowMajor,CblasNoTrans,numentities_n,numentities_n,1.0,entityinfo->a,numentities_n,kx_tilde,1,1.0,Y,1);
         			sum = cblas_ddot(numentities_n,kx_tilde,1,Y,1);	//kx_tilde' * v * kx_tilde	 
					
					f_exh = (1-(2*sum1/numentities_n)+konesum/(numentities_n*numentities_n)-sum); 
					
					if (f_exh > alphastar)
					{  
						if(jstar>1){
							for (l= 0;l<numentities_n;++l)
								for (j=0;j<numattributes_m;++j) 
									zmiusxalpha[l*numattributes_m+j] = zmiusx[l*numattributes_m+j]+alpha[jstar-1]*xbzx[j];

							for (j =0,sum1=0.0;j<numentities_n;++j) { 
									kx[j] =1.0; 
								for (l=0;l<numattributes_m;++l) kx[j] *= pow(p,zmiusxalpha[j*numattributes_m+l]*zmiusxalpha[j*numattributes_m+l]); //kx
									sum1 += kx[j];
							}

							for (j=0;j<numentities_n;++j)	kx_tilde[j]=kx[j] - (kone[j]-sum1)/numentities_n+konesum/(numentities_n*numentities_n); //kx_tilde 

							for (m=0;m<numattributes_m;++m)  { 		
								jkxcolumsum[m]=0.0;		
								for (j = 0;j<numentities_n;++j)	 {  
									jkx[j][m] =   kx[j] *2*log(p)*(zmiusxalpha[j*numattributes_m+m]); //Jkx
									jkxcolumsum[m] += jkx[j][m];
								} 
								for (j=0;j<numentities_n;++j) jt[j*numattributes_m+m]= jkx[j][m]-jkxcolumsum[m]/numentities_n;//jkx_tilde
							}
							
							memcpy(gprime,jkxcolumsum,sizeof(double)*numattributes_m);	
							memset(C,0,numattributes_m*numentities_n*sizeof(double));
							cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,numattributes_m,numentities_n,numentities_n,1.0,jt,numattributes_m,entityinfo->a,numentities_n,0.0,C,numentities_n);
							cblas_dgemv(CblasRowMajor,CblasNoTrans,numattributes_m,numentities_n,-2.0,C,numentities_n,kx_tilde,1,-2.0/numentities_n,gprime,1); //xbzx 
							gp_exh = cblas_ddot(numattributes_m,gprime,1,xbzx,1); 
							jstar =(gp_exh > 0.0)? jstar-1:jstar; 		
							break;					
						}   	
						
					} 
					alphastar = f_exh;
				}
				
	  										
				if (jstar > 1) { 
					for (l=0;l<numattributes_m;++l) xhat[l] = points_XT[z*numattributes_m+l]+0.5*(alpha[jstar-1]+alpha[jstar])*xbzx[l];
				}
				else  {  					
					for (l=0;l<numattributes_m;++l) xhat[l] =points_XT[z*numattributes_m+l]+(alpha[jstar])*xbzx[l];	
				}	

			}  
			 
			else 
			{  
				for (l=0;l<numattributes_m;++l) xhat[l] = points_XT[z*numattributes_m+l];
			}
			
		 	for (l=0;l<numattributes_m;++l) fprintf(myoutput,"%g \t",xhat[l]);
			fprintf(myoutput,"\n");  
		
		}		
	
 


	fclose(myoutput); 
	
	free(Y);
	free(jt);
	free(C);
	free(zmiusx);
	free(jkxcolumsum);
	free(maxA);
	free(zmiusxalpha);
	free(points_XT);
	free(kx);
	free(kone);
	free(kx_tilde);
	free(xbzx);
	free(jkx);
	free(jkx_tilde);
	free(gprime); 
	return 0;
}


































 








