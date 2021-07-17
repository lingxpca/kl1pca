#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "type.h" 
 
int pfkpca (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);  

int pfkpca (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo)
{
	int i,j,k,l,m,z; 
	double rho = entityinfo->rho; 
	int numattributes_m = entityinfo->numattributes_m;
	int numentities_n = entityinfo->numentities_n;
	double *points_XT = entityinfo->points_XT; 


	double *xbzx = entityinfo->XBZX; 
	double sum,sum1,sum2,sum3,sum4,sum5, A;
	double* alpha; 
	double I = 500.0;
	double h = 1.2;
	double alphastar; 
	double xhat[numattributes_m];
	int jstar=0;
	double prod;  
	double kone[numentities_n];
	double *kx;
	double *kx_tilde; 
	double **jkx;
	double **jkx_tilde;
	double f,p;
	
	FILE *myoutput=fopen("output.txt","w");
 			   				  
	 
  	kx=(double *)malloc(numentities_n*sizeof(double));
  	kx_tilde = (double *)malloc(numentities_n*sizeof(double));
 
  
  	jkx=(double **)calloc(numentities_n,sizeof(double *));
  	jkx_tilde=(double **)calloc(numentities_n,sizeof(double *));
	for (i=0;i<numentities_n;++i)  {
		jkx[i]=(double *)calloc(numattributes_m,sizeof(double));
		jkx_tilde[i]=(double *)calloc(numattributes_m,sizeof(double));
	}
	alpha =(double *)malloc(10000*sizeof(double)); 
 
  
   if (strcmp(entityinfo->kernel,"pol")==0) { 
		for (z=0;z<numentities_n;z++){  

			for (i=0;i<numentities_n;++i)	{ //******S1****** 
				sum1 = 0.0; 
				for (j =0;j<numentities_n;++j) 	{ 
					kone[j] = 0.0;
					kx[j] = 1.0;
					for (k=0;k<numentities_n;++k) kone[j] += entityinfo->KK[j][k];
					for (l=0;l<numattributes_m;++l) kx[j] += points_XT[z*numattributes_m+l]*points_XT[j*numattributes_m+l]; //kx
					sum1 += kx[j];
				} 
				sum2=0.0;
				for (k=0;k<numentities_n;++k) sum2 += kone[k];
				kx_tilde[i]=kx[i] - (kone[i]+sum1)/numentities_n+sum2/(numentities_n*numentities_n); //centerlized kx refer to Nojun 2013
			}

			for (i=0;i<numattributes_m;++i) { 
				sum1=0.0; 
				for (j = 0;j<numentities_n;++j)	 {  
					jkx[j][i] =  points_XT[j*numattributes_m+i]; //Jkx
					sum1 += jkx[j][i]/numentities_n;
				}
				for (j=0;j<numentities_n;++j) jkx_tilde[j][i]= jkx[j][i]-sum1;//centralized jkx  
  			}

  															//calculate f
			sum =0.0; 
			for(j=0;j<numentities_n;++j) {
				sum4 = 0.0;
				for (l=0;l<numentities_n;++l) sum4 += kx_tilde[l]*entityinfo->a[l*numentities_n+j];
				sum += sum4*kx_tilde[j];		
			} 
			
			sum5 = 0.0; 
			sum1 = 0.0;
			for (l=0;l<numentities_n;++l) {
				sum5 += kx[l];  
				sum1 += pow(points_XT[z*numattributes_m+l],2);  
			}
			f = 1+sum1-(2*sum5/numentities_n)+sum2/(numentities_n*numentities_n)-sum;   

 			if (f > pow(10,-9)) {	  //calculate steepest descent    
  				sum4=0.0;
				for (i=0;i<numattributes_m;++i){
					xbzx[i]=0.0;
					sum1 =sum3 = 0.0;  
					for(j=0;j<numentities_n;++j) {
						sum1 +=  jkx[j][i];
						sum2 = 0.0;
						for (l=0;l<numentities_n;++l) sum2 += jkx_tilde[l][i]*entityinfo->a[l*numentities_n+j];
						sum3 += sum2*kx_tilde[j];
					} 
					xbzx[i]= 2*sum1/numentities_n+2*sum3;	 
					sum4 += pow(2*sum1/numentities_n+2*sum3,2);	 
				}	
            	for (i = 0;i<numattributes_m;++i) xbzx[i] = xbzx[i]/pow(sum4,.5);  																																		
		 																			//******S2****** 	
				A=0.0;
				for (i=0;i<numentities_n;++i) {	
					if (i != z) {
					sum3=0.0;
					for (j=0;j<numattributes_m;++j) sum3 += (points_XT[i*numattributes_m+j]-points_XT[z*numattributes_m+j])*xbzx[j];
					if (sum3 > A) A = sum3;	}
				}  
	
				int r = 0;
				alpha[0] =0.0;
				while(alpha[r]<A) {
					alpha[r+1] = alpha[r] + A*pow(h,r)/(I-1.0); 
					r = r + 1;
				}; 
				alpha[r]=A; 
	
	 													//******S3*****calculate g(alpha);  																																															
  				alphastar = pow(2,36);
				for(i=1;i<r+1;++i) 	{	 
					for (m=0;m<numentities_n;++m)	{
						sum1 = 0.0; 
						for (j =0;j<numentities_n;++j)	{
							kx[j] = 1.0;
							kone[j] = 0.0; 
							for (k=0;k<numentities_n;++k) kone[j] += entityinfo->KK[j][k];
							for (l = 0;l<numattributes_m;++l) kx[j] += (points_XT[z*numattributes_m+l]+alpha[i]*xbzx[l])*points_XT[j*numattributes_m+l]; //kx
 							sum1 += kx[j]; 
						}
						sum2 = 0.0;
						for (k =0;k<numentities_n;++k) sum2 += kone[k];
						kx_tilde[m]=kx[m] - (kone[m]+sum1)/numentities_n+sum2/(numentities_n*numentities_n); //kx_tilde
					}
					
					sum3 = 0.0;
					for (l=0;l<numentities_n;++l) sum3 += kone[l]/(numentities_n*numentities_n);
					
					sum =0.0; 
					for(j=0;j<numentities_n;++j)	{
						sum4 = 0.0;
						for (l=0;l<numentities_n;++l) sum4 += kx_tilde[l]*entityinfo->a[l*numentities_n+j];
						sum += sum4*kx_tilde[j];		//kx_tilde' * v * kx_tilde	
					}
							
				sum5 =0.0;
				sum1 =0.0;
				for (l=0;l<numentities_n;++l) {
					sum5 += kx[l];
					sum1 += pow(points_XT[z*numattributes_m+l]+alpha[i]*xbzx[l],2);		
				}		
					if (1+sum1-(2*sum5/numentities_n)+sum2/pow(numentities_n,2)-sum < alphastar) alphastar = 1+sum1-(2*sum5/numentities_n)+sum2/pow(numentities_n,2)-sum;
					else {jstar = i-1; break;};
				}	
	
								/////////////////////////calculate g'(alphastar1)		
																														
			for (i=0;i<numentities_n;++i) {
				sum1=0.0;
				for (j =0;j<numentities_n;++j) 	{
					kone[j]=0.0;
					kx[j] =1.0;
					for (k=0;k<numentities_n;++k) kone[j] += entityinfo->KK[j][k];
					for (l=0;l<numattributes_m;++l) kx[j] += (points_XT[z*numattributes_m+l]+alpha[jstar]*xbzx[l])*points_XT[j*numattributes_m+l]; //kx
					sum1 += kx[j];
				}
				sum2=0.0;
				for (k=0;k<numentities_n;++k) sum2 += kone[k];
				kx_tilde[i]=kx[i] - (kone[i]+sum1)/numentities_n+sum2/(numentities_n*numentities_n); //kx_tilde
			}																													
																																
																														
			for (i=0;i<numattributes_m;++i)	{ 
	 			sum1=0.0; 
    			for (j = 0;j<numentities_n;++j)	  {  
					jkx[j][i] = points_XT[j*numattributes_m+i]; //Jkx
					sum1 += jkx[j][i]/numentities_n;
    			}
   			for (j=0;j<numentities_n;++j) jkx_tilde[j][i]= jkx[j][i]-sum1;//jkx_tilde
  			}
  
  			sum4=0.0;
			for (i=0;i<numattributes_m;++i)  {
				sum1 =0.0; 
				sum3 =0.0;
				for(j=0;j<numentities_n;++j)	{
					sum1 +=  jkx[j][i];
					sum2 = 0.0;
					for (l=0;l<numentities_n;++l) sum2 += jkx_tilde[l][i]*entityinfo->a[l*numentities_n+j];
					sum3 += sum2*kx_tilde[j];
	  			} 
				sum4 += 2*(sum1/numentities_n+sum3)*xbzx[i];
			}
			 		
			if ( sum4 < 0 ) {	
				for (l=0;l<numattributes_m;++l) xhat[l] = points_XT[z*numattributes_m+l]+0.5*(alpha[jstar]+alpha[jstar+1])*xbzx[l];
			}
			else if (sum4 == 0) {
				for (l=0;l<numattributes_m;++l) xhat[l] =points_XT[z*numattributes_m+l]+(alpha[jstar])*xbzx[l];
			}
			else if (sum4> 0) { 
				for (l=0;l<numattributes_m;++l) xhat[l] = points_XT[z*numattributes_m+l]+0.5*(alpha[jstar-1]+alpha[jstar])*xbzx[l];
			} 
			}
	
			else {
				for (l=0;l<numattributes_m;++l) xhat[l] = points_XT[z*numattributes_m+l];
			}

			for (l=0;l<numattributes_m;++l) fprintf(myoutput,"%g \t",xhat[l]);
			fprintf(myoutput,"\n"); 
	  }
	}


									//rbf
	else if (strcmp(entityinfo->kernel,"rbf")==0) {  
	p = exp(-1.0/(entityinfo->vox*rho));	
		for (z=0;z<numentities_n;z++){  	 //for each point																	 
											
			for (i=0;i<numentities_n;++i)	{ //******S1****** 
				sum1=0.0;
				for (j =0;j<numentities_n;++j) 	{
					kone[j]=0.0;
					kx[j] =1.0;
					for (k=0;k<numentities_n;++k) kone[j] += entityinfo->KK[j][k];
					for (l=0;l<numattributes_m;++l) kx[j] *= pow(p,pow(points_XT[z*numattributes_m+l]-points_XT[j*numattributes_m+l],2)); //kx
					sum1 += kx[j]; 
				}
				sum2=0.0;
				for (k=0;k<numentities_n;++k) sum2 += kone[k];
				kx_tilde[i]=kx[i] - (kone[i]+sum1)/numentities_n+sum2/(numentities_n*numentities_n); //kx_tilde  
			}

	 
			
			for (i=0;i<numattributes_m;++i) { 
				sum1=0.0; 
				for (j = 0;j<numentities_n;++j)	 { 
					jkx[j][i]=0.0;
					prod = 1.0;
					for (l=0;l<numattributes_m;++l) prod *= pow(p,pow(points_XT[z*numattributes_m+l]-points_XT[j*numattributes_m+l],2));//kx
					jkx[j][i] = prod *2*log(p)*(points_XT[z*numattributes_m+i]-points_XT[j*numattributes_m+i]); //Jkx
					sum1 += jkx[j][i]/numentities_n;
				}
				for (j=0;j<numentities_n;++j) jkx_tilde[j][i]= jkx[j][i]-sum1;//jkx_tilde  
			}

						//calculate f
			sum =0.0; 
			for(j=0;j<numentities_n;++j)	{
				sum4 = 0.0;
				for (l=0;l<numentities_n;++l) sum4 += kx_tilde[l]*entityinfo->a[l*numentities_n+j];
				sum += sum4*kx_tilde[j];		//kx_tilde' * v * kx_tilde	
			} 
				
			sum5 =0.0;
			for (l=0;l<numentities_n;++l) sum5 += kx[l];	
			f = 1-(2*sum5/numentities_n)+sum2/(numentities_n*numentities_n)-sum; 

 

			if (f > pow(10,-9)) { //calculate steepest descent 
				sum4=0.0;
				for (i=0;i<numattributes_m;++i){
					xbzx[i]= 0.0;
					sum1 =sum3 = 0.0;  
					for(j=0;j<numentities_n;++j)	{
						sum1 +=  jkx[j][i];
						sum2 = 0.0;
						for (l=0;l<numentities_n;++l) sum2 += jkx_tilde[l][i]*entityinfo->a[l*numentities_n+j];
						sum3 += sum2*kx_tilde[j]; 
					} 
					xbzx[i]= 2*sum1/numentities_n+2*sum3;	 
					sum4 += pow(2*sum1/numentities_n+2*sum3,2);	 
				}	
				for (i = 0;i<numattributes_m;++i) xbzx[i] = xbzx[i]/pow(sum4,.5);  
											//******S2****** 	
			
	 			A=0.0;
				for (i=0;i<numentities_n;++i) {	
					if (i != z) {
					sum3=0.0;
					for (j=0;j<numattributes_m;++j) sum3 += (points_XT[i*numattributes_m+j]-points_XT[z*numattributes_m+j])*xbzx[j];
 
				   if (sum3 > A) A = sum3;}
				}  
	
				int r = 0;
				alpha[0] =0.0;
				while(alpha[r]<A) {
					alpha[r+1] = alpha[r] + A*pow(h,r)/(I-1.0); 
					r = r + 1;
				}; 
				alpha[r]=A;  				 	
				                   //******S3******  	//calculate g(alpha);  		 														
  				alphastar = pow(2,36);
				for(i=1;i<r+1;++i) 	{	 
					for (m=0;m<numentities_n;++m)	{
						sum1=0.0;
						for (j =0;j<numentities_n;++j)	{
							kx[j] =  1.0;
							kone[j] = 0.0;
							for (k=0;k<numentities_n;++k) kone[j] += entityinfo->KK[j][k];
							for (l = 0;l<numattributes_m;++l) kx[j] *= pow(p,pow(points_XT[z*numattributes_m+l]+alpha[i]*xbzx[l]-points_XT[j*numattributes_m+l],2)); //kx
							sum1 += kx[j];
						}
						sum2 = 0.0;
						for (k =0;k<numentities_n;++k) sum2 += kone[k];
						kx_tilde[m]=kx[m] - (kone[m]-sum1)/numentities_n+sum2/(numentities_n*numentities_n); //kx_tilde
					}
					
					sum3 = 0.0;
					for (l=0;l<numentities_n;++l) sum3 += kone[l]/(numentities_n*numentities_n);
					
					sum =0.0; 
					for(j=0;j<numentities_n;++j)	{
						sum4 = 0.0;
						for (l=0;l<numentities_n;++l) sum4 += kx_tilde[l]*entityinfo->a[l*numentities_n+j];
						sum += sum4*kx_tilde[j];		//kx_tilde' * v * kx_tilde	
					}
							
					sum5 =0.0;
					for (l=0;l<numentities_n;++l) sum5 += kx[l];	
					
					if (1-(2*sum5/numentities_n)+sum2/pow(numentities_n,2)-sum < alphastar) alphastar = 1-(2*sum5/numentities_n)+sum2/pow(numentities_n,2)-sum;
					else {jstar = i-1; break;};
				}	
	
	  									//calculate g'(alphastar1)		
																														
				for (i=0;i<numentities_n;++i) 	{
					sum1=0.0;
					for (j =0;j<numentities_n;++j)	{
						kone[j]=0.0;
						kx[j] =1.0;
						for (k=0;k<numentities_n;++k) kone[j] += entityinfo->KK[j][k];
						for (l=0;l<numattributes_m;++l) kx[j] *= pow(p,pow(points_XT[z*numattributes_m+l]+alpha[jstar]*xbzx[l]-points_XT[j*numattributes_m+l],2)); //kx
						sum1 += kx[j];
					}
					sum2=0.0;
					for (k=0;k<numentities_n;++k) sum2 += kone[k];
					kx_tilde[i]=kx[i] - (kone[i]-sum1)/numentities_n+sum2/(numentities_n*numentities_n); //kx_tilde
				}																													
																																		
																														
				for (i=0;i<numattributes_m;++i)  { 
					sum1=0.0; 
					for (j = 0;j<numentities_n;++j)	   { 
						jkx[j][i]=0.0;
						prod = 1.0;
						for (l=0;l<numattributes_m;++l) prod *= pow(p,pow(points_XT[z*numattributes_m+l]+alpha[jstar]*xbzx[l]-points_XT[j*numattributes_m+l],2));//kx
						jkx[j][i] = prod *2*log(p)*(points_XT[z*numattributes_m+i]+alpha[jstar]*xbzx[i]-points_XT[j*numattributes_m+i]); //Jkx
						sum1 += jkx[j][i]/numentities_n;
					}
					for (j=0;j<numentities_n;++j) jkx_tilde[j][i]= jkx[j][i]-sum1;//jkx_tilde
				}
  
  				sum4=0.0;
				for (i=0;i<numattributes_m;++i) {
					sum1 = sum3 = 0.0; 
					for(j=0;j<numentities_n;++j) {
						sum1 +=  jkx[j][i];
						sum2 = 0.0;
						for (l=0;l<numentities_n;++l) sum2 += jkx_tilde[l][i]*entityinfo->a[l*numentities_n+j];
						sum3 += sum2*kx_tilde[j];
					}	 
					sum4 += -2*(sum1/numentities_n+sum3)*xbzx[i];
				}
							
				if ( sum4 < 0 ) {	
					for (l=0;l<numattributes_m;++l) xhat[l] = points_XT[z*numattributes_m+l]+0.5*(alpha[jstar]+alpha[jstar+1])*xbzx[l];} 
				else if (sum4 == 0) {
					for (l=0;l<numattributes_m;++l) xhat[l] =points_XT[z*numattributes_m+l]+(alpha[jstar])*xbzx[l];	}
				else if (sum4> 0) { 
					for (l=0;l<numattributes_m;++l) xhat[l] = points_XT[z*numattributes_m+l]+0.5*(alpha[jstar-1]+alpha[jstar])*xbzx[l]; } 
			}

			else {
				for (l=0;l<numattributes_m;++l) xhat[l] = points_XT[z*numattributes_m+l];
			}
 		 
			for (l=0;l<numattributes_m;++l) fprintf(myoutput,"%g \t",xhat[l]);
			fprintf(myoutput,"\n"); 
		}
	}
	fclose(myoutput);
	return 0;
}


































 








