#include <stdlib.h>
#include <stdio.h> 
#include <ctype.h>
#include "type.h"

int readdata(IOINFOptr ioinfo, ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);

void rc_find(FILE *fp,int* rows,int* cols)
{    
    *rows = 0;
    int i,j;
    *cols = 0; 
     
    while((i=fgetc(fp))!=EOF)
    {
            if (i == ' ') {
                ++j;
            } 
            else if (i == '\n') {
                (*rows)++; 
                *cols=j+1;
                j = 0;   
            }
    } 
    fclose(fp);
}

int readdata(IOINFOptr ioinfo, ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo)
{
  FILE *datafile = ioinfo->datafile;
  int status = ioinfo->status; 
  int i = probleminfo->i;
  int j = probleminfo->j;
  FILE *getrc = fopen(ioinfo->problem,"r");   
  rc_find(getrc,&(entityinfo->numentities_n),&(entityinfo->numattributes_m)); 
	
  status = 0;
  /* read number of observations */
  fscanf  (datafile, "%d %d", &(entityinfo->numentities_n), &(entityinfo->numattributes_m));
  entityinfo->points_XT = (double *) malloc (entityinfo->numentities_n*entityinfo->numattributes_m*sizeof (double)); 
   
  /* read in data */
  for (i = 0; i < entityinfo->numentities_n; ++i) {
    for (j = 0; j < entityinfo->numattributes_m; ++j) {
      status = fscanf (datafile, "%lf ", &(entityinfo->points_XT[entityinfo->numattributes_m*i+j]));
      if (status == EOF) {
        fprintf (stderr, "Data format incorrect.  Terminating ...\n");
	return 1;
      }
      else {
	status = 0;
      }
    }
  }
  

  fclose(datafile);

  return 0;
}
