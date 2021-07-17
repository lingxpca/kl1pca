#include <stdlib.h>
#include <stdio.h>
#include "type.h"

int readdata(IOINFOptr ioinfo, ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);

int readdata(IOINFOptr ioinfo, ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo)
{
  FILE *datafile = ioinfo->datafile;
  int status     = ioinfo->status;

  int i = probleminfo->i;
  int j = probleminfo->j;

  status = 0;
  /* read number of observations */
  fscanf  (datafile, "%d %d", &(entityinfo->numentities_n), &(entityinfo->numattributes_m));

  /* allocate memory for points at each step of projection, entities are in columns */
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
