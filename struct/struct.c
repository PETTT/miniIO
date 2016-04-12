#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <mpi.h>
#include "open-simplex-noise.h"

/* #include <limits.h> */
/* #include <assert.h> */

void print_usage(int rank, const char *errstr)
{
  if(rank != 0)  return;

  if(errstr)
    fprintf(stderr, "%s\n\n", errstr);
    fprintf(stderr, "Usage: mpi_launcher [-n|-np NPROCS] ./struct "
	    "\n");

    /*## Add Output Modules' Usage String ##*/

    /*## End of Output Module Usage Strings ##*/
}


int main(int argc, char **argv)
{

  /* MPI vars */
  int rank, nprocs; 



  /* Init MPI  */
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  printf("Hi: %d of %d\n", rank, nprocs);

  MPI_Finalize();

  return 0;
}
