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
  double noisespacefreq = 10;    /* Spatial frequency of noise */
  double noisetimefreq = 0.25;    /* Temporal frequency of noise */
  int ni, nj, nk;
  int i, numtask, point_id, tmp_id;
  int x_id, y_id, z_id;
  int xy_dim,  x_dim;
  float deltax, deltay, deltaz;
  int numPoints;
  float data;
  float height;
  int height_id;
  int mask;
  float t=0.0;
  float mask_thres;
  int mask_thres_id;
  struct osn_context *simpnoise;    /* Open simplex noise context */
  
  /* MPI vars */
  int rank, nprocs;
  int cni, cnj, cnk;   /* Points in this task */
  int is, js, ks;     /* Global index starting points */
  float xs, ys, zs;    /* Global coordinate starting points */  /* Init MPI  */  /* Init MPI  */

  ni = nk = nj = 2;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  deltax = 1.f/(ni-1);
  deltay = 1.f/(nj-1);
  deltaz = 1.f/(nk-1);

  
  numPoints = ni*nj*nk;
  numtask =  numPoints/nprocs;


  xy_dim = ni * nj;
  x_dim = ni;

  mask_thres_id = (int)  mask_thres/deltaz;
   
  /* Set up osn */
  open_simplex_noise(12345, &simpnoise);   /* Fixed seed, for now */

  
  for (i=0; i<numtask; i++)
    {
      point_id = (rank*numtask)+i;

      /* Calculate point ids */
      z_id =  point_id / xy_dim;
      tmp_id =  point_id % xy_dim;
      y_id = tmp_id/x_dim ;
      x_id = tmp_id%x_dim ;

      /* calculate grid positions */
      xs = x_id * deltax;
      ys = y_id * deltay;
      zs = z_id * deltaz;

      /* later implement smarter way that requires only the first xy to calulate height something like the following  */
      /* if ( z_id == 0) */
      /* { */
      /*   /\* calculate height filed *\/ */
      /*   height =  (float)open_simplex_noise2(simpnoise, xs*noisespacefreq, ys*noisespacefreq); */
      /*   printf("Height: %f\n", height); */

      /* } */

      height =  (float)open_simplex_noise2(simpnoise, xs*noisespacefreq, ys*noisespacefreq);
      data = (float)open_simplex_noise4(simpnoise, xs*noisespacefreq, ys*noisespacefreq, zs*noisespacefreq, t*noisetimefreq);

      height_id = (int) height/deltaz;
      
      /* Calculate mask values */
      if (height_id == z_id )
      	{
	   if (z_id > mask_thres_id)
	    {
	      
	      mask = 4;  /* land surface above sea level*/
	    }
	   else if (z_id == mask_thres_id)
	    {
	      
	      mask = 3;  /* land surface at sea level*/
	    }
	    else if (z_id < mask_thres_id)
	    {
	      
	      mask = 2;  /* land surface below sea level*/
	    }
	}
      else if (height_id > z_id )
      	{
      	  mask = 1;  /* below land surface */
      	}
       else if (height_id < z_id )
      	{
	   if (z_id > mask_thres_id)
	    {
	      
	      mask = 4;  /* Atmosphere*/
	    }
	   else if (z_id == mask_thres_id)
	    {
	      
	      mask = 0;  /* ocean at sea level*/
	    }
	    else if (z_id < mask_thres_id)
	    {
	      
	      mask = -1;  /* ocean below sea level*/
	    }
      	}

      printf("++++++++++++++++++++++++++++++++++++++++++++\n");      
      printf("rank=%d: %d of %d:  point_id=%d\n", rank, rank+1, nprocs, point_id+1);
      printf("Point_id: (%d, %d, %d) = %f\n", x_id, y_id, z_id, data);
      printf("Point_pos: (%f, %f, %f) = %f\n", xs, ys, zs, data);
      printf("Height: %f mask=%d\n", height, mask);
    }
  
  open_simplex_noise_free(simpnoise);
  MPI_Finalize();

  return 0;
}
