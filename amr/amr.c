#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <mpi.h>
#include "open-simplex-noise.h"

/* #include <limits.h> */
/* #include <assert.h> */

void refine(int t, float rpId, float thres, float value, int level, float x, float y, float z, float dx, float dy, float dz, struct osn_context *osn);
void print_usage(int rank, const char *errstr);


int main(int argc, char **argv)
{

  int i, numtask, block_id, realblock_id, tmp_id;
  int x_id, y_id, z_id;
  int xy_dim,  x_dim;
  float deltax, deltay, deltaz;
  float deltax_center, deltay_center, deltaz_center;

  float threshold=0.0;
  
  int ni = 0;      /* Global grid size */
  int nj = 0;
  int nk = 0;

  int unit_val = 2;
  
  int nx = unit_val;
  int ny = unit_val;
  int nz = unit_val;
  
  int numBlocks;
  struct osn_context *simpnoise;    /* Open simplex noise context */
  
  /* MPI vars */
  int rank, nprocs; 
  int cni, cnj, cnk;   /* Points in this task */
  int is, js, ks;     /* Global index starting points */
  float xs, ys, zs;    /* Global coordinate starting points */  /* Init MPI  */
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  ni = nx + 1;
  nj = ny + 1;
  nk = nz + 1;
  
  numBlocks = nx*ny*nz;
  numtask =  numBlocks/nprocs;
  
  
  deltax = 1.f/(ni-1);
  deltay = 1.f/(nj-1);
  deltaz = 1.f/(nk-1);

  deltax_center = deltax/2;
  deltaz_center = deltay/2;
  deltaz_center = deltaz/2;

  /* cni = ni / inp; */
  /* cnj = nj / jnp; */
  /* cnk = nk / knp; */

  xy_dim = nx * ny;
  x_dim = nx;
  printf("Hi: xy_dim=%d xy_dim=%d\n", xy_dim, x_dim);

  /* Set up osn */
  open_simplex_noise(12345, &simpnoise);   /* Fixed seed, for now */
  for (i=0; i<numtask; i++)
  {
    block_id = (rank*numtask)+i;
    /* realblock_id = block_id - 1; */

    /* Calculate id for lower-left-front corner of blck */
    /* z_id =  realblock_id / xy_dim; */
    z_id =  block_id / xy_dim;
    /* tmp_id = realblock_id - (k_id * xy_dim); */
    /* tmp_id =  realblock_id % xy_dim; */
    tmp_id =  block_id % xy_dim;
    y_id = tmp_id/x_dim ;
    x_id = tmp_id%x_dim ;

    xs = x_id * deltax;
    ys = y_id * deltay;
    zs = z_id * deltaz;
    
    printf("Hi: rank=%d: %d of %d:  block_id=%d\n", rank, rank+1, nprocs, block_id+1);
    printf("Point_id: (%d, %d, %d) \n", x_id, y_id, z_id);
    refine(0, (block_id+1)/10.0, threshold, .1 + threshold, 0, xs, ys, zs, deltax, deltay, deltaz, simpnoise);
    
  }

  open_simplex_noise_free(simpnoise);
  MPI_Finalize();

  return 0;
}

void refine(int t, float rpId, float thres, float value, int level, float x, float y, float z, float dx, float dy, float dz, struct osn_context *osn)
{
  double noisespacefreq = 10;    /* Spatial frequency of noise */
  double noisetimefreq = 0.25;    /* Temporal frequency of noise */
  float x_center, y_center, z_center, center_val;
  float xpts[8], ypts[8], zpts[8];
  float data[8];
  float conns[12];
  int i;
  int below = 0;
  int above = 0;
  int split = 0;
  int exact = 0;
  int maxLevel = 4;
  float refinePathID = 0;
  float new_refinePathID = 0;
  
  /* given a unit cube (8 points) */
  /* - calulate poistion x,y,z for all points */
  xpts[0] = x;
  ypts[0] = y;
  zpts[0] = z;
  xpts[1] = x + dx;
  ypts[1] = y;
  zpts[1] = z;
  xpts[2] = x;
  ypts[2] = y + dy;
  zpts[2] = z;
  xpts[3] = x + dx;
  ypts[3] = y + dy;
  zpts[3] = z;
  xpts[4] = x;
  ypts[4] = y;
  zpts[4] = z + dz;
  xpts[5] = x + dx;
  ypts[5] = y;
  zpts[5] = z + dz;
  xpts[6] = x;
  ypts[6] = y + dy;
  zpts[6] = z + dz;
  xpts[7] = x + dx;
  ypts[7] = y + dy;
  zpts[7] = z + dz;

  /* refinePathID = pow(10,level) * rpId; */
  refinePathID = 10 * rpId;

  printf("Block refinePathID = %d\n", (int) refinePathID);
  printf("Point_1: (%f %f, %f) \n", xpts[0], ypts[0], zpts[0]);
  printf("Point_8: (%f %f, %f) \n", xpts[7], ypts[7],  zpts[7]);

  /* Get block center value */
  x_center = x + dx/2;
  y_center = y + dy/2;
  z_center = z + dz/2;
  center_val = (float)open_simplex_noise4(osn, x_center*noisespacefreq,
					  y_center*noisespacefreq, z_center*noisespacefreq, t*noisetimefreq);
  
  printf("Block_Center: (%f %f, %f) with val=%f\n", x_center, y_center, z_center, center_val);

  /* - calulate value using open_simplex_noise4 */
  for (i=0; i<8; i++)
  {
    data[i] = (float)open_simplex_noise4(osn, xpts[i]*noisespacefreq,
					 ypts[i]*noisespacefreq, zpts[i]*noisespacefreq, t*noisetimefreq);

    if (data[i] == thres)
      exact = 1;
    
    if (data[i] > thres)
    	above = 1;
    
    if (data[i] < thres)
    	 below = 1;
    
    printf("Point_%d: (%f %f, %f) = %f \n", i, xpts[i], ypts[i], zpts[i], data[i]);
  }

  /* if (center_val > (thres + (2*.02))) */
  /*   split = 1; */

  /* printf("%d!= %d\n", (int) (value*10000), (int) (center_val*10000)); */



  /* - check refinement criteria */
  if (below && above  && !exact && (level < maxLevel))
    split = 1;

  
  /*   Refine if criteria satisfied */
  /* if (split && ( (int) (value*10000) != (int) (center_val*10000)) && !exact) */
  if (split)
  {
    xpts[0] = x;
    ypts[0] = y;
    zpts[0] = z;
    xpts[1] = x + dx/2;
    ypts[1] = y;
    zpts[1] = z;
    xpts[2] = x;
    ypts[2] = y + dy/2;
    zpts[2] = z;
    xpts[3] = x + dx/2;
    ypts[3] = y + dy/2;
    zpts[3] = z;
    xpts[4] = x;
    ypts[4] = y;
    zpts[4] = z + dz/2;
    xpts[5] = x + dx/2;
    ypts[5] = y;
    zpts[5] = z + dz/2;
    xpts[6] = x;
    ypts[6] = y + dy/2;
    zpts[6] = z + dz/2;
    xpts[7] = x + dx/2;
    ypts[7] = y + dy/2;
    zpts[7] = z + dz/2;

    printf("Start to refine block: %f ==============================\n", refinePathID);
    level++;
    for (i=0; i<8; i++)
    {
      new_refinePathID = refinePathID + (i+1)/10.0;
      printf("Refine on cell %d: (%f %f, %f) = %f ... refinePathID=%f\n", i, xpts[i], ypts[i], zpts[i], data[i], new_refinePathID);
      refine(t, new_refinePathID, thres, center_val, level, xpts[i], ypts[i], zpts[i], dx/2, dy/2, dz/2, osn);
    }
  }
  else
    printf("Block %f connot be Refined ... refinement level=%d +++++++++++++++\n", refinePathID, level);
 
  /* - print values and connectivity */
    
}


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
