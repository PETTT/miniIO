#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <mpi.h>
#include "open-simplex-noise.h"
#include "cubes.h"

/* #include <limits.h> */
/* #include <assert.h> */
#ifdef HAS_ADIOS
#  include "adiosamr.h"
#endif


static const int fnstrmax = 4095;

void print_usage(int rank, const char *errstr);


int main(int argc, char **argv) {
  int debug=0;
  int i, a, numtask, block_id, realblock_id, tmp_id , t;                  /* loop indices */
  int tt;                       /* Actual time step from tstart */
  double noisespacefreq = 10;   /* Spatial frequency of noise */
  double noisetimefreq = 0.25;  /* Temporal frequency of noise */
  int tstart = 0;
  int nt = 50;                  /* Number of time steps */
  int x_id, y_id, z_id;
  int xy_dim,  x_dim;
  float deltax, deltay, deltaz;
  float deltax_center, deltay_center, deltaz_center;
  int maxLevel = 8;
  float threshold=0.0;
  uint64_t numOcts=0;        /* Number of cubes */
  uint64_t numPoints=0;    /* Number of points */
  uint64_t numCoords=0;    /* Number of coords */
  uint64_t numDarrays=0;    /* Number of data arrays */
  int ni = 0;      /* Global grid size */
  int nj = 0;
  int nk = 0;
  
  int nx = 0;
  int ny = 0;
  int nz = 0;
  cubeInfo cubedata;

  int numBlocks;
  struct osn_context *simpnoise;    /* Open simplex noise context */
  
  /* MPI vars */
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, nprocs; 
  float xs, ys, zs;    /* Global coordinate starting points */ 

#ifdef HAS_ADIOS
  char      *adios_groupname="amr";
  char      *adios_method="MPI";
  struct adiosamrinfo adiosamr_nfo;
#endif
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);
  
  /* Parse command line */
  for(a = 1; a < argc; a++) {
    if(!strcasecmp(argv[a], "--blockdim")) {
      nx = ny = nz = atoi(argv[++a]);
    }else if(!strcasecmp(argv[a], "--threshold")) {
      threshold = strtof(argv[++a], NULL);
    } else if(!strcasecmp(argv[a], "--levels")) {
      maxLevel = atoi(argv[++a]);
    } else if(!strcasecmp(argv[a], "--noisespacefreq")) {
      noisespacefreq = strtod(argv[++a], NULL);
    } else if(!strcasecmp(argv[a], "--noisetimefreq")) {
      noisetimefreq = strtod(argv[++a], NULL);
    } else if(!strcasecmp(argv[a], "--tsteps")) {
      nt = atoi(argv[++a]);
    } else if(!strcasecmp(argv[a], "--tstart")) {
      tstart = atoi(argv[++a]);
    }else if(!strcasecmp(argv[a], "--debug")) {
      debug = 1;
    } else {
      if(rank == 0)   fprintf(stderr, "Option not recognized: %s\n\n", argv[a]);
      print_usage(rank, NULL);
      MPI_Abort(comm, 1);
    }
  }

  numBlocks = nx*ny*nz;
  
  /* Check arguments & proc counts */
  if(nx < 1 || ny < 1 || nz < 1) {
    print_usage(rank, "Error: blockdim not specified or incorrect");
    MPI_Abort(comm, 1);
  }
  
  if( numBlocks % nprocs) {
    print_usage(rank, "Error: number of grid points is not evenly divisible "
                "by number of processors.\n   This is required for proper load balancing.");
    MPI_Abort(comm, 1);
  }

  if(maxLevel < 0 ) {
    print_usage(rank, "Error: number of levels not specified or incorrect");
    MPI_Abort(comm, 1);
  }
 
  numtask =  numBlocks/nprocs;
  
  ni = nx + 1;
  nj = ny + 1;
  nk = nz + 1;
  
  deltax = 1.f/(ni-1);
  deltay = 1.f/(nj-1);
  deltaz = 1.f/(nk-1);

  deltax_center = deltax/2;
  deltaz_center = deltay/2;
  deltaz_center = deltaz/2;

  xy_dim = nx * ny;
  x_dim = nx;

  /* Set up osn */
  open_simplex_noise(12345, &simpnoise);   /* Fixed seed, for now */

  /* Allocate arrays */
  cubesinit(&cubedata, numtask, maxLevel, debug);
  
  /* init ADIOS */
#ifdef HAS_ADIOS
  adiosamr_init(&adiosamr_nfo, adios_method, adios_groupname, comm, rank, nprocs, nt);
  adiosamr_addxvar(&adiosamr_nfo, "data");
#endif
  
  for(t = 0, tt = tstart; t < nt; t++, tt++) {

    if (debug) {
	printf("Hi: rank=%d: %d of %d:  timestep=%d\n", rank, rank+1, nprocs, tt);
    }
    
    for (i=0; i<numtask; i++) {
      block_id = (rank*numtask)+i;
      /* realblock_id = block_id - 1; */

      /* Calculate id for lower-left-front corner of block */
      z_id =  block_id / xy_dim;
      tmp_id =  block_id % xy_dim;
      y_id = tmp_id/x_dim ;
      x_id = tmp_id%x_dim ;

      xs = x_id * deltax;
      ys = y_id * deltay;
      zs = z_id * deltaz;

      if (debug) {
	printf("Block_id=%d\n", block_id+1);
	printf("Point_id: (%d, %d, %d) \n", x_id, y_id, z_id);
      }
      refine(&cubedata, tt, (block_id+1), threshold, .1 + threshold, 0, xs, ys, zs, deltax, deltay, deltaz, simpnoise, maxLevel);
    }

    /* print out data */
    if (debug)  {
      cubeprint(&cubedata);
    }
    
#ifdef HAS_ADIOS
    adiosamr_write(&adiosamr_nfo, tt, cubedata.npoints, cubedata.points, &cubedata.data);
#endif
    
  }

  if (debug) 
    printf("Finalizing:rank %d \n", rank);
  
  /* finalize ADIOS */
#ifdef HAS_ADIOS
  adiosamr_finalize(&adiosamr_nfo);
#endif

  open_simplex_noise_free(simpnoise);
  cubesfree(&cubedata);
  MPI_Finalize();

  return 0;
}



void print_usage(int rank, const char *errstr)
{
  if(rank != 0)  return;
  if(errstr)
    fprintf(stderr, "%s\n\n", errstr);
  fprintf(stderr,
	  "Usage: mpi_launcher [-n|-np NPROCS] ./arm --blockdim N [options]\n"
	  "    NPROCS : # of tasks launched by MPI; may or may not be implied or required by system\n\n"
	  "  Required:\n"
	  "    --blockdim N : Specifies the initial number of blocks along each axes\n"
	  "      N : Number of blocks along the I,J,K axes respectively\n"
	  "      valid values are > 0\n\n" 
	  "  Optional:\n"
	  "    --debug: Turns on debugging statements \n"
	  "    --threshold T : Mask theshold; valid values are floats between -1.0 and 1.0 \n"
	  "      T : threshold value; Default: 0.0\n"
	  "    --levels L : Maximum levels of refinement; valid values are >= 0 \n"
	  "      L : max refinmentment levels value; Default: 8\n"
	  "    --noisespacefreq FNS : Spatial frequency of noise function\n"
	  "      FNS : space frequency value; Default: 10.0\n"
	  "    --noisetimefreq FNT : Temporal frequency of noise function\n"
	  "      FNT : time frequency value;  Default: 0.25\n"
	  "    --tsteps NT : Number of time steps; valid values are > 0;  Default:  50)\n"
	  "    --tstart TS : Starting time step; valid values are > 0;  Default: 0\n"
	  );

  /*## End of Output Module Usage Strings ##*/
}
