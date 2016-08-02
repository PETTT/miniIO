#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <mpi.h>
#include "open-simplex-noise.h"


/* #include <limits.h> */
/* #include <assert.h> */

#ifdef HAS_ADIOS
#  include "adiosstruct.h"
#endif

static const int fnstrmax = 4095;

void print_usage(int rank, const char *errstr);

int main(int argc, char **argv)
{
  int debug=0;
  int a, i, t;                  /* loop indices */
  int tt;                       /* Actual time step from tstart */
  double noisespacefreq = 10;   /* Spatial frequency of noise */
  double noisetimefreq = 0.25;  /* Temporal frequency of noise */
  int tstart = 0;
  int nt = 50;                  /* Number of time steps */
  int ni = 0;                   /* Global grid size */
  int nj = 0;
  int nk = 0;
  int numtask;                  /* task per processor */
  int point_id, tmp_id;         /* grid point id */
  int x_id, y_id, z_id;         /* point id along each axis */
  int xy_dim,  x_dim;
  float deltax, deltay, deltaz; 
  int numPoints;
  float *data;
  float *height;
  int height_id;
  uint64_t *ola_mask;
  float mask_thres=0.0;             /* mask threshold */
  int mask_thres_id;
  struct osn_context *simpnoise;    /* Open simplex noise context */
  
  /* MPI vars */
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, nprocs;
  int cni, cnj, cnk;   /* Points in this task */
  int is, js, ks;     /* Global index starting points */
  float xs, ys, zs;    /* Global coordinate starting points */  

  /* ADIOS vars */
  uint64_t cstart=0;
  uint64_t cnpoints=0;
  uint64_t npoints=0;

#ifdef HAS_ADIOS
  char      *adios_groupname="struct";
  char      *adios_method="MPI";
  struct adiosstructinfo adiosstruct_nfo;
#endif

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  /* Parse command line */
  for(a = 1; a < argc; a++) {
    if(!strcasecmp(argv[a], "--size")) {
      ni = atoi(argv[++a]);
      nj = atoi(argv[++a]);
      nk = atoi(argv[++a]);
    } else if(!strcasecmp(argv[a], "--maskthreshold")) {
      mask_thres = strtof(argv[++a], NULL);
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

  numPoints = ni*nj*nk;
  npoints  = numPoints;
 
  /* Check arguments & proc counts */
  if(ni < 2 || nj < 2 || nk < 2) {
    print_usage(rank, "Error: size not specified or incorrect");
    MPI_Abort(comm, 1);
  }
  
  if( numPoints % nprocs) {
    print_usage(rank, "Error: number of grid points is not evenly divisible "
                "by number of processors.\n   This is required for proper load balancing.");
    MPI_Abort(comm, 1);
  }

  if(nt < 1 ) {
    print_usage(rank, "Error: number of timesteps not specified or incorrect");
    MPI_Abort(comm, 1);
  }
  
  cni = ni / nprocs;
  cnj = nj / nprocs;
  cnk = nk / nprocs;
 
  numtask =  numPoints/nprocs;
  cnpoints  = numtask;
  
  xy_dim = ni * nj;
  x_dim = ni;
  
  deltax = 1.f/(ni-1);
  deltay = 1.f/(nj-1);
  deltaz = 1.f/(nk-1);
  
  mask_thres_id = (int) mask_thres/deltaz;

  if (debug) {
    printf("Grid size= (%d x %d x %d) = %llu points\n", ni, nj, nk, npoints);
    printf("Mask theshold = %f   mask threshold id = %d\n", mask_thres, mask_thres_id);
  }
  
  /* Set up osn */
  open_simplex_noise(12345, &simpnoise);   /* Fixed seed, for now */

  /* Allocate arrays */
  data = (float *) malloc((size_t)numtask*sizeof(float));
  height = (float *) malloc((size_t)numtask*sizeof(float));
  ola_mask = (uint64_t *) malloc((size_t)numtask*sizeof(uint64_t));
  /* data_vars = 3; */

  /* init ADIOS */
#ifdef HAS_ADIOS

  adiosstruct_init(&adiosstruct_nfo, adios_method, adios_groupname, comm, rank, nprocs, nt);
  adiosstruct_addxvar(&adiosstruct_nfo, "data", data);
  adiosstruct_addxvar(&adiosstruct_nfo, "height", height);
  
  /* adios_define_var(adios_groupid, "ola_mask", "", adios_integer, "numtask", */
  /* 		   "numPoints", "var_start"); */
#endif
  
  for(t = 0, tt = tstart; t < nt; t++, tt++) {
    cstart = rank*numtask;
    for (i=0; i<numtask; i++) {
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
      /*   height =  (float)open_simplex_noise2(simpnoise, xs*noisespacefreq, ys*noisespacefreq); */
      /* } */

      height[i] =  (float)open_simplex_noise2(simpnoise, xs*noisespacefreq, ys*noisespacefreq);
      data[i] = (float)open_simplex_noise4(simpnoise, xs*noisespacefreq, ys*noisespacefreq, zs*noisespacefreq, tt*noisetimefreq);

      height_id = (int) height[i]/deltaz;

      /* /\* other way of calculating ola_mask values *\/ */
      /* if (height_id == z_id ) { */
      /* 	if (z_id > mask_thres_id) {	     */
      /* 	  ola_mask[i] = 2;  /\* land surface above *\/ */
      /* 	} else if (z_id == mask_thres_id) {	     */
      /* 	  ola_mask[i] = 2;  /\* land surface at sea level*\/ */
      /* 	} else if (z_id < mask_thres_id) {	     */
      /* 	  ola_mask[i] = 2;  /\* land surface below sea level*\/ */
      /* 	} */
      /* } else if ( z_id < height_id ) {     	 */
      /* 	ola_mask[i] = 1;  /\* below land surface *\/ */
      /* } else if ( z_id  > height_id) {	 */
      /* 	if (z_id > mask_thres_id) {	     */
      /* 	  ola_mask[i] = 3;  /\* Atmosphere*\/ */
      /* 	} else if (z_id == mask_thres_id) {	     */
      /* 	  ola_mask[i] = 0;  /\* ocean at sea level*\/ */
      /* 	} else if (z_id < mask_thres_id) {	     */
      /* 	  ola_mask[i] = 0;  /\* ocean below sea level*\/ */
      /* 	} */
      /* } */

      
      /* Calculate ola_mask values */
      if (height_id == z_id ) {
	if (z_id > mask_thres_id) {	    
	  ola_mask[i] = 2;  /* land surface above */
	} else if (z_id == mask_thres_id) {	    
	  ola_mask[i] = 2;  /* land surface at sea level*/
	} else if (z_id < mask_thres_id) {	    
	  ola_mask[i] = 2;  /* land surface below sea level*/
	}
      } else if (z_id == mask_thres_id) {	    
	ola_mask[i] = 0;  /* ocean at sea level*/
      } else if (z_id < mask_thres_id) {
	ola_mask[i] = 0;  /* ocean below sea level*/
      } else if (z_id > mask_thres_id) {	    
	  ola_mask[i] = 3;  /* Atmosphere*/
      }
      

      if (debug) {
	printf("++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("timestep=%d rank=%d: %d of %d:  point_id=%d\n", tt, rank, rank+1, nprocs, point_id+1);
	printf("timestep=%d Point_id: (%d, %d, %d) = %f\n", tt, x_id, y_id, z_id, data[i]);
	printf("timestep=%d Point_pos: (%f, %f, %f) = %f\n", tt, xs, ys, zs, data[i]);
	printf("timestep=%d Height: %f Height id: %d mask=%llu\n", tt, height[i], height_id, ola_mask[i]);

	/*  /\*  write out land surface data -- later can be changed *\/ */
	/* if  (ola_mask[i] == 2) { */
	/* 	printf("++++++++++++++++++++++++++++++++++++++++++++\n");       */
	/* 	printf("timestep=%d rank=%d: %d of %d:  point_id=%d\n", tt, rank, rank+1, nprocs, point_id+1); */
	/* 	printf("timestep=%d Point_id: (%d, %d, %d) = %f\n", tt, x_id, y_id, z_id, data[i]); */
	/* 	printf("timestep=%d Point_pos: (%f, %f, %f) = %f\n", tt, xs, ys, zs, data[i]); */
	/* 	printf("timestep=%d Height: %f Height id: %d mask=%d\n", tt, height[i], height_id, ola_mask[i]); */
	/* } */
      }
    }
    
#ifdef HAS_ADIOS

    adiosstruct_write(&adiosstruct_nfo, tt, cnpoints, npoints, cstart, ola_mask);
    /* adios_write(adios_handle, "ola_mask", ola_mask); */
    
#endif 

  }


    /* finalize ADIOS */
#ifdef HAS_ADIOS
    /* adios_finalize(rank); */
    adiosstruct_finalize(&adiosstruct_nfo);
#endif

  open_simplex_noise_free(simpnoise);
  free(data);
  free(height);
  free(ola_mask);

  MPI_Finalize();

  return 0;
}


void print_usage(int rank, const char *errstr)
{
  if(rank != 0)  return;
  if(errstr)
    fprintf(stderr, "%s\n\n", errstr);
  fprintf(stderr,
	  "Usage: mpi_launcher [-n|-np NPROCS] ./struct --size NI NJ NK [options]\n"
	  "    NPROCS : # of tasks launched by MPI; may or may not be implied or required by system\n\n"
	  "  Required:\n"
	  "    --size NI NJ NK : Specifies the size of the grid\n"
	  "      NI, NJ, NK : Number of grid points along the I,J,K axes respectively\n"
	  "      valid values are > 1\n\n" 
	  "  Optional:\n"
	  "    --debug: Turns on debugging statements \n"
	  "    --maskthreshold MT : Mask theshold; valid values are floats between -1.0 and 1.0 \n"
	  "      MT : mask threshold value; Default: 0.0\n"
	  "    --noisespacefreq FNS : Spatial frequency of noise function\n"
	  "      FNS : space frequency value; Default: 10.0\n"
	  "    --noisetimefreq FNT : Temporal frequency of noise function\n"
	  "      FNT : time frequency value;  Default: 0.25\n"
	  "    --tsteps NT : Number of time steps; valid values are > 0 (Default value 50)\n"
	  "    --tstart TS : Starting time step; valid values are >= 0  (Default value 0)\n"
	  );

  /*## End of Output Module Usage Strings ##*/
}
