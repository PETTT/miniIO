/*
* Copyright (c) DoD HPCMP PETTT.  All rights reserved.
* See LICENSE file for details.
*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <mpi.h>
#include "open-simplex-noise.h"
#include "timer.h"


/* #include <limits.h> */
/* #include <assert.h> */

static const float FILLVALUE = -999;

#ifdef HAS_ADIOS
#  include "adiosstruct.h"
#endif

#ifdef HAS_HDF5
#  include "hdf5struct.h"
#endif

void print_usage(int rank, const char *errstr);

int main(int argc, char **argv)
{
  int debug=0;                  /* Flag to generate debug prints statements */
  int debugIO=0;                /* Flag to generate debug IO*/  
  int a, i, j, k, t;            /* loop indices */
  int tt;                       /* Actual time step from tstart */
  float x, y, z;
  double noisespacefreq = 10.0; /* Spatial frequency of noise */
  double noisetimefreq = 0.25;  /* Temporal frequency of noise */
  int tstart = 0;
  int nt = 10;                  /* Number of time steps */
  int ni = 0;                   /* Global grid size */
  int nj = 0;
  int nk = 0;
  int inp = 0;      /* Number of tasks in i */
  int jnp = 0;      /* Number of tasks in j */
  int knp = 1;      /* Number of tasks in k */
  int numtask;                   /* task per processor */
  int point_id, tmp_id;          /* grid point id */
  int x_index, y_index, z_index; /* point index along each axis */
  int xy_dims,  x_dims;
  float deltax, deltay, deltaz; 
  int numPoints;
  float *data;
  float *height;
  int height_index;
  int hindex;
  int maskTindex;
  int *ola_mask;
  int *ol_mask;
  float mask_thres=0.0;             /* mask threshold */
  int mask_thres_index;
  struct osn_context *simpnoise;    /* Open simplex noise context */
  double heighttime, computetime, outtime;   /* Timers */
  
  const int num_varnames=4;
  char *varnames[num_varnames];
  

  /* MPI vars */
  MPI_Comm comm = MPI_COMM_WORLD;
  int cprocs[3], cpers[3], crnk[3];  /* MPI Cartesian info */
  int rank, nprocs;
  int cni, cnj, cnk;   /* Points in this task */
  int is, js, ks;     /* Global index starting points */
  float xs, ys, zs;    /* Global coordinate starting points */  

  /* ADIOS vars */
  uint64_t cstart=0;
  uint64_t cnpoints=0;
  uint64_t npoints=0;

#ifdef HAS_HDF5
  int hdf5out = 0;
#endif

#ifdef HAS_ADIOS
  char      *adios_groupname="struct";
  char      *adios_method="MPI";
  struct adiosstructinfo adiosstruct_nfo;
#endif

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    

  /* Parse command line */
  for(a = 1; a < argc; a++) {


    if(!strcasecmp(argv[a], "--tasks")) {
            inp = atoi(argv[++a]);
            jnp = atoi(argv[++a]);
    } else if(!strcasecmp(argv[a], "--size")) {
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
    }else if(!strcasecmp(argv[a], "--debugIO")) {
      debugIO = 1; 
    }else if(!strcasecmp(argv[a], "--hdf5")) {
#ifdef HAS_HDF5
      hdf5out = 1;
#else
      if(rank == 0)   fprintf(stderr, "HDF5 option not available: %s\n\n", argv[a]);
      print_usage(rank, NULL);
      MPI_Abort(MPI_COMM_WORLD, 1); 
#endif
    } else {
      if(rank == 0)   fprintf(stderr, "Option not recognized: %s\n\n", argv[a]);
      print_usage(rank, NULL);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  numPoints = ni*nj*nk;
  npoints  = numPoints;
 
  /* Check arguments & proc counts */
  if(inp < 1 || jnp < 1 ) {
    print_usage(rank, "Error: tasks not specified or incorrect");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if(ni < 1 || nj < 1 || nk < 1) {
    print_usage(rank, "Error: size not specified or incorrect");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if(inp*jnp != nprocs) {
    print_usage(rank, "Error: product of tasks does not equal total MPI tasks");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if(ni % inp || nj % jnp || nk % knp) {
    print_usage(rank, "Error: number of points on an axis is not evenly divisible "
                "by axis tasks.\n   This is required for proper load balancing.");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if(nt < 1 ) {
    print_usage(rank, "Error: number of timesteps not specified or incorrect");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

   /* Set up Cartesian communicator */
  cprocs[0] = inp;  cprocs[1] = jnp;  cprocs[2] = knp;
  cpers[0] = 0;  cpers[1] = 0;  cpers[2] = 0;    /* No periodicity */
  MPI_Cart_create(MPI_COMM_WORLD, 3, cprocs, cpers, 1, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Cart_coords(comm, rank, 3, crnk);

  deltax = 1.f/(ni-1);
  deltay = 1.f/(nj-1);
  deltaz = 1.f/(nk-1);
  cni = ni / inp;
  cnj = nj / jnp;
  cnk = nk / knp;
  is = crnk[0] * cni;
  js = crnk[1] * cnj;
  ks = crnk[2] * cnk;
  xs = is * deltax;
  ys = js * deltay;
  zs = ks * deltaz;

  xy_dims = ni * nj;
  x_dims = ni;

  /* mask_thres_index = (int) mask_thres/deltaz; */
  mask_thres_index = (int) ( ((mask_thres+1)/2) * (nk-1));
  maskTindex = nk-1;
  
  if (debug) {
    /* printf("rank_cord(%d,%d,%d) rank=%d: %d of %d\n", crnk[0], crnk[1], crnk[2] , rank, rank+1, nprocs); */
    printf("Grid size= (%d x %d x %d) = %llu points\n", ni, nj, nk, npoints);
    /* printf("Grid deltas= (%f x %f x %f)\n", deltax, deltay, deltaz); */
    printf("Mask theshold = %f   mask threshold index = %d\n", mask_thres, mask_thres_index);
  }
  /* Set up osn */
  open_simplex_noise(12345, &simpnoise);   /* Fixed seed, for now */

  /* Allocate arrays */
  data = (float *) malloc((size_t)cni*cnj*cnk*sizeof(float));
  height = (float *) malloc((size_t)cni*cnj*cnk*sizeof(float));
  ola_mask = (int *) malloc((size_t)cni*cnj*cnk*sizeof(int));
  ol_mask = (int *) malloc((size_t)cni*cnj*cnk*sizeof(int));

  varnames[0] = "data";
  varnames[1] = "height";
  varnames[2] = "ola_mask";
  varnames[3] = "ol_mask";

  /* init ADIOS */
#ifdef HAS_ADIOS

  adiosstruct_init(&adiosstruct_nfo, adios_method, adios_groupname, comm, rank, nprocs, nt,
		   ni, nj, nk, is, cni, js, cnj, ks, cnk, deltax, deltay, deltaz, FILLVALUE);
  adiosstruct_addrealxvar(&adiosstruct_nfo, varnames[0], data);

  if (debugIO) {
    adiosstruct_addrealxvar(&adiosstruct_nfo, varnames[1], height);
    adiosstruct_addintxvar(&adiosstruct_nfo, varnames[2], ola_mask);
    adiosstruct_addintxvar(&adiosstruct_nfo, varnames[3], ol_mask);
  }
#endif

  /* generate masked grid */
  /* Spatial loops */
  size_t ii;     /* data index */

  timer_tick(&heighttime, comm, 1);
  z = zs;
  for(k = 0, ii = 0; k < cnk; k++) {
    y = ys;
    for(j = 0; j < cnj; j++) {
      x = xs;
      for(i = 0; i < cni; i++, ii++) {
	x_index = (int) (x/deltax);
	y_index = (int) (y/deltay);
	z_index = (int) (z/deltaz);

	/* calculate point index */
	point_id = (z_index * xy_dims) + (y_index * x_dims) + x_index;
	    
	height[ii] =  (float)open_simplex_noise2(simpnoise, x*noisespacefreq, y*noisespacefreq);

	/* height_index = (int) height[ii]/deltaz; */
	height_index = (int) (((height[ii]+1)/2) * (nk-1));
	  
	/* Calculate ola_mask values */
	if (z_index > mask_thres_index  && z_index > height_index) {
	  ola_mask[ii] = 2;  /* Atmosphere */
	}
	else if (z_index < mask_thres_index  && z_index > height_index) {
	  ola_mask[ii] = 0;  /* ocean */
	}
	else if (z_index <= height_index) {
	  if (height[ii] >= mask_thres  || z_index < height_index) {
	    ola_mask[ii] = 1;  /* land */
	  }
	  else {
	    ola_mask[ii] = 0;  /* ocean */
	  }
	}
	else if (z_index == mask_thres_index  && height[ii] <= mask_thres) {
	  ola_mask[ii] = 0;  /* ocean */
	}
	else {
	  printf("WARNING: ola_mask condition not considered for Point_index: (%d,%d,%d)\n"
		 "Point_id: %d  Height: %f HeightID: %d  mask_thres_index=%d\n",
		 x_index, y_index, z_index, point_id+1, height[ii], height_index, mask_thres_index); 
	}

	
	hindex = (int) ( (((height[ii]+1)/2) * (nk-1)) / (mask_thres_index+1)) * (nk-1);

	hindex = (int) ((height[ii]+1) / (mask_thres+1) * (nk-1));
	
	if (hindex > maskTindex) {
	  hindex = maskTindex;
	}

	/* Calculate ol_mask values */
	if (z_index < maskTindex  && z_index > hindex) {
	  ol_mask[ii] = 0;  /* ocean */
	}
	else if (z_index <= hindex) {
	  if (height[ii] >= mask_thres  || z_index < hindex) {
	    ol_mask[ii] = 1;  /* land */
	  }
	  else {
	    ol_mask[ii] = 0;  /* ocean */
	  }
	}
	else if (z_index == maskTindex  && height[ii] <= mask_thres) {
	  ol_mask[ii] = 0;  /* ocean */
	}
	else {
	  printf("WARNING: ol_mask condition not considered for Point_index: (%d,%d,%d)\n"
		 "Point_id: %d  Height: %f HeightID: %d maskTindex=%d\n",
		 x_index, y_index, z_index, point_id+1, height[ii], hindex, maskTindex); 
	}

	if (debug) {
	  printf("++++++++++++++++++++++++++++++++++++++++++++\n");
	  printf("rank_cord(%d,%d,%d) rank=%d: %d of %d\n", crnk[0], crnk[1], crnk[2] , rank, rank+1, nprocs);
	  printf("LDims: (%d,%d,%d)\n", cni, cnj, cnk); 
	  printf("GDims: (%d,%d,%d)\n", ni, nj, nk);
	  printf("SDims: (%d,%d,%d)\n", is, js, ks); 
	  printf("Point_index: (%d,%d,%d), Point_id:  %d\n", x_index, y_index, z_index,  point_id+1);
	  printf("Point_pos: (%f, %f, %f)  mask_thres_index %d -> %d\n", x, y, z, mask_thres_index, maskTindex);
	  printf("Height: %f HeightID: %d -> %d ola_mask=%d -> %d\n", height[ii], height_index,  hindex, ola_mask[ii], ol_mask[ii]);
	}	      
	
	x += deltax;
      }
      y += deltay;
    }
    z += deltaz;
  }	
  timer_tock(&heighttime);
 
  /* generate ocean land data */
  for(t = 0, tt = tstart; t < nt; t++, tt++) {
    /* Spatial loops */

    timer_tick(&computetime, comm, 1);
    z = zs;
    for(k = 0, ii = 0; k < cnk; k++) {
      y = ys;
      for(j = 0; j < cnj; j++) {
	x = xs;
	for(i = 0; i < cni; i++, ii++) {

	  if ( ol_mask[ii] == 0) {
	    /* if  ( ola_mask[ii] == 0) { */
	    data[ii] = (float)open_simplex_noise4(simpnoise, x*noisespacefreq, y*noisespacefreq, z*noisespacefreq, tt*noisetimefreq);
	  }
	  else {
	    data[ii] = FILLVALUE;
	  }
       
	  x += deltax;
	}
	y += deltay;
      }
      z += deltaz;
    }	
    timer_tock(&computetime);


    timer_tick(&outtime, comm, 1);
#ifdef HAS_ADIOS

    adiosstruct_write(&adiosstruct_nfo, tt);
    
#endif
 
#ifdef HAS_HDF5
    if(hdf5out) {
      if(rank == 0) {
	printf("      Writing hdf5...\n");   fflush(stdout);
      }
      writehdf5(num_varnames, varnames, comm, rank, nprocs, tt,
		is, js, ks,
		ni, nj, nk, cni, cnj, cnk,  
		deltax, deltay, deltaz,
		data, height, ola_mask, ol_mask);
    }
#endif

    timer_tock(&outtime);
    timer_collectprintstats(computetime, comm, 0, "   Compute");
    timer_collectprintstats(outtime, comm, 0, "   Output");
    timer_collectprintstats(heighttime, comm, 0, "   Height");

  }

    /* finalize ADIOS */
#ifdef HAS_ADIOS
    adiosstruct_finalize(&adiosstruct_nfo);
#endif

  open_simplex_noise_free(simpnoise);
  free(data);
  free(height);
  free(ola_mask);
  free(ol_mask);

  MPI_Finalize();

  return 0;
}


void print_usage(int rank, const char *errstr)
{
  if(rank != 0)  return;
  if(errstr)
    fprintf(stderr, "%s\n\n", errstr);
  fprintf(stderr,
	  "Usage: mpi_launcher [-n|-np NPROCS] ./struct --tasks INP JNP --size NI NJ NK [options]\n"
	  "    NPROCS : # of tasks launched by MPI; may or may not be implied or required by system\n\n"
	  "  Required:\n"
	  "    --tasks INP JNP: Specifies the parallel decomposition of tasks\n"
	  "      INP : # of tasks along the I (X) axis\n"
	  "      JNP : # of tasks along the J (Y) axis\n"
	  "        NOTE that INP * JNP * KNP == NPROCS is required!\n"
	  "    --size NI NJ NK : Specifies the size of the grid\n"
	  "      NI, NJ, NK : Number of grid points along the I,J,K axes respectively\n"
	  "      valid values are > 1\n\n" 
	  "  Optional:\n"
	  "    --debug: Turns on debugging print statements \n"
	  "    --debugIO: Turns on debugging IO (corrently only works with ADIOS IO) \n"
	  "    --maskthreshold MT : Mask theshold; valid values are floats between -1.0 and 1.0 \n"
	  "      MT : mask threshold value; Default: 0.0\n"
	  "    --noisespacefreq FNS : Spatial frequency of noise function\n"
	  "      FNS : space frequency value; Default: 10.0\n"
	  "    --noisetimefreq FNT : Temporal frequency of noise function\n"
	  "      FNT : time frequency value;  Default: 0.25\n"
	  "    --tsteps NT : Number of time steps; valid values are > 0 (Default value 10)\n"
	  "    --tstart TS : Starting time step; valid values are >= 0  (Default value 0)\n"
#ifdef HAS_HDF5
	  "    --hdf5 : Enable HDF5 output (i.e. XDMF)\n"
#endif
	  );
}
