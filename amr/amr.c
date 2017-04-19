/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <mpi.h>
#include "open-simplex-noise.h"
#include "cubes.h"
#include "timer.h"

#ifdef HAS_VTKOUT
#include "vtkout.h"
#endif

#ifdef HAS_ADIOS
#  include "adiosamr.h"
#endif

#ifdef HAS_NC
#  include "ncamr.h"
#endif

#ifdef HAS_HDF5
#  include "hdf5amr.h"
#endif


void print_usage(int rank, const char *errstr);


int main(int argc, char **argv) {
  int debug=0;
  int i, j, k, a, block_id, t;  /* loop indices */
  int tt;                        /* Actual time step from tstart */
  float x, y, z;
  double noisespacefreq = 10.0; /* Spatial frequency of noise */
  double noisetimefreq = 0.25;  /* Temporal frequency of noise */
  int tstart = 0;
  int nt = 50;                   /* Number of time steps */
  float deltax, deltay, deltaz;
  int maxLevel = 4;
  float threshold=0.0;
  int inp = 0;              /* Number of tasks in i */
  int jnp = 0;              /* Number of tasks in j */
  int knp = 0;              /* Number of tasks in k */
  int ni = 0;               /* Global grid size */
  int nj = 0;
  int nk = 0;
  cubeInfo cubedata;
  struct osn_context *simpnoise;    /* Open simplex noise context */
  double computetime, outtime;   /* Timers */

  /* MPI vars */
  MPI_Comm comm = MPI_COMM_WORLD;
  int cprocs[3], cpers[3], crnk[3];  /* MPI Cartesian info */
  int rank, nprocs;
  int cni, cnj, cnk;   /* Points in this task */
  int is, js, ks;      /* Global index starting points */
  float xs, ys, zs;    /* Global coordinate starting points */

#ifdef HAS_VTKOUT
    int vtkout = 0;
#endif

#ifdef HAS_ADIOS
  char      *adios_groupname="amr";
  char      *adios_method=NULL;   /* POSIX|MPI|MPI_LUSTRE|MPI_AGGREGATE|PHDF5   */
  struct adiosamrinfo adiosamr_nfo;
#endif

#ifdef HAS_NC
  char      *nc_groupname="amr";
  struct ncamrinfo ncamr_nfo;
  int ncout = 0;
#endif

#ifdef HAS_HDF5
  char      *hdf5_groupname="amr";
  struct hdf5amrinfo hdf5amr_nfo;
  int hdf5out = 0;
#endif

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  /* Parse command line */
  for(a = 1; a < argc; a++) {
    if(!strcasecmp(argv[a], "--tasks")) {
      inp = atoi(argv[++a]);
      jnp = atoi(argv[++a]);
      knp = atoi(argv[++a]);
    } else if(!strcasecmp(argv[a], "--size")) {
      ni = atoi(argv[++a]);
      nj = atoi(argv[++a]);
      nk = atoi(argv[++a]);
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
    }

#ifdef HAS_ADIOS
    else if(!strcasecmp(argv[a], "--adios")) {
      adios_method = argv[++a];
    }
#endif

    else if(!strcasecmp(argv[a], "--hdf5")) {
#ifdef HAS_HDF5
      hdf5out = 1;
#else
      if(rank == 0)   fprintf(stderr, "HDF5 option not available: %s\n\n", argv[a]);
      print_usage(rank, NULL);
      MPI_Abort(comm, 1);
#endif
    }

#ifdef HAS_NC
    else if(!strcasecmp(argv[a], "--nc")) {
      ncout = 1;
    }
#endif

#ifdef HAS_VTKOUT
    else if(!strcasecmp(argv[a], "--vtkout")) {
      vtkout = 1;
    }
#endif


    else {
      if(rank == 0)   fprintf(stderr, "Option not recognized: %s\n\n", argv[a]);
      print_usage(rank, NULL);
      MPI_Abort(comm, 1);
    }
  }

  /* Check arguments & proc counts */
  if(inp < 0 || jnp < 0 || knp < 0) {
    print_usage(rank, "Error: tasks not specified or incorrect");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if(ni < 2 || nj < 2 || nk <  2) {
    print_usage(rank, "Error: size not specified or incorrect");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if(inp*jnp*knp != nprocs) {
    print_usage(rank, "Error: product of tasks does not equal total MPI tasks");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if((ni-1) % inp || (nj-1) % jnp || (nk-1) % knp) {
    print_usage(rank, "Error: number of points-1 on an axis is not evenly divisible "
		"by axis tasks.\n   This is required for proper load balancing.");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }


  if(maxLevel < 0 ) {
    print_usage(rank, "Error: number of levels not specified or incorrect");
    MPI_Abort(comm, 1);
  }

  if (nt < 1) {
    print_usage(rank, "Error: number of timesteps not specified or incorrect");
    MPI_Abort(comm, 1);
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
  cni = (ni-1) / inp;
  cnj = (nj-1) / jnp;
  cnk = (nk-1) / knp;
  is = crnk[0] * cni;
  js = crnk[1] * cnj;
  ks = crnk[2] * cnk;
  xs = is * deltax;
  ys = js * deltay;
  zs = ks * deltaz;



  /* Set up osn */
  open_simplex_noise(12345, &simpnoise);   /* Fixed seed, for now */

  /* Allocate arrays */
  cubesinit(&cubedata, cni*cnj*cnk, maxLevel, debug);

  /* init ADIOS */
#ifdef HAS_ADIOS
  if (adios_method) {
    adiosamr_init(&adiosamr_nfo, adios_method, adios_groupname, comm, rank, nprocs, nt);
    adiosamr_addxvar(&adiosamr_nfo, "data");
  }
#endif

#ifdef HAS_HDF5
  if(hdf5out) {
    hdf5_init(&hdf5amr_nfo, hdf5_groupname, comm, rank, nprocs, nt);
    hdf5_addxvar(&hdf5amr_nfo, "data");
  }
#endif

#ifdef HAS_NC
  if(ncout) {

  }
#endif

  if (debug) {
    printf("(cni=%d, cnj=%d, cnk=%d) \n", cni, cnj, cnk);
  }

  for(t = 0, tt = tstart; t < nt; t++, tt++) {
    size_t ii;     /* data index */

    cubedata.npoints = 0;
    cubedata.ncubes = 0;


    if (debug) {
	printf("Hi: rank=%d: %d of %d:  timestep=%d\n", rank, rank+1, nprocs, tt);
    }

    timer_tick(&computetime, comm, 1);

    z = zs;
    for(k = 0, ii = 0; k < cnk; k++) {
      y = ys;
      for(j = 0; j < cnj; j++) {
	x = xs;
	for(i = 0; i < cni; i++, ii++) {

	  /* calculate block_id */
	  block_id = ii;

	  if (debug) {
	    printf("Start from main Block_id=%d\n", block_id+1);
	  }
	  refine(&cubedata, tt, (block_id+1), threshold, 0, x, y, z, deltax, deltay, deltaz, simpnoise, maxLevel, noisespacefreq, noisetimefreq);

	  x += deltax;
	}
	y += deltay;
      }
      z += deltaz;
    }

    timer_tock(&computetime);

    /* print out data */
    if (debug)  {
      cubeprint(&cubedata);
    }


#ifdef HAS_VTKOUT
    if(vtkout) {
      if(rank == 0) {
	printf("      Writing VTK ...\n");   fflush(stdout);
      }
      writevtk("amr.out", comm, rank, nprocs, tt, cubedata.npoints, cubedata.ncubes,
	       cubedata.points, cubedata.data, "data", debug);
    }
#endif

    timer_tick(&outtime, comm, 1);

#ifdef HAS_ADIOS
    if (adios_method) {
      if(rank == 0) {
	printf("      Writing ADIOS ...\n");   fflush(stdout);
      }
      adiosamr_write(&adiosamr_nfo, tt, cubedata.npoints, cubedata.points, &cubedata.data);
    }
#endif

#ifdef HAS_HDF5
    if(hdf5out) {
      if(rank == 0) {
	printf("      Writing HDF5 ...\n");   fflush(stdout);
      }
      hdf5_write(&hdf5amr_nfo, tt, cubedata.npoints, cubedata.points, &cubedata.data);
    }
#endif

#ifdef HAS_NC
    if(ncout) {
      if(rank == 0) {
        printf("     Writing netCDF ...\n"); fflush(stdout);
      }
      nc_write(&ncamr_nfo, tt, cubedata.npoints, cubedata.points, &cubedata.data);
    }
#endif

    timer_tock(&outtime);
    timer_collectprintstats(computetime, comm, 0, "   Compute");
    timer_collectprintstats(outtime, comm, 0, "   Output");
  }


  if (debug)
    printf("Finalizing:rank %d \n", rank);

  /* finalize ADIOS */
#ifdef HAS_ADIOS
  if (adios_method) adiosamr_finalize(&adiosamr_nfo);
#endif

#ifdef HAS_HDF5
  if(hdf5out) {
    hdf5_finalize(&hdf5amr_nfo);
  }
#endif

#ifdef HAS_NC
  if(ncout) {
    nc_finalize(&ncamr_nfo);
  }
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
	  "Usage: mpi_launcher [-n|-np NPROCS] ./arm --tasks INP JNP KNP --size NI NJ NK [options]\n"
	  "    NPROCS : # of tasks launched by MPI; may or may not be implied or required by system\n\n"
	  "  Required:\n"
	  "    --tasks INP JNP KNP : Specifies the parallel decomposition of tasks\n"
	  "      INP : # of tasks along the I (X) axis\n"
	  "      JNP : # of tasks along the J (Y) axis\n"
	  "      KNP : # of tasks along the K (Z) axis\n"
	  "        NOTE that INP * JNP * KNP == NPROCS is required!\n"
	  "    --size NI NJ NK : Specifies the size of the grid\n"
	  "      NI, NJ, NK : Number of grid points along the I,J,K axes respectively\n\n"
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

#ifdef HAS_NC
          "\n   --nc: Enable netCDF output\n"
#endif

#ifdef HAS_ADIOS
	  "    --adios [POSIX|MPI|MPI_LUSTRE|MPI_AGGREGATE|PHDF5]: Enable ADIOS output\n"
#endif
	  );

#ifdef HAS_VTKOUT
    fprintf(stderr, "    --vtkout : Enable VTK output.\n");
#endif

  /*## End of Output Module Usage Strings ##*/
}
