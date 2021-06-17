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
#include "args.h"
#include "splitspace.h"

/* #include <limits.h> */
/* #include <assert.h> */

static const float FILLVALUE = -999;

/*## Add Output Modules' Includes Here ##*/

#ifdef HAS_ADIOS
#  include "adiosstruct.h"
#endif

#ifdef HAS_NC
#  include "netcdfstruct.h"
#endif

#ifdef HAS_HDF5
#  include "hdf5struct.h"
#endif

/*## End of Output Module Includes ##*/


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
	  "        NOTE that INP * JNP == NPROCS is required!\n"
	  "  Recommended:\n"
	  "    --size NI[x] NJ[x] NK : Specifies the size of the grid\n"
	  "      NI, NJ, NK : Number of grid points along the I,J,K axes respectively\n"
	  "      Put an x after NI or NJ to indicate that the value scales with INP/JNP\n"
	  "      valid values are > 1 with or without x \n"
	  "      Default: 128x 128x 128\n\n"
	  "  Optional:\n"
	  "    --debug : Turns on debugging print statements \n"
	  "    --debugBal : Turns on debugging print statements for balancing\n"
	  "    --maskthreshold MT : Mask theshold; valid values are floats between -1.0 and 1.0 \n"
	  "      MT : mask threshold value; Default: 0.0\n"
	  "    --noisespacefreqmask FNSi[x] FNSj[x] : Spatial frequency of noise function \n"
	  "                                           for land mask\n"
	  "      FNS[x] : space frequency value; Default: 3.5 3.5\n"
	  "    --noisespacefreq FNSi[x] FNSj[x] FNSk : Spatial frequency of noise function \n"
	  "                                            for data\n"
	  "      FNS[x] : space frequency value; Default: 10.0x 10.0x 10.0\n"
	  "    --noisetimefreq FNT : Temporal frequency of noise function\n"
	  "      FNT : time frequency value;  Default: 0.25\n"
	  "    --tsteps NT : Number of time steps; valid values are > 0 (Default value 10)\n"
	  "    --tstart TS : Starting time step; valid values are >= 0  (Default value 0)\n"
	  "    --balance : Turns on computational load balancing \n"

    /*## Add Output Modules' Usage String ##*/

#ifdef HAS_ADIOS
	  "    --adios [POSIX|MPI|MPI_LUSTRE|MPI_AGGREGATE|PHDF5]: Enable ADIOS output\n"
      "    --adiosopts OPTS : Pass options to ADIOS\n"
      "    --adios_transform TRANSFORM : Pass a transform, e.g., compression, to ADIOS\n"
	  "    --debugIO : Turns on debugging IO (corrently only works with ADIOS IO) \n"
#endif


#ifdef HAS_HDF5
	  "    --hdf5 : Enable HDF5 output (i.e. XDMF)\n"
	  "    --hdf5_chunk y z : Chunk Size y z \n"
	  "      valid values are  NJ/JNP/y,NK/z\n"
          "    --hdf5_compress : enable compression. Valid value is a comma seperate (no spaces) list:  \n"
          "        <compression type: gzip or szip>,<compression parameter(s) corresponding to HDF5 compression API>  \n"
          "        gzip,<value is level (see H5Pset_deflate)> \n"
          "        szip,<value is <options_mask>,<pixels_per_block (see H5Pset_szip)> \n"
          "             For example, --hdf5_compress szip,H5_SZIP_NN_OPTION_MASK,8 \n" 
          "        NOTE: compression requires chunked datasets \n"
#endif

#ifdef HAS_NC
      "    --nc : Enable NetCDF Output\n"
#endif
      
    /*## End of Output Module Usage Strings ##*/

	  );
}

int main(int argc, char **argv)
{
  int debug=0;                  /* Flag to generate debug prints statements */
  int debugIO=0;                /* Flag to generate debug IO*/
  int debugBal=0;               /* Flag to generate debug for balance */
  int balance=0;                /* Flag to enable/disable load balancing */
  int a, i, j, k, t;            /* loop indices */
  int tt;                       /* Actual time step from tstart */
  float x, y, z;
  double noisefreqmask_i = 3.5; /* Spatial frequency of noise for mask */
  double noisefreqmask_j = 3.5; /* Spatial frequency of noise for mask */
  int noisefreqmask_iscale = 0;  /* Scale mask spatial frequency with i tasks */
  int noisefreqmask_jscale = 0;  /* Scale mask spatial frequency with j tasks */
  double noisefreq_i = 10.0;    /* Spatial frequency of noise for data */
  double noisefreq_j = 10.0;    /* Spatial frequency of noise for data */
  double noisefreq_k = 10.0;    /* Spatial frequency of noise for data */
  int noisefreq_iscale = 1;    /* Scale spatial frequency with i tasks */
  int noisefreq_jscale = 1;    /* Scale spatial frequency with j tasks */
  double noisefreq_t = 0.25;    /* Temporal frequency of noise */
  int tstart = 0;
  int nt = 10;                  /* Number of time steps */
  int ni = 128;                   /* Global grid size */
  int nj = 128;
  int nk = 128;
  int niscale = 1;              /* Scale global grid size with tasks */
  int njscale = 1;
  int inp = 0;      /* Number of tasks in i */
  int jnp = 0;      /* Number of tasks in j */
  int knp = 1;      /* Number of tasks in k */
  int numtask;                   /* task per processor */
  uint64_t point_id;          /* grid point id */
  int x_index, y_index, z_index; /* point index along each axis */
  int xy_dims,  x_dims;
  float deltax, deltay, deltaz;
  float *data;
  float *height;
  int height_index;
  int hindex;
  int maskTindex;
  int *ola_mask;
  int *ol_mask;
  float mask_thres=0.5;     /* upper mask threshold  (range -1 to 1) */
  float bot_mask_thres=0.5; /* bottom mask threshold (range 0.0 to (mask_thres+1)/2 ) */
  int mask_thres_index;
  struct osn_context *simpnoise;    /* Open simplex noise context */
  double heighttime, computetime, outtime;   /* Timers */
  double initheighttime, loadbaltime;    /* Initialization timers */

  const int num_varnames=1;
  char *varnames[num_varnames];

  /* balance spacesplit vars */
  struct factorstruct factors;
  recList rects_list;
  int *ol_mask_2dsum;
  int tmp_mask;
  int sum_id;
  float h_tmp;
  
  /* MPI vars */
  MPI_Comm comm = MPI_COMM_WORLD;
  int cprocs[3], cpers[3], crnk[3];  /* MPI Cartesian info */
  int rank, nprocs;
  int cni, cnj, cnk;   /* Points in this task */
  int is, js, ks;      /* Global index starting points */
  float xs, ys, zs;    /* Global coordinate starting points */

  /* ADIOS vars */
  uint64_t cstart=0;
  uint64_t cnpoints=0;
  uint64_t npoints=0;

  /*## Add Output Modules' Variables Here ##*/
  
#ifdef HAS_NC
    int ncout = 0;
#endif

#ifdef HAS_HDF5

#define STR_MAX  64

  int hdf5out = 0;
  hsize_t *hdf5_chunk=NULL;
  char hdf5_compress[STR_MAX];
  unsigned int compress_par[10];
#endif

#ifdef HAS_ADIOS
  char      *adios_groupname="struct";
  char      *adios_method=NULL;   /* POSIX|MPI|MPI_LUSTRE|MPI_AGGREGATE|PHDF5   */
  char      *adios_transform=NULL;    /* e.g., compression */
  struct adiosstructinfo adiosstruct_nfo;
  char      *adiosopts = NULL;
#endif

  /*## End of Output Module Variables ##*/
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /* Parse command line */
  for(a = 1; a < argc; a++) {

    if(!strcasecmp(argv[a], "--tasks")) {
      inp = atoi(argv[++a]);
      jnp = atoi(argv[++a]);
    } else if(!strcasecmp(argv[a], "--size")) {
      atoix(argv[++a], &ni, &niscale);
      atoix(argv[++a], &nj, &njscale);
      nk = atoi(argv[++a]);
    } else if(!strcasecmp(argv[a], "--maskthreshold")) {
      mask_thres = strtof(argv[++a], NULL);
      bot_mask_thres = (mask_thres+1.0)/2;
    } else if(!strcasecmp(argv[a], "--noisespacefreqmask")) {
      atodx(argv[++a], &noisefreqmask_i, &noisefreqmask_iscale);
      atodx(argv[++a], &noisefreqmask_j, &noisefreqmask_jscale);
    } else if(!strcasecmp(argv[a], "--noisespacefreq")) {
      atodx(argv[++a], &noisefreq_i, &noisefreq_iscale);
      atodx(argv[++a], &noisefreq_j, &noisefreq_jscale);
      noisefreq_k = strtod(argv[++a], NULL);
    } else if(!strcasecmp(argv[a], "--noisetimefreq")) {
      noisefreq_t = strtod(argv[++a], NULL);
    } else if(!strcasecmp(argv[a], "--tsteps")) {
      nt = atoi(argv[++a]);
    } else if(!strcasecmp(argv[a], "--tstart")) {
      tstart = atoi(argv[++a]);
    } else if(!strcasecmp(argv[a], "--debug")) {
      debug = 1;
    } else if(!strcasecmp(argv[a], "--debugBal")) {
      debugBal = 1;
    }
    else if(!strcasecmp(argv[a], "--balance")) {
      balance = 1;
    }

    /*## Add Output Modules' Command Line Arguments Here ##*/
    
#ifdef HAS_ADIOS
    else if(!strcasecmp(argv[a], "--adios")) {
      adios_method = argv[++a];
    }
    else if(!strcasecmp(argv[a], "--adios_transform")) {
      adios_transform = argv[++a];
    }
    else if(!strcasecmp(argv[a], "--adiosopts")) {
      adiosopts = argv[++a];
    }

    else if(!strcasecmp(argv[a], "--debugIO")) {
      debugIO = 1;
    }
#endif

    else if( !strcasecmp(argv[a], "--hdf5") )
      {
#ifdef HAS_HDF5
	hdf5out = 1;
      }
    else if(!strcasecmp(argv[a], "--hdf5_chunk")) {
      hdf5_chunk = malloc(2 * sizeof(hsize_t));
      hdf5_chunk[1] = (hsize_t)strtoul(argv[++a], NULL, 0);
      hdf5_chunk[0] = (hsize_t)strtoul(argv[++a], NULL, 0);
      if (hdf5_chunk[0] <= 0 || hdf5_chunk[1] <= 0 ) {
	print_usage(rank, "Error: Illegal chunk dim sizes");
	MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    else if(!strcasecmp(argv[a], "--hdf5_compress")) {
      char tmp[STR_MAX];
      strncpy(tmp, argv[++a], STR_MAX-1);
      char * pch;
      pch = strtok (tmp, " ,");
      strncpy( hdf5_compress, pch, strlen(pch) + 1 );
      pch = strtok (NULL, " ,");
      if ( strcmp(hdf5_compress,"szip") == 0 ) {
        int icnt = 0;
        while (pch != NULL) {
          if(icnt == 0) {
            if ( strcmp(pch,"H5_SZIP_EC_OPTION_MASK") == 0 ) {
              compress_par[icnt] = H5_SZIP_EC_OPTION_MASK;
            } else if ( strcmp(pch,"H5_SZIP_NN_OPTION_MASK") == 0 ) {
              compress_par[icnt] = H5_SZIP_NN_OPTION_MASK;
            } else {
              if(rank == 0) fprintf(stderr, "szip option not recognized: %s\n\n",pch);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
          } else if(icnt == 1) {
            compress_par[icnt] = (unsigned int)strtoul(pch, NULL, 0);
            /* pixels_per_block and must be even and not greater than 32 */
            if(compress_par[icnt] % 2 != 0 || compress_par[icnt] > 32 || compress_par[icnt] < 2 ) {
              if(rank == 0) fprintf(stderr, "szip pixels_per_block and must be even and not greater than 32: %d\n\n",compress_par[icnt]);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
          }
          pch = strtok (NULL, " ,");
          icnt++;
        }
      } else if ( strcmp(hdf5_compress,"gzip") == 0 || strcmp(hdf5_compress,"shuffle+gzip") == 0 ) {
        compress_par[0] = (unsigned int)strtoul(pch, NULL, 0);
      }
    }
#else
      if(rank == 0)   fprintf(stderr, "HDF5 option not available: %s\n\n", argv[a]);
      print_usage(rank, NULL);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
#endif

  else if (!strcasecmp(argv[a],"--nc") )
    {
#ifdef HAS_NC
      ncout = 1;
    }
#else
    if(rank == 0)  fprintf(stderr, "NC Option not available: %s\n\n", argv[a]);
    print_usage(rank, NULL);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
#endif

    /*## End of Output Module Command Line Arguments ##*/

    else {
      if(rank == 0)   fprintf(stderr, "Option not recognized: %s\n\n", argv[a]);
      print_usage(rank, NULL);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
   }

  /* Scale arguments specified to scale by ranks */
  if(niscale)   ni *= inp;
  if(njscale)   nj *= jnp;
  if(noisefreqmask_iscale)  noisefreqmask_i *= inp;
  if(noisefreqmask_jscale)  noisefreqmask_j *= jnp;
  if(noisefreq_iscale)  noisefreq_i *= inp;
  if(noisefreq_jscale)  noisefreq_j *= jnp;

  /* Check arguments & proc counts */
  if(inp < 1 || jnp < 1 ) {
    print_usage(rank, "Error: tasks not specified or incorrect");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if(ni <= 1 || nj <= 1 || nk <= 1) {
    print_usage(rank, "Error: size not specified or incorrect dim must be > 1");
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

  npoints  = (uint64_t)ni*nj*nk;

  /* Print options */
  if(rank == 0) {
    printf("struct options:\n");
    printf("      tasks: %d %d\n", inp, jnp);
    printf("      size: %d %d %d\n", ni, nj, nk);
    printf("      datasize: %llu bytes\n", npoints*4);
    printf("      maskthreshold: %f\n", mask_thres);
    printf("      noisespacefreqmask: %f %f\n", noisefreqmask_i, noisefreqmask_i);
    printf("      noisespacefreq: %f %f %f\n", noisefreq_i, noisefreq_j, noisefreq_k);
    printf("      noisetimefreq: %f\n", noisefreq_t);
    printf("      tsteps: %d\n", nt);
    printf("      tstart: %d\n", tstart);
    printf("      balance: %d\n", balance);
  }

  cpers[0] = 0;  cpers[1] = 0;  cpers[2] = 0;    /* No periodicity */
  crnk[0] = 0; crnk[1] = 0; crnk[2] = 0;

  /* if not balancing use  Cartesian communicator for proc decompostion */
  if (!balance) {
    /* Set up Cartesian communicator */
    cprocs[0] = inp;  cprocs[1] = jnp;  cprocs[2] = knp;
    MPI_Cart_create(MPI_COMM_WORLD, 3, cprocs, cpers, 1, &comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Cart_coords(comm, rank, 3, crnk);
  }

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

  /* adjust mask threshold  to compensate by bottom threshold */
  mask_thres = mask_thres - bot_mask_thres;
  mask_thres_index = (int) ( ((mask_thres+1)/2) * (nk-1));
  maskTindex = nk-1;

  /* Set up osn */
  open_simplex_noise(12345, &simpnoise);   /* Fixed seed, for now */

  /* Allocate arrays */
  if (balance) {
    timer_tick(&initheighttime, comm, 1);
    ol_mask_2dsum =  (int *) malloc((size_t)xy_dims*sizeof(int));

    /* initalize sums to zero */
    for(j = 0; j < nj; j++) {
      for(i = 0; i < ni; i++) {
	/* calculate index */
	x_index = i;
	y_index = j;
	z_index = k;
	sum_id = (y_index * x_dims) + x_index;
	ol_mask_2dsum[sum_id] = 0;
      }
    }
    sum_id = 0;

    /* Compute 2d sums   */
    y = ys;
    for(j = 0; j < nj; j++) {
      x = xs;
      for(i = 0; i < ni; i++) {

        /* calculate index */
        x_index = i;
        y_index = j;
        sum_id = (y_index * x_dims) + x_index;

        /* Get height and subtract bottom threshold */
        h_tmp =  (float)open_simplex_noise2(simpnoise, x*noisefreqmask_i, y*noisefreqmask_j)  - bot_mask_thres;

        hindex = (int) ((h_tmp+1) / (mask_thres+1) * (nk-1));

        if (hindex > maskTindex) {
          hindex = maskTindex;
        }

        for(k = 0; k < nk; k++) {
          z_index = k;

          /* Calculate ol_mask values */
          if (z_index < maskTindex  && z_index > hindex) {
            tmp_mask = 1;  /* ocean */
          }
          else if (z_index <= hindex) {
            if (h_tmp >= mask_thres  || z_index < hindex) {
              tmp_mask = 0;  /* land */
            }
            else {
              tmp_mask = 1;  /* ocean */
            }
          }
          else if (z_index == maskTindex  && h_tmp <= mask_thres) {
            tmp_mask = 0;  /* ocean */
          }
          else {
            printf("WARNING: ol_mask condition not considered for Point_index: (%d,%d,%d)\n"
        	   "Point_id: %d  Height: %f HeightID: %d maskTindex=%d\n",
        	   x_index, y_index, z_index, sum_id+1, h_tmp, hindex, maskTindex);
          }
          ol_mask_2dsum[sum_id] += tmp_mask;
	}
        
        x += deltax;
      }
      y += deltay;
    }
    timer_tock(&initheighttime);
    timer_collectprintstats(initheighttime, comm, 0, "Initial Height Field");

    if (rank == 0 && debug) {
      
      printf("Struct sum of layers %d to %d\n", 0, nk-1);
      
      for(i = 0; i < ni; i++) {
	for(j = 0; j < nj; j++) {
	  /* calculate index */
	  x_index = i;
	  y_index = j;
	  z_index = k;
	  sum_id = (y_index * x_dims) + x_index;
	  
	  if( j == 0) {
	    printf("[ %d ",  ol_mask_2dsum[sum_id]);
	  }
	  else  {
	    printf("%d ",  ol_mask_2dsum[sum_id]);
	  }
	}
	printf("]\n");
      }
      printf("\n");
    }

    /*if (rank == 0) {
      FILE *f;
      f = fopen("struct.olmask.dat", "w");
      fwrite(ol_mask_2dsum, sizeof(int), ni*nj, f);
      fclose(f);
    } */

    /* Use splitspace for proc decompostion */   
    timer_tick(&loadbaltime, comm, 1);
    init_splitspace(&factors, &rects_list, nprocs);
    splitspace(&rects_list, ni, nj, &factors, ol_mask_2dsum);
    timer_tock(&loadbaltime);
    timer_collectprintstats(loadbaltime, comm, 0, "Load Balancing");

    if (rank == 0 && debugBal) {
      print_rects(&rects_list);
    }

    /* set data bound based on rank */
    if (rank == 0 && debugBal) {
      printf("Proc %d has rect: (%d,%d:%d,%d)\n", rank, rects_list.recs[rank].x0, rects_list.recs[rank].y0,
	     rects_list.recs[rank].x1, rects_list.recs[rank].y1);
    }
    free (ol_mask_2dsum);
    
    is = rects_list.recs[rank].x0;
    js = rects_list.recs[rank].y0;
    ks = 0;
    cni = rects_list.recs[rank].x1+1 - rects_list.recs[rank].x0;
    cnj = rects_list.recs[rank].y1+1 - rects_list.recs[rank].y0;
    cnk = nk; 
    xs = is * deltax;
    ys = js * deltay;
    zs = ks * deltaz;

    free_splitspace (&factors, &rects_list);
  }

  /* Allocate arrays */
  data = (float *) malloc((size_t)cni*cnj*cnk*sizeof(float));
  height = (float *) malloc((size_t)cni*cnj*cnk*sizeof(float));
  ola_mask = (int *) malloc((size_t)cni*cnj*cnk*sizeof(int));
  ol_mask = (int *) malloc((size_t)cni*cnj*cnk*sizeof(int));

  varnames[0] = "data";
  /*varnames[1] = "height";
  varnames[2] = "ola_mask";
  varnames[3] = "ol_mask"; */

  /*## Add Output Modules' Initialization Here ##*/
  
 /* init ADIOS */
#ifdef HAS_ADIOS
  if (adios_method) {
    adiosstruct_init(&adiosstruct_nfo, adios_method, adios_transform, adios_groupname, 
                     adiosopts, comm, rank, nprocs, nt, ni, nj, nk, is, cni, js, cnj, ks, cnk, 
                     deltax, deltay, deltaz, FILLVALUE);
    adiosstruct_addrealxvar(&adiosstruct_nfo, varnames[0], data);

    if (debugIO) {
      adiosstruct_addrealxvar(&adiosstruct_nfo, "height", height);
      adiosstruct_addintxvar(&adiosstruct_nfo, "ola_mask", ola_mask);
      adiosstruct_addintxvar(&adiosstruct_nfo, "ol_mask", ol_mask);
    }
  }
#endif

  /*## End of Output Module Initialization ##*/
 
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
	point_id = ((uint64_t)z_index * xy_dims) + (y_index * x_dims) + x_index;

	/* Get height and subtract bottom threshold */
	height[ii] =  (float)open_simplex_noise2(simpnoise, x*noisefreqmask_i, y*noisefreqmask_j)  - bot_mask_thres;

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
		 "Point_id: %llu  Height: %f HeightID: %d  mask_thres_index=%d\n",
		 x_index, y_index, z_index, point_id+1, height[ii], height_index, mask_thres_index);
	}

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
		 "Point_id: %llu  Height: %f HeightID: %d maskTindex=%d\n",
		 x_index, y_index, z_index, point_id+1, height[ii], hindex, maskTindex);
	}

	if (debug) {
	  printf("++++++++++++++++++++++++++++++++++++++++++++\n");
	  printf("rank_cord(%d,%d,%d) rank=%d: %d of %d\n", crnk[0], crnk[1], crnk[2] , rank, rank+1, nprocs);
	  printf("LDims: (%d,%d,%d)\n", cni, cnj, cnk);
	  printf("GDims: (%d,%d,%d)\n", ni, nj, nk);
	  printf("SDims: (%d,%d,%d)\n", is, js, ks);
	  printf("Point_index: (%d,%d,%d), Point_id:  %llu\n", x_index, y_index, z_index,  point_id+1);
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
	    data[ii] = (float)open_simplex_noise4(simpnoise, x*noisefreq_i, y*noisefreq_j, z*noisefreq_k, tt*noisefreq_t);
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

    /*## Add OUTPUT Modules' Function Calls Per Timestep Here ##*/
    
#ifdef HAS_ADIOS
    if (adios_method) {
       if(rank == 0) {
	printf("      Writing Adios...\n");   fflush(stdout);
       }
       adiosstruct_write(&adiosstruct_nfo, tt);
    }
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
		data, hdf5_chunk, hdf5_compress, compress_par);
    }
#endif

#ifdef HAS_NC
    if(ncout) {
      if(rank == 0) {
        printf("     Writing netCDF...\n"); fflush(stdout);
      }
      writenc(num_varnames, varnames, comm, rank, nprocs, tt,
              is, js, ks,
              ni, nj, nk, cni, cnj, cnk,
              data);
    }
#endif

    /*## End of OUTPUT Module Function Calls Per Timestep ##*/

    timer_tock(&outtime);
    timer_collectprintstats(computetime, comm, 0, "   Compute");
    timer_collectprintstats(outtime, comm, 0, "   Output");
    timer_collectprintstats(heighttime, comm, 0, "   Height");

  }

  /*## Add Output Modules' Cleanup Here ##*/
    
  /* finalize ADIOS */
#ifdef HAS_ADIOS
  if (adios_method) adiosstruct_finalize(&adiosstruct_nfo);
#endif

#ifdef HAS_HDF5
  if(hdf5_chunk)
    free(hdf5_chunk);
#endif

  /*## End of Output Module Cleanup ##*/
  
  open_simplex_noise_free(simpnoise);
  free(data);
  free(height);
  free(ola_mask);
  free(ol_mask);

  MPI_Finalize();

  return 0;
}


