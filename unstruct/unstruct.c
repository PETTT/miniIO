#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <strings.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <mpi.h>
#include <inttypes.h>
#include "timer.h"
#include "open-simplex-noise.h"

/*## Add Output Modules' Includes Here ##*/

#ifdef HAS_PRZM
#  include "przm.h"
#endif

#ifdef HAS_ADIOS
#  include "adiosunstruct.h"
#endif

#ifdef HAS_NC
#  include "ncunstruct.h"
#endif

#ifdef HAS_HDF5
#  include "hdf5unstruct.h"

#  ifdef H5Z_ZFP_USE_PLUGIN
#     include "H5Zzfp_plugin.h"
#  endif


#define CD_VALUES_MAX 10

    size_t cd_nelmts=CD_VALUES_MAX;
    unsigned int cd_values[CD_VALUES_MAX];

#endif

/*## End of Output Module Includes ##*/

void print_usage(int rank, const char *errstr)
{
    if(rank != 0)  return;
    if(errstr)
        fprintf(stderr, "%s\n\n", errstr);
    fprintf(stderr,
"Usage: mpi_launcher [-n|-np NPROCS] ./unstruct --points PTS [options]\n"
"     OR\n"
"       mpi_launcher [-n|-np NPROCS] ./unstruct --pointspertask PTST [options]\n"
"    NPROCS : # of tasks launched by MPI; may or may not be implied or required by system\n\n"
"  Required (one or the other, but not both):\n"
"    --points PTS : Specifies total number of points for all tasks\n"
"    --pointspertask PTST : Specified number of points per single task\n\n"
"  Optional:\n"
"    --roundness UR VR : Shape of superquadric for base grid\n"
"       0 0 is a cube, 1 1 is a sphere, 2 2 is an octahedron, >2 >2 is increasingly concave\n"
"       Defaults: 0.3 0.3\n"
"    --animroundness UR VR : Shape of superquadric at final time step\n"
"       If specified, linearly interpolates over time from starting roundness\n"
"       Defaults: no change from initial roundness\n"
"    --tsteps NT : Total number of time steps\n"
"       Default: 50\n"
"    --noisespacefreq FNS : Spatial frequency of noise function\n"
"      FNS : space frequency value; Default: 10.0\n"
"    --noisetimefreq FNT : Temporal frequency of noise function\n"
"      FNT : time frequency value;  Default: 0.25\n"
    );

    /*## Add Output Modules' Usage String ##*/

#ifdef HAS_PRZM
    fprintf(stderr, "   --przm : Enable PRZM output.\n");
#endif
#ifdef HAS_HDF5
    fprintf(stderr, "   --hdf5 : Enable HDF5 output.\n"
                    "   --hdf5_chunk x \n"
                    "      values of chunk size; x are pointspertask/x, must be divisible\n"
                    "      setting x to zero disables chunking\n"
                    "    --hdf5_compress : enable compression. Valid value is a comma seperated (no spaces) list:  \n"
                    "        <compression type: gzip or szip>,<compression parameter(s) corresponding to HDF5 compression API>  \n"
                    "        gzip,<value is level (see H5Pset_deflate)> \n"
                    "        szip,<value is <options_mask>,<pixels_per_block (see H5Pset_szip)> \n"
                    "             For example, --hdf5_compress szip,H5_SZIP_NN_OPTION_MASK,8 \n"
                    "        zfp,<value is zfpmode,[zfp dependant parameter(s)]> \n"
                    "        NOTE: compression requires chunked datasets \n"
            
	    );
#endif
#ifdef HAS_NC
    fprintf(stderr, "    --nc : Enable netCDF Output \n");
    fprintf(stderr, "    --nclog [log level] : Enable netCDF log level\n");
#endif

    /*## End of Output Module Usage Strings ##*/
}

/* Split n into prime factors then recombine into factors n1 & n2 */

void prime_split(int n, int *n1, int *n2)
{
	int i, sq=(int)sqrt(n);
    int vpos = 0;
    int *v;
    int p1 = 1, p2 = n;
    int orign = n;

    /* Trivial cases */
    if(n < 1) {
        *n1 = *n2 = 0;
        return;
    } else if(n == 1) {
        *n1 = *n2 = 1;
    }

    v = (int *) malloc((sq+2)*sizeof(int));   /* List of prime factors */

    v[vpos++] = 1;   /* Better for cases of prime numbers, incl 2 */

    /* Prime seive */
	while(n%2 == 0) {
		v[vpos++] = 2;
        n /= 2;
    }
	for(i = 3; i <= sq; i += 2)	{
		while(n%i==0) {
            v[vpos++] = i;
			n /= i;
        }
	}
	if(n > 1) {   /* For prime numbers */
        v[vpos++] = n;
    }

    /* Bisect the list of primes iteratively, compute 2 factors at the split:
     *    n1*n2 is current bisection, p1*p2 is previous */
    for(i = 1; i < vpos; ++i) {   /* Try the i-th bisection */
        int j;
        *n1 = 1;  *n2 = 1;
        for(j = 0; j < i; ++j)  *n1 *= v[j];
        for(j = i; j < vpos; ++j) *n2 *= v[j];
        if(*n1 >= *n2)  break;    /* Stop when 2nd factor is greater or equal */
        p1 = *n1;  p2 = *n2;   /* Otherwise, set previous factors */
    }
    if(*n1 - *n2 > p2 - p1) {   /* Whichever set of factors is closest together wins */
        *n1 = p1;  *n2 = p2;
    }

    free(v);

    assert(orign == *n1 * *n2);

    if(*n2 > *n1) {   /* Make sure n1 >= n2, else swap */
        p1 = *n1;   *n1 = *n2;   *n2 = p1;
    }
}

/* Signum function */
double sgn(double x)
{
    return copysign(1.0, x);
}

/* Superquadric ellipsoidal "c" function */
float sqc(double w, double m)
{
    return (float)(sgn(cos(w)) * pow(fabs(cos(w)), m));
}

/* Superquadric ellipsoidal "s" function */
float sqs(double w, double m)
{
    return (float)(sgn(sin(w)) * pow(fabs(sin(w)), m));
}


int main(int argc, char **argv)
{
    int a;
    uint64_t npoints = 0;    /* Total number of points */
    uint64_t nptstask = 0;     /* Number of point per task */
    int nu, nv, nlyr;    /* Points per task per u,v,lyr spherical coord axis */
    float du, dv;          /* delta's along u & v points */
    float u0, u1, v0, v1;      /* Starting/ending points along u & v */
    int i, j, k, t;
    int nt = 50;               /* Number of time steps */
    uint64_t ii;
    float *xpts, *ypts, *zpts;    /* Grid points */
    uint64_t nelems2, *conns2;      /* Number of grid triangles & connection array in 2D */
    uint64_t nelems3, *conns3;      /* Number of triangular prisms & connection array */
    uint64_t datasize;              /* Size of all grid and variable expected to be output */
    float *data;                  /* Data array */
    float uround0 = 0.3f;        /* Superquadric roundness u parameter, starting */
    float vround0 = 0.3f;        /* Superquadric roundness v parameter, starting */
    float uround1 = -1.f;       /* Superquadric roundness u, ending over time */
    float vround1 = -1.f;       /* Superquadric roundness v, ending over time */
    float uround, vround;       /* Current superquadric roundness */
    double noisespacefreq = 10;    /* Spatial frequency of noise */
    double noisetimefreq = 0.25;    /* Temporal frequency of noise */
    struct osn_context *osn;    /* Open simplex noise context */
    double gridtime, computetime, outtime;    /* Timers */

    /* MPI vars */
    int rank, nprocs;
    int uprocs, vprocs;      /* Number of tasks on u & v spherical axes respectively */
    int urank, vrank;        /* 2D rank along u & v */

    /*## Add Output Modules' Variables Here ##*/

#ifdef HAS_PRZM
    int przmout = 0;
#endif

#ifdef HAS_ADIOS
    char *adiosmethod = NULL;
    struct adiosinfo adiosnfo;
#endif

#ifdef HAS_NC
    int ncloglevel = 0;
    int ncout = 0;
#endif

#ifdef HAS_HDF5

#define STR_MAX 128
#define CD_VALUES_MAX 10

    int hdf5out = 0;
    hsize_t *hdf5_chunk=NULL;
    char hdf5_compress[STR_MAX];
    hdf5_compress[0] = '\0';

    int filter_id = 0;
    size_t cd_nelmts=CD_VALUES_MAX;
    unsigned int compress_par[CD_VALUES_MAX];
#endif

    /*## End of Output Module Variables ##*/

    /* Init MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Parse command line */
    for(a = 1; a < argc; a++) {
        if(!strcasecmp(argv[a], "--points")) {
            npoints = strtoull(argv[++a], NULL, 10);
            nptstask = 0;
        } else if(!strcasecmp(argv[a], "--pointspertask")) {
            nptstask = strtoull(argv[++a], NULL, 10);
            npoints = 0;
        } else if(!strcasecmp(argv[a], "--roundness")) {
            uround0= atof(argv[++a]);
            vround0= atof(argv[++a]);
        } else if(!strcasecmp(argv[a], "--animroundness")) {
            uround1= atof(argv[++a]);
            vround1= atof(argv[++a]);
        } else if(!strcasecmp(argv[a], "--tsteps")) {
            nt = atoi(argv[++a]);
        } else if(!strcasecmp(argv[a], "--noisespacefreq")) {
            noisespacefreq = strtod(argv[++a], NULL);
        } else if(!strcasecmp(argv[a], "--noisetimefreq")) {
            noisetimefreq = strtod(argv[++a], NULL);
        }

        /*## Add Output Modules' Command Line Arguments Here ##*/

#ifdef HAS_PRZM
        else if(!strcasecmp(argv[a], "--przm")) {
            przmout = 1;
        }
#endif

#ifdef HAS_ADIOS
        else if(!strcasecmp(argv[a], "--adios")) {
            adiosmethod = argv[++a];
        }
#endif

#ifdef HAS_HDF5
        else if(!strcasecmp(argv[a], "--hdf5")) {
            hdf5out = 1;
        }
        else if(!strcasecmp(argv[a], "--hdf5_chunk")) {
	  hdf5_chunk = malloc(1 * sizeof(hsize_t));
	  hdf5_chunk[0] = (hsize_t)strtoul(argv[++a], NULL, 0);
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
                  cd_values[icnt] = H5_SZIP_EC_OPTION_MASK;
                } else if ( strcmp(pch,"H5_SZIP_NN_OPTION_MASK") == 0 ) {
                  cd_values[icnt] = H5_SZIP_NN_OPTION_MASK;
                } else {
                  if(rank == 0) fprintf(stderr, "szip option not recognized: %s\n\n",pch);
                  MPI_Abort(MPI_COMM_WORLD, 1);
                }
              } else if(icnt == 1) {
                cd_values[icnt] = (unsigned int)strtoul(pch, NULL, 0);
                /* pixels_per_block and must be even and not greater than 32 */
                if(cd_values[icnt] % 2 != 0 || cd_values[icnt] > 32 || cd_values[icnt] < 2 ) {
                  if(rank == 0) fprintf(stderr, "szip pixels_per_block and must be even and not greater than 32: %d\n\n",cd_values[icnt]);
                  MPI_Abort(MPI_COMM_WORLD, 1);
                }
              }
              pch = strtok (NULL, " ,");
              icnt++;
            }
            filter_id = H5Z_FILTER_SZIP;
            cd_nelmts = 2;
          } else if ( strcmp(hdf5_compress,"gzip") == 0 ) {
            filter_id = H5Z_FILTER_DEFLATE;
            cd_nelmts = 1;
            cd_values[0] = (unsigned int)strtoul(pch, NULL, 0);
          }
#ifdef H5Z_ZFP_USE_PLUGIN
          else if ( strcmp(hdf5_compress,"zfp") == 0 ) {

            filter_id = H5Z_FILTER_ZFP;

            /* setup zfp filter via generic (cd_values) interface */
            /*
              cd_vals    0       1        2         3         4         5
              ----------------------------------------------------------------
              rate:      1    unused    rateA     rateB     unused    unused
              precision: 2    unused    prec      unused    unused    unused
              accuracy:  3    unused    accA      accB      unused    unused
              expert:    4    unused    minbits   maxbits   maxprec   minexp
              reversubke 5    unused    unused    unused    unused    unused
            */

            int zfpmode=0;

            if (strcmp(pch,"H5Z_ZFP_MODE_RATE") == 0 ) {
              zfpmode = 1;
            } else if (strcmp(pch,"H5Z_ZFP_MODE_PRECISION") == 0 ) {
              zfpmode = 2;
            } else if (strcmp(pch,"H5Z_ZFP_MODE_ACCURACY") == 0 ) {
              zfpmode = 3;
            } else if (strcmp(pch,"H5Z_ZFP_MODE_EXPERT") == 0 ) {
              zfpmode = 4;
            } else if (strcmp(pch,"H5Z_ZFP_MODE_REVERSIBLE") == 0 ) {
              zfpmode = 5;
            }

            pch = strtok (NULL, " ,");

            if (zfpmode == H5Z_ZFP_MODE_RATE) {
              double rate = (double)strtoul(pch, NULL, 0);
              H5Pset_zfp_rate_cdata(rate, cd_nelmts, cd_values);
            } else if (zfpmode == H5Z_ZFP_MODE_PRECISION) {
              uint prec = (uint)strtoul(pch, NULL, 0);
              H5Pset_zfp_precision_cdata(prec, cd_nelmts, cd_values);
            } else if (zfpmode == H5Z_ZFP_MODE_ACCURACY) {
              double acc = (double)strtoul(pch, NULL, 0);
              H5Pset_zfp_accuracy_cdata(acc, cd_nelmts, cd_values);
            } else if (zfpmode == H5Z_ZFP_MODE_EXPERT) {
              uint minbits = (uint)strtoul(pch, NULL, 0);
              pch = strtok (NULL, " ,");
              uint maxbits = (uint)strtoul(pch, NULL, 0);
              pch = strtok (NULL, " ,");
              uint maxprec = (uint)strtoul(pch, NULL, 0);
              pch = strtok (NULL, " ,");
              int minexp = (int)strtoul(pch, NULL, 0);
              H5Pset_zfp_expert_cdata(minbits, maxbits, maxprec, minexp, cd_nelmts, cd_values);
            } else if (zfpmode == H5Z_ZFP_MODE_REVERSIBLE) {
              H5Pset_zfp_reversible_cdata(cd_nelmts, cd_values);
            } else {
              cd_nelmts = 0; /* causes default behavior of ZFP library */
            }
            if(rank == 0) {
              printf("%d cd_values= ",cd_nelmts);
              for (int i = 0; i < cd_nelmts; i++)
                printf("%u,", cd_values[i]);
              printf("\n");
            }
          }

          if(filter_id != 0) {
            /* Check if filter is available */
            htri_t avail = H5Zfilter_avail(filter_id);
            if (avail <= 0) {
              printf("error: %s compression is not available\n", hdf5_compress);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
          }
#endif
	}
#endif

#ifdef HAS_NC
        else if(!strcasecmp(argv[a],"--nc")) {
          ncout = 1;
        }
        else if(!strcasecmp(argv[a],"--nclog")) {
          ncloglevel = (int)strtoimax(argv[++a], NULL, 0);
        }

#endif

        /*## End of Output Module Command Line Arguments ##*/

        else {
            if(rank == 0)  fprintf(stderr, "Option not recognized: %s\n\n", argv[a]);
            print_usage(rank, NULL);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    /* Check arguments */
    if(nprocs % 2 > 0 && nprocs > 1) {
        print_usage(rank, "Error: NPROCS is not even; that is required for load balancing");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    if(npoints == 0 && nptstask == 0) {
        print_usage(rank, "Error: neither points or pointspertask specified, or there was\n"
                          "       an error parsing them");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if(npoints % nprocs > 0) {
        print_usage(rank, "Error: points must be evenly divisible by NPROCS for best\n"
                          "       load balancing");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if(uround1 == -1.f)   uround1 = uround0;
    if(vround1 == -1.f)   vround1 = vround0;

    /* Determine a volumetric spherical topology that meets # of points requested */
    prime_split(nprocs, &uprocs, &vprocs);   /* Split tasks across spherical axes */
    if(npoints)   /* split npoints among tasks */
        nptstask = npoints / nprocs;
    nu = (int) ceil(cbrt(nptstask * uprocs * 2 / vprocs));
                  /* need uprocs/vprocs proportional; want 2x nu's as layers */
    nv = nu * vprocs / uprocs;
    nlyr = nu / 2;
    nptstask = (uint64_t)nu * nv * nlyr;   /* nptstask won't be exactly as requested */
    npoints = nptstask * nprocs;
    if(rank == 0)
        printf("Actual points: %"PRIu64" , points per task: %"PRIu64"\n"
                "uprocs: %d , vprocs: %d\n"
                "nu: %d , nv: %d , nlyr: %d\n",
                npoints, nptstask, uprocs, vprocs, nu, nv, nlyr);

    /* Divide spherical topology tiles along u,v */
    du = 2 * M_PI / uprocs / (nu-1);
    dv = M_PI / vprocs / (nv-1);
    urank = rank % uprocs;
    vrank = rank / uprocs;
    u0 = urank * 2 * M_PI / uprocs - M_PI;
    u1 = (urank+1) * 2 * M_PI / uprocs - M_PI;   /* #points to be continuous across task */
    v0 = vrank * M_PI / vprocs - M_PI/2;
    v1 = (vrank+1) * M_PI / vprocs - M_PI/2;
    /*DBG printf("%d: %f-%f, %f-%f, %f, %f\n", rank, u0, u1, v0, v1, du, dv); */

    /* Allocate grid points */
    xpts = (float *) malloc(nptstask*sizeof(float));
    ypts = (float *) malloc(nptstask*sizeof(float));
    zpts = (float *) malloc(nptstask*sizeof(float));
    datasize = nptstask*sizeof(float)*3;

    /* Add surface connections, all are triangles */
    nelems2 = (nu-1) * (nv-1) * 2;
    conns2 = (uint64_t *) malloc(nelems2*3*sizeof(uint64_t));
    datasize += nelems2*3*sizeof(uint64_t);
    for(i = 0, ii = 0; i < nu-1; i++) {
        for(j = 0; j < nv-1; j++) {
            uint64_t ndx = rank * nptstask + i * nv + j;   /* Global point index */
            conns2[ii++] = ndx;
            conns2[ii++] = ndx + nv + 1;
            conns2[ii++] = ndx + nv;
            conns2[ii++] = ndx;
            conns2[ii++] = ndx + 1;
            conns2[ii++] = ndx + nv + 1;
        }
    }

    /* Add volume connections, all are triangular prisms extended from base 2D grid */
    nelems3 = nelems2 * (nlyr-1);
    conns3 = (uint64_t *) malloc(nelems3*6*sizeof(uint64_t));
    datasize += nelems3*6*sizeof(uint64_t);
    for(k = 0, ii = 0; k < nlyr-1; k++) {
        uint64_t euv, iuv, nuv = nu * nv * k, nuv2 = nu * nv * (k+1);
        for(euv = 0, iuv = 0; euv < nelems2; euv++) {
            conns3[ii++] = conns2[iuv] + nuv;
            conns3[ii++] = conns2[iuv+1] + nuv;
            conns3[ii++] = conns2[iuv+2] + nuv;
            conns3[ii++] = conns2[iuv++] + nuv2;
            conns3[ii++] = conns2[iuv++] + nuv2;
	    conns3[ii++] = conns2[iuv++] + nuv2;
        }
    }

    /* Set up osn */
    open_simplex_noise(12345, &osn);   /* Fixed seed, for now */

    /* Allocate data */
    data = (float *) malloc(nptstask*sizeof(float));
    datasize += nptstask*sizeof(float);

    if(rank == 0) {
        printf("nelems2: %"PRIu64" , nelems3: %"PRIu64"\n"
               "data size per task = %"PRIu64" , all tasks = %"PRIu64"\n",
               nelems2, nelems3, datasize, datasize*nprocs);
#ifdef HAS_HDF5
	printf("Enable HDF5: ");
	if(hdf5out)
	  printf("yes\n");
	else
	  printf("no\n");

	if(hdf5_chunk)
	  printf("Enable HDF5 chunking: %lld \n", hdf5_chunk[0]);
	
	printf("Enable HDF5 compression: ");
	if(strcasecmp(hdf5_compress, "\0") != 0)
	  printf("yes, %s\n", hdf5_compress);
	else
	  printf("no\n");
#endif
    }
    /*## Add Output Modules' Initialization Here ##*/

#ifdef HAS_ADIOS
    if(adiosmethod) {
        adiosunstruct_init(&adiosnfo, adiosmethod, "unstruct.out", MPI_COMM_WORLD,
                           rank, nprocs, nt, nptstask, nelems3, nelems2);
        adiosunstruct_addvar(&adiosnfo, "noise");
    }
#endif

    /*## End of Output Module Initialization ##*/

    /* Main loops */
    for(t = 0; t < nt; t++) {
        float tpar = (float)t / (nt-1);    /* Time anim. parameter [0,1] */
        if(rank == 0) {
            printf("t = %d\n", t);   fflush(stdout);
        }

        timer_tick(&gridtime, MPI_COMM_WORLD, 1);

        /* Generate grid points with superquadric */
        uround = (1-tpar)*uround0 + tpar*uround1;
        vround = (1-tpar)*vround0 + tpar*vround1;
        if( uround0 != uround1 || t == 0 ) {   /* Grid only updated if animating or time 0 */
            for(k = 0, ii = 0; k < nlyr; k++) {
                /* layer w in [1,3] by squares: */
                float w = 1.f + powf((float)k/(nlyr-1), 2.f) * 2.f;
                for(i = 0; i < nu; i++) {
                    float u = u0 + du*i;
                    for(j = 0; j < nv; j++, ii++) {
                        float v = v0 + dv*j;
                        xpts[ii] = w * sqc(v, vround) * sqc(u, uround);
                        ypts[ii] = w * sqc(v, vround) * sqs(u, uround);
                        zpts[ii] = w * sqs(v, vround);
                    }
                }
            }
        }

        timer_tock(&gridtime);
        timer_tick(&computetime, MPI_COMM_WORLD, 1);

        /* Set data */
        for(i = 0; i < nptstask; i++) {
            data[i] = (float)open_simplex_noise4(osn, xpts[i]*noisespacefreq,
                    ypts[i]*noisespacefreq, zpts[i]*noisespacefreq, t*noisetimefreq);
        }

        timer_tock(&computetime);
        timer_tick(&outtime, MPI_COMM_WORLD, 1);

        /*## Add Output Modules' Function Calls Per Timestep Here ##*/

#ifdef HAS_PRZM
        if(przmout) {
            if(rank == 0) {
                printf("   Writing przm...\n");   fflush(stdout);
            }
            writeprzm("unstruct", MPI_COMM_WORLD, t, nptstask, xpts, ypts, zpts,
                      nelems3, conns3, nelems2, conns2, "noise", data);
        }
#endif

#ifdef HAS_ADIOS
        if(adiosmethod) {
            if(rank == 0) {
                printf("   Writing adios...\n");    fflush(stdout);
            }
            adiosunstruct_write(&adiosnfo, t, xpts, ypts, zpts, conns3, conns2, &data);
        }
#endif

#ifdef HAS_HDF5
        if(hdf5out) {
            if(rank == 0) {
                printf("   Writing hdf5...\n");   fflush(stdout);
            }
            writehdf5("unstruct", MPI_COMM_WORLD, t, npoints, nptstask, xpts, ypts, zpts,
                      nelems3, conns3, nelems2, conns2, "noise", data, hdf5_chunk, filter_id, cd_nelmts, cd_values);


        }
#endif

#ifdef HAS_NC
        nc_set_log_level(ncloglevel);

        if(ncout) {
          if(rank == 0) {
            printf("    Writing netCDF...\n"); fflush(stdout);
          }

          writenc("unstruct", MPI_COMM_WORLD, t, npoints, nptstask, xpts, ypts, zpts,
                  nelems3, conns3, nelems2, conns2, "noise", data);
        }
#endif

        /*## End of Output Module Functions Calls Per Timestep ##*/

        timer_tock(&outtime);
        timer_collectprintstats(gridtime, MPI_COMM_WORLD, 0, "   Grid");
        timer_collectprintstats(computetime, MPI_COMM_WORLD, 0, "   Compute");
        timer_collectprintstats(outtime, MPI_COMM_WORLD, 0, "   Output");
    }

    /*## Add Output Modules' Cleanup Here ##*/

#ifdef HAS_ADIOS
    if(adiosmethod)
        adiosunstruct_finalize(&adiosnfo);
#endif

    /*## End of Output Module Cleanup ##*/

    /* Cleanup */
    open_simplex_noise_free(osn);
    free(xpts);  free(ypts);  free(zpts);
    free(conns2);  free(conns3);
    free(data);
    MPI_Finalize();

    return 0;
}
