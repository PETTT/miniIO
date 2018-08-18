/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <stdint.h>
#include <math.h>
#include <mpi.h>
#include "sfc.h"
#include "iso.h"
#include "timer.h"
#include "open-simplex-noise.h"

/*## Add Output Modules' Includes Here ##*/

#ifdef HAS_PVTI
#  include "pvti.h"
#endif

#ifdef HAS_PVTP
#  include "pvtp.h"
#endif

#ifdef HAS_ADIOS
#  include "adiosfull.h"
#  include "adiosiso.h"
#endif

#ifdef HAS_NC
#  include "nccartiso.h"
#endif

#ifdef HAS_HDF5
#  include "hdf5cartiso.h"
#endif

/*## End of Output Module Includes ##*/


void print_usage(int rank, const char *errstr)
{
    if(rank != 0)  return;
    if(errstr)
        fprintf(stderr, "%s\n\n", errstr);
    fprintf(stderr,
"Usage: mpi_launcher [-n|-np NPROCS] ./cartiso --tasks INP JNP KNP --size NI NJ NK [options]\n"
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
"    --sigma SX SY SZ : Specifies the (initial) standard deviations of the Gaussian\n"
"                       function\n"
"      SX, SY, SZ : Sigma values for X, Y, Z axes respectively\n"
"      Defaults: SX = SY = SZ = 0.2\n"
"      For gaussmove and sin2gauss modes, this is the sigma value\n"
"      For gaussresize mode, this is the initial sigma value\n"
"    --sigmaend SX1 SY1 SZ1 : Specifies the final standard deviation of the Gaussian\n"
"                             for gaussresize\n"
"      SX1, SY1, SZ1 : Final sigma values for X, Y, Z axes respectively\n"
"      Defaults: SX1 = SY1 = SZ1 = 0.4\n"
"    --center X0 Y0 Z0 : Specifies the center of the Gaussian function for sin2gauss\n"
"                        or gaussresize\n"
"      X0, Y0, Z0 : Coordinates of center of Gaussian for X, Y, Z axes\n"
"      Valid ranges are from 0.0 to 1.0\n"
"      Is overridden by --centertask or --gaussmove\n"
"    --centertask : Centers the Gaussian function in the middle of the 1st task's\n"
"                   domain\n"
"      Is overridden by --center (whichever is last), or forced by --gaussmove\n"
"    --amp AMP : Amplitude of the sinusoid function from the top (1.0)\n"
"      AMP : Sinusoid tops out at 1.0, so AMP is the distance it goes down from 1.0\n"
"      Valid ranges are from 0.0 to 1.0\n"
"      Default: near the value of gaussian of 1*sigma, i.e. e^(-0.5)\n"
"    --freq FX FY FZ : Frequency of sinusoid function\n"
"      FX, FY, FZ : Indicates the number of sinusoid periods across the domain\n"
"    --noisespacefreq FNX FNY FNZ : Spatial frequency of noise function\n"
"      FNS : space frequency value; Default: 10.0\n"
"    --noisetimefreq FNT : Temporal frequency of noise function\n"
"      FNT : time frequency value;  Default: 0.25\n"
"    --noiserange NMIN NMAX : Specify min & max values of noise; Default: 0 1\n"
"    --fulloutputnoiseonly : Only output the noise field for full output\n"
"    --tsteps NT : Number of time steps; valid values are > 0\n"
"    --tstart TS : Starting time step; valid values are > 0\n"
"    --sin2gauss : Mode that 'morphs' between a sinusoid and a Gaussian over all\n"
"                  time steps.  This is the default if no mode is specified\n"
"    --gaussmove : Mode that moves a Gaussian through the spatial domain of all tasks\n"
"    --gaussresize : Mode that resizes a Gaussian from start to end sigma sizes\n"
"    --backward : Reverses the direction and starting point of gaussmove mode\n"
    );

    /*## Add Output Modules' Usage String ##*/

#ifdef HAS_PVTI
    fprintf(stderr, "    --pvti : Enable PVTI full output.\n");
#endif

#ifdef HAS_PVTP
    fprintf(stderr, "    --pvtp : Enable PVTP isosurface output.\n");
#endif

#ifdef HAS_ADIOS
    fprintf(stderr, 
       "    --adiosfull METHOD : Enable ADIOS full output with corresponding method.\n"
       "    --adiosiso METHOD : Enable ADIOS isosurface output with corresponding method.\n"
       "    --adiosopts OPTS : Pass options to ADIOS.\n"
       "    --adios_transform TRANSFORM : Pass a transform, e.g., compression, to ADIOS\n"
    );
#endif

#ifdef HAS_HDF5
    fprintf(stderr, "    --hdf5i : Enable HDF5 full output.\n"
                    "    --hdf5p : Enable HDF5 isosurface output.\n"
                    "    --hdf5i_chunk CI CJ CK : Chunk Size (CI,CJ,CK); full output.\n"
                    "      valid values are CI <= NI, CJ <= NJ, CK <= NK\n"
                    "    --hdf5p_chunk CI : Integer percentage of triangles (CI); isosurface output.\n"
                    "      valid values are 2 <= CI <= 100 \n"
                    "    --hdf5_compress : enable compression \n"
    );
#endif

#ifdef HAS_NC
    fprintf(stderr, "    --nci : Enable NC full output.\n"
            "    --ncp : Enable NC isosurface output.\n");
#endif

    /*## End of Output Module Usage Strings ##*/
}

void print_stats(MPI_Comm comm, int rank, int nprocs, uint64_t ntris)
{
    uint64_t *rntris=NULL;   /* All triangle counts */
    int i;
    double Lmax = -1., Lbar = 0.;
    double Ltot = 0., Limb, Lstd = 0.;

    if(rank == 0)
        rntris = (uint64_t *) malloc(nprocs*sizeof(uint64_t));
    MPI_Gather(&ntris, 1, MPI_UINT64_T, rntris, 1, MPI_UINT64_T, 0, comm);

    if(rank == 0) {
        for(i = 0; i < nprocs; i++) {
            if(rntris[i] > Lmax)   Lmax = rntris[i];
            Lbar += rntris[i];
        }
        Ltot = Lbar;
        Lbar /= nprocs;
        Limb = (Lmax / Lbar - 1);
        for(i = 0; i < nprocs; i++)
            Lstd += (rntris[i]-Lbar)*(rntris[i]-Lbar);
        Lstd = sqrt(Lstd/nprocs);
        printf("      Total tris = %.0f, Mean = %.1f, Std = %.1f, Load imbalance = %0.2f\n",
                Ltot, Lbar, Lstd, Limb);
        free(rntris);
    }
}


int main(int argc, char **argv)
{
    int i, j, k, t;    /* Indices */
    int tt;         /* Actual time step from tstart */
    int a;
    float x, y, z;
    float deltax, deltay, deltaz;
    float *data, *xdata;
    int inp = 0;      /* Number of tasks in i */
    int jnp = 0;      /* Number of tasks in j */
    int knp = 0;      /* Number of tasks in k */
    int ni = 0;      /* Global grid size */
    int nj = 0;
    int nk = 0;
    float sigma0x = 0.2f;   /* Gaussian sigma initial */
    float sigma0y = 0.2f;
    float sigma0z = 0.2f;
    float sigma1x = 0.4f;   /* Gaussian sigma final, for gaussresize mode */
    float sigma1y = 0.4f;
    float sigma1z = 0.4f;
    float sigmax = sigma0x;
    float sigmay = sigma0y;
    float sigmaz = sigma0z;
    float x0 = 0.5f;       /* Gaussian center */
    float y0 = 0.5f;
    float z0 = 0.5f;
    int centertask = 1;    /* Set Gaussian center at first task */
    float A = -1.f;   /* Ampl. of sinusoid from top relative when Gaussian present, -1 invalid */
    float fx = 15.f;      /* Freq. of sinusoid */
    float fy = 15.f;
    float fz = 15.f;
    float omegax;      /* Angular freq. */
    float omegay;
    float omegaz;
    double noisespacefreqx = 10.0;   /* Spatial frequency of noise */
    double noisespacefreqy = 10.0;   /* Spatial frequency of noise */
    double noisespacefreqz = 10.0;   /* Spatial frequency of noise */
    double noisetimefreq = 0.25;    /* Temporal frequency of noise */
    float noisemin = -1.f, noisemax = 1.f;   /* Range of noise */
    float noisescale, noiseshift;   /* Range of noise as scale + shift */
    int fulloutputnoiseonly = 0;    /* Boolean: output only noise in full output */
    int tstart = 0;
    int nt = 11;  /* Number of time steps */
    typedef enum { sin2gauss, gaussmove, gaussresize } modetype;
    modetype mode = sin2gauss;       /* Time animation mode */
    int gaussmovebackward = 0;       /* Whether gaussmove goes backward */
    struct sfc3_ctx sfc;     /* Space filling curve for gaussmove */
    int sfc0i=0, sfc0j=0, sfc0k=0;    /* Space filling curve 1st of 2 points */
    struct isoinfo iso;       /* Isosurface context */
    struct osn_context *osn;    /* Open simplex noise context */
    float isothresh = -1.f;    /* Threshold of isosurface, -1 invalid */
    double computetime, fullouttime, isotime, isoouttime;   /* Timers */

    /* MPI vars */
    int rank, nprocs;
    int cprocs[3], cpers[3], crnk[3];  /* MPI Cartesian info */
    MPI_Comm comm;    /* MPI Cartesian communicator */
    int cni, cnj, cnk;   /* Points in this task */
    int is, js, ks;     /* Global index starting points */
    float xs, ys, zs;    /* Global coordinate starting points */

    /*## Add Output Modules' Variables Here ##*/

#ifdef HAS_PVTI
    int pvtiout = 0;
#endif

#ifdef HAS_PVTP
    int pvtpout = 0;
#endif

#ifdef HAS_ADIOS
    char *adiosfullmethod = NULL;
    struct adiosfullinfo adiosfull_nfo;
    char *adiosisomethod = NULL;
    struct adiosisoinfo adiosiso_nfo;
    char *adios_transform=NULL;    /* e.g., compression */
    char *adiosopts = NULL;
#endif

#ifdef HAS_NC
    int nciout = 0;
    int ncpout = 0;
#endif

#ifdef HAS_HDF5
    int hdf5iout = 0;
    int hdf5pout = 0;
    hsize_t *hdf5i_chunk=NULL;
    hsize_t *hdf5p_chunk=NULL;
    int hdf5_compress = 0;
#endif

    /*## End of Output Module Variables ##*/

    /* Init MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

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
        } else if(!strcasecmp(argv[a], "--sigma")) {
            sigmax = sigma0x = strtof(argv[++a], NULL);
            sigmay = sigma0y = strtof(argv[++a], NULL);
            sigmaz = sigma0z = strtof(argv[++a], NULL);
        } else if(!strcasecmp(argv[a], "--sigmaend")) {
            sigma1x = strtof(argv[++a], NULL);
            sigma1y = strtof(argv[++a], NULL);
            sigma1z = strtof(argv[++a], NULL);
        } else if(!strcasecmp(argv[a], "--center")) {
            centertask = 0;
            x0 = strtof(argv[++a], NULL);
            y0 = strtof(argv[++a], NULL);
            z0 = strtof(argv[++a], NULL);
        } else if(!strcasecmp(argv[a], "--centertask")) {
            centertask = 1;
        } else if(!strcasecmp(argv[a], "--amp")) {
            A = strtof(argv[++a], NULL);
        } else if(!strcasecmp(argv[a], "--freq")) {
            fx = strtof(argv[++a], NULL);
            fy = strtof(argv[++a], NULL);
            fz = strtof(argv[++a], NULL);
        } else if(!strcasecmp(argv[a], "--noisespacefreq")) {
            noisespacefreqx = strtod(argv[++a], NULL);
            noisespacefreqy = strtod(argv[++a], NULL);
            noisespacefreqz = strtod(argv[++a], NULL);
        } else if(!strcasecmp(argv[a], "--noisetimefreq")) {
            noisetimefreq = strtod(argv[++a], NULL);
        } else if(!strcasecmp(argv[a], "--noiserange")) {
            noisemin = strtof(argv[++a], NULL);
            noisemax = strtof(argv[++a], NULL);
        } else if(!strcasecmp(argv[a], "--fulloutputnoiseonly")) {
            fulloutputnoiseonly = 1;
        } else if(!strcasecmp(argv[a], "--tsteps")) {
            nt = atoi(argv[++a]);
        } else if(!strcasecmp(argv[a], "--tstart")) {
            tstart = atoi(argv[++a]);
        } else if(!strcasecmp(argv[a], "--sin2gauss")) {
            mode = sin2gauss;
        } else if(!strcasecmp(argv[a], "--gaussmove")) {
            mode = gaussmove;
        } else if(!strcasecmp(argv[a], "--gaussresize")) {
            mode = gaussresize;
        } else if(!strcasecmp(argv[a], "--backward")) {
            gaussmovebackward = 1;
        }

        /*## Add Output Modules' Command Line Arguments Here ##*/

#ifdef HAS_PVTI
        else if(!strcasecmp(argv[a], "--pvti")) {
            pvtiout = 1;
        }
#endif

#ifdef HAS_PVTP
        else if(!strcasecmp(argv[a], "--pvtp")) {
            pvtpout = 1;
        }
#endif

#ifdef HAS_ADIOS
        else if(!strcasecmp(argv[a], "--adiosfull")) {
            adiosfullmethod = argv[++a];
        }
        else if(!strcasecmp(argv[a], "--adiosiso")) {
            adiosisomethod = argv[++a];
        }
        else if(!strcasecmp(argv[a], "--adios_transform")) {
            adios_transform = argv[++a];
        }
        else if(!strcasecmp(argv[a], "--adiosopts")) {
            adiosopts = argv[++a];
        }
#endif

#ifdef HAS_NC
        else if(!strcasecmp(argv[a], "--nci")) {
          nciout = 1;
        } else if(!strcasecmp(argv[a], "--ncp")) {
          ncpout = 1;
        }
#endif

#ifdef HAS_HDF5
        else if(!strcasecmp(argv[a], "--hdf5i")) {
	    hdf5iout = 1;
        }
        else if(!strcasecmp(argv[a], "--hdf5p")) {
	    hdf5pout = 1;
        }
        else if(!strcasecmp(argv[a], "--hdf5i_chunk")) {
            hdf5i_chunk = malloc(3 * sizeof(hsize_t));
            hdf5i_chunk[2] = (hsize_t)strtoul(argv[++a], NULL, 0);
            hdf5i_chunk[1] = (hsize_t)strtoul(argv[++a], NULL, 0);
            hdf5i_chunk[0] = (hsize_t)strtoul(argv[++a], NULL, 0);
        }
        else if(!strcasecmp(argv[a], "--hdf5p_chunk")) {
            hdf5p_chunk = malloc(1 * sizeof(hsize_t));
            hdf5p_chunk[0] = (hsize_t)strtoul(argv[++a], NULL, 0);
        }
        else if(!strcasecmp(argv[a], "--hdf5_compress")) {
            hdf5_compress=1;
        }
#endif

        /*## End of Output Module Command Line Arguments ##*/

        else {
            if(rank == 0)   fprintf(stderr, "Option not recognized: %s\n\n", argv[a]);
            print_usage(rank, NULL);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    /* Check arguments & proc counts */
    if(inp < 1 || jnp < 1 || knp < 1) {
        print_usage(rank, "Error: tasks not specified or incorrect");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if(ni < 1 || nj < 1 || nk < 1) {
        print_usage(rank, "Error: size not specified or incorrect");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if(nt < 1) {
        print_usage(rank, "Error: tsteps incorrect, must be > 0");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if(inp*jnp*knp != nprocs) {
        print_usage(rank, "Error: product of tasks does not equal total MPI tasks");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if(ni % inp || nj % jnp || nk % knp) {
        print_usage(rank, "Error: number of points on an axis is not evenly divisible "
                "by axis tasks.\n   This is required for proper load balancing.");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Set up Cartesian communicator */
    cprocs[0] = inp;  cprocs[1] = jnp;  cprocs[2] = knp;
    cpers[0] = 0;  cpers[1] = 0;  cpers[2] = 0;    /* No periodicity */
    MPI_Cart_create(MPI_COMM_WORLD, 3, cprocs, cpers, 1, &comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Cart_coords(comm, rank, 3, crnk);

    /* Assign default arguments for A & isothresh */
    if(isothresh == -1.f) {
        /*isothresh = exp(-0.5);   * Default is at the gaussian of 1*sigma */
        isothresh = 0.68;   /* Needs to be > 2/3 to extract the "top" half of sin topology */
    }
    if(A == -1.f)
        A = (1 - exp(-0.5));   /* Default is at gaussian of 1*sigma, but from bottom */
    if(rank==0)  printf("isothresh = %f, A = %f\n", isothresh, A);

    /* Data inits */
    omegax = fx * 2 * M_PI;
    omegay = fy * 2 * M_PI;
    omegaz = fz * 2 * M_PI;
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
    /* add a ghost point to the far side of each axis for pvti & iso continuity */
    /* note: this seems to confuse readability of formats that don't expect continuity,
     *       like ADIOS and ?, but it does not hurt for mere benchmarking of those */
    if(1 && crnk[0] < inp-1)   cni++;
    if(1 && crnk[1] < jnp-1)   cnj++;
    if(1 && crnk[2] < knp-1)   cnk++;

    /* Set up space filling curve for gaussmove
          The curve moves through each parallel task */
    if(mode == gaussmove || centertask) {
        sfc3_serpentine_init(&sfc, inp, jnp, knp, gaussmovebackward);
        sfc0i = sfc.i;  sfc0j = sfc.j;  sfc0k = sfc.k;
        sfc3_serpentine_next(&sfc);
    }
    /* Set up center coords if centertask is on */
    if(centertask) {
        x0 = (sfc0i + 1.f) / inp - .5f / inp;
        y0 = (sfc0j + 1.f) / jnp - .5f / jnp;
        z0 = (sfc0k + 1.f) / knp - .5f / knp;
        if(rank==0)  printf("x0=%f, y0=%f, z0=%f\n", x0, y0, z0);
    }
    /* Set up isosurfacing structure */
    isoinit(&iso, xs, ys, zs, deltax, deltay, deltaz, cni, cnj, cnk, 1);
    /* Set up osn */
    open_simplex_noise(12345, &osn);   /* Fixed seed, for now */

    /* Allocate arrays */
    data = (float *) malloc((size_t)cni*cnj*cnk*sizeof(float));
    xdata = (float *) malloc((size_t)cni*cnj*cnk*sizeof(float));

    /*## Add Output Modules' Initialization Here ##*/

#ifdef HAS_ADIOS
    if(adiosfullmethod) {
        adiosfull_init(&adiosfull_nfo, adiosfullmethod, adios_transform, "cartiso.full", comm,
                       rank, nprocs, nt, ni, nj, nk, is, cni, js, cnj, ks, cnk, deltax, 
                       deltay, deltaz, adiosopts);
        if(!fulloutputnoiseonly)  adiosfull_addvar(&adiosfull_nfo, "value", data);
        adiosfull_addvar(&adiosfull_nfo, "noise", xdata);
    }
    if(adiosisomethod) {
        adiosiso_init(&adiosiso_nfo, adiosisomethod, adios_transform, "cartiso.iso", comm, 
                       rank, nprocs, nt, ni, nj, nk, cni, cnj, cnk, adiosopts);
        adiosiso_addxvar(&adiosiso_nfo, "noise");
    }
#endif

    /*## End of Output Module Initialization ##*/

    /* Main loops */

    for(t = 0, tt = tstart; t < nt; t++, tt++) {
        float tpar = (float)t;    /* Time anim. parameter [0,1] */
        float alpha = 1.f;        /* Time parameter: changes for sin2gauss */
        float sinshift, sinscale, sigmax2, sigmay2, sigmaz2;
        float tpar_sfc;    /* Time parameter between two sfc points */
        size_t ii;     /* data index */

        if(nt > 1)      /* Rescale tpar [0,1] if num timesteps > 1 */
           tpar = (float)t / (nt-1);

        if(rank == 0) {
            printf("t = %d\n", tt);
            fflush(stdout);
        }

        timer_tick(&computetime, comm, 1);

        switch(mode) {
            case sin2gauss:
                /* Easy, this just sets the alpha morphing parameter */
                alpha = tpar;
                break;
            case gaussresize:
                /* Interpolate sigma with the time parameter */
                sigmax = (1.f-tpar)*sigma0x + tpar*sigma1x;
                sigmay = (1.f-tpar)*sigma0y + tpar*sigma1y;
                sigmaz = (1.f-tpar)*sigma0z + tpar*sigma1z;
                break;
            case gaussmove:
                /* Interpolate between two sfc points with time parameter */
                tpar_sfc = tpar*(nprocs-1) - sfc.ndx + 1;  /* Param between points */
                while(tpar_sfc > 1.f) {        /* while param>1, it's time to advance sfc */
                    sfc0i = sfc.i;  sfc0j = sfc.j;  sfc0k = sfc.k;
                    sfc3_serpentine_next(&sfc);
                    tpar_sfc -= 1.f;
                }
                x0 = (1-tpar_sfc) * ((sfc0i + 1.f) / inp - .5f / inp) +
                       (tpar_sfc) * ((sfc.i + 1.f) / inp - .5f / inp);
                y0 = (1-tpar_sfc) * ((sfc0j + 1.f) / jnp - .5f / jnp) +
                       (tpar_sfc) * ((sfc.j + 1.f) / jnp - .5f / jnp);
                z0 = (1-tpar_sfc) * ((sfc0k + 1.f) / knp - .5f / knp) +
                       (tpar_sfc) * ((sfc.k + 1.f) / knp - .5f / knp);
                break;
        }

        /* Intermediate terms */
        sinshift = (2/(alpha*A-alpha+1)-1)*3;   /* Shift sine up to top 1, with alpha */
        sinscale = (alpha*A-alpha+1)/2/3;   /* Scale sine down to A ampl, with alpha */
        sigmax2 = 2*sigmax*sigmax;    /* 2*sigma^2 */
        sigmay2 = 2*sigmay*sigmay;
        sigmaz2 = 2*sigmaz*sigmaz;

        /* Spatial loops */
        z = zs;
        for(k = 0, ii = 0; k < cnk; k++) {
            y = ys;
            for(j = 0; j < cnj; j++) {
                x = xs;
                for(i = 0; i < cni; i++, ii++) {
                    float sinusoid = ( sin(omegax*x) + sin(omegay*y) + \
                                       cos(omegaz*z) + sinshift ) * sinscale;
                    data[ii] = exp( -alpha*( (x-x0)*(x-x0)/sigmax2 + \
                                             (y-y0)*(y-y0)/sigmay2 +	\
                                             (z-z0)*(z-z0)/sigmaz2 ) ) * sinusoid;
                    xdata[ii] = (float)open_simplex_noise4(osn, x * noisespacefreqx,
                                  y * noisespacefreqy, z * noisespacefreqz, tt*noisetimefreq);
                    x += deltax;
                }
                y += deltay;
            }
            z += deltaz;
        }

        /* Fully normalize noise noisemin to noisemax */
        if(1) {       /* for now, pretend that this will be optional */
            size_t all = (size_t)cni * cnj * cnk;
            float min = xdata[0], max = xdata[0];
            float rmin, rmax, rrange;
            for(ii = 1; ii < all; ii++) {
                if(xdata[ii] < min)   min = xdata[ii];
                if(xdata[ii] > max)   max = xdata[ii];
            }
            MPI_Allreduce(&min, &rmin, 1, MPI_FLOAT, MPI_MIN, comm);
            MPI_Allreduce(&max, &rmax, 1, MPI_FLOAT, MPI_MAX, comm);
            /* Noise scale & shift */
            rrange = rmax - rmin;
            noisescale = noisemax/rrange - noisemin/rrange;
            noiseshift = noisemin*rmin/rrange - noisemax*rmin/rrange + noisemin;
            for(ii = 0; ii < all; ii++)
                xdata[ii] = xdata[ii] * noisescale + noiseshift;
        }

        timer_tock(&computetime);

        if(rank == 0) {
            printf("   Full Output...\n");   fflush(stdout);
        }
        timer_tick(&fullouttime, comm, 1);

        /*## Add FULL OUTPUT Modules' Function Calls Per Timestep Here ##*/

#ifdef HAS_PVTI
        if(pvtiout) {
            if(rank == 0) {
                printf("      Writing pvti...\n");   fflush(stdout);
            }
            if(!fulloutputnoiseonly) {
                writepvti("cartiso", "value", comm, rank, nprocs, tt, ni, nj, nk,
                      is, is+cni-1, js, js+cnj-1, ks, ks+cnk-1,
                      deltax, deltay, deltaz, data);
            }
            writepvti("cartiso", "noise", comm, rank, nprocs, tt, ni, nj, nk,
                      is, is+cni-1, js, js+cnj-1, ks, ks+cnk-1,
                      deltax, deltay, deltaz, xdata);
        }
#endif

#ifdef HAS_ADIOS
        if(adiosfullmethod) {
            if(rank == 0) {
                printf("      Writing adios full...\n");   fflush(stdout);
            }
            adiosfull_write(&adiosfull_nfo, tt);
	}
#endif

#ifdef HAS_NC
        if(nciout) {
          if(rank == 0) {
            printf("      Writing nci...\n"); fflush(stdout);
          }

          if(!fulloutputnoiseonly) {
              writenci("cartiso", "value", comm, rank, nprocs, tt, ni, nj, nk,
                   is, is+cni-1, js, js+cnj-1, ks, ks+cnk-1,
                   deltax, deltay, deltaz, cni, cnj, cnk, data);
          }
          writenci("cartiso", "noise", comm, rank, nprocs, tt, ni, nj, nk,
                   is, is+cni-1, js, js+cnj-1, ks, ks+cnk-1,
                   deltax, deltay, deltaz, cni, cnj, cnk, xdata);
        }
#endif

#ifdef HAS_HDF5
        if(hdf5iout) {
            if(rank == 0) {
                printf("      Writing hdf5i...\n");   fflush(stdout);
            }
            if(!fulloutputnoiseonly) {
                writehdf5i("cartiso", "value", comm, rank, nprocs, tt, ni, nj, nk,
                    is, is+cni-1, js, js+cnj-1, ks, ks+cnk-1,
		            deltax, deltay, deltaz, cni, cnj, cnk, data, hdf5i_chunk, hdf5_compress);
            }
            writehdf5i("cartiso", "noise", comm, rank, nprocs, tt, ni, nj, nk,
                   is, is+cni-1, js, js+cnj-1, ks, ks+cnk-1,
		           deltax, deltay, deltaz, cni, cnj, cnk, xdata, hdf5i_chunk, hdf5_compress);
        }
#endif

        /*## End of FULL OUTPUT Module Function Calls Per Timestep ##*/

        timer_tock(&fullouttime);

        /* Generate isosurface */
        if(rank == 0) {
            printf("   Isosurface...\n");   fflush(stdout);
        }
        timer_tick(&isotime, comm, 1);
        isosurf(&iso, isothresh, data, xdata);
        timer_tock(&isotime);
        /*printf("      %d tris = %llu\n", rank, iso.ntris);*/
        print_stats(comm, rank, nprocs, iso.ntris);

        if(rank == 0) {
            printf("   Isosurface Output...\n");   fflush(stdout);
        }
        timer_tick(&isoouttime, comm, 1);

        /*## Add ISOSURFACE OUTPUT Modules' Function Calls Per Timestep Here ##*/

#ifdef HAS_PVTP
        if(pvtpout) {
            if(rank == 0) {
                printf("      Writing pvtp...\n");   fflush(stdout);
            }
            writepvtp("cartiso", "iso", comm, rank, nprocs, tt, iso.ntris,
                      iso.points, iso.norms, iso.xvals, "noise");
        }
#endif

#ifdef HAS_ADIOS
        if(adiosisomethod) {
            if(rank == 0) {
                printf("      Writing adios iso...\n");   fflush(stdout);
            }
            adiosiso_write(&adiosiso_nfo, tt, iso.ntris, iso.points, iso.norms,
                           &iso.xvals);
	}
#endif

#ifdef HAS_NC
        if(ncpout) {
          if(rank == 0) {
            printf("     Writing ncp...\n"); fflush(stdout);
          }
          writencp("cartiso", "iso", comm, rank, nprocs, tt, iso.ntris,
                     iso.points, iso.norms, iso.xvals, "noise");
        }
#endif

#ifdef HAS_HDF5
        if(hdf5pout) {
            if(rank == 0) {
                printf("      Writing hdf5p...\n");   fflush(stdout);
            }
            writehdf5p("cartiso", "iso", comm, rank, nprocs, tt, iso.ntris,
		       iso.points, iso.norms, iso.xvals, "noise", hdf5p_chunk, hdf5_compress);
        }
#endif

        /*## End of ISOSURFACE OUTPUT Module Function Calls Per Timestep ##*/

        timer_tock(&isoouttime);
        timer_collectprintstats(computetime, comm, 0, "   Compute");
        timer_collectprintstats(fullouttime, comm, 0, "   FullOutput");
        timer_collectprintstats(isotime, comm, 0, "   Isosurface");
        timer_collectprintstats(isoouttime, comm, 0, "   IsoOutput");
    }

    /*## Add Output Modules' Cleanup Here ##*/

#ifdef HAS_ADIOS
    if(adiosfullmethod)  adiosfull_finalize(&adiosfull_nfo);
    if(adiosisomethod) adiosiso_finalize(&adiosiso_nfo);
#endif

    /*## End of Output Module Cleanup ##*/

    open_simplex_noise_free(osn);
    isofree(&iso);
    free(data);
    free(xdata);

#ifdef HAS_HDF5
    free(hdf5i_chunk);
    free(hdf5p_chunk);
#endif

    MPI_Finalize();

    return 0;
}
