#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <strings.h>
#include <assert.h>
#include <stdint.h>
#include <mpi.h>

/*## Add Output Modules' Includes Here ##*/

/*## End of Output Module Includes ##*/

void print_usage(int rank, const char *errstr)
{
    if(rank != 0)  return;
    if(errstr)
        fprintf(stderr, "%s\n\n", errstr);
    fprintf(stderr,
"Usage: mpi_launcher [-n|-np NPROCS] ./unstruct "
"\n");

    /*## Add Output Modules' Usage String ##*/

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
float sgn(float x)
{
    return copysignf(1.f, x);
}

/* Superquadric ellipsoidal "c" function */
float sqc(float w, float m)
{
    return sgn(cosf(w)) * powf(fabsf(cosf(w)), m);
}

/* Superquadric ellipsoidal "s" function */
float sqs(float w, float m)
{
    return sgn(sinf(w)) * powf(fabsf(sinf(w)), m);
}


int main(int argc, char **argv)
{
    int a;
    uint64_t npoints = 0;    /* Total number of points */
    uint64_t nptstask = 0;     /* Number of point per task */
    int nu, nv, nlyr;    /* Points per task per u,v,lyr spherical coord axis */
    float du, dv;          /* delta's along u & v points */
    float u0, u1, v0, v1;      /* Starting/ending points along u & v */
    int i, j, k;
    uint64_t ii;
    float *xpts, *ypts, *zpts;    /* Grid points */
    uint64_t nelems2, *conns2;      /* Number of grid triangles & connection array in 2D */
    uint64_t nelems3, *conns3;      /* Number of triangular prisms & connection array */
    float uround = 0.3f;        /* Superquadric roundness u parameter */
    float vround = 0.3f;        /* Superquadric roundness v parameter */

    /* MPI vars */
    int rank, nprocs;
    int uprocs, vprocs;      /* Number of tasks on u & v spherical axes respectively */
    int urank, vrank;        /* 2D rank along u & v */

    /*## Add Output Modules' Variables Here ##*/

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
            uround = atof(argv[++a]);
            vround = atof(argv[++a]);
        }

        /*## Add Output Modules' Command Line Arguments Here ##*/

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
        MPI_Abort(MPI_COMM_WORLD, 1);
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
        printf("Actual points: %llu, points per task: %llu\n"
                "uprocs: %d, vprocs: %d\n"
                "nu: %d, nv: %d, nlyr: %d\n", 
                npoints, nptstask, uprocs, vprocs, nu, nv, nlyr); 

    /* Divide spherical topology tiles along u,v */
    du = 2 * M_PI / uprocs / nu;
    dv = M_PI / vprocs / nv;
    urank = rank % uprocs;
    vrank = rank / uprocs;
    u0 = urank * du * nu - M_PI;   u1 = u0 + du * nu;   /* #points to be continuous across task */
    v0 = vrank * dv * nv - M_PI/2;   v1 = v0 + dv * nv;
    /*DBG*/printf("%d: %f-%f, %f-%f, %f, %f\n", rank, u0, u1, v0, v1, du, dv);
        
    /* Generate grid points with superquadric */
    xpts = (float *) malloc(nptstask*sizeof(float));
    ypts = (float *) malloc(nptstask*sizeof(float));
    zpts = (float *) malloc(nptstask*sizeof(float));
    for(k = 0, ii = 0; k < nlyr; k++) {
        float w = 1.f + powf((float)k/(nlyr-1), 2.f) * 2.f;  /* layer w in [1,3] by squares */
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

    /* Add surface connections, all are triangles */
    nelems2 = (nu-1) * (nv-1) * 2;
    conns2 = (uint64_t *) malloc(nelems2*3*sizeof(uint64_t));
    for(i = 0, ii = 0; i < nu-1; i++) {
        for(j = 0; j < nv-1; j++) {
            uint64_t ndx = i * nv + j;
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
    for(k = 0, ii = 0; k < nlyr-1; k++) {
        uint64_t euv, iuv, nuv = nu * nv * k;
        for(euv = 0, iuv = 0; euv < nelems2; euv++) {
            conns3[ii++] = conns2[iuv++];
            conns3[ii++] = conns2[iuv++];
            conns3[ii++] = conns2[iuv++];
            conns3[ii++] = conns2[iuv++] + nuv;
            conns3[ii++] = conns2[iuv++] + nuv;
            conns3[ii++] = conns2[iuv++] + nuv;
        }
    }

    /* DBG: write csv file in alternating series */
    {
        FILE *f;
        int r;
        for(r = 0; r < nprocs; ++r) {
            if(r == rank) {
                f = fopen("tstunstruct.csv", r==0 ? "wb" : "ab");
                for(k = 0, ii = 0; k < nlyr; k++)
                    for(i = 0; i < nu; i++) 
                        for(j = 0; j < nv; j++, ii++) 
                            fprintf(f, "%f,%f,%f,%d\n", xpts[ii], ypts[ii], zpts[ii], r);
                fclose(f);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    /*## Add Output Modules' Initialization Here ##*/

    /*## End of Output Module Initialization ##*/

    /* Main loops */

        /*## Add Output Modules' Function Calls Per Timestep Here ##*/

        /*## End of Output Module Functions Calls Per Timestep ##*/

    /*## Add Output Modules' Cleanup Here ##*/

    /*## End of Output Module Cleanup ##*/

    /* Cleanup */

    free(xpts);  free(ypts);  free(zpts);
    free(conns2);  free(conns3);
    MPI_Finalize();

    return 0;
}

