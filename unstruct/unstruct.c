#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <strings.h>
#include <assert.h>
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
        if(*n1 > *n2)  break;    /* Stop when 2nd factor is greater */
        p1 = *n1;  p2 = *n2;   /* Otherwise, set previous factors */
    }
    if(*n1 - *n2 > p2 - p1) {   /* Whichever set of factors is closest together wins */
        *n1 = p1;  *n2 = p2;
    }

    free(v);

    assert(orign == *n1 * *n2);
}


int main(int argc, char **argv)
{
    int i, a;
    uint64_t npoints = 0;    /* Total number of points */
    uint64_t nptstask = 0;     /* Number of point per task */
    int ntheta, nphi, nlayer;    /* Points per task per spherical coord axis */
    int thetaprocs, phiprocs;             /* Number of tasks on theta & phi respectively */

    /* MPI vars */
    int rank, nprocs;

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
    prime_split(nprocs, &thetaprocs, &phiprocs);   /* Split tasks across spherical axes */
    if(npoints)   /* split npoints among tasks */
        nptstask = npoints / nprocs;
    ntheta = (int) ceil(cbrt(nptstask * thetaprocs * 2 / phiprocs));
                  /* need thetaprocs/phiprocs proportional; want 2x ntheta's as layers */
    nphi = ntheta * phiprocs / thetaprocs;
    nlayer = ntheta / 2;
    nptstask = (uint64_t)ntheta * nphi * nlayer;   /* nptstask won't be exactly as requested */
    npoints = nptstask * nprocs;
    if(rank == 0)  
        printf("Actual points: %llu, points per task: %llu\n"
                "thetaprocs: %d, phiprocs: %d\n"
                "ntheta: %d, nphi: %d, nlayer: %d\n", 
                npoints, nptstask, thetaprocs, phiprocs, ntheta, nphi, nlayer); 

    /* Generate grid */
    

    /*## Add Output Modules' Initialization Here ##*/

    /*## End of Output Module Initialization ##*/

    /* Main loops */

        /*## Add Output Modules' Function Calls Per Timestep Here ##*/

        /*## End of Output Module Functions Calls Per Timestep ##*/

    /*## Add Output Modules' Cleanup Here ##*/

    /*## End of Output Module Cleanup ##*/

    /* Cleanup */

    MPI_Finalize();

    return 0;
}

