#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
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
    );

    /*## Add Output Modules' Usage String ##*/

    /*## End of Output Module Usage Strings ##*/
}

int main(int argc, char **argv)
{
    int i;
    uint64_t npoints = 0;
    uint64_t nptstask = 0;

    /* MPI vars */
    int rank, nprocs;

    /*## Add Output Modules' Variables Here ##*/

    /*## End of Output Module Variables ##*/

    /* Init MPI */
    MPI_Init(&argc, &argv);
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
    if(npoints == 0 && nptstask == 0) {
        print_usage(rank, "Error: neither points or pointspertask specified, or there was\n"
                          "       an error parsing them");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Determine a topology that meets # of points requested */
    
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

