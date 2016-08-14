/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#include <mpi.h>
#include <adios.h>

struct adiosinfo {
    char *name;
    MPI_Comm comm;
    int rank;
    int nprocs;
    int tsteps;
    uint64_t cnpoints;   /* number of points in local task */
    uint64_t cnelems3;   /* local number of triangular elements in 2D base grid */
    uint64_t cnelems2;   /* local number of prism elements in 3D grid */

    int nvars;
    int maxvars;
    char **varnames;

    int bufallocsize;

    int64_t gid;
};

void adiosunstruct_init(struct adiosinfo *nfo, char *method, char *name, MPI_Comm comm, 
        int rank, int nprocs, int tsteps, uint64_t nptstask, 
        uint64_t nelems3, uint64_t nelems2);

void adiosunstruct_addvar(struct adiosinfo *nfo, char *varname);

void adiosunstruct_write(struct adiosinfo *nfo, int tstep, float *xpts, float *ypts,
                         float *zpts, uint64_t *conns3, uint64_t *conns2, float **vars);

void adiosunstruct_finalize(struct adiosinfo *nfo);

