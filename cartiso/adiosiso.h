/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#include <mpi.h>
#include <adios.h>

struct adiosisoinfo {
    char *name;
    MPI_Comm comm;
    int rank;
    int nprocs;
    int tsteps;
    int ni, nj, nk;
    int cni, cnj, cnk;

    int numxvars;
    int maxxvars;
    char **xvarnames;

    int bufallocsize;

    int64_t gid;
};

void adiosiso_init(struct adiosisoinfo *nfo, char *method, char *name,
        MPI_Comm comm, int rank, int nprocs, int tsteps, int ni, int nj, int nk,
        int cni, int cnj, int cnk, char *adiosopts);

void adiosiso_addxvar(struct adiosisoinfo *nfo, char *varname);

void adiosiso_write(struct adiosisoinfo *nfo, int tstep, uint64_t ntris, float *points,
        float *norms, float **xvals);

void adiosiso_finalize(struct adiosisoinfo *nfo);

