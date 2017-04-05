/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#include <mpi.h>
#include <adios.h>

struct adiosfullinfo {
    char *name;
    MPI_Comm comm;
    int rank;
    int nprocs;
    int tsteps;
    int ni, nj, nk;
    int is, cni, js, cnj, ks, cnk;
    float deltax, deltay, deltaz;

    int nvars;
    int maxvars;
    char **varnames;
    float **datas;

    int bufallocsize;

    int64_t gid;
};

void adiosfull_init(struct adiosfullinfo *nfo, char *method,
               char *name, MPI_Comm comm, int rank, int nprocs,
               int tsteps, int ni, int nj, int nk, int is, int cni, int js, int cnj,
               int ks, int cnk, float deltax, float deltay, float deltaz, char *adiosopts);

void adiosfull_addvar(struct adiosfullinfo *nfo, char *varname, float *data);

void adiosfull_write(struct adiosfullinfo *nfo, int tstep);

void adiosfull_finalize(struct adiosfullinfo *nfo);

