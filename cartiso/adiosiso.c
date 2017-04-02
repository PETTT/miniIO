/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#include <stdio.h>
#include <stdlib.h>

#include "adiosiso.h"

static const int fnstrmax = 4095;


void adiosiso_init(struct adiosisoinfo *nfo, char *method, char *name,
        MPI_Comm comm, int rank, int nprocs, int tsteps, int ni, int nj, int nk,
        int cni, int cnj, int cnk, char *adiosopts)
{
    static char emptystr[] = "";

    /* Set up struct, haven't decided if using all of these yet */
    nfo->name = name;
    nfo->comm = comm;
    nfo->rank = rank;
    nfo->nprocs = nprocs;
    nfo->tsteps = tsteps;
    nfo->ni = ni;
    nfo->nj = nj;
    nfo->nk = nk;
    nfo->cni = cni;
    nfo->cnj = cnj;
    nfo->cnk = cnk;
    nfo->numxvars = 0;
    nfo->maxxvars = 1000;
    nfo->xvarnames = (char **) malloc(nfo->maxxvars * sizeof(char *));
    nfo->bufallocsize = 0;

    if(!adiosopts)  adiosopts = emptystr;

    /* Set up ADIOS */
    adios_init_noxml(comm);    /* Not sure if this would conflict with adiosfull */
    adios_declare_group(&nfo->gid, nfo->name, "", adios_flag_no);
    adios_select_method(nfo->gid, method, adiosopts, "");

    /* Define output variables */
    adios_define_var(nfo->gid, "rank", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "tstep", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "npoints", "", adios_unsigned_long, "", "", "");
    adios_define_var(nfo->gid, "cnpoints", "", adios_unsigned_long, "", "", "");
    adios_define_var(nfo->gid, "cstart", "", adios_unsigned_long, "", "", "");
    adios_define_var(nfo->gid, "npoints3", "", adios_unsigned_long, "", "", "");
    adios_define_var(nfo->gid, "cnpoints3", "", adios_unsigned_long, "", "", "");
    adios_define_var(nfo->gid, "cstart3", "", adios_unsigned_long, "", "", "");
    adios_define_var(nfo->gid, "coords", "", adios_real, 
                     "cnpoints3", "npoints3", "cstart3");
    adios_define_var(nfo->gid, "norms", "", adios_real, 
                     "cnpoints3", "npoints3", "cstart3");
    adios_define_var(nfo->gid, "connections", "", adios_unsigned_long, 
                     "cnpoints", "npoints", "cstart");
}    

void adiosiso_addxvar(struct adiosisoinfo *nfo, char *varname)
{
    if(nfo->numxvars > nfo->maxxvars)
        return;   /* Just ignore too many variables, for now */
    nfo->xvarnames[nfo->numxvars] = varname;
    nfo->numxvars++;
    adios_define_var(nfo->gid, varname, "", adios_real, "cnpoints", "npoints", "cstart");
}

void adiosiso_write(struct adiosisoinfo *nfo, int tstep, uint64_t ntris, float *points,
        float *norms, float **xvals)
{
    char fname[fnstrmax+1];
    int timedigits = 4;
    uint64_t groupsize, totalsize;
    int64_t handle;
    int ret;
    int bufneeded;
    uint64_t i, npoints, cnpoints, cstart, npoints3, cnpoints3, cstart3, *npts_all;
    uint64_t *conns;

    /* Set local sizes */
    cnpoints = ntris * 3;
    cnpoints3 = cnpoints * 3;
    groupsize = sizeof(int) /*rank*/ + sizeof(int) /*tstep*/ +
                sizeof(uint64_t)*3 /*npoints-cstart*/ +
                sizeof(uint64_t)*3 /*npoints3-cstart3*/ +
                sizeof(float)*cnpoints3 /*coords*/ +
                sizeof(float)*cnpoints3 /*norms*/ +
                sizeof(uint64_t)*cnpoints /*connections*/ +
                sizeof(float)*cnpoints*nfo->numxvars; /*xvars*/ 

    /* Allocate buffer large enough for all data to write, if not done already */
    bufneeded = (int)(groupsize/(1024*1024));
    bufneeded += bufneeded/10 + 5;   /* Add an extra 10% & 5MB to be sure */
    if(nfo->bufallocsize < bufneeded) {
#       if ADIOS_VERSION_GE(1,10,0)
        adios_set_max_buffer_size(bufneeded);
#       else
        adios_allocate_buffer(ADIOS_BUFFER_ALLOC_NOW, bufneeded);
#       endif
        nfo->bufallocsize = bufneeded;
    }

    /* Determine global sizes */
    npts_all = (uint64_t *) malloc(nfo->nprocs * sizeof(uint64_t));
    MPI_Allgather(&cnpoints, 1, MPI_UNSIGNED_LONG_LONG,
                  npts_all, 1, MPI_UNSIGNED_LONG_LONG, nfo->comm);
    for(npoints = 0, i = 0; i < nfo->nprocs; ++i) 
        npoints += npts_all[i];
    for(cstart = 0, i = 0; i < nfo->rank; ++i)
        cstart += npts_all[i];
    npoints3 = npoints * 3;
    cstart3 = cstart *3;
    free(npts_all);

    /* Create connections */
    conns = (uint64_t *) malloc(cnpoints * sizeof(uint64_t));
    for(i = 0; i < cnpoints; i++)
        conns[i] = i;

    /* Set filename */
    snprintf(fname, fnstrmax, "%s.%0*d.bp", nfo->name, timedigits, tstep);

    /* Open & Write */
    ret = adios_open(&handle, nfo->name, fname, "w", nfo->comm);
    if(ret) {
        fprintf(stderr, "Error opening ADIOS file: %s\n", fname);
        return;
    }
#   if ADIOS_VERSION_LE(1,9,0)
    adios_group_size(handle, groupsize, &totalsize);
#   endif

    adios_write(handle, "rank", &nfo->rank);
    adios_write(handle, "tstep", &tstep);
    adios_write(handle, "npoints", &npoints);
    adios_write(handle, "cnpoints", &cnpoints);
    adios_write(handle, "cstart", &cstart);
    adios_write(handle, "npoints3", &npoints3);
    adios_write(handle, "cnpoints3", &cnpoints3);
    adios_write(handle, "cstart3", &cstart3);
    adios_write(handle, "coords", points);
    adios_write(handle, "norms", norms);
    adios_write(handle, "connections", conns);
    for(i = 0; i < nfo->numxvars; i++)
        adios_write(handle, nfo->xvarnames[i], xvals[i]);

    adios_close(handle);
    free(conns);
}

void adiosiso_finalize(struct adiosisoinfo *nfo)
{
    free(nfo->xvarnames);
    adios_finalize(nfo->rank);
}

