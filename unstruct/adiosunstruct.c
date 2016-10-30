/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#include <stdio.h>
#include <stdlib.h>

#include "adiosunstruct.h"

static const int fnstrmax = 4095;

void adiosunstruct_init(struct adiosinfo *nfo, char *method, char *name, MPI_Comm comm, 
        int rank, int nprocs, int tsteps, uint64_t nptstask, 
        uint64_t nelems3, uint64_t nelems2)
{
    /* Set up struct */
    nfo->name = name;
    nfo->comm = comm;
    nfo->rank = rank;
    nfo->nprocs = nprocs;
    nfo->tsteps = tsteps;
    nfo->cnpoints = nptstask;
    nfo->cnelems3 = nelems3;
    nfo->cnelems2 = nelems2;
    nfo->nvars = 0;
    nfo->maxvars = 1000;
    nfo->varnames = (char **) malloc(nfo->maxvars * sizeof(char *));
    nfo->bufallocsize = 0;
    
    /* Set up ADIOS */ 
    adios_init_noxml(comm);
    adios_declare_group(&nfo->gid, nfo->name, "", adios_flag_no);
    /*chkdir1task(dirname, comm);*/
    adios_select_method(nfo->gid, method, "", /*dirname*/"");

    /* Define output variables */
    adios_define_var(nfo->gid, "rank", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "tstep", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "npoints", "", adios_unsigned_long, "", "", "");
    adios_define_var(nfo->gid, "cnpoints", "", adios_unsigned_long, "", "", "");
    adios_define_var(nfo->gid, "cspoints", "", adios_unsigned_long, "", "", "");
    adios_define_var(nfo->gid, "nconns3", "", adios_unsigned_long, "", "", "");    
    adios_define_var(nfo->gid, "cnconns3", "", adios_unsigned_long, "", "", "");    
    adios_define_var(nfo->gid, "csconns3", "", adios_unsigned_long, "", "", "");    
    adios_define_var(nfo->gid, "nconns2", "", adios_unsigned_long, "", "", "");
    adios_define_var(nfo->gid, "cnconns2", "", adios_unsigned_long, "", "", "");
    adios_define_var(nfo->gid, "csconns2", "", adios_unsigned_long, "", "", "");    
    adios_define_var(nfo->gid, "xpts", "", adios_real,
                     "cnpoints", "npoints", "cspoints");
    adios_define_var(nfo->gid, "ypts", "", adios_real,
                     "cnpoints", "npoints", "cspoints");
    adios_define_var(nfo->gid, "zpts", "", adios_real,
                     "cnpoints", "npoints", "cspoints");
    adios_define_var(nfo->gid, "conns3", "", adios_unsigned_long,
                     "cnconns3", "nconns3", "csconns3");
    adios_define_var(nfo->gid, "conns2", "", adios_unsigned_long,
                     "cnconns2", "nconns2", "csconns2");
}

void adiosunstruct_addvar(struct adiosinfo *nfo, char *varname)
{
    if(nfo->nvars > nfo->maxvars)
        return;   /* Just ignore too many variables, for now */
    nfo->varnames[nfo->nvars] = varname;
    nfo->nvars++;
    adios_define_var(nfo->gid, varname, "", adios_real, 
                     "cnpoints", "npoints", "cspoints");
}

void adiosunstruct_write(struct adiosinfo *nfo, int tstep, float *xpts, float *ypts,
                         float *zpts, uint64_t *conns3, uint64_t *conns2, float **vars)
{
    char fname[fnstrmax+1];
    int timedigits = 4;
    uint64_t groupsize, totalsize;
    int64_t handle;
    int ret;
    int bufneeded;
    uint64_t i;
    uint64_t npoints, cspoints; 
    uint64_t cnconns3, nconns3, csconns3;
    uint64_t cnconns2, nconns2, csconns2;

    /* Set global sizes: we assume all tasks have the same size unstructured data!!! */
    npoints = nfo->cnpoints * nfo->nprocs;    /* global # of points */
    cspoints = nfo->cnpoints * nfo->rank;     /* local starting points */
    cnconns3 = nfo->cnelems3 * 6;           /* 6 connections per 3D prism element */
    nconns3 = cnconns3 * nfo->nprocs;    /* global # of 3D connections */
    csconns3 = cnconns3 * nfo->rank;     /* local starting 3D connections */
    cnconns2 = nfo->cnelems2 * 3;           /* 3 connections per 2D triangle element */
    nconns2 = cnconns2 * nfo->nprocs;    /* global # of 2D connections */
    csconns2 = cnconns2 * nfo->rank;     /* local starting 2D connections */
    
    /* ADIOS group size */
    groupsize = sizeof(int) /*rank*/ + sizeof(int) /*tstep*/ +
                sizeof(uint64_t)*9 /*npoints-cselems2*/ +
                sizeof(float)*3*nfo->cnpoints /*xpts-zpts*/ +
                sizeof(uint64_t)*cnconns3 /*conns3*/ +
                sizeof(uint64_t)*cnconns2 /*conns2*/ +
                sizeof(float)*nfo->cnpoints*nfo->nvars; /*vars*/

    /* Allocate buffer large enough for all data to write, if not done already */
    bufneeded = (int)(groupsize/(1024*1024));
    bufneeded += bufneeded/8 + 5;   /* Add an extra 12.5% & 5MB to be sure */
    if(nfo->bufallocsize < bufneeded) {
#       if ADIOS_VERSION_GE(1,10,0)
        adios_set_max_buffer_size(bufneeded);
#       else
        adios_allocate_buffer(ADIOS_BUFFER_ALLOC_NOW, bufneeded);
#       endif
        nfo->bufallocsize = bufneeded;
    }
    
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
    adios_write(handle, "cnpoints", &nfo->cnpoints);
    adios_write(handle, "cspoints", &cspoints);
    adios_write(handle, "nconns3", &nconns3);
    adios_write(handle, "cnconns3", &cnconns3);
    adios_write(handle, "csconns3", &csconns3);
    adios_write(handle, "nconns2", &nconns2);
    adios_write(handle, "cnconns2", &cnconns2);
    adios_write(handle, "csconns2", &csconns2);
    adios_write(handle, "xpts", xpts);
    adios_write(handle, "ypts", ypts);
    adios_write(handle, "zpts", zpts);
    adios_write(handle, "conns3", conns3);
    adios_write(handle, "conns2", conns2);
    for(i = 0; i < nfo->nvars; i++)
        adios_write(handle, nfo->varnames[i], vars[i]);

    adios_close(handle);
}

void adiosunstruct_finalize(struct adiosinfo *nfo)
{
    free(nfo->varnames);
    adios_finalize(nfo->rank);
}

