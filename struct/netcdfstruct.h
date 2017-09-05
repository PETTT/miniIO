/*
* Copyright (c) DoD HPCMP PETTT.  All rights reserved.
* See LICENSE file for details.
*/

#ifdef HAS_NC
#include "netcdf.h"
#endif

void writenc(const int num_varnames, char **varnames, MPI_Comm comm, int rank, int nprocs, int tstep,
	       int is, int js, int ks,
             int ni, int nj, int nk, int cni, int cnj, int cnk,
             float *data);
