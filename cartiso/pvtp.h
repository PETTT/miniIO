/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

void writepvtp(char *name, char *varname, MPI_Comm comm, int rank, int nprocs, 
               int tstep, uint64_t ntris, float *points, float *norms, 
               float *xvals, char *xname);
