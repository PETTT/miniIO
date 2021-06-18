/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#include <hdf5.h>

void writehdf5(char *name, MPI_Comm comm, int tstep, uint64_t npoints, uint64_t nptstask,
               float *xpts, float *ypts, float *zpts, uint64_t nelems3, uint64_t *conns3,
               uint64_t nelems2, uint64_t *conns2, char *varname, float *data, hsize_t *h5_chunk, char *hdf5_compress, unsigned int *compress_par);
