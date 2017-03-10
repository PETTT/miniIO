#ifdef HAS_HDF5
#  include "hdf5.h"
#endif

void writehdf5p(char *name, char *varname, MPI_Comm comm, int rank, int nprocs, 
               int tstep, uint64_t ntris, float *points, float *norms, 
		float *xvals, char *xname, hsize_t *h5_chunk, int hdf5_compress);

void writehdf5i(char *name, char *varname, MPI_Comm comm, int rank, int nprocs, 
               int tstep, int ni, int nj, int nk, int is, int ie, int js, int je,
		int ks, int ke, float deltax, float deltay, float deltaz, int nci, int ncj, int nck, float *data, hsize_t *h5_chunk, int hdf5_compress);

