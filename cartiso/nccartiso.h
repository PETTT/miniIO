#ifdef HAS_NC
#include "netcdf.h"
#endif

void writencp(char *name, char *varname, MPI_Comm comm, int rank, int nprocs,
               int tstep, uint64_t ntris, float *points, float *norms,
              float *xvals, char *xname);

void writenci(char *name, char *varname, MPI_Comm comm, int rank, int nprocs,
              int tstep, int ni, int nj, int nk, int is, int ie, int js, int je,
              int ks, int ke, float deltax, float deltay, float deltaz, int nci, int ncj, int nck,
              float *data);
