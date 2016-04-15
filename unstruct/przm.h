
void writeprzm(char *name, MPI_Comm comm, int tstep, uint64_t npoints,
               float *xpts, float *ypts, float *zpts, uint64_t nelems3, uint64_t *conns3,
               uint64_t nelems2, uint64_t *conns2, char *varname, float *data);

