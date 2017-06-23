/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

//void timer_tick(double *, MPI_Comm, int);
//void timer_tock(double *);
//void timer_collectstats(double,MPI_Comm,int,struct timer_statinfo*);

void writenc(char *name, MPI_Comm comm, int tstep, uint64_t npoints, uint64_t nptstask,
             float *xpts, float *ypts, float *zpts, uint64_t nelems3, uint64_t *conns3,
             uint64_t nelems2, uint64_t *conns2, char *varname, float *data);
