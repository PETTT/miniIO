/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include <pdirs.h>

#include "netcdf.h"
#include "netcdf_par.h"
#include "ncunstruct.h"

//#define TIMEIO

static const int fnstrmax = 4095;

#define NCERR {if(err != NC_NOERR) {printf("(rank %d) Error at line %d: %s\n",rank,__LINE__,nc_strerror(err)); fflush(stdout); MPI_Abort(comm,1);}}

uint64_t nelems_in[2];
uint64_t nelems_out[2];

void writenc(char *name, MPI_Comm comm, int tstep, uint64_t npoints,
             uint64_t nptstask, float *xpts, float *ypts, float *zpts,
             uint64_t nelems3, uint64_t *conns3, uint64_t nelems2,
             uint64_t *conns2, char *varname, float *data) {

  char dirname[fnstrmax+1];
  char fname[fnstrmax+1];
  char rel_fname[fnstrmax+1];
  char fname_xdmf[fnstrmax+1];
  int rank, nprocs;
  int timedigits = 4;
  MPI_Info info = MPI_INFO_NULL;

  int conns3id;
  int conns2id;
  int ncid;
  int grid_id;
  int varsid;

  int phony_dim_0_id;
  int phony_dim_1_id;
  int phony_dim_2_id;
  int phony_dim_3_id;

  int xyz_ids[3] = {0,0,0};

  int err = 0;

  size_t start[1] = {0};
  size_t count[1] = {0};

  size_t dims[1] = {0};


  int chunk_pid;
  size_t chunk;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /* Make dir for all output and subdir for timestep */
  snprintf(dirname, fnstrmax, "%s.nc.d", name);
  mkdir1task(dirname, comm);
  snprintf(dirname, fnstrmax, "%s.nc.d/t%0*d.d", name, timedigits, tstep);
  mkdir1task(dirname, comm);

  /* Set up MPI info */
  MPI_Info_create(&info);
  /*     MPI_Info_set(info, "striping_factor", "1");  */

  chkdir1task(dirname, comm);

  snprintf(fname, fnstrmax, "unstruct.nc.d/t%0*d.d/r.nc", timedigits, tstep);
  snprintf(rel_fname, fnstrmax, "t%0*d.d/r.nc", timedigits, tstep);
  snprintf(fname_xdmf, fnstrmax, "unstruct.nc.d/t%0*d.d.xmf", timedigits, tstep);

  nelems_in[0] = nelems3 ;
  nelems_in[1] = nelems2 ;

  MPI_Reduce( nelems_in, nelems_out, 2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD );

#ifdef TIMEIO
  double createfile, prewrite, write, postwrite;   /* Timers */
  timer_tick(&createfile, comm, 0);
#endif

  //if(rank==0) {


    err = nc_create_par(fname,NC_NETCDF4|NC_MPIIO,comm,info,&ncid); NCERR;


    /* Optional grid points */
    if(xpts && ypts && zpts) {
      /* Create the dataspace for the dataset. */
      dims[0] = (size_t)npoints;
      count[0] =(size_t)nptstask;
      //filespace = H5Screate_simple(1, dims, NULL);

      /* Create Grid Group */
      err = nc_def_grp(ncid,"grid points",&grid_id); NCERR;
      err = nc_def_dim(grid_id,"phony_dim_0",dims[0],&phony_dim_0_id); NCERR;
      err = nc_def_var(grid_id,"x",NC_FLOAT,1,&phony_dim_0_id,&xyz_ids[0]); NCERR;
      err = nc_def_var(grid_id,"y",NC_FLOAT,1,&phony_dim_0_id,&xyz_ids[1]); NCERR;
      err = nc_def_var(grid_id,"z",NC_FLOAT,1,&phony_dim_0_id,&xyz_ids[2]); NCERR;

      err = nc_var_par_access(grid_id,xyz_ids[0],NC_COLLECTIVE); NCERR;
      err = nc_var_par_access(grid_id,xyz_ids[1],NC_COLLECTIVE); NCERR;
      err = nc_var_par_access(grid_id,xyz_ids[2],NC_COLLECTIVE); NCERR;

    } /* if x,y,z */


    /* MSB is it possible that some processors have 0? */

    /* Optional grid connections, writes a 64-bit 0 if no connections */
    if(conns3 && nelems3) {

      /* Create the dataspace for the dataset. */
      dims[0] = (size_t)nelems_out[0]*6;
      err = nc_def_dim(ncid,"phony_dim_3",dims[0],&phony_dim_3_id); NCERR;

      err = nc_def_var(ncid,"conns3",NC_UINT64,1,&phony_dim_3_id,&conns3id); NCERR;
      err = nc_var_par_access(ncid,conns3id,NC_COLLECTIVE); NCERR;

    } /* if conns & nelems3 */

    /* Optional 2D surface triangle connections, writes a 64-bit 0 if none */
    if(conns2 && nelems2) {

      /* Create the dataspace for the dataset. */
      dims[0] = (size_t)nelems_out[1]*3;

      err = nc_def_dim(ncid,"phony_dim_1",dims[0],&phony_dim_1_id); NCERR;
      err = nc_def_var(ncid,"conns2",NC_UINT64,1,&phony_dim_1_id,&conns2id); NCERR;
      err = nc_var_par_access(ncid,conns2id,NC_COLLECTIVE); NCERR;

    } /* if conns2 & nelems2 */

    /* Optional variable data, starting with number of variables */
    if(data && varname) {
      /* Create the dataspace for the dataset. */
      dims[0] = (size_t)npoints;
      err = nc_def_dim(ncid,"phony_dim_2",dims[0],&phony_dim_2_id); NCERR;
      err = nc_def_var(ncid,"vars",NC_FLOAT,1,&phony_dim_2_id,&varsid);
      err = nc_var_par_access(ncid,varsid,NC_COLLECTIVE); NCERR;
    } /* if data & varname */

    //} /* Rank == 0 */

    err = nc_enddef(ncid); NCERR;
    MPI_Barrier(MPI_COMM_WORLD);



#ifdef TIMEIO
  timer_tock(&createfile);
  timer_collectprintstats(createfile, comm, 0, "CreateFile");
  timer_tick(&prewrite, comm, 0);
#endif

#ifdef TIMEIO
  timer_tock(&prewrite);
  timer_collectprintstats(prewrite, comm, 0, "PreWrite");

#endif
#ifdef TIMEIO
  timer_tick(&write, comm, 0);
#endif
  /* Optional grid points */
  if(xpts && ypts && zpts) {

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    start[0] =(size_t)(nptstask*rank);
    count[0] =(size_t)nptstask;

    err = nc_put_vara_float(grid_id,xyz_ids[0],start,count,xpts); NCERR;
    err = nc_put_vara_float(grid_id,xyz_ids[1],start,count,ypts); NCERR;
    err = nc_put_vara_float(grid_id,xyz_ids[2],start,count,zpts); NCERR;

  } /* if x,y,z pts */


  /* MSB is it possible that some processors have 0? */

  /* Optional grid connections, writes a 64-bit 0 if no connections */
  if(conns3 && nelems3) {

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    start[0] =(size_t)(nelems3*6*rank);
    count[0] =(size_t)nelems3*6;

    err = nc_put_vara(ncid,conns3id,start,count,conns3); NCERR;

  }

  /*     Optional 2D surface triangle connections, writes a 64-bit 0 if none */
  if(conns2 && nelems2) {
    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    start[0] =(size_t)(nelems2*3*rank);
    count[0] =(size_t)nelems2*3;

    err = nc_put_vara(ncid,conns2id,start,count,conns3); NCERR;
  }

  /* Optional variable data, starting with number of variables */
  if(data && varname) {

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    start[0] =(size_t)(nptstask*rank);
    count[0] =(size_t)nptstask;

    err = nc_put_vara_float(ncid,varsid,start,count,data); NCERR;
  }

#ifdef TIMEIO
  timer_tock(&write);
  timer_collectprintstats(write, comm, 0, "write");
#endif

#ifdef TIMEIO
  timer_tick(&postwrite, comm, 0);
#endif

#ifdef TIMEIO
  timer_tock(&postwrite);
  timer_collectprintstats(postwrite, comm, 0, "PostWrite");
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  err = nc_close(ncid);

}
