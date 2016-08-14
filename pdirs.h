/*
 * Functions to create and check directories from an MPI code
 * on a parallel file system
 *
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#include <sys/stat.h>
#include <errno.h>
#include <unistd.h>
#include <mpi.h>

/* Create a directory from one rank (rank=0)
 *    dirname: name of directory
 *    comm: MPI communicator 
 * note: this can be called from any/all ranks as long as rank 0 calls it
 * note: aborts for any error except when the directory already exists */

static void mkdir1task(char *dirname, MPI_Comm comm)
{
    int rank, ret;

    MPI_Comm_rank(comm, &rank);
    if(rank == 0) {
        ret = mkdir(dirname, 0777);
        if( ret && errno != EEXIST ) {  /* An error for anything but pre-existing */
            fprintf(stderr, "pvtk error: Could not create directory.\n");
            MPI_Abort(comm, 1);
        }
    }
}

/* Make sure a directory exists from one rank (rank=size-1)
 *    dirname: name of directory
 *    comm: MPI communicator
 * note: intended to be used in conjunction with mkdir1task to insure that
 *       a created directory shows up on other ranks
 * note: this can be called from any/all ranks as long as rank size-1 calls it
 * note: tries several times and aborts if it doesn't exist
 */

static void chkdir1task(char *dirname, MPI_Comm comm)
{
    int rank, nprocs, ret;
    struct stat fsbuf;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    /* Have the last rank ensure that the directory is available */
    if(rank == nprocs-1) { 
        do {
            int retry = 0;
            usleep(200000);
            ret = stat(dirname, &fsbuf);
            if(ret && retry > 50) {   /* 50 retries * 200000 us = 10 s timeout */
                fprintf(stderr, "pvtk error: timeout waiting for directory %s\n",
                        dirname);
                MPI_Abort(comm, 1);
            }
        } while(ret);
    }
    MPI_Barrier(comm);   /* Everyone is safe to use directory after barrier */
}

