#!/usr/bin/env python

# Process an output log from unstruct and produce a results table
#   + Determines output method from log contents
#   + Can also check filename for method
#     - Filename must contain .przm., .posix., .mpi., etc (see fname2method)
#     - If filename & log contents don't agree, will abort with error
#     - If you didn't name file properly (or at all), just gives a warning
#   + Each column is a time step
#   + Table rows: method, cores, grid points (Million), full file size, time, rate (GB/s)
#   + Multiple methods/types can be pasted into Excel or similar for plotting/analysis

import sys
import re

fname2method = { ".przm.": "MPI-Indiv", ".posix.": "ADIOS-POSIX", 
        ".mpi.": "ADIOS-MPI", ".mpilus.": "ADIOS-MPI_Lustre", ".mpiag.": "ADIOS-MPI_Aggr",
        ".phdf5.": "ADIOS-PHDF5", ".hdf5.": "HDF5", }
lname2method = { "przm": "MPI-Indiv", "POSIX": "ADIOS-POSIX", 
        "MPI ": "ADIOS-MPI", "MPI_LUSTRE": "ADIOS-MPI_Lustre",
        "MPI_AGGREGATE": "ADIOS-MPI_Aggr", "PHDF5": "ADIOS-PHDF5", 
        "hdf5": "HDF5" }

v = []
totalrate = -1.0
a = 1
if sys.argv[a] == "-a":   # produce average rate
    totalrate = 0.0
    a += 1

for fname in sys.argv[a:]:
    if totalrate > 0.0: totalrate = 0.0

    # Determine method from file name, if it's there
    method1key = next((x for x in fname2method.keys() if x in fname), False)
    if not method1key: print "WARNING: method not specified in filename"

    f = open(fname)
    for l in f:
        ls = l.split()
        if(len(ls) < 1): continue
        if(len(ls) < 2): ls.append("")   # Ensure all outer if's can handle ls

        if ls[0][0:4] == "tsk=":
            if ls[0][4] == "\"": continue    # Skip LSF copy of batch script

            cores = int(ls[0][4:])

            # Determine method from line; make sure matches one from file name
            method2key = next((x for x in lname2method.keys() if re.search(x, l)), False)
            if method2key:
                if method1key:
                    if fname2method[method1key] != lname2method[method2key]:
                        print "Error: method in filename doesn't match method in file"
                        sys.exit(1)
                    else:
                        method = fname2method[method1key]
                else:
                    method = lname2method[method2key]
            else:
                print "WARNING: method not specified in log file"
                if method1key:
                    method = fname2method[method1key]
                else:
                    print "ERROR: not method specified, cannot proceed"
                    sys.exit(1)

        elif ls[0] == "Actual" and ls[1] == "points:":
            pointsM = long(ls[2]) / 1000000
        elif ls[0] == "data" and ls[1] == "size":
            sizegb = float(ls[10]) / 1024**3
        elif ls[0] == "Output":
            outtime = float(ls[11][:-1])
            rate = sizegb/outtime
            if totalrate >= 0.0: totalrate += rate
            v.append([method, cores, pointsM, sizegb, outtime, rate])

print "%-16s %5s %8s %7s %6s %6s" % ("method", "cores", "points*M", "f_GB", "ftime", "fGB/s")
for i in v:
    print "%-16s %5d %8s %7.1f %6.2f %6.2f" % tuple(i)
if totalrate >= 0.0: print "Avg rate =", totalrate/len(v)

