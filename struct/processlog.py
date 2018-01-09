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

totalrate = -1.0
a = 1
printFName = False

for arg in sys.argv[1:2]:
    if arg == "-a":   # produce average rate
        totalrate = 0.0
        a += 1

if len(sys.argv[a:]) > 1: printFName = True

for fname in sys.argv[a:]:

    v = []
    if totalrate > 0.0: totalrate = 0.0     # reset for new file
    bal = False

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
            
            print ls
            cores = int(ls[0][4:]) *  int(ls[1])
            #size = ls[2][3:] + 'x' +ls[3][3:] + 'x' + ls[4][3:]
            #sizeprod = float(ls[2][3:]) * float(ls[3][3:]) * float(ls[4][3:])
            #sizegb = sizeprod*4/1024**3
            
            # Determine method from line; make sure matches one from file name
            method2key = next((x for x in lname2method.keys() if re.search(x, l)), False)
            if(re.search("balance", l)):
               bal= True

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

            #if bal: 
            #   method = method + "_Bal"
            #   bal = False
        elif ls[0] == "size:":
            size = ls[1] + "x" + ls[2] + "x" + ls[3]
        
        elif ls[0] == "datasize:":
            sizegb = float(ls[1]) / 1024**3

        elif ls[0] == "Output":
            outtime = float(ls[11][:-1])
            rate = sizegb/outtime
            balstr = "Yes" if bal else "No"
            if totalrate >= 0.0: totalrate += rate
            v.append([method, balstr, cores, size, sizegb, outtime, rate])

    if printFName: print fname
    print "%-16s %3s %5s %-16s %7s %6s %6s" % ("method", "bal", "cores", "size", "f_GB", "ftime", "fGB/s")
    for i in v:
        print "%-16s %3s %5d %-16s %7.1f %6.2f %6.2f" % tuple(i)
    if totalrate >= 0.0: print "Avg rate =", totalrate/len(v)

