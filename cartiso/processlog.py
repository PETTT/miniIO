#!/usr/bin/env python

"""
Process an output log from cartiso and produce a results table
  + Determines output type (full or iso) and method from log contents
  + Can also check filename for type & method
    - Filename must contain .pvti., .pvtp., .posix., .mpi., etc (see fname2method)
    - If filename & log contents don't agree, will abort with error
    - If you didn't name file properly (or at all), just gives a warning
  + Each column is a time step
  + Table rows(Full): method, cores, grid size, full file size, time, rate (GB/s)
  + Table rows(Iso): method, cores, grid size, full file size, virtual time, virtual rate
                     iso phase, load imbalance, imbalance %, triangles, iso size, 
                     iso time, iso rate
    - "virtual" time is iso computation time + iso output time
  + Multiple methods/types can be pasted into Excel or similar for plotting/analysis
  + Give -a option to compute an average rate at the end
  + Options:
      -a => print average rate(s) after each file
      -n => don't print the column headings (good for auto-copy-paste to excel)
"""

import sys
import re

if len(sys.argv) <= 1:
    print __doc__

phases = ["sin2gauss", "gaussmove", "gaussresize", "gaussmove_r"]
name2otype = { "pvti": "Full", "pvtp": "Iso", "adiosfull": "Full", "adiosiso": "Iso",
        "hdf5i": "Full", "hdf5p": "Iso" }
fname2method = { ".pvti.": "MPI-Indiv", ".pvtp.": "MPI-Indiv", ".posix.": "ADIOS-POSIX", 
        ".mpi.": "ADIOS-MPI", ".mpilus.": "ADIOS-MPI_Lustre", ".mpiag.": "ADIOS-MPI_Aggr",
        ".phdf5.": "ADIOS-PHDF5", ".hdf5i.": "HDF5", ".hdf5p.": "HDF5" }
lname2method = { "pvti": "MPI-Indiv", "pvtp": "MPI-Indiv", "POSIX": "ADIOS-POSIX", 
        "MPI ": "ADIOS-MPI", "MPI$": "ADIOS-MPI", "MPI_LUSTRE": "ADIOS-MPI_Lustre",
        "MPI_AGGREGATE": "ADIOS-MPI_Aggr", "PHDF5": "ADIOS-PHDF5", 
        "hdf5i": "HDF5", "hdf5p": "HDF5" }

totalrate = -1.0
totalvrate = -1.0
header = True
a = 1
for arg in sys.argv[1:3]:
    if arg == "-a":   # produce average rate
        totalrate = 0.0
        totalvrate = 0.0
        a += 1
    if arg == "-n":   # no header
        header = False
        a += 1

# Iter over all files
for fname in sys.argv[a:]:
    v = []
    if totalrate > 0.0: totalrate = 0.0
    if totalvrate > 0.0: totalvrate = 0.0

    # Determine type and method from file name, if it's there
    otype1key = next((x for x in name2otype.keys() if x in fname), False)
    if not otype1key: print "WARNING: type not specified in filename"
    method1key = next((x for x in fname2method.keys() if x in fname), False)
    if not method1key: print "WARNING: method not specified in filename"

    f = open(fname)
    phase = -1
    for l in f:
        ls = l.split()
        if(len(ls) < 1): continue
        if(len(ls) < 2): ls.append("")   # Ensure all outer if's can handle ls

	if ls[0][0:4] == "tsk=":
            if ls[0][4] == "\"": continue    # Skip LSF copy of batch script

	    cores = int(ls[0][4:]) * int(ls[1]) * int(ls[2])
            if ls[3][:3] == "nc=":    # Added nc= field to SGI output, and size not cubed
                size = ls[7][2:] + 'x' + ls[8] + 'x' + ls[9]
                sizeprod = float(ls[7][2:]) * float(ls[8]) * float(ls[9])
                fullgb = sizeprod*4*2/1024**3
            else:
                size = ls[6][2:]
                sizec = int(ls[6][2:])
                fullgb = float(sizec)**3*4*2/1024**3

            # Determine type & method from line; make sure matches one from file name
            otype2key = next((x for x in name2otype.keys() if x in l), False)
            if otype2key:
                if otype1key: 
                    if otype1key != otype2key:
                        print "ERROR: type in filename doesn't match type in file"
                        sys.exit(1)
                    else:
                        otype = name2otype[otype1key]
                else:
                    otype = name2otype[otype2key]
            else:
                print "WARNING: type not specified in log file"
                if otype1key:
                    otype = name2otype[otype1key]
                else:
                    print "ERROR: no type specified, cannot proceed"
                    sys.exit(1)
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

	elif ls[0] == "isothresh":
	    phase += 1
	elif ls[0] == "Total" and ls[1] == "tris":
	    tris = float(ls[3][:-1])
            # Size of output depends upon method
            if method == "MPI-Indiv":
                # ... (coords + var + norms + conns + offsets) ...
                isogb = tris*(3*3*4 + 3*4 + 3*3*4 + 3*8 + 8) / 1024**3
            else:
                # ... (coords + var + norms + conns) ...
                isogb = tris*(3*3*4 + 3*4 + 3*3*3 + 3*8) / 1024**3
	    limb = float(ls[13])
	elif ls[0] == "FullOutput":
	    fulltime = float(ls[11][:-1])
	elif (ls[0],ls[1]) == ("Isosurface","timer"):
	    isoctime = float(ls[11][:-1])
	elif ls[0] == "IsoOutput":
	    isotime = float(ls[11][:-1])
            if otype == "Full":
                rate = fullgb/fulltime
                if totalrate >= 0.0: totalrate += rate
                v.append([method, cores, size, fullgb, fulltime, rate])
            elif otype == "Iso":
                vtime = isoctime+isotime  #equivlent time if iso+isoout replaces full output
                vrate = fullgb/vtime
                rate = isogb/isotime
                if totalrate >= 0.0: totalrate += rate
                if totalvrate >= 0.0: totalvrate += vrate
                v.append([method, cores, size, fullgb, vtime, vrate, phases[phase],
                          limb, limb/cores, tris, isogb, isotime, rate])

    if otype == "Full":
        if header: print "%-16s %5s %15s %5s %5s %5s" % ("method", "cores", "size", "f_GB", "ftime", "fGB/s")
        for i in v:
            print "%-16s %5d %15s %5.0f %5.2f %5.2f" % tuple(i)
        if totalrate >= 0.0: print "Avg rate =", totalrate/len(v)
    elif otype == "Iso":
        if header: print "%-16s %5s %15s %5s %5s %6s %11s %8s %5s %11s %7s %4s %6s" % ("method", "cores", 
            "size", "f_GB", "vtime", "vGB/s", "phase", "loadimb", "limb%", "tris", "i_GB", "itim", "i_GB/s")
        for i in v:
            print "%-16s %5d %15s %5.0f %5.2f %6.2f %11s %8.2f %5.3f %11.0f %7.3f %4.2f %6.3f" % tuple(i)
        if totalrate >= 0.0: print "Avg vrate =", totalvrate/len(v), ", rate =", totalrate/len(v)
    
