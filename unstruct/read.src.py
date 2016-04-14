# ParaView Programmable Source: vtkUnstructuredGrid
# Read a przm output collection
# For now, reads parallel output serially
from paraview import vtk
import numpy as np
from glob import glob
#filename = "/Users/sean/sync/wrk/ace4/miniIO/unstruct/unstruct.przm/t0000.d/r00.dat"
przmname = "/Users/sean/sync/wrk/ace4/miniIO/unstruct/unstruct.przm"
tstep = 0

tfiles = glob(przmname+"/t*.d")
rfiles = glob(tfiles[tstep]+"/r*.dat")
nr = len(rfiles)

output = self.GetOutput()
xpts = [None]*nr; ypts = [None]*nr; zpts = [None]*nr; conns = [None]*nr;
hasGrid = False;  numCells = 0;  numConns = 0;

for r,rfile in enumerate(rfiles):
    f = open(rfile)

    numPoints = np.fromfile(f, np.uint64, 1)[0]
    hasGrid = np.fromfile(f, np.uint32, 1)[0]
    if hasGrid:
        xpts[r] = np.fromfile(f, np.float32, numPoints)
        ypts[r] = np.fromfile(f, np.float32, numPoints)
        zpts[r] = np.fromfile(f, np.float32, numPoints)

    numCells = np.fromfile(f, np.uint64, 1)[0] 
    numConns = numCells * np.uint64(6)
    if numCells:
        conns[r] = np.fromfile(f, np.uint64, numConns)
    f.close()

if hasGrid:
    pts = vtk.vtkPoints()
    for r in range(nr):
            for i in range(len(xpts[r])):
                pts.InsertNextPoint(xpts[r][i], ypts[r][i], zpts[r][i])
    output.SetPoints(pts)

if numCells:
    output.Allocate(numCells * np.uint64(nr), 1000)
    for r in range(nr):
        for i in range(numCells):
            pointIds = vtk.vtkIdList()
            for j in range(6):
                pointIds.InsertId(j, conns[r][i*6+j])
            output.InsertNextCell(13, pointIds)

