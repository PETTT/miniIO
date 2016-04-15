# ParaView Programmable Source: vtkUnstructuredGrid
# Read a przm output collection
# For now, reads parallel output serially
from paraview import vtk
import numpy as np
from glob import glob
przmname = "/Users/sean/sync/wrk/ace4/miniIO/unstruct/unstruct.przm"
tstep = 0

tfiles = glob(przmname+"/t*.d")
rfiles = glob(tfiles[tstep]+"/r*.dat")
nr = len(rfiles)

f = [None]*nr;   # One file handle per rank file
hasGrid = False
pts = vtk.vtkPoints()

for r,rfile in enumerate(rfiles):
    f[r] = open(rfile)

    numPoints = np.fromfile(f[r], np.uint64, 1)[0]
    hasGrid = np.fromfile(f[r], np.uint32, 1)[0]
    if hasGrid:
        xpts = np.fromfile(f[r], np.float32, numPoints)
        ypts = np.fromfile(f[r], np.float32, numPoints)
        zpts = np.fromfile(f[r], np.float32, numPoints)
        for i in range(len(xpts)):
            pts.InsertNextPoint(xpts[i], ypts[i], zpts[i])
if hasGrid: output.SetPoints(pts)

needalloc = True

for r in range(nr):
    numCells = np.fromfile(f[r], np.uint64, 1)[0] 
    numConns = numCells * np.uint64(6)
    if numCells:
        if needalloc:
            output.Allocate(numCells * np.uint64(nr), 1000)
            needalloc = False
        conns = np.fromfile(f[r], np.uint64, numConns)
        for i in range(numCells):
            pointIds = vtk.vtkIdList()
            for j in range(6):
                pointIds.InsertId(j, conns[i*6+j])
            output.InsertNextCell(13, pointIds)

tris = []
numPoints = numPoints * np.uint64(nr)

for r in range(nr):
    numTris = int(np.fromfile(f[r], np.uint64, 1)[0])
    if numTris:
        tris = np.concatenate((tris, np.fromfile(f[r], np.uint64, numTris)))
        print tris.shape
if type(tris) is not list:
    sdata = np.zeros(numPoints, dtype=np.uint8)
    sdata[np.uint32(tris)] = 1
    output.PointData.append(sdata, "surfflag")
        
for r in range(nr):
    f[r].close()

