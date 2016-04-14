# ParaView Programmable Source: vtkUnstructuredGrid
from paraview import vtk
import numpy as np
filename = "/Users/sean/sync/wrk/ace4/miniIO/unstruct/unstruct.przm/t0000.d/r00.dat"

f = open(filename)
output = self.GetOutput()

numPoints = np.fromfile(f, np.uint64, 1)[0]
hasGrid = np.fromfile(f, np.uint32, 1)[0]
if hasGrid:
    xpts = np.fromfile(f, np.float32, numPoints)
    ypts = np.fromfile(f, np.float32, numPoints)
    zpts = np.fromfile(f, np.float32, numPoints)
    pts = vtk.vtkPoints()
    for i in range(numPoints):
        pts.InsertNextPoint(xpts[i], ypts[i], zpts[i])
    output.SetPoints(pts)

numCells = int( np.fromfile(f, np.uint64, 1)[0] )
if numCells:
    conns = np.fromfile(f, np.uint64, numCells*6)
    output.Allocate(numCells, 1000)
    for i in range(numCells):
        pointIds = vtk.vtkIdList()
        for j in range(6):
            pointIds.InsertId(j, int(conns[i*6+j]))
        output.InsertNextCell(13, pointIds)
f.close()

