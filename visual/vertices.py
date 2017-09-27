import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

import topo


class VertexLoader(object):
    def __init__(self, fileName):
        self.fileLocation = fileName

    def loadVertex(self):
        rawData = np.loadtxt(self.fileLocation)
        size = rawData.shape

        (self.numSnapshots, self.numVertices) = \
            (size[0] // 3, 1) if len(size) < 2 else (size[0], size[1] // 3)
        self.data = rawData.reshape(self.numSnapshots, self.numVertices, 3)

    def loadTopo(self):
        topoLoader = topo.MeshworkLoader(r"C:\Users\drels\OneDrive\Documents\Source\Repos\MembraneSimulation\MeshGeneration\neighbors.txt")
        self.meshwork = topo.Meshwork()
        topoLoader.loadTo(self.meshwork)
    
    def loadAll(self):
        self.loadVertex()
        self.loadTopo()


class Plottor(object):
    def __init__(self):
        self.fig = plt.figure()
        self.fig.patch.set_alpha(0.0)
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.patch.set_alpha(0.0)

    
        self.ax.axis('off')

    def plotVertices(self, loader, snapshotNum, sort=None):
        """
        sort, if specified, is a function taking vertex info (numVertices by 3),
        and give a boolean list of length of numVertices
        """
        toBePlotted = loader.data[snapshotNum, :, :]
        if sort != None:
            toBePlotted = toBePlotted[sort(toBePlotted)]
        
        self.ax.scatter(
            toBePlotted[:, 1],
            toBePlotted[:, 2],
            toBePlotted[:, 0]
        )
    
    def plotTriangles(self, loader, snapshotNum, sort=None):
        toBePlotted = loader.data[snapshotNum, :, :]

        areSortedIndices = sort(toBePlotted) if sort else [True] * loader.numVertices

        newToBePlotted = toBePlotted[areSortedIndices]

        # newIndexMap is something like [N, 0, N, 1, 2, N, N, 3, ...]
        indexIter = iter(range(len(newToBePlotted)))
        newIndexMap = [(next(indexIter) if areSortedIndices[x] else None) for x in range(len(toBePlotted))]

        newTri = [[newIndexMap[eachIndex] for eachIndex in eachTri.vIndices] for eachTri in loader.meshwork.facets if np.all([(newIndexMap[eachIndex] is not None) for eachIndex in eachTri.vIndices])]

        tri = matplotlib.tri.Triangulation(
            newToBePlotted[:, 1], # x (actually y)
            newToBePlotted[:, 2], # y (actually z)
            newTri,
            mask=None
        )
        
        self.ax.plot_trisurf(tri, newToBePlotted[:, 0], alpha=0.3, linewidth=1, edgecolor='k')

if __name__ == '__main__':
    loader = VertexLoader(r'C:\Users\drels\OneDrive\Documents\Source\Repos\MembraneSimulation\MembraneSimulation\p_out.SimOut')
    loader.loadAll()

    plottor = Plottor()
    sortVertices = lambda x: x[:, 0] > 0.97e-6
    wantedSnapshot = 6
    filamentTipX = wantedSnapshot * 0.005e-6 + 0.99
    #plottor.plotVertices(loader, wantedSnapshot, sort=sortVertices)
    plottor.plotTriangles(loader, wantedSnapshot, sort=sortVertices)
    plottor.ax.set_xlim(-0.25e-6, 0.25e-6)
    plottor.ax.set_ylim(-0.25e-6, 0.25e-6)
    plottor.ax.set_zlim(0.97e-6, 1.10e-6)

    plt.show()
