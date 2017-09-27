import numpy as np

class Vertex(object):
    def __init__(self):
        self.nIndices = None # Indices, in counter-clockwise direction

class Facet(object):
    def __init__(self, v1, v2, v3):
        self.vIndices = (v1, v2, v3) # Indices of vertices in counter-clockwise direction
    def __eq__(self, obj):
        # Given that all three indices are different
        objIndex = next((x for x in range(3) if obj.vIndices[x] == self.vIndices[0]), None)
        if objIndex is None:
            return False
        return all(self.vIndices[x] == obj.vIndices[(x + objIndex) % 3] for x in range(3))

class Meshwork(object):
    def __init__(self):
        self.vertices = []
        self.facets = []

class MeshworkLoader(object):
    def __init__(self, fileName):
        self.fileName = fileName
    
    def loadTo(self, meshworkObj):
        print("Loading topology...")
        
        vertexIndex = 0
        with open(self.fileName) as f:
            for eachLine in f:
                neighborInfo = np.fromstring(eachLine, dtype=int, sep='\t')
                self._parseNeighborInfo(vertexIndex, neighborInfo, meshworkObj)
                vertexIndex += 1

                if(vertexIndex % 50 == 0):
                    print("Processing vertices: %d..." % vertexIndex)
        
        print("Total facets: %d" % len(meshworkObj.facets))
    
    def _parseNeighborInfo(self, vertexIndex, neighborInfo, meshworkObj):
        numNeighbors = len(neighborInfo)

        # Vertices
        newVertex = Vertex()
        newVertex.nIndices = neighborInfo
        meshworkObj.vertices.append(newVertex)

        # Facets
        for neighborInfoIdx in range(numNeighbors):
            newFacet = Facet(vertexIndex, neighborInfo[neighborInfoIdx], neighborInfo[(neighborInfoIdx + 1) % numNeighbors])
            if next((x for x in meshworkObj.facets if x == newFacet), None) is None:
                meshworkObj.facets.append(newFacet)


        