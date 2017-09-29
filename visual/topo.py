import numpy as np

class Vertex(object):
    def __init__(self):
        self.nIndices = None # Neighbor indices, in counter-clockwise direction

class Facet(object):
    def __init__(self, vIndices):
        self.vIndices = vIndices # Indices of vertices in counter-clockwise direction (must have length 3)
    def __eq__(self, obj):
        """Returns whether two facets are identical.

        Two triangles are considered identical if they contain same indices
        and are in the same loop sequence. For example (1, 4, 5) == (4, 5, 1),
        but (1, 4, 5) != (1, 5, 4)
        """
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
    def __init__(self, neighborFileName, triangleFileName):
        self.neighborFileName = neighborFileName
        self.triangleFileName = triangleFileName
    
    def loadTo(self, meshworkObj):
        print("Loading topology...")
        
        # Vertices
        try:
            with open(self.neighborFileName) as f:
                for eachLine in f:
                    neighborInfo = np.fromstring(eachLine, dtype=int, sep='\t')
                    self._parseNeighborInfo(neighborInfo, meshworkObj)

        except IOError as e:
            print("Cannot open the file containing neighboring vertices.")
            raise
        print("Total vertices: %d" % len(meshworkObj.vertices))
        
        # Triangles
        try:
            with open(self.triangleFileName) as f:
                for eachLine in f:
                    triangleInfo = np.fromstring(eachLine, dtype=int, sep='\t')
                    self._parseTriangleInfo(triangleInfo, meshworkObj)
        except IOError as e:
            print("Cannot open the file containing triangles.")
            raise
        print("Total triangles: %d" % len(meshworkObj.facets))
    
    def _parseNeighborInfo(self, neighborInfo, meshworkObj):
        newVertex = Vertex()
        newVertex.nIndices = neighborInfo
        meshworkObj.vertices.append(newVertex)

    def _parseTriangleInfo(self, triangleInfo, meshworkObj):
        meshworkObj.facets.append(Facet(triangleInfo))
        