import networkx as nx
import numpy as np

class PrimaryTriangle :
    """
    A Primary Triangle of the Sierpinski Gasket is defined as a triangle that
    will be subdivided into three smaller triangles at a successive iteration.
    """
    def __init__ (self,vertices) :
        self.vertices = vertices
    def edges (self) :
        v = self.vertices
        return [(v[i-1],v[i]) for i in range(3)]

def iterate_triangle (triangle,num_vertices) :
    newv = [num_vertices+i for i in range(3)]
    v = triangle.vertices
    new_v = [[v[0],newv[0],newv[2]],[newv[0],v[1],newv[1]],[newv[2],newv[1],v[2]]]
    return [PrimaryTriangle(new_vertex) for new_vertex in new_v]

def allprimarytriangles (trianglelist) :
    pt = []
    for triangle in trianglelist :
        if isinstance (triangle,list) :
            pt.extend (allprimarytriangles(triangle))
        else :
            pt.append (triangle)
    return pt

class SierpinskiGasket :
    """
    This class instantiates an object that contains a nested list of lists, at
    the bottom of which there are PrimaryTriangle objects. The class can be used
    to draw the triangle (in its network representation) or get the adjacency
    matrix of the network.
    """
    def __init__ (self,generation=1) :
        firstvertices = [0,1,2]
        self.triangles = [PrimaryTriangle (firstvertices)]
        self.num_vertices = 3
        if generation != 1 :
            for i in range(1,generation) :
                self.iterate ()
    def iterate_trianglelist (self,triangle) :
        if isinstance (triangle,list) :
            nested_list = [self.iterate_trianglelist(t) for t in triangle]
            return nested_list
        else :
            new_triangles = iterate_triangle (triangle,self.num_vertices)
            self.num_vertices += 3
            return new_triangles
    def iterate (self) :
        self.triangles = self.iterate_trianglelist (self.triangles)
    def get_graph (self) :
        G = nx.Graph()
        for i in range (self.num_vertices) :
            G.add_node (i)
        for t in allprimarytriangles(self.triangles) :
            for edge in t.edges() :
                G.add_edge (edge[0],edge[1])
        return G
    def get_adjacency_matrix (self) :
        return np.asarray (nx.to_numpy_matrix(self.get_graph()))
    def draw (self) :
        nx.draw_spectral (self.get_graph(),with_labels=True)
