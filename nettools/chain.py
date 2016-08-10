import networkx as nx
from .lattice import Lattice

class Chain (Lattice) :
    def __init__ (self,N) :
        self.N = N
    def get_graph (self) :
        G = nx.Graph()
        for i in range (self.N) :
            G.add_node(i)
        for i in range (1,self.N) :
            G.add_edge(i-1,i)
        return G
