import networkx as nx
from .lattice import Lattice

def site2D_id (i,j,n) :
    """
    A function that returns the index of the site on the square lattice
    """
    return i*n + j

class Site2D :
    """
    The lattice is made of simple two-index sites
    """
    def __init__ (self,i,j) :
        self.x = (i,j)

class SquareLattice (Lattice) :
    """
    The square lattice class. It is instantiated using the number of sites on
    each edge.
    """
    def __init__ (self,n) :
        self.sites = []
        self.n = n
        for i in range (n) :
            for j in range (n) :
                self.sites.append (Site2D(i,j))
    def get_graph (self) :
        G = nx.Graph()
        # first add all sites, then all edges
        for site in self.sites :
            G.add_node (site)
        for site in self.sites :
            i,j = site.x
            if i>0 :
                site_id = site2D_id(i-1,j,self.n)
                G.add_edge (site,self.sites[site_id])
            if j>0 :
                site_id = site2D_id(i,j-1,self.n)
                G.add_edge (site,self.sites[site_id])
        return G
