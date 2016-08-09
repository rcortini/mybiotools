import networkx as nx

def site3d_id (i,j,k,n) :
    """
    A function that returns the index of the site on the cubic lattice
    """
    return i*n**2 + j*n + k

class Site3D :
    """
    The lattice is made of simple three-index sites
    """
    def __init__ (self,i,j,k) :
        self.x = (i,j,k)

class CubicLattice :
    """
    The cubic lattice class. It is instantiated using the number of sites on
    each edge.
    """
    def __init__ (self,n) :
        self.sites = []
        self.n = n
        for i in range (n) :
            for j in range (n) :
                for k in range (n) :
                    self.sites.append (Site3D(i,j,k))
    def get_graph (self) :
        G = nx.Graph()
        # first add all sites, then all edges
        for site in self.sites :
            G.add_node (site)
        for site in self.sites :
            i,j,k = site.x
            if i>0 :
                site_id = site3d_id (i-1,j,k,self.n)
                G.add_edge (site,self.sites[site_id])
            if j>0 :
                site_id = site3d_id (i,j-1,k,self.n)
                G.add_edge (site,self.sites[site_id])
            if k>0 :
                site_id = site3d_id (i,j,k-1,self.n)
                G.add_edge (site,self.sites[site_id])
        return G
