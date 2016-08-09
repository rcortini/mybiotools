import networkx as nx
import numpy as np

class Lattice :
    def get_adjacency_matrix (self) :
        return np.asarray (nx.to_numpy_matrix(self.get_graph()))
    def draw (self) :
        nx.draw_spectral (self.get_graph(),with_labels=True)
