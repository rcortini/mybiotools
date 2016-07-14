import MDAnalysis as mda
from MDAnalysis.analysis.distances import contact_matrix, distance_array

class hoomdsim :
    def __init__ (self,xml,dcd=None) :
        if dcd is not None :
            u = mda.Universe (xml,dcd)
        else :
            u = mda.Universe (xml)
        self.u = u

    def calculate_hic (self,polymer_text,teq,tsample,threshold=2.5) :
        """
        Calculates the 'Hi-C' matrix of the polymer, given the polymer_text
        variable that expresses how the code should select the particles that
        belongs to the polymer. User should provide the 'teq' and 'tsample'
        variables, which express, respectively, the number of frames to exclude
        at the start of the simulation, and the frequency at which the contact
        matrix should be calculated.
        
        Optional 'threshold' parameter for the
        thresholding of the contacts.
        """
        u = self.u
        polymer = u.select_atoms (polymer_text)
        N = polymer.n_atoms
        hic = np.zeros((N,N),dtype=int)
        for ts in u.trajectory[teq::tsample] :
            hic += contact_matrix(polymer.positions,
                                  cutoff=threshold,
                                  box=ts.dimensions).astype(int)
        self.hic = hic
