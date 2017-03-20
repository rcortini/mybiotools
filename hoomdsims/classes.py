import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import contact_matrix, distance_array
from . import simanalysis
from mybiotools.moremath import KL_divergence

class hoomdsim :
    def __init__ (self,xml,dcd=None) :
        if dcd is not None :
            u = mda.Universe (xml,dcd,format='DCD')
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

    def calculate_chipseq (self,polymer_text,tracer_text,teq,tsample,
                           threshold=2.5,aggregate=True) :
        """
        Calculates the 'ChIP-seq' profile of the tracers' contacts with the
        polymer. Users should provide the text to select the particles that
        belong to the polymer (polymer_text) and the one to select the tracers
        (tracer_text). User should provide the 'teq' and 'tsample' variables,
        which express, respectively, the number of frames to exclude at the
        start of the simulation, and the frequency at which the contact matrix
        should be calculated.
        
        Optional 'threshold' parameter for the thresholding of the contacts.
        """
        u = self.u
        tracers = u.select_atoms (tracer_text)
        polymer = u.select_atoms (polymer_text)
        ntracers = tracers.n_atoms
        N = polymer.n_atoms
        chip_seq = np.zeros ((N,ntracers),dtype=int)
        for ts in u.trajectory[teq::tsample] :
            d = distance_array (polymer.positions,tracers.positions,
                                box=ts.dimensions)
            chip_seq += d<threshold
        if aggregate :
            self.chip_seq = np.sum (chip_seq,axis=1)
        else :
            self.chip_seq = chip_seq

    def calculate_contact_trace (self,polymer_text,tracer_text,teq,tsample,
                                 threshold=2.5,nbins=50) :
        """
        Calculates the indices of the monomers with which each tracer was in
        contact in the trajectory. The return data is under the form of a list
        of lists, in which each element corresponds to each tracer. User should
        supply the text to select the polymer particles (polymer_text), and
        tracer particles (tracer_text). User should also provide the 'teq' and
        'tsample' variables, which express, respectively, the number of frames
        to exclude at the start of the simulation, and the frequency at which
        the contact matrix should be calculated. The jump size distribution is
        calculated based on these values.

        Optional 'threshold' parameter for the thresholding of the contacts.
        The 'nbins' optional parameter specifies the number of bins to use to
        calculate the distribution of jump sizes.
        """
        u = self.u
        polymer = u.select_atoms (polymer_text)
        tracers = u.select_atoms (tracer_text)
        ntracers = len (tracers)
        contact_trace = [[] for i in range (ntracers)]
        for ts in u.trajectory [teq::tsample] :
            d = distance_array (polymer.positions,tracers.positions,
                                box=ts.dimensions)
            c = np.where (d<threshold)
            if c :    
                for i in range(len(c[1])) :
                    tracer_i = c[1][i]
                    contact_i = c[0][i]
                    contact_trace[tracer_i].append (contact_i)
        self.contact_trace = contact_trace
        d = np.empty(0,dtype=int)
        for contact in contact_trace :
            d = np.concatenate ((d,np.abs(np.ediff1d (contact))))
        self.d_dist = np.histogram (d,bins=nbins,normed=True)

    def calculate_tracer_msd (self,tracer_text,teq,tsample) :
        self.msd_t = simanalysis.msd_t (self,tracer_text,teq,tsample)

    def calculate_r2_KLdiv_t (self,polymer_text,tracer_text,teq,tsample,threshold) :
        """
        Calculate r2 and KLdiv for the simulation C and R vectors (virtual
        ChIP-seq and row sum of virtual Hi-C, as a function of time
        """
        # get number of slices in the simulation and init KLdiv_t and r2_t
        u = self.u
        nframes = u.trajectory.n_frames
        nslice = simanalysis.traj_nslice (u,teq,tsample)
        self.KLdiv_t = np.zeros (nslice)
        self.r2_t = np.zeros (nslice)
        # get polymer and tracers
        polymer = u.select_atoms(polymer_text)
        tracers = u.select_atoms(tracer_text)
        npolymer = polymer.n_atoms
        ntracers = tracers.n_atoms
        # init C, R, and H vectors
        C = np.zeros (npolymer)
        R = np.zeros (npolymer)
        H = np.zeros ((npolymer,npolymer))
        # iterate on all frames in slice
        for i,ts in enumerate(u.trajectory[teq::tsample]) :
            # ChIP-seq
            d = distance_array (polymer.positions,tracers.positions,box=ts.dimensions)
            C += np.sum(d<threshold,axis=1)
            # Hi-C
            d = distance_array (polymer.positions,polymer.positions,box=ts.dimensions)
            H[d<threshold] += 1.0
            R = H.sum(axis=1)
            # now calculate KLdiv and r2 at this time point
            self.KLdiv_t[i] = KL_divergence(C,R)
            self.r2_t[i] = np.corrcoef(C,R)[0,1]**2
