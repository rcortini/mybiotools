import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import contact_matrix, distance_array

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
                           threshold=2.5) :
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
        chip_seq = np.zeros (N,dtype=int)
        for ts in u.trajectory[teq::tsample] :
            d = distance_array (polymer.positions,tracers.positions,
                                box=ts.dimensions)
            chip_seq += np.sum (d<threshold,axis=1)
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
        # initialize universe, tracers selection
        u = self.u
        tracers = u.select_atoms (tracer_text)
        ntracers = tracers.n_atoms
        # get the number of frames in the slice (http://stackoverflow.com/a/7223557)
        traj_slice = u.trajectory[teq::tsample]
        nslice = sum(1 for _ in traj_slice)
        # initialize the matrix containing all the positions
        # of the tracers at all the sampling frames
        tracers_pos = np.zeros ((nslice,ntracers,3))
        for i,ts in enumerate(u.trajectory[teq::tsample]) :
            tracers_pos[i,:,:] = tracers.positions
        # now initialize the Delta matrix, which contains the
        # squared differences between the tracers' positions
        # at different time delays
        Nt = int(nslice/2)
        Delta = np.zeros((ntracers,Nt,Nt))
        for delay in xrange(1,Nt+1) :
            for t0 in xrange (Nt) :
                t1 = t0 + delay
                pos1 = tracers_pos[t1,:,:]
                pos0 = tracers_pos[t0,:,:]
                Delta[:,delay-1,t0] = np.sum((pos1-pos0)**2,axis=1)
        # first calculate the mean over the tracers
        msd_average_tracers = np.mean(Delta,axis=0)
        # then calculate the mean over the starting points
        self.msd_t = np.mean(msd_average_tracers,axis=1)
