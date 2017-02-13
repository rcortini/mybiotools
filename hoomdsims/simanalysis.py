import numpy as np
from MDAnalysis.analysis.distances import distance_array
import mybiotools as mbt

def traj_nslice (u,teq,tsample) :
    """
    Returns the number of frames in the trajectory in universe u, using teq as
    equilibration time and tsample as sampling time
    """
    # get the number of frames in the slice (http://stackoverflow.com/a/7223557)
    traj_slice = u.trajectory[teq::tsample]
    return sum(1 for _ in traj_slice)

def hic_chipseq_r2 (hic, chipseq) :
    """
    Calculate the Pearson correlation coefficient between the row sum of the
    given Hi-C matrix and the given ChIP-seq profile.
    """
    hic_rowsum = np.sum(hic,axis=1)/float(np.sum(hic))
    return np.corrcoef(hic_rowsum,chipseq)[0,1]**2

def ps (H) :
    """
    Calculate the normalized probability of contact between a monomer and all
    others as a function of the linear distance s.
    """
    p = np.array ([np.mean (np.diagonal (H, offset=k))
                   for k in range (H.shape[0])])
    return p/np.sum(p)

def contacts_with (sim,polymer_text,tracers_text,bindingsites_text,teq,tsample,threshold) :
    """
    Calculate the relative proportion of contacts of the tracers with binding
    sites compared with non-binding sites. As usual user should supply
    equilibration time, sampling time, and contact threshold value.
    """
    # select polymer, tracers, and binding sites
    polymer = sim.u.select_atoms (polymer_text)
    tracers = sim.u.select_atoms (tracers_text)
    bss = sim.u.select_atoms (bindingsites_text)
    # select binding site indices
    bs_n = bss.n_atoms
    bs_idx = bss.indices
    # select non-binding site indices
    polymer_idx = polymer.indices
    nbs_idx = np.setdiff1d (polymer_idx,bs_idx)
    nbs_n = nbs_idx.size
    # evaluate contacts with binding sites and non-binding sites for each
    # independent simulation snapshot
    c = []
    for i,ts in enumerate(sim.u.trajectory[teq::tsample]) :
        d = distance_array (polymer.positions,tracers.positions,
                            box=ts.dimensions)
        contacts = d<threshold
        cB = np.sum (contacts[bs_idx]).astype('float')
        cA = np.sum (contacts[nbs_idx]).astype('float')
        if cA != 0 :
            c.append ((cB/cA) / (float(bs_n)/nbs_n))
    return np.mean(np.array(c))

def fit_msd (msd,cutoff,delta_t,scale_l) :
    """
    Perform a simple fit of the supplied time-dependent MSD, using a linear
    regression of the logarithms of the values. User must supply the conversion
    factor from time to real time and from length to real length. Also, user
    must supply the cutoff value: from there on the values will be considered.
    This is because the long-time behaviour is generally what matters really.
    """
    # prepare the values to fit: exclude the first value because it is zero
    t = np.arange(msd.size)*delta_t
    x = np.log(t[cutoff:])
    y = np.log(msd[cutoff:]*scale_l**2)
    # perform fit to y = ax + b with their errors
    b,a,db,da = mbt.linear_regression (x,y,0.99)
    # now convert the value of b into a diffusion coefficient
    D = np.exp(b)/6.0
    dD = np.exp(db)/6.0
    return a,da,D,dD

def msd_t (sim,particles_text,teq,tsample) :
    """
    Calculate the mean square displacement of the particles defined by
    'particles_text' in simulation sim, using sampling tsample and equilibration
    time teq. Returns the matrix corresponding to the mean square displacement
    of each particle, along with a matrix corresponding to the variance in the
    estimate of this quantity.
    """
    u = sim.u
    particles = u.select_atoms (particles_text)
    nparticles = particles.n_atoms
    nslice = traj_nslice (u,teq,tsample)
    # initialize the matrix containing all the positions
    # of the particles at all the sampling frames
    particles_pos = np.zeros ((nslice,nparticles,3))
    for i,ts in enumerate(u.trajectory[teq::tsample]) :
        particles_pos[i,:,:] = particles.positions
    # now initialize the Delta matrix, which contains the
    # squared differences between the particles' positions
    # at different time delays
    Nt = int(nslice/2)
    Delta = np.zeros((nparticles,Nt,Nt))
    for delay in xrange(1,Nt+1) :
        for t0 in xrange (Nt) :
            t1 = t0 + delay
            pos1 = particles_pos[t1,:,:]
            pos0 = particles_pos[t0,:,:]
            Delta[:,delay-1,t0] = np.sum((pos1-pos0)**2,axis=1)
    # return the matrices of MSD and its variance
    return np.mean(Delta,axis=2),np.var(Delta,axis=2)

def dmin_sel (sim,sel1_text,sel2_text,teq,tsample) :
    """
    Calculate the minimum distance between the atoms defined in sel1 and the
    atoms defined in sel2, as a function of time. Returns a matrix that contains
    the minimum distance for each atom defined in sel1. As usual user should
    supply equilibration time, sampling time, and contact threshold value.
    """
    # define atom selections
    sel1 = sim.u.select_atoms (sel1_text)
    sel2 = sim.u.select_atoms (sel2_text)
    # get number of atoms in selection 1
    natoms = sel1.n_atoms
    nslice = traj_nslice (sim.u,teq,tsample)
    dmin = np.zeros((natoms,nslice))
    for i,ts in enumerate(sim.u.trajectory[teq::tsample]) :
        d = distance_array (sel1.positions,sel2.positions,
                            box=ts.dimensions)
        dmin[:,i] = d.min(axis=1)
    return dmin

def particle_images (sim,frame_id) :
    """
    Get the image index of all particles in simulation, at the frame 'frame_id'
    """
    # get positions of all particles: define first the atom selection, then jump to
    # the user-requested trajectory frame, get the box dimensions (currently works
    # only for orthorhombic boxes, then calculate the image indices
    atoms = sim.u.select_atoms ('all')
    ts = sim.u.trajectory[frame_id]
    L = ts.dimensions[:3]
    pos = atoms.positions + L/2.
    return pos//L

def jumping_matrix (sim,polymer_text,tracer_text,teq,tsample,threshold) :
    """
    Calculate the matrix that represents the number of times that the tracers
    (defined by 'tracer_text') jump from one site to another site of the polymer
    (defined by 'polymer_text'). The simulation 'sim' is sampled at 'tsample',
    excluding the first 'teq' time frames. Contact between a tracer and the
    polymer is defined by the distance being smaller than 'threshold'.
    """
    # define polymer and tracers
    u = sim.u
    polymer = u.select_atoms(polymer_text)
    tracers = u.select_atoms(tracer_text)
    n_polymer = polymer.n_atoms
    n_tracers = tracers.n_atoms
    # initialize jumping matrix and first distance matrix d_prev
    J = np.zeros ((n_polymer,n_polymer),dtype=np.int32)
    ts = u.trajectory [teq]
    d_prev = distance_array (polymer.positions,tracers.positions,
                            box=ts.dimensions)
    D_prev = d_prev<threshold
    for ts in u.trajectory [teq::tsample] :
        # get distance matrix at current time step
        d_next = distance_array (polymer.positions,tracers.positions,
                            box=ts.dimensions)
        D_next = d_next<threshold
        # get jumps of all tracers and add it to the jumping matrix
        for i in xrange (n_tracers) :
            t_prev = D_prev [:,i]
            t_next = D_next [:,i].reshape ((n_polymer,1))
            t = t_prev * t_next
            J += t
        D_prev = D_next.copy()
    return J
