import numpy as np
from MDAnalysis.analysis.distances import distance_array
import mybiotools as mbt

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

def fit_msd (msd,delta_t,scale_l) :
    """
    Perform a simple fit of the supplied time-dependent MSD, using a linear
    regression of the logarithms of the values. User must supply the conversion
    factor from time to real time and from length to real length
    """
    # prepare the values to fit: exclude the first value because it is zero
    t = np.arange(msd.size)*delta_t
    x = np.log(t[1:])
    y = np.log(msd[1:]*scale_l**2)
    # perform fit to y = ax + b with their errors
    b,a,db,da = mbt.linear_regression (x,y,0.99)
    # now convert the value of b into a diffusion coefficient
    D = np.exp(b)/6.0
    dD = np.exp(db)/6.0
    return a,da,D,dD
