import numpy as np
import os
import mybiotools as mbt

# general variables that are used throughout the calculations
sim_root_dir = os.getenv('HOME')+'/work/data/hoomd/sbs_tracers/production/2016-10-24'
sim_base_name = 'sbs_tracers'

# simulation variables
nsims = 10
N = 1024
nframes = 10000
phivals = np.loadtxt('%s/phi_vals'%(sim_root_dir))
evals = np.loadtxt('%s/e_vals'%(sim_root_dir))

# units of measure
scale_l = 1.5e-6                                  # cm
scale_m = 5.73e-18                                # g
scale_e = 4.1e-14                                 # erg
scale_t = scale_l * np.sqrt(scale_m/scale_e)      # s

def sim_name (phi,e,n) :
    """
    Name of the simulation and basename of the simulation
    """
    basename = 'sbs_tracers-phi-%.2f-e-%.1f'%(phi,e)
    simname = '%s-%d'%(basename,n)
    return basename, simname

def sim_directory(sim_root_dir,phi,e,n) :
    """
    Returns the directory containing the simulation data
    """
    basename, simname = sim_name(phi,e,n)
    return '%s/%s/%s'%(sim_root_dir,basename,simname)

class sbs_tracers_sim :
    """
    This class allows to load the data of the simulations. It is conceived so
    that it loads data that is _already_ stored into files on the hard drive. To
    perform the calculations directly on the trajectory, one should use the
    'hoomdsim' class, which allows to open the trajectory information.
    """
    def __init__ (self,sim_root_dir,phi,e,n) :
        self.phi = phi
        self.e = e
        self.n = n
        basename,simname = sim_name(phi,e,n)
        simdir = sim_directory(sim_root_dir,phi,e,n)
        self.chipseq = None
        self.hic = None
        chipseq_file = '%s/chipseq.npy'%simdir
        hic_file = '%s/hic.npy'%simdir
        if os.path.exists(chipseq_file) and os.path.exists(hic_file) :
            self.chipseq = np.load('%s/chipseq.npy'%simdir)
            self.hic = np.load('%s/hic.npy'%simdir)

def load_all_sims () :
    """
    Returns a dictionary containing all the simulations, in a way that one can
    access simulation data using sims[(phi,e,n)].
    """
    sims = {}
    for phi in phivals :
        for e in evals :
            mbt.log_message('load_all_sims','Loading (%.2f,%.1f)'%(phi,e))
            for sim in xrange(1,nsims+1) :
                sims[(phi,e,sim)] = sbs_tracers_sim(sim_root_dir,phi,e,sim)
    return sims
