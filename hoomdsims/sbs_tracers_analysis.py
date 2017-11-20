import numpy as np
import os
import mybiotools as mbt

# general variables that are used throughout the calculations
root_dir = os.getenv('HOME')+'/work/data/hoomd/sbs_tracers'
production_dir = '%s/production'%root_dir
sim_base_name = 'sbs_tracers'

# simulation variables
N = 1024
ntracers = 10
nframes = 10000
phivals = np.loadtxt('%s/phi_vals'%(production_dir))
evals = np.loadtxt('%s/e_vals'%(production_dir))
nsims = int (np.loadtxt('%s/n_sims'%(production_dir)))
dt = 0.005
dcd_freq = 10000

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

def sim_root_directory(run_id) :
    """
    Root directory of the production run
    """
    return '%s/%s'%(production_dir,run_id)

def sim_directory(run_id,phi,e,n) :
    """
    Returns the directory containing the simulation data
    """
    sim_root_dir = sim_root_directory(run_id)
    basename, simname = sim_name(phi,e,n)
    return '%s/%s/%s'%(sim_root_dir,basename,simname)

def load_sim(run_id,phi,e,n) :
    """
    Returns a loaded simulation corresponding to the given simulation
    parameters.
    """
    basename, simname = sim_name(phi,e,n)
    simdir = sim_directory(run_id,phi,e,n)
    xml = '%s/%s.xml'%(simdir,simname)
    dcd = '%s/%s.dcd'%(simdir,simname)
    gsd = '%s/%s.gsd'%(simdir,simname)
    if os.path.exists(xml) and os.path.exists(dcd) :
        sim = mbt.hoomdsim (xml,dcd)
    elif os.path.exists(gsd) :
        sim = mbt.hoomdsim(gsd)
    sim.phi = phi
    sim.e = e
    sim.n = n
    sim.basename = basename
    sim.simname = simname
    sim.simdir = simdir
    return sim

def data_directory(run_id) :
    """
    Returns the directory containing the analysis data of the simulations
    identified by 'run_id'
    """
    return '%s/data'%(sim_root_directory(run_id))

def figure_directory(run_id) :
    """
    Returns the directory containing the figures corresponding to the
    simulations identified by 'run_id'
    """
    return '%s/figures'%(sim_root_directory(run_id))
