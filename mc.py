import numpy as np

def metropolis (beta,H_initial,H_final) :
    """
    Metropolis Monte Carlo basic step. User provides 'beta', the inverse
    temperature, and the two energies 'H_initial' and 'H_final'. Returns
    True or False depending on whether the move was accepted or rejected.
    """
    if H_final <= H_initial :
        return True
    else :
        if (np.random.rand() < np.exp (-beta*(H_final-H_initial))) :
            return True
        else :
            return False
