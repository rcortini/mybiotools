import numpy as np
from scipy.linalg import eig

def jump_to (p):
    """
    Select a site according to the probability vector p
    """
    value = np.random.random ()
    return np.argmax (value-p<0.)

def row_normalize_matrix (M) :
    """
    From the matrix M, return the matrix M', defined as
    M'_ij = M_ij/sum_k(M_ik)
    """
    n = np.sum (M,axis=1)
    N = M.shape [0]
    Mnorm = M.copy ()
    for i in range (N) :
        if n[i]!=0. :
            Mnorm[i] /= n[i]
    return Mnorm

def random_walk (startsite,P,nsteps) :
    """
    Perform a random walk on the graph described by the cumulative sum
    row-normalized matrix P. The random walk starts at startsite, and lasts
    nsteps. Returns the sequence of the sites that were visited.
    """
    rw = np.zeros (nsteps,dtype=np.int32)
    value = np.random.random ()
    rw[0] = startsite
    site = startsite
    for i in xrange (1,nsteps) :
        site = jump_to(P[site])
        rw[i] = site
    return rw

def occupancy (rw,N) :
    """
    Given a random walk rw performed on a graph with N nodes, return the
    observed occupancy of each node.
    """
    c = np.zeros (N,dtype=np.int32)
    for step in rw :
        c[step] += 1
    return c

def occupancy_theory (A) :
    """
    Given the adjacency matrix A, return the equilibrium population associated
    to the graph represented by A.
    """
    P = row_normalize_matrix (A.astype('f'))
    # get left eigenspace of the row-normalized matrix
    eigs = eig(P,left=True)
    # get normalized eigenvector associated to largest eigenvalue, which is
    # supposed to be 1.
    imax = np.argmax (eigs[0])
    w = eigs[1][:,imax].real
    return w/w.sum()

def adjacency (rw,N) :
    """
    Given a random walk rw performed on a graph with N nodes, return the
    observed adjacency matrix.
    """
    A = np.zeros ((N,N),dtype=np.int32)
    step_prev = rw[0]
    for step in rw[1:] :
        A[step_prev,step] += 1
        step_prev = step
    return A
