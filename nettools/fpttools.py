import numpy as np

def jump_to (p):
    """
    Select a site according to the probability vector p
    """
    value = np.random.random ()
    return np.argmax (value-p<0.)

def FPT (P,startsite,endsite) :
    """
    Returns the first passage time of a single search for 'endsite', starting
    from 'startsite', on the row-normalized cumulative sum probability matrix
    P. Caution: no check is performed on the sanity of P.
    """
    site = startsite
    t = 0
    while site!=endsite :
        site = jump_to(P[site])
        t+=1
    return t

def FPT_distribution (P,startsite,endsite,
                      ntrials=None,bins=100) :
    """
    For the row-normalized cumulative sum probability matrix P, return the
    first passage time distribution for random walks starting at 'startsite' and
    ending at 'endsite'.

    Optional arguments:
        - ntrials: number of FPTs to extract (default, N*10, where N is the
        dimension of P
        - bins: number of bins to use for the histogram containing the
        distribution (default: 100)
    """
    if ntrials is None :
        # number of nodes in the network
        N = P.shape[0]
        ntrials = N*10
    fpt = np.zeros (ntrials)
    for i in range (ntrials) :
        fpt[i] = FPT (P,startsite,endsite)
    return np.histogram (fpt,bins=bins,normed=True)

def GFPT (Q,target,ntrials=None,bins=100) :
    """
    Given the adjacency matrix Q, compute the global first passage time
    distribution, that is, the first passage time distribution averaged over the
    starting sites, with a weight that corresponds to the stationary
    distribution.
    """
    N = Q.shape[0]
    P = np.cumsum (row_normalize_matrix (Q),axis=1)
    fpt_startsite = np.zeros ((bins,N))
    binsvals = np.zeros (N+1)
    for i in xrange(N) :
        if i == target : continue
        fpt_startsite[:,i],binsvals = FPT_distribution(P,startsite,target,ntrials,bins)
    W = np.sum (Q,axis=1)
    W /= np.sum(W)
    gfpt = np.mean (W*fpt_startsite,axis=1)
    return gfpt,binsvals,fpt_startsite
