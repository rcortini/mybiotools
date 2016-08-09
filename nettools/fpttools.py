import numpy as np
from pathos.multiprocessing import ProcessingPool as Pool

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

def FPT_distribution (P,startsite,endsite,bins,
                      ntrials=None) :
    """
    For the row-normalized cumulative sum probability matrix P, return the
    first passage time distribution for random walks starting at 'startsite' and
    ending at 'endsite'. Note that the bins of the distribution need to be
    computed beforehand, and passed to the function.

    Optional arguments:
        - ntrials: number of FPTs to extract (default, N*10, where N is the
        dimension of P
    """
    if ntrials is None :
        # number of nodes in the network
        N = P.shape[0]
        ntrials = N*10
    fpt = np.zeros (ntrials)
    for i in range (ntrials) :
        fpt[i] = FPT (P,startsite,endsite)
    return np.histogram (fpt,bins=bins,density=True)[0]*np.ediff1d(bins)

def GFPT (Q,target,bins,ntrials=None,nthreads=1) :
    """
    Given the adjacency matrix Q, compute the global first passage time
    distribution, that is, the first passage time distribution averaged over the
    starting sites, with a weight that corresponds to the stationary
    distribution. The bins of the distribution need to be supplied to the
    function.
    """
    N = Q.shape[0]
    P = np.cumsum (row_normalize_matrix (Q),axis=1)
    nbins = bins.shape[0]-1
    fpt_startsite = np.zeros ((nbins,N))
    if nthreads == 1 :
        for i in xrange(N) :
            if i == target : continue
            fpt_startsite[:,i] = FPT_distribution(P,i,target,bins,ntrials)
    else :
        pool = Pool (nthreads)
        def FPT_partial (i) :
            return FPT_distribution(P,i,target,bins,ntrials)
        fpt_map = pool.map (FPT_partial, range(N))
        for i in xrange (N) :
            fpt_startsite[:,i] = fpt_map [i]
    W = np.sum (Q,axis=1)
    W /= np.sum(W)
    gfpt = np.sum (W*fpt_startsite,axis=1)
    return gfpt,fpt_startsite

def MFPT (gfpt,bins) :
    """
    Given the global mean first passage time distribution, together with the
    bins for which it was calculated, return the mean. Note that we need this
    function because the np.histogram function returns a probability density
    which is not a mass function (i.e. its sum is not one), so that we need to
    evaluate the spacing between the bins.
    """
    x = np.array([0.5*(bins[i-1]+bins[i]) for i in range (1,len(bins))])
    return np.sum(x*gfpt*np.ediff1d(bins))

def GFPT_theory (T,nu) :
    """
    This function returns the theoretical GFPT distribution. Taken from
    Bénichou2011, equation 3. User should supply the values of the rescaled
    times to compute, and the "nu" parameter, which is the ratio between the
    fractal dimension and the walk dimension.
    """
    if nu>=1 :
        # non-compact case
        return np.exp(-T)
    else :
        # compact case
        nterms = 100
        A = 2.0*(1-nu**2)/nu
        ak = np.array([besselzeros.n(n,-nu) for n in xrange(1,nterms)])
        jnu = jv(nu,ak)
        j1_nu = jv(1.0-nu,ak)
        gt = np.zeros_like(T)
        for i,t in enumerate(T) :
            sum_terms = np.power(ak,1.0-2.0*nu) * jnu/j1_nu * np.exp (-ak**2/A * t)
            gt[i] = 2.0**(2.0*nu+1)/A * gamma(1.0+nu)/gamma(1.0-nu) * np.sum (sum_terms)
        return gt
