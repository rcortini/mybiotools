import numpy as np
from pathos.multiprocessing import ProcessingPool as Pool
from scipy.special import gamma, jv
from scipy.linalg import eig
from .random_walks import jump_to, row_normalize_matrix
from mybiotools import error_message

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
    return np.histogram (fpt,bins=bins,density=True)[0]

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
    return gfpt

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
    Benichou2011, equation 3. User should supply the values of the rescaled
    times to compute, and the "nu" parameter, which is the ratio between the
    fractal dimension and the walk dimension.
    """
    if nu>=1 :
        # non-compact case
        return np.exp(-T)
    else :
        try :
            import besselzeros
        except ImportError :
            error_message ("GFPT_theory","Could not import besselzeros module")
            return np.zeros_like (T)
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

def GMFPT_theory (A,weighted=True) :
    """
    According to the theory of Lin et al., 2012, the global mean first passage
    time can be calculated by finding the eigenspectrum of the Laplacian matrix
    of the graph. This function calculates the GMFPT from their formula, for the
    graph described by the adjacency matrix A, to all sites. Optional parameter
    'weighted' allows for the choice of having the same quantity but weighted
    with the stationary distribution.
    """
    N = A.shape[0]
    d = np.sum(A,axis=1)
    E = np.sum(d)/2.
    L = np.diag(d) - A
    L_eigs = eig(L)
    sortidx = np.argsort(L_eigs[0])
    l = np.array([L_eigs[0][i].real for i in sortidx])
    v = np.array([L_eigs[1][:,i].real for i in sortidx])
    T = np.zeros(N)
    dv = np.dot (v,d)
    if not weighted :
        for j in xrange(N) :
            for i in range(1,N) :
                T[j] += 1.0/l[i] * (2*E*v[i,j]**2 - v[i,j]*dv[i])
        return float(N)/(N-1.0) * T
    else :
        for j in xrange(N) :
            for i in range(1,N) :
                dvi = v[i,j]*dv[i]
                T[j] += 1.0/l[i]*((2*E)**2*v[i,j]**2 - 2*v[i,j]*2*E*dvi - dvi**2)
        return T/(2*E)

def extend_adjacency_matrix (A0,p_void) :
    """
    This function takes the adjacency matrix 'A0' and adds a node to it. The
    node represents a state that is equally probably reachable from any other
    node, with probability 'p_void'.
    """
    N = A0.shape[0]
    A = np.zeros((N+1,N+1))
    A[:N,:N] = A0
    d_j = np.sum(A0,axis=1)
    lambda_j = p_void * d_j / (1-p_void)
    A[N,:-1] = lambda_j
    A[:-1,N] = lambda_j
    return A
