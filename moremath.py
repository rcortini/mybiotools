import numpy as np
from scipy import stats

def autocorrelation (x) :
    """
    Compute the autocorrelation of the signal, based on the properties of the
    power spectral density of the signal.
    """
    xp = x-np.mean(x)
    f = np.fft.fft(xp)
    p = np.array([np.real(v)**2+np.imag(v)**2 for v in f])
    pi = np.fft.ifft(p)
    return np.real(pi)[:x.size/2]/np.sum(xp**2)

def linear_fit (x,y) :
    """
    Fit (x,y) to a linear function, using unweighted least-square minimization
    procedure from numpy linear algebra. Returns the coefficients a and b of y =
    a + b x
    """
    A = np.vstack([x, np.ones(x.size)]).T
    return np.linalg.lstsq(A,y)[0]

def linear_regression (x,y,prob) :
    """
    Fit (x,y) to a linear function, using maximum likelihood estimation of the
    confidence intervals on the coefficients, given the user-supplied
    probability *prob*
    """
    n = len(x)
    xy = x*y
    xx = x*x
    # estimates
    xmean = x.mean()
    ymean = y.mean()
    xxmean = xx.mean()
    xymean = xy.mean()
    b1 = (xymean-xmean*ymean) / (xxmean-xmean**2)
    b0 = ymean-b1*xmean
    s2 = 1./n * sum([(y[i] - b0 - b1 * x[i])**2 for i in xrange(n)])
    #confidence intervals
    alpha = 1 - prob
    c1 = stats.chi2.ppf(alpha/2.,n-2)
    c2 = stats.chi2.ppf(1-alpha/2.,n-2)
    # get results and return
    c = -1 * stats.t.ppf(alpha/2.,n-2)
    bb1 = c * (s2 / ((n-2) * (xxmean - (xmean)**2)))**.5
    bb0 = c * ((s2 / (n-2)) * (1 + (xmean)**2 / (xxmean - xmean**2)))**.5
    return b0,b1,bb0,bb1

def wlinear_fit (x,y,w) :
    """
    Fit (x,y,w) to a linear function, using exact formulae for weighted linear
    regression. This code was translated from the GNU Scientific Library (GSL),
    it is an exact copy of the function gsl_fit_wlinear.
    """
    # compute the weighted means and weighted deviations from the means
    # wm denotes a "weighted mean", wm(f) = (sum_i w_i f_i) / (sum_i w_i) */
    mask = [w>0]
    W = np.sum(w[mask])
    wm_x = np.average(x,weights=w)
    wm_y = np.average(y,weights=w)
    dx = x-wm_x
    dy = y-wm_y
    wm_dx2 = np.average(dx**2,weights=w)
    wm_dxdy = np.average(dx*dy,weights=w)
    # In terms of y = a + b x
    d2 = 0.0
    b = wm_dxdy / wm_dx2
    a = wm_y - wm_x*b
    cov_00 = (1.0/W) * (1.0 + wm_x**2/wm_dx2)
    cov_11 = 1.0 / (W*wm_dx2)
    cov_01 = -wm_x / (W*wm_dx2)
    # Compute chi^2 = \sum w_i (y_i - (a + b * x_i))^2
    chi2 = np.sum (w * (y-(a+b*x))**2)
    return a,b,cov_00,cov_11,cov_01,chi2

def KL_divergence (P,Q) :
    return stats.entropy (P,Q)

def LJ_potential (r,sigma,epsilon) :
    """Lennard-Jones potential"""
    r6 = (sigma/r)**6
    r12 = r6*r6
    return 4.0 * epsilon * (r12 - r6)

def new_average(N,old_average,new_datapoint) :
    """
    Given N observations that have an average 'old_average', returns the new
    average given a new observation at 'new_datapoint'
    """
    if N==0 :
        return new_datapoint
    else :
        return 1.0/(N+1) * (old_average * N + new_datapoint)

def fit_powerlaw(x,y) :
    """
    Given a set of observation y = f(x), fit the data to a power law, and return
    the amplitude and exponent of the fit.
    """
    mask = np.logical_and(x>0,y>0)
    xfit = np.log(x[mask])
    yfit = np.log(y[mask])
    res = mbt.linear_fit(xfit,yfit)
    return np.exp(res[1]),res[0]
