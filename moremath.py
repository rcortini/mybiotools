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
