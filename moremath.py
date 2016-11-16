import numpy as np

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
