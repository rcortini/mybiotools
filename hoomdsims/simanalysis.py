import numpy as np

def hic_chipseq_r2 (hic, chipseq) :
    """
    Calculate the Pearson correlation coefficient between the row sum of the
    given Hi-C matrix and the given ChIP-seq profile.
    """
    hic_rowsum = np.sum(hic,axis=1)/float(np.sum(hic))
    return np.corrcoef(hic_rowsum,chipseq)[0,1]**2

def ps (H) :
    """
    Calculate the normalized probability of contact between a monomer and all
    others as a function of the linear distance s.
    """
    p = np.array ([np.mean (np.diagonal (H, offset=k))
                   for k in range (H.shape[0])])
    return p/np.sum(p)
