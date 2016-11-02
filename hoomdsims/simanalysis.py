import numpy as np

def hic_chipseq_r2 (hic, chipseq) :
    """
    Calculate the Pearson correlation coefficient between the row sum of the
    given Hi-C matrix and the given ChIP-seq profile.
    """
    hic_rowsum = np.sum(H,axis=1)/float(np.sum(H))
    return np.corrcoef(hic_rowsum,chipseq)[0,1]**2
