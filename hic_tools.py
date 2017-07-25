import numpy as np

def counts_hic (H,threshold) :
    """
    Given the Hi-C matrix H, calculate the approximate number of contacts that
    the genomic sites make with the others. First, remove the rows that contain
    less than 'threshold' counts, then do the exponential of the entropy.
    """
    counts = []
    mask = H.sum(axis=1)>=threshold
    for row in H[mask] :
        cleanrow = row[mask]
        counts.append(np.exp(scipy.stats.entropy(cleanrow)))
    return np.array(counts),mask


