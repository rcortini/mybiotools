import numpy as np

def make_zerone_output_index(a,chromosome_list=None) :
    # if the chromosome names are not given, get them
    if chromosome_list is None :
        chromosome_list = np.unique(a['chr'])
    # init the dictionary and init the iteration
    c_idx = {}
    c_start = 0
    prev_c = a[0]['chr']
    # we don't sort the input array: Zerone already does this by default
    for i,b in enumerate(a) :
        this_c = b['chr']
        if this_c != prev_c :
            c_end = i-1
            c_idx[prev_c] = (c_start,c_end)
            c_start = i
            prev_c = this_c
    # the last chromosome needs to be manually added
    c_idx[this_c] = (c_start,i)
    return c_idx

def parse_zerone_output(fname,chromosome_list=None) :
    """
    Parses a Zerone output and returns a numpy array. The values of the numpy array
    are: chromosome, start, end, enrichment, read_1, read_2, ..., read_n, p.
    The number of `read_i` columns depends on the invocation of Zerone and cannot
    be known beforehand.
    """
    # first, we start by reading the first non-comment line in the Zerone file, to
    # determine the number of `read` columns in the file
    comment = True
    with open(fname,'r') as f :
        for line in f :
            if not line.startswith('#') :
                break
    n_readcols = len(line.split())-6
    zerone_dtype = [('chr','S256'),
                    ('start',np.int64),
                    ('end',np.int64),
                    ('enrichment',np.int32),
                    ('control',np.int64)]
    for i in range(n_readcols) :
        zerone_dtype.append(('read_%d'%(i),np.int64))
    zerone_dtype.append(('p',float))
    # now we parse the file using the `genfromtxt` function from numpy
    a = np.genfromtxt(fname,dtype=np.dtype(zerone_dtype))
    # next, we exclude the values of the array that pertain to chromosomes that are not
    # included in the chromosome list that was passed by the user (if any)
    if chromosome_list is not None :
        a = np.array([s for s in a if s['chr'] in chromosome_list])
    # now pass the array to the index maker, and return the array along with the index
    c_idx = make_zerone_output_index(a,chromosome_list)
    return a,c_idx

def find_zerone_peak(a,c_idx,peak,bin_size=300) :
    """
    Returns the values of the `a` array corresponding to the genomic coordinates
    of the `peak`. Uses the `c_idx` dictionary to rapidly calculate which are the indices
    of the `a` array that correspond to the peak
    """
    c_start,c_end = c_idx[peak['chr']]
    peak_idx_start = peak['start']//bin_size
    peak_idx_end = peak['end']//bin_size
    if peak_idx_start == peak_idx_end :
        return [a[c_start+peak_idx_start]]
    else :
        return a[c_start+peak_idx_start:c_start+peak_idx_end+1]
