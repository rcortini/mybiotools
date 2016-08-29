import numpy as np
import os
import pysam
import pandas as pd

def parse_sam (samfilename, mapq_threshold=20) :
    """
    Parses a .sam file and produces an numpy array with the
    'chr', 'start', 'end', 'name', 'mapq', 'strand' fields.
    For the moment can handle only very simple cases of flags.
    """
    sam_dtype = np.dtype([
                    ('chr','S254'),
                    ('start',np.int64),
                    ('end',np.int64),
                    ('name','S2048'),
                    ('mapq','i8'),
                    ('strand','S2')
                   ])
    samfile = pysam.AlignmentFile(samfilename, "r")
    reads = []
    for read in samfile.fetch(until_eof=True) :
        if read.mapq > mapq_threshold :
            c = samfile.get_reference_name(read.rname)
            mapq = read.mapq
            start = read.positions[0]
            end = read.positions[-1]
            if read.flag == 0 :
                strand = '+'
            elif read.flag == 16 :
                strand = '-'
            else :
                raise ValueError
            name = read.qname
            r = (c,start,end,name,mapq,strand)
            reads.append (r)
    return np.array (reads, dtype=sam_dtype)

def parse_broadpeak (fname) :
    """
    Parses a broadpeak file type, and produces a numpy array with the fields of
    'chr', 'start', 'end', 'name', 'score', 'strand', 'val', 'p', 'q'.
    """
    broadpeak_dtype = np.dtype([
                                ('chr','S10'),
                                ('start',np.int64),
                                ('end',np.int64),
                                ('name','S10'),
                                ('score','i8'),
                                ('strand','S10'),
                                ('val','f'),
                                ('p','f'),
                                ('q','f'),
                               ])
    return np.genfromtxt (fname,dtype=broadpeak_dtype)

def parse_narrowpeak (fname) :
    """
    Parses a broadpeak file type, and produces a numpy array with the fields of
    'chr', 'start', 'end', 'name', 'score', 'strand', 'val', 'p', 'q', 'peak'.
    """
    narrowpeak_dtype = np.dtype([
                                ('chr','S10'),
                                ('start',np.int64),
                                ('end',np.int64),
                                ('name','S10'),
                                ('score','i8'),
                                ('strand','S10'),
                                ('val','f'),
                                ('p','f'),
                                ('q','f'),
                                ('peak',np.int64),
                               ])
    return np.genfromtxt (fname,dtype=narrowpeak_dtype)

def res_string (res) :
    """
    Converts the integer resolution 'res' into a string of type
    5kb, 10kb, 1mb...
    """
    m = res/1000
    if m>=1000 :
        s = 'm'
        m = m/1000
    else :
        s = 'k'
    return '%d%sb'%(m,s)

def load_hic_Rao (hic_res,name,normed=True) :
    """
    Load the Hi-C matrices from the experiments of Rao et al, 2014, for the
    lymphoblastoid cell line GM12878. User must specify the resolution, the name
    of the chromosome, and whether or not to apply the normalization suggested
    in the paper
    """
    Rao_datadir = '/mnt/ant-login/rcortini/work/data/GM12878_replicate/'
    hic_res_string = res_string (hic_res)
    # build the directory name that contains the data that we want to analyze
    d = '%s/%s_resolution_intrachromosomal/chr%s/MAPQGE30'%(Rao_datadir,
                                                            hic_res_string,
                                                            name)
    fname = '%s/chr%s_%s.RAWobserved'%(d,name,hic_res_string)
    normname = '%s/chr%s_%s.KRnorm'%(d,name,hic_res_string)
    if not os.path.exists (fname) or not os.path.exists (fname) :
        raise ValueError('Data for chromosome %s at resolution %d does not exist'
                         %(name,hic_res))
    norm = np.loadtxt (normname)
    N = norm.shape[0]
    H = np.zeros ((N,N))
    with open(fname,'r') as f :
        for line in f :
            c = line.strip('\n').split()
            i = int(c[0])/hic_res
            j = int(c[1])/hic_res
            if np.isnan (norm[i]) or np.isnan (norm[j]) :
                M = float(c[2])
            else :
                if normed :
                    M = float(c[2])/(norm[i]*norm[j])
                else :
                    M = float(c[2])
            H[i,j] = H[j,i] = M
    return H

def chromosome_size (name) :
    """
    Returns the size in base pairs of a given chromosome, according to the
    genome version h19
    """
    fname = '/mnt/ant-login/rcortini/work/data/human/genome_size.dat'
    with open (fname,'r') as f :
        for line in f :
            curatedline = line.strip('\n').split()
            if curatedline[0] == 'chr%s'%(name) :
                return int (curatedline[1])
    return 0

def parse_hic (name) :
    """
    Parses a Hi-C file. Supposed to be a xls Excel file. (TODO: extend)
    """
    return pd.read_excel(f).as_matrix()
