import numpy as np
import os
import pysam
import pandas as pd
import pysam

def parse_sam (samfilename, mapq_threshold=20) :
    """
    Parses a .sam file and produces an numpy array with the
    'chr', 'start', 'end', 'name', 'mapq', 'strand' fields.
    For the moment can handle only very simple cases of flags.
    """
    sam_dtype = np.dtype([
                    ('chr','S256'),
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
                                ('name','S2048'),
                                ('score','i8'),
                                ('strand','S2'),
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
                                ('name','S256'),
                                ('score','i8'),
                                ('strand','S2'),
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

def chromosome_size (name) :
    """
    Returns the size in base pairs of a given chromosome, according to the
    genome version h19
    """
    fname = '%s/work/data/human/genome_size.dat'%(os.getenv('HOME'))
    with open (fname,'r') as f :
        for line in f :
            curatedline = line.strip('\n').split()
            if curatedline[0] == '%s'%(name) :
                return int (curatedline[1])
    return 0

def load_hic_Rao (hic_res,name,normed=True,
                  Rao_datadir = '/mnt/ant-login/rcortini/work/data/GM12878_replicate/') :
    """
    Load the Hi-C matrices from the experiments of Rao et al, 2014, for the
    lymphoblastoid cell line GM12878. User must specify the resolution, the name
    of the chromosome, and whether or not to apply the normalization suggested
    in the paper
    """
    hic_res_string = res_string (hic_res)
    # build the directory name that contains the data that we want to analyze
    d = '%s/%s_resolution_intrachromosomal/%s/MAPQGE30'%(Rao_datadir,
                                                            hic_res_string,
                                                            name)
    fname = '%s/%s_%s.RAWobserved'%(d,name,hic_res_string)
    normname = '%s/%s_%s.KRnorm'%(d,name,hic_res_string)
    if not os.path.exists (fname) or not os.path.exists (fname) :
        raise ValueError('Data for chromosome %s at resolution %d does not exist'
                         %(name,hic_res))
    if normed :
        norm = np.loadtxt (normname)
    N = chromosome_size(name)/hic_res + 1
    H = np.zeros ((N,N))
    with open(fname,'r') as f :
        for line in f :
            c = line.strip('\n').split()
            i = int(c[0])/hic_res
            j = int(c[1])/hic_res
            if normed :
                if np.isnan (norm[i]) or np.isnan (norm[j]) :
                    M = float(c[2])
                else :
                    M = float(c[2])/(norm[i]*norm[j])
            else :
                M = float(c[2])
            H[i,j] = H[j,i] = M
    return H

def parse_hic (name) :
    """
    Parses a Hi-C file.
    """
    if name.endswith ('.xls') or name.endswith('.xlsx') :
        return pd.read_excel(name).as_matrix()
    elif name.endswith ('.tsv.gz') or name.endswith ('.tsv') :
        if 'raw' in name :
            hic_dtype=np.dtype({'names':['chr','i','j','val'],
                                'formats':['S12',np.int64,np.int64,np.int64]})
        else :
            hic_dtype=np.dtype({'names':['chr','i','j','val'],
                                'formats':['S12',np.int64,np.int64,np.float64]})
        return np.genfromtxt(name,dtype=hic_dtype)

def parse_kallisto_rnaseq (name) :
    """
    Parses an annotated tsv file that was the output of Kallisto.
    """
    kallisto_dtype = np.dtype([
            ('target_id','S2048'),
            ('length',np.int64),
            ('eff_length',np.float64),
            ('est_counts',np.float64),
            ('tpm',np.float64)
        ])
    return np.genfromtxt(name,dtype=kallisto_dtype,skip_header=1)

def bam_to_matrix(bam,chromosome,start,end,resolution,
                 flag=1807) :
    """
    Using pysam (that is, samtools), take a bam file and extract the reads of a
    hi-c file corresponding to the given region (chromosome, start, end) and
    having a flag that is less than the given one.
    Returns a matrix of counts at the given resolution.
    """
    # the 'b' flag here indicates that we are dealing with a bam file
    with pysam.AlignmentFile(bam,'rb') as samfile :
        # init the matrix
        N = (end-start)/resolution
        H = np.zeros((N,N),dtype=np.int32)
        # init the iterator
        samiter = samfile.fetch(chromosome,start,end)
        # iterate
        for read in samiter :
            if read.flag <= flag :
                i = read.pos
                j = read.mpos
                # we don't check 'i' here: this is already guaranteed by the
                # invocation to the "fetch" in samtools
                if j>=start and j<=end :
                    i_mat = (i-start)/resolution
                    j_mat = (j-start)/resolution
                    H[i_mat,j_mat] += 1
                    H[j_mat,i_mat] += 1
    return H

def parse_simple_bed (fname) :
    """
    A parser for a BED file that contains only the indication of chromosome,
    region start and region end.
    """
    simple_bed_dtype = np.dtype([
                                ('chr','S10'),
                                ('start',np.int64),
                                ('end',np.int64)
                               ])
    return np.genfromtxt (fname,dtype=simple_bed_dtype)


def counts_to_hic (counts,start,end,resolution) :
    """
    Returns a complete filled matrix given the 'counts' array, by taking
    for granted that the counts correspond to a given chromosome.
    """
    N = (end-start)/resolution + 1
    H = np.zeros((N,N),dtype=counts['val'].dtype)
    for h in counts :
        i = (h['i']-start)/resolution
        j = (h['j']-start)/resolution
        H[i,j] = H[j,i] = h['val']
    return H
