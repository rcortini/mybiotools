import numpy as np
import pysam

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
