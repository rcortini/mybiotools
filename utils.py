from __future__ import print_function
import time, sys, errno, os
import numpy as np

def time_string () :
    return time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime ())

def error_message (program_name, message) :
    full_message = "%s %s: ERROR: %s"%(time_string (), program_name, message)
    print (full_message, file=sys.stderr)

def log_message (program_name, message) :
    full_message = "%s %s: INFO: %s"%(time_string (), program_name, message)
    print (full_message)

def warn_message (program_name, message) :
    full_message = "%s %s: WARNING: %s"%(time_string (), program_name, message)
    print (full_message)

def consecutive_true(data):
    """
    A function that computes the number of consecutive 'True' values in the
    array 'data'. Taken from http://stackoverflow.com/a/24343375/2312821
    """
    return np.diff(np.where(np.concatenate(([data[0]],
                                     data[:-1] != data[1:],
                                     [True])))[0])[::2]

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def map_seq_to_int(seq) :
    """
    Map a sequence of DNA to an integer
    """
    base = 4
    n = 0
    base_map = {'A':0,'T':1,'C':2,'G':3}
    for i,c in enumerate(seq) :
        n += base_map[c] * base**i
    return n

def map_int_to_seq(num) :
    """
    Map an integer to a DNA sequence
    """
    i = 0
    n = 0
    base_map = {0:'A',1:'T',2:'C',3:'G'}
    seq = ''
    while n!=num :
        f = base**(i+1)
        digit = int( (num%base**(i+1)) / base**i )
        n += digit*base**i
        i += 1
        seq += base_map[digit]
    return seq
