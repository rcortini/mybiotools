import numpy as np
import pandas as pd
import os
from .utils import log_message, warn_message
from .parsers import res_string, parse_kallisto_rnaseq

def load_beato_metadata (
    metadata_file='/home/rcortini/work/data/beato_lab_metadata.xlsx') :
    return pd.read_excel(metadata_file)

def load_hic_metadata (
    hic_metadata_file='/home/rcortini/work/data/beato_lab_hic_metadata.xlsx') :
    return pd.read_excel(hic_metadata_file)

def track_location (typ,sample_id,xavi_datadir='/mnt/xavi/data') :
    # build the directory name where the files are
    d = "%s/%s/samples/%s/peaks"%(xavi_datadir,
                                  typ,
                                  sample_id)
    # select all files that end with ".narrowPeak" in the directory, and
    # then prefer to read the one that is in the directory that has
    # "with_control"
    peakfiles = []
    for root,sub,files in os.walk(d) :
        for f in files :
            if f.endswith (".narrowPeak") :
                peakfiles.append('%s/%s'%(root,f))
    fin = None
    for peakfile in peakfiles :
        if 'with_control' in peakfile :
            fin = peakfile
            break
        else :
            fin = peakfile
    if fin is None :
        warn_message('cell_load_tracks','Data not found for %s'%sample_id)
    return fin

def cell_load_tracks (cell,tracks,resolution=10000,xavi_datadir='/mnt/xavi/data') :
    for i,track in tracks.iterrows() :
        metadata = track.to_dict()
        sample_id = track['SAMPLE_ID']
        metadata['type'] = sample_id.split('_')[-1]
        metadata['resolution'] = resolution
        fin = track_location(metadata['type'],sample_id)
        if fin is not None :
            metadata['fname'] = fin
            cell.load_track(metadata)

def hic_location(sample_id,resolution,datatype_string='raw',xavi_hic_datadir='/mnt/hic') :
    d = "%s/samples/%s/downstream"%(xavi_hic_datadir,sample_id)
    fname = None
    res = res_string (resolution)
    if os.path.exists(d) :
        for root,subdirs,files in os.walk(d) :
            for f in files :
                if datatype_string in f and f.endswith ('%s.tsv.gz'%res) :
                    fname = '%s/%s'%(root,f)
    if fname is None :
        warn_message('cell_load_hic','Data not found for %s'%sample_id)
    return fname

def cell_load_hic (cell,tracks,resolution,
                   datatype_string='raw',
                   xavi_hic_datadir='/mnt/hic') :
    for i,track in tracks.iterrows() :
        metadata = track.to_dict()
        sample_id = track['SAMPLE_ID']
        metadata['type'] = 'hic'
        metadata['resolution'] = resolution
        metadata['fname'] = hic_location (sample_id,resolution)
        if metadata['fname'] is not None :
            log_message('cell_load_hic','Loading %s'%(sample_id))
            cell.load_track(metadata)

def load_rnaseq (sample_id,xavi_datadir='/mnt/xavi/data') :
    # build directory name
    rnaseq_datadir = '%s/rnaseq/samples'%xavi_datadir
    this_datadir = '%s/%s/quantifications/kallisto/'%(rnaseq_datadir,sample_id)
    # search for our file in the directory
    fname = None
    for root,sub,files in os.walk(this_datadir) :
        for f in files :
            if f=='abundance.tsv' :
                fname = '%s/abundance.tsv'%root
    # check that the file was found
    if fname is not None :
        ref_genome = fname.replace(this_datadir,'').replace('/paired_end/abundance.tsv','')
    return parse_kallisto_rnaseq(fname), ref_genome

def hic_bam_location(sample_id,hic_bam_datadir='/mnt/hic_bam') :
    """
    Returns the location of the BAM file of the corresponding 'sample_id'
    """
    d = "%s/samples/%s"%(hic_bam_datadir,sample_id)
    fname = None
    if os.path.exists(d) :
        for root,subdirs,files in os.walk(d) :
            for f in files :
                if f.endswith('.bam') :
                    fname = '%s/%s'%(root,f)
    if fname is None :
        mbt.warn_message('hic_bam_location','Data not found for %s'%sample_id)
    return fname
