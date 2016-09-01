import numpy as np
import pandas as pd
import os
from .utils import warn_message
from .parsers import res_string

def load_beato_metadata (
    metadata_file='/home/rcortini/work/data/beato_lab_metadata.xlsx') :
    return pd.read_excel(metadata_file)

def load_hic_metadata (
    hic_metadata_file='/home/rcortini/work/data/beato_lab_hic_metadata.xlsx') :
    return pd.read_excel(hic_metadata_file)

def cell_load_tracks (cell,tracks,resolution=10000,xavi_datadir='/mnt/xavi/data') :
    for i,track in tracks.iterrows() :
        metadata = track.to_dict()
        sample_id = track['SAMPLE_ID']
        metadata['type'] = sample_id.split('_')[-1]
        metadata['resolution'] = resolution
        # build the directory name where the files are
        d = "%s/%s/samples/%s/peaks"%(xavi_datadir,
                                      metadata['type'],
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
        if fin is not None :
            metadata['fname'] = fin
            cell.load_track(metadata)
        else :
            warn_message('cell_load_tracks','Data not found for %s'%sample_id)

def cell_load_hic (cell,tracks,resolution,
                   datatype_string='raw',
                   xavi_hic_datadir='/mnt/hic') :
    res = res_string (resolution)
    for i,track in tracks.iterrows() :
        metadata = track.to_dict()
        sample_id = track['SAMPLE_ID']
        metadata['type'] = 'hic'
        metadata['resolution'] = resolution
        d = "%s/samples/%s/downstream"%(xavi_hic_datadir,sample_id)
        metadata['fname'] = None
        if os.path.exists(d) :
            for root,subdirs,files in os.walk(d) :
                for f in files :
                    if datatype_string in f and f.endswith ('%s.tsv.gz'%res) :
                        fin = '%s/%s'%(root,f)
                        metadata['fname'] = fin
        if metadata['fname'] is not None :
            cell.load_track(metadata)
        else :
            warn_message('cell_load_hic','Data not found for %s'%sample_id)
