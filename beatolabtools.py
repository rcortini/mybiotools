import numpy as np
import pandas as pd
import os
from .utils import warn_message

def load_beato_metadata (
    metadata_file='/home/rcortini/work/data/beato_lab_metadata.xlsx') :
    return pd.read_excel(metadata_file)

def cell_load_tracks (cell,tracks,xavi_datadir='/mnt/xavi/data') :
    for i,track in tracks.iterrows() :
        metadata = track.to_dict()
        sample_id = track['SAMPLE_ID']
        metadata['type'] = sample_id.split('_')[-1]
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
