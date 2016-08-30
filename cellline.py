import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from .parsers import parse_hic, parse_narrowpeak
from .vistools import line_plot

def select_tracks (tracks,conditions) :
    selected_tracks = []
    if not conditions :
        return tracks
    for track in tracks :
        for key,val in conditions.iteritems() :
            try :
                if track[key] != val :
                    break
            except KeyError :
                break
            # the iteration gets here only if all condition keys exist, and
            # if the track satisfies all the desired conditions.
            selected_tracks.append(track)
    return selected_tracks

class CellLine :
    def __init__ (self,name) :
        self.name = name
        self._data = []
        # these are the parsers that are currently available
        self.parsers = {'hic'         :           parse_hic,
                        'chipseq'     :           parse_narrowpeak}
    def load_track (self,metadata) :
        """
        Loads a single track of data to the CellLine's data track list. To do
        so, it uses the user-supplied 'metadata' dictionary to evince two
        things:
            - the file name, as given through the 'fname' key
            - the data type, as given through the 'type' key
        The function then attempts to read the fname using the corresponding
        parser. If successful, then a track is added to the records, which can
        then be accessed by the 'get_data' function. All the remaining keys and
        values in the metadata dictionary are copied to the record.
        """ 
        # sanity check on the metadata: has fname and type keys
        try :
            fname = metadata['fname']
            datatype = metadata['type']
        except KeyError as e:
            failedkey = e.args[0]
            raise KeyError ("User must supply the '%s' key to the metadata"%
                            failedkey)
        # sanity check on the metadata: parser exist for given datatype
        if datatype not in self.parsers.keys() :
            raise KeyError ("No parser exists for data type %s"%datatype)
        # sanity check on the metadata: fname exists
        if not os.path.exists (fname) :
            raise IOError ("%s does not exist"%fname)
        # everything's fine: load the data and append the record to the records.
        datadict = metadata.copy()
        datadict['data'] = self.parsers[datatype](fname)
        self._data.append (datadict)
    def get_data (self,conditions) :
        """
        Returns a list of tracks loaded into the CellLine by specifying one or
        more 'conditions'. The conditions are to be specified by a dictionary in
        which the keys correspond to metadata fields, and the values to the
        desired value of that field.
        """
        return select_tracks (self.data,conditions)

    @property
    def data (self) :
        return self._data
    @data.setter
    def data(self,track) :
        # TODO: sanity checks
        self._data.append (track)

class Region :
    def __init__ (self,chromosome,start,end,resolution) :
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.resolution = resolution
        self._data = []
        self.xvals = np.arange (start,end,resolution)
        self.N = self.xvals.shape[0]
    def set_chipseq_track (self,track,q_threshold=None) :
        """
        This function passes the 'track' dictionary to the Region class and
        coarse-grains its values to the Region's xvals. The 'track' must have
        the structure as given by the CellLine class. The optional argument
        'q_threshold' can be passed to select the peaks that meet a certain
        quality threshold.
        """
        mytrack = {}
        for key,val in track.iteritems () :
            if key != 'data' :
                mytrack[key] = val
        mytrack['track'] = np.zeros_like(self.xvals,dtype=np.float64)
        for peak in track['data'] :
            if (peak['chr'] == self.chromosome) and\
               (peak['start'] > self.start)     and\
               (peak['start'] < self.end) :
                start = (peak['start']-self.start)/self.resolution
                if q_threshold is not None :
                    if peak['q'] < q_threshold :
                        continue
                mytrack['track'][start] = peak['val']
        self._data.append (mytrack)
    def set_hic (self,hic) :
        # sanity check on the consistency of the provided Hi-C data with the
        # data within the Region
        if np.sum(hic['data']['start']%self.resolution) != 0 and\
           np.sum(hic['data']['end']%self.resolution) != 0 :
            raise ValueError ("Data resolution does not match Region's resolution")
        mytrack = {}
        for key,val in hic.iteritems () :
            if key != 'data' :
                mytrack[key] = val
        # now get all the values that correspond to the Region's chromosome and
        # extension
        mask = np.logical_and (hic['data']['chr']==self.chromosome,
                               np.logical_and(hic['data']['start']>self.start,
                                              hic['data']['end']<self.end))
        rawH = hic['data'][mask]
        # set the values of the matrix
        H = np.zeros((self.N,self.N),dtype=rawH['val'].dtype)
        for h in rawH :
            i = (h['start']-self.start)/self.resolution
            j = (h['end']-self.start)/self.resolution
            H[i,j] = H[j,i] = h['val']
        # and finally update the Region's data records
        mytrack['track'] = H
        self._data.append (mytrack)
    def get_tracks (self,conditions) :
        return select_tracks (self._data,conditions)
    def plot_tracks (self,conditions) :
        selected_tracks = select_tracks (self._data,conditions)
        n = len (selected_tracks)
        fig, axes = plt.subplots (n,1,figsize=(10,n*3))
        for i,track in enumerate(selected_tracks) :
            if i==n-1 :
                line_plot (axes[i],self.xvals,track['track'],show_xaxis=True)
            else :
                line_plot (axes[i],self.xvals,track['track'])
        return fig