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
        is_valid = True
        for key,val in conditions.iteritems() :
            try :
                if track[key] != val :
                    is_valid = False
                    break
            except KeyError :
                is_valid = False
                break
        if is_valid :
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
        self._data = []
    def set_chipseq_track (self,track,q_threshold=None) :
        """
        This function passes the 'track' dictionary to the Region class and
        coarse-grains its values to the xvals desired. The 'track' must have
        the structure as given by the CellLine class, plus an additional
        parameter 'resolution', which is the coarse-graining parameter. The
        optional argument 'q_threshold' can be passed to select the peaks that
        meet a certain quality threshold.
        """
        mytrack = {}
        for key,val in track.iteritems () :
            if key != 'data' :
                mytrack[key] = val
        resolution = mytrack['resolution']
        mytrack['xvals'] = np.arange (self.start,self.end,resolution)
        mytrack['track'] = np.zeros_like (mytrack['xvals'],dtype=np.float64)
        for peak in track['data'] :
            if (peak['chr'] == self.chromosome) and\
               (peak['start'] > self.start)     and\
               (peak['start'] < self.end) :
                start = (peak['start']-self.start)/resolution
                if q_threshold is not None :
                    if peak['q'] < q_threshold :
                        continue
                mytrack['track'][start] = peak['val']
        self._data.append (mytrack)
    def set_hic (self,hic) :
        mytrack = {}
        for key,val in hic.iteritems () :
            if key != 'data' :
                mytrack[key] = val
        # now get all the values that correspond to the Region's chromosome and
        # extension
        mask = np.logical_and (hic['data']['chr']==self.chromosome,
                               np.logical_and(hic['data']['start']>=self.start,
                                              hic['data']['end']<self.end))
        rawH = hic['data'][mask]
        # set the values of the matrix
        N = (self.end-self.start)/mytrack['resolution']
        H = np.zeros((N,N),dtype=rawH['val'].dtype)
        for h in rawH :
            i = (h['start']-self.start)/mytrack['resolution']
            j = (h['end']-self.start)/mytrack['resolution']
            H[i,j] = H[j,i] = h['val']
        # and finally update the Region's data records
        mytrack['track'] = H
        self._data.append (mytrack)
    def set_data (self,cell,conditions,q_threshold=None) :
        """
        This is a convenient way to interface a CellLine instance with a Region
        instance. You pass the CellLine object 'cell' to this method, along with
        the given conditions, and the data tracks are added to the Regions's
        records.
        """
        tracks = cell.get_data (conditions)
        for track in tracks :
            if track['type'] == 'chipseq' :
                self.set_chipseq_track(track,q_threshold=q_threshold)
            elif track['type'] == 'hic' :
                self.set_hic(track)
            else :
                raise ValueError ("Unsupported data type %s"%track['type'])
    def get_data (self,conditions) :
        return select_tracks (self._data,conditions)
    def plot_tracks (self,conditions) :
        selected_tracks = select_tracks (self._data,conditions)
        n = len (selected_tracks)
        fig, axes = plt.subplots (n,1,figsize=(10,n*3))
        for i,track in enumerate(selected_tracks) :
            if i==n-1 :
                line_plot (axes[i],track['xvals'],track['track'],show_xaxis=True)
            else :
                line_plot (axes[i],track['xvals'],track['track'])
        return fig
