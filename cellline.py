import numpy as np
import pandas as pd
import os
from .parsers import parse_hic, parse_narrowpeak

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
        tracks = []
        for track in self._data :
            for key,val in conditions.iteritems() :
                try :
                    if track[key] != val :
                        break
                except KeyError :
                    break
                # the iteration gets here only if all condition keys exist, and
                # if the track satisfies all the desired conditions.
                tracks.append(track)
        return tracks
