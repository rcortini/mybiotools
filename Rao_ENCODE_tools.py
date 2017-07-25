from .parsers import res_string

def datadir (chromosome_name,resolution,threshold,global_datadir) :
    return "%s/%s/res_%s/threshold_%d"%(global_datadir,
                                        chromosome_name,
                                        res_string(resolution),
                                        threshold)
