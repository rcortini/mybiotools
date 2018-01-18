import numpy as np
import re

def restore_images (images_file,system) :
    # get information on the image indices of each particle from the images file
    images = np.loadtxt (images_file)
    snapshot = system.take_snapshot()
    for i,im in enumerate (images) :
        snapshot.particles.image[i] = list (im)
    system.restore_snapshot (snapshot)

def getpar(fin,parname) :
    with open(fin,'r') as f :
        for line in f :
            if re.search('%s = '%(parname),line) :
                break
    return line.split('=')[1].split('#')[0].replace(' ','')
