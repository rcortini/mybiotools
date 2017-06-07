import numpy as np

def restore_images (images_file,system) :
    # get information on the image indices of each particle from the images file
    images = np.loadtxt (images_file)
    snapshot = system.take_snapshot()
    for i,im in enumerate (images) :
        snapshot.particles.image[i] = list (im)
    system.restore_snapshot (snapshot)
