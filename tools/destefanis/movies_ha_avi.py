def movies_avi(path,name):
    import movies as pf
    from astropy.io import fits
    import os
    import numpy as np
    list = os.listdir(path)
    list.sort() 
    list2 = []
    path = path + "/" 
    for file in list: list2.append(path + file)
    img1 = fits.open(list2[0])
    tam = img1[1].data.shape
    cubo = np.zeros((len(list),tam[0],tam[1]))
    for i,file in enumerate(list2) : cubo[i,:,:] = fits.open(file)[1].data
    pf.save_movie(name,cubo,fps=5)
    return
