from plot_hmi import save_img

infile="/gehme/data/sdo/hmi/20101111/hmi.m_45s.2010.11.11_06_21_45_TAI.magnetogram.fits"
opath="/gehme/scratch/destefanis"

save_img(infile,opath,point=[100,100], point_px=[30,50])