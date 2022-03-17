from tools.destefanis.plot_hmi import save_img

infile="/gehme/data/sdo/hmi/20101111/hmi.m_45s.2010.11.11_06_21_45_TAI.magnetogram.fits"
opath="/gehme/scratch/destefanis"
p2plot=[100,100]
save_img(infile,opath,point=p2plot)