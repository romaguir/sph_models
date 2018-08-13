import numpy as np
from sph_models.sph import write_epix

sph_file='/home/romaguir/Documents/sph_models/sph_models/models/S40RTS.sph'

dz = 25.0
dstart = 0.0 + dz
depths = np.arange(dstart,2900.0+dz,dz)
depths[0] = 24.4
depths[-1] = 2890.0

for i in range(0,len(depths)-1):
    fname='s40rts_{:1.1f}_{:1.1f}.epix'.format(depths[i],depths[i+1])
    midpoint = (depths[i] + depths[i+1]) / 2.0
    print 'writing epix layer for depth ',midpoint
    write_epix(sph_file=sph_file,depth=midpoint,pixel_width=1.0,
               out_file='epix_files/'+fname)
