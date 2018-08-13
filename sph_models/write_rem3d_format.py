import glob
import numpy as np
import fortranformat as ff
import matplotlib.pyplot as plt

fname='S40RTS_pixel.r3d'
kernel_set='BOX25+I1D'
epix_dir='epix_files'
ref_model='PREM'
r_desc='(SH+SV)*0.5'

dz = 25.0
dstart = 0.0 + dz
depths = np.arange(dstart,2900.0+dz,dz)
depths[0] = 24.4
depths[-1] = 2890.0

epixfiles = glob.glob(epix_dir+'/*')

#----------------------------------------------------------------------
#write header
#----------------------------------------------------------------------
def write_header(rem3d_file,ref_model,kernel_set):
    rem3d_file.write('REFERENCE MODEL: {} \n'.format(ref_model))
    rem3d_file.write('KERNEL SET: {}\n'.format(kernel_set))

#----------------------------------------------------------------------
#write radial kernels
#----------------------------------------------------------------------
def write_radial_kernels(rem3d_file,epix_dir,r_desc):

    n_rad_kernels = len(depths)-1
    rem3d_file.write('RADIAL STRUCTURE KERNELS: {}\n'.format(n_rad_kernels))
    for i in range(0,n_rad_kernels):
        rem3d_file.write('DESC  {:3.0f}: {}, {:1.1f} - {:1.1f} km\n'.format(i+1,
                         r_desc,depths[i],depths[i+1]))

#----------------------------------------------------------------------
#write horizontal parameterization(s)
#----------------------------------------------------------------------
def write_horizontal_param(rem3d_file,n_hpar=1,pixel_width=1.0):
    rem3d_file.write('HORIZONTAL PARAMETERIZATIONS: {}\n'.format(n_hpar))

    for i in range(0,n_hpar):
        rem3d_file.write('HPAR   {}: PIXELS,  {:1.1f} x {:1.1f}\n'.format(
                         i+1,pixel_width,pixel_width))

        shape = (int(180.0/pixel_width),int(360.0/pixel_width))
        epixarr0 = np.loadtxt(epixfiles[0])
        lats = epixarr0[:,0]
        lons = epixarr0[:,1]
        lats = np.reshape(lats,shape,order='F')
        lons = np.reshape(lons,shape,order='F')
        lats_pts = lats.flatten()
        lons_pts = lons.flatten()

        for j in range(0,len(lons_pts)):
            lon_here = lons_pts[j]
            lat_here = lats_pts[j]
            if lon_here > 180.0:
               lon_here -= 360.0
            rem3d_file.write('{:5.1f} {:5.1f} {:5.1f}\n'.format(lon_here,
                             lat_here, pixel_width))

#----------------------------------------------------------------------
#write model coefficients
#----------------------------------------------------------------------
def write_model_coefs(rem3d_file,pixel_width=1.0):

    line = ff.FortranRecordWriter('(6E12.4)')

    shape = (int(180.0/pixel_width),int(360.0/pixel_width))
    for i in range(0,len(depths)-1):
        print 'writing coefficients for layer ',i+1
        epixname='s40rts_{:1.1f}_{:1.1f}.epix'.format(depths[i],depths[i+1])
        epixarr = np.loadtxt(epix_dir+'/'+epixname)
        coefs = epixarr[:,3]
        coefs = np.reshape(coefs,shape,order='F')
        coefs_pts = coefs.flatten()

        rem3d_file.write('STRU  {:3.0f}:  {:1.0f}\n'.format(i+1,pixel_width))
        rem3d_file.write(line.write(coefs_pts)+'\n')

#----------------------------------------------------------------------
#main program
#----------------------------------------------------------------------
def main(fname,ref_model,kernel_set,epix_dir,r_desc):
    rem3d_file = open(fname,'w')
    write_header(rem3d_file=rem3d_file,ref_model=ref_model,
                 kernel_set=kernel_set)
    write_radial_kernels(rem3d_file=rem3d_file,epix_dir=epix_dir,
                         r_desc=r_desc)
    write_horizontal_param(rem3d_file=rem3d_file)
    write_model_coefs(rem3d_file=rem3d_file)

#run main
main(fname=fname,ref_model=ref_model,kernel_set=kernel_set,
     epix_dir=epix_dir,r_desc=r_desc)
