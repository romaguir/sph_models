#!/home/romaguir/anaconda2/bin/python
from sys import argv
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import pyshtools
from scipy.interpolate import interp1d
from copy import deepcopy

n_splines = 21

def read_splines(splines_dir='./splines'):

   splines = []
   for i in range(1,n_splines+1):
      g = np.loadtxt(splines_dir+'/spline{}.dat'.format(i))
      splines.append(g)

   return splines

def read_sph(sph_file,lmin=0,lmax=40):

   #Read header
   f = open(sph_file)
   lines = f.readlines()
   sph_degree = int(lines[0].strip().split()[0])
   f.close()
   
   #Read spherical harmonic coefficients
   vals = []
   for line in lines[1:]:
      vals_list = line.strip().split()
      for val in vals_list:
         vals.append(np.float(val))

   vals = np.array(vals)
   sph_splines = []
   count = 0

   for k in range(0,n_splines):
      coeffs = np.zeros((2,sph_degree+1,sph_degree+1))
      for l in range(0,sph_degree+1):
         ind = 0
         nm = (l*2) + 1
         c_row = []
         n_m = 1. 
   
         for m in range(0,(2*l)+1):
            if m == 0:
               order = 0
               coeffs[0,l,0] = vals[count]
               count += 1
            elif np.mod(m,2) == 1:
               order = int(m * (n_m)/m)
               coeffs[0,l,order] = vals[count]
               count += 1
            elif np.mod(m,2) == 0:
               order = int(m * (n_m)/m)
               coeffs[1,l,order] = vals[count]
               count += 1
               n_m += 1

      #filter out unwanted degrees
      if lmin != 0:
         coeffs[0,0:lmin,:] = 0.0
         coeffs[1,0:lmin,:] = 0.0
      if lmax < 40:
         coeffs[0,lmax+1:,:] = 0.0
         coeffs[1,lmax+1:,:] = 0.0

      clm = pyshtools.SHCoeffs.from_array(coeffs,normalization='ortho',csphase=-1)
      sph_splines.append(deepcopy(clm))

   return sph_splines


def find_spl_vals(depth):

   splines = read_splines()
   spl_vals = []

   if depth < 0 or depth > 2891:
      raise ValueError("depth should be in range 0 - 2891 km")

   for spline in splines:
      z = spline[:,0]
      f_z = spline[:,1]
      get_val = interp1d(z,f_z)
      spl_val = get_val(depth)
      spl_vals.append(spl_val)

   return spl_vals 

def plot_tomo_map(sph_file,depth,vmin=-2.0,vmax=2.0,lmin=0,lmax=40):
   sph_splines = read_sph(sph_file,lmin,lmax)
   spl_vals = find_spl_vals(depth)
   map_dv = 0.0

   for i,sph_spline in enumerate(sph_splines):
      grid = sph_spline.expand()
      map_dv += spl_vals[i] * grid.data

   lons,lats = np.meshgrid(grid.lons(),grid.lats())
   tomo_map = plt.axes(projection=ccrs.Mollweide(180))
   tomo_map.pcolor(lons,lats,map_dv*100.0,
                     transform=ccrs.PlateCarree(),
                     cmap='jet_r',vmin=vmin,vmax=vmax)
   tomo_map.coastlines()
   plt.show()

def extract_dep_map(sph_file,depth,lmin=0,lmax=40):
   sph_splines = read_sph(sph_file,lmin,lmax)
   spl_vals = find_spl_vals(depth)
   map_dv = 0.0

   for i,sph_spline in enumerate(sph_splines):
      grid = sph_spline.expand()
      map_dv += spl_vals[i] * grid.data

   return map_dv

def radial_correlation_function(sph_file,lmin,lmax):

   depths = np.arange(50,2800,20)
   maps_list = []
   for depth in depths:
      map_dv = extract_dep_map(sph_file,depth,lmin,lmax)
      maps_list.append(map_dv)

   rad_corr_func = np.zeros((len(depths),len(depths)))
   for i in range(0,len(depths)):
      for j in range(0,len(depths)):
         corr_coef = np.corrcoef(np.ravel(maps_list[i]),np.ravel(maps_list[j]))
         rad_corr_func[i,j] = corr_coef[0,-1]

   plt.imshow(rad_corr_func,cmap='jet',vmin=-1.0,vmax=1.0)
   plt.colorbar()
   plt.show()

def radial_rms_function(sph_file,lmin,lmax):

   depths = np.arange(50,2800,10)
   maps_list = []
   for depth in depths:
      map_dv = extract_dep_map(sph_file,depth,lmin,lmax)
      maps_list.append(map_dv)

   rad_rms_func = np.zeros((len(depths),len(depths)))
   for i in range(0,len(depths)):
      for j in range(0,len(depths)):
         a = np.ravel(maps_list[i]) / np.abs(np.max(maps_list[i]))
         b = np.ravel(maps_list[j]) / np.abs(np.max(maps_list[j]))
         #resid = np.ravel(maps_list[i]) - np.ravel(maps_list[j])
         resid = a - b
         rms = np.sqrt(np.mean(resid**2))
         rad_rms_func[i,j] = rms

   plt.imshow(rad_rms_func,cmap='jet')
   plt.colorbar()
   plt.show()
      

#radial_correlation_function('models/S40RTS.sph')
#radial_correlation_function('models/S40RTS.sph',lmin=0,lmax=4)
#radial_rms_function('models/S40RTS.sph',lmin=0,lmax=4)

plot_tomo_map(sph_file='models/S40RTS.sph',depth=2800.0,vmin=-2.0,vmax=2.0,lmin=0,lmax=40)
