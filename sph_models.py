#!/home/romaguir/anaconda2/bin/python
from sys import argv
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import pyshtools
from scipy.interpolate import interp1d,RegularGridInterpolator
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
      grid = sph_spline.expand('DH2')
      print grid.info
      map_dv += spl_vals[i] * grid.data

   return map_dv

def get_epixarr(sph_file,depth,pixel_width=1.0):
   lon_start = pixel_width / 2.0
   lat_start = -90 + (pixel_width / 2.0)
   lons = np.arange(lon_start,360.0,pixel_width/2.0)
   lats = np.arange(lat_start,90.0,pixel_width/2.0)
   print lons,lats


def radial_correlation_function(sph_file,lmin,lmax):

   depths = np.arange(50,2800,10)
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

def plot_slice(sph_file,lon0,lat0,lon1,lat1,save_pts=False):
   #sph_splines = read_sph(sph_file,lmin,lmax)
   #spl_vals = find_spl_vals(depth)
   dep_maps = []

   depths = np.arange(50,2800,50)
   lats = np.linspace(-90,90,82)
   lons = np.linspace(0,360,164)

   for depth in depths:
      print depth
      dep_maps.append(extract_dep_map(sph_file,depth))

   points = (depths,lats,lons)
   values = np.array(dep_maps)
   rgi = RegularGridInterpolator(points=points,values=values)

   lons_slice = np.linspace(lon0,lon1,100)
   lats_slice = np.linspace(lat1,lat1,100)
   sample_pts = []

   for i in range(0,len(lons_slice)):
      for j in range(0,len(lats_slice)):
         for k in range(0,len(depths)):
            sample_pts.append([depths[k],lats_slice[j],lons_slice[i]])

   sample_vals = rgi(sample_pts)
   print sample_vals
   x = []
   y = []
   z = []
   for i in range(0,len(sample_vals)):
      ####THIS IS A TEST FOR A LONGITUDINAL SLICE
      x.append(sample_pts[i][2])
      y.append(sample_pts[i][0]) 
      z.append(sample_vals[i])

   #NOTE, ONLY VALID FOR LONGITUDINAL SLICES!
   if save_pts:
      pts_file = open('tomo_pts.dat','w')
      for i in range(0,len(sample_vals)):
         z_here = sample_pts[i][0]
         rad_here = 6371.0 - z_here
         lon_here = sample_pts[i][2]
         val_here = sample_vals[i]*100.0
         pts_file.write('{} {} {}\n'.format(lon_here,rad_here,val_here))

   #plt.scatter(x,y,c=z,edgecolor='none',cmap='jet_r',vmin=-0.02,vmax=0.02)
   #plt.gca().invert_yaxis()
   #plt.show()

   a = np.reshape(z,((len(lons_slice)*len(lats_slice),len(depths)) ))
   plt.imshow(a,aspect='auto',vmin=-0.02,vmax=0.02,cmap='jet_r')
   plt.colorbar()
   plt.show()

def find_min_max(sph_file):
   dep_maps = []

   depths = np.arange(50,2800,50)
   lats = np.linspace(-90,90,82)
   lons = np.linspace(0,360,164)

   for depth in depths:
      dep_maps.append(extract_dep_map(sph_file,depth))
   model_vals = np.array(dep_maps)
   return np.min(model_vals), np.max(model_vals)
