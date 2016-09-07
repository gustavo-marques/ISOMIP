from netCDF4 import Dataset
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt

# use the IC file from MOM to adust the thickness of the ice shelf
# this script should be used after ice_sheet_setup.py and then running the model

# ocean thickness
h = Dataset('../ncfiles/ISOMIP_IC_OCEAN1.nc').variables['h'][0,:]
# ssh, proxy for base of ice shelf
ssh = -Dataset('../ncfiles/ISOMIP_IC_OCEAN1.nc').variables['ave_ssh'][0,:]
# bedrockTopography
topo = Dataset('../ncfiles/ocean_geometry.nc').variables['D'][:] 
# read 3D file
file3D = netCDF4.Dataset('../ncfiles/Ocean1_3D.nc','r+')
area = file3D.variables['area'][:]
thick = file3D.variables['thick'][:]

h_min = 20.0 # min thickness, in m
jm,im = area.shape
for j in range(jm):
	for i in range(im):
		if area[j,i]>0: # ice shelf
			if ((ssh[j,i]+0.1)>topo[j,i]):
				print 'Grounded at j,i',j,i
			else:
				if (h[:,j,i].sum() < h_min/2.): # add ice
					thick[j,i] = thick[j,i] + (h_min - h[:,j,i].sum())
					print 'Added at j,i',j,i
				elif (h[:,j,i].sum() < h_min): # remove ice
					thick[j,i] = thick[j,i] - h_min/2.
					print 'Removed at j,i',j,i
				else:
					print 'Open at j,i',j,i


h_all = h.sum(axis=0)
h_all[h_all<1.0e-5]=-9999
for j in range(jm):
        for i in range(im):
            if (h_all[j,i]>0): 
               if (h_all[j,i]<h_min):
                  if (h_all[j,i]<h_min/2.):
                     thick[j,i] = thick[j,i] + h_all[j,i]
                  else:
                     thick[j,i] = thick[j,i] - h_all[j,i]

area[thick==0]=0.0

# 3D
file3D.variables['area'][:,:] = area[:,:]
file3D.variables['thick'][:,:] = thick[:,:]
file3D.close()

