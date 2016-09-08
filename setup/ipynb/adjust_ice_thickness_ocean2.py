from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

# use the IC file from MOM to adust the thickness of the ice shelf
# this script should be used after ice_sheet_setup.py and then running the model

# ocean thickness
h = Dataset('../ncfiles/ISOMIP_IC_OCEAN2.nc').variables['h'][0,:]
h_all = h.sum(axis=0); h_all[h_all<1.0e-5] = np.nan

# ssh, proxy for base of ice shelf
ssh = -Dataset('../ncfiles/ISOMIP_IC_OCEAN2.nc').variables['ave_ssh'][0,:]
# bedrockTopography
topo = Dataset('../ncfiles/ocean_geometry.nc').variables['D'][:] 
# read 3D file
file3D = Dataset('../ncfiles/Ocean2_3D.nc','r+')
area = file3D.variables['area'][:]
thick = file3D.variables['thick'][:]

h_min = 20.0 # min thickness, in m
jm,im = area.shape
for j in range(jm):
	for i in range(im):
            if not np.isnan(h_all[j,i]):
		if (h_all[j,i]<h_min):
			print 'h_min not met at j,i',j,i
			if (h_all[j,i] < h_min/2.): # add ice
				thick[j,i] = thick[j,i] + (h_min - h_all[j,i])
				print 'Added at j,i',j,i
		        else: # remove ice
			        thick[j,i] = thick[j,i] - h_min/2.
				print 'Removed at j,i',j,i


#for j in range(jm):
#   for i in range(50,70):
#       if np.isnan(h_all[j,i]):
#          thick[j,i] = thick[j,i] - 2*h_min
#          print 'Removed at j,i',j,i
#
#       if (h_all[j,i] < 2*h_min):
#          thick[j,i] = thick[j,i] - (2*h_min - h_all[j,i])
#          print 'Removed at j,i',j,i
       

area[thick==0]=0.0
# 3D
file3D.variables['area'][:,:] = area[:,:]
file3D.variables['thick'][:,:] = thick[:,:]
file3D.close()

