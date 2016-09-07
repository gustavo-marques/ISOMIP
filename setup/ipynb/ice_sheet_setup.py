import netCDF4
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
from scipy import interpolate

#topography
# define parameters (in m)
Bmax=720.0 # max depth of bedrock topography
B0=-150.0  # bedrock topography at x =0
B2=-728.8  # Second bedrock topography coeff
B4=343.91  # Third bedrock topography coeff
B6=-50.57  # Forth bedrock topography coeff
x_bar=300.0e3 # Characteristic along-flow lenght scale of the bedrock
dc=500.0  # depth of the trough compared with side walls
fc=4.0e3  # Characteristic width of the side walls of the channel
wc=24.0e3 # half-width of the trough
Ly=80.0e3 # domain width (across ice flow)

# domain
dx=2.0e3 # resolution
x=np.arange(320.0e3 + dx*0.5,800.0e3,dx)
y=np.arange(0.0+dx*0.5,80.0e3,dx)

# eq. 3
x_til=x/x_bar
# eq. 2
Bx=B0+B2*x_til**2 + B4*x_til**4 + B6*x_til**6
# eq 4
By=(dc/(1+np.exp(-2*(y- Ly/2. - wc)/fc))) + (dc/(1+np.exp(2*(y- Ly/2. + wc)/fc)))
# eq 1
B=np.zeros((len(y),len(x)))
for i in range(len(x)):
    for j in range(len(y)):
        B[j,i]=np.max([Bx[i]+By[j],-Bmax])

[X,Y]=np.meshgrid(x,y)
plt.figure()
plt.contourf(X/1.0e3,Y/1.0e3,B)
plt.colorbar()
plt.xlabel('x (km)')
plt.ylabel('y (km)')
#plt.show()

# Initial T/S
def get_rho(S,T):
    rho0=1027.51; Sref=34.2; Tref=-1.0
    alpha = 3.733e-5; beta = 7.843e-4
    rho = rho0 * (1 - alpha * (T-Tref) + beta * (S-Sref))
    return rho

# initial T and S eqs(18 and 19)
z = np.linspace(0,-B.min(),1000)
# COLD
T0=-1.9; Tb=-1.9
S0=33.8; Sb=34.55
temp_cold = -1.9
salt_cold = S0 + (Sb - S0)* z/Bmax
rho_cold = get_rho(salt_cold,temp_cold)

# WARM
T0=-1.9; Tb=1.0
S0=33.8; Sb=34.7
temp_warm = T0 + (Tb - T0)* z/Bmax
salt_warm = S0 + (Sb - S0)* z/Bmax
rho_warm = get_rho(salt_warm,temp_warm)

plt.figure()
plt.plot(rho_cold-1000,z,'r',rho_warm-1000,z,'b')
plt.title('Cold (r), Warm (b)')
plt.xlabel('Density (kg/m^3)')
plt.ylabel('Depth (m)')
plt.grid()
#plt.show()
plt.close('all')

# ice shelf
# read provided netcdf
ncpath='../ncfiles/Ocean1_input_geom_v1.01.nc'
x = netCDF4.Dataset(ncpath).variables['x'][:]
y = netCDF4.Dataset(ncpath).variables['y'][:]
upperSurface = netCDF4.Dataset(ncpath).variables['upperSurface'][:]
lowerSurface = netCDF4.Dataset(ncpath).variables['lowerSurface'][:]
[X,Y]=np.meshgrid(x/1.0e3,y/1.0e3)

# read 3D file
file3D = netCDF4.Dataset('../ncfiles/Ocean1_3D.nc','r+')
# read 2D file
file2D = netCDF4.Dataset('../ncfiles/Ocean1_2D.nc','r+')

thick = upperSurface - lowerSurface
# calve when thick<100
#upperSurface[thick<100]=0.0
#lowerSurface[thick<100]=0.0
#thick = upperSurface - lowerSurface

# original data
plt.figure(figsize=(12,6))
plt.contourf(X,Y,upperSurface,50)
plt.colorbar()
plt.title('upperSurface [m]')
plt.figure(figsize=(12,6))
plt.contourf(X,Y,thick,50)
plt.colorbar()
plt.title('Thickness [m]')
plt.figure(figsize=(12,6))
plt.contourf(X,Y,lowerSurface, 50)
plt.colorbar()
plt.title('lowerSurface [m]')
plt.xlabel('x [km]')
plt.ylabel('y [km]')
#plt.show()

# find ice shelf front, it does not depend on y
front=np.nonzero(lowerSurface[0,:]==0)[0][0]
#front=np.nonzero(lowerSurface[0,:]==0)[0][-1]

# smoothing
lowerSurface_smoth = np.zeros(lowerSurface.shape)
upperSurface_smoth = np.zeros(upperSurface.shape)

sigma = [2,4] # (y,x) the standard deviation of the distribution
#lowerSurface_smoth[:,0:front] = gaussian_filter(lowerSurface[:,0:front],sigma)
lowerSurface_smoth = gaussian_filter(lowerSurface,sigma)
#upperSurface_smoth[:,0:front] = gaussian_filter(upperSurface[:,0:front],sigma)
upperSurface_smoth = gaussian_filter(upperSurface,sigma)

thick_smoth = upperSurface_smoth - lowerSurface_smoth

# apply calving
#upperSurface_smoth[thick_smoth<100]=0.0
#lowerSurface_smoth[thick_smoth<100]=0.0
#thick_smoth = upperSurface_smoth - lowerSurface_smoth

# plot after first smoothing
plt.figure(figsize=(12,6))
plt.contourf(X,Y,upperSurface_smoth - upperSurface,50)
plt.colorbar()
plt.title('upperSurface (differece)')
plt.figure(figsize=(12,6))
plt.contourf(X,Y,thick_smoth - thick,50)
plt.colorbar()
plt.title('Thickness (differece)')
plt.figure(figsize=(12,6))
plt.contourf(X,Y,lowerSurface_smoth - lowerSurface, 50)
plt.colorbar()
plt.title('lowerSurface (difference)')
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.close('all')

# smooth the shelf fron, just for layer model
#sigma = 1 # this seems to be reasonable
#front10=front-10
#for j in range(lowerSurface.shape[0]):
#    lowerSurface_smoth[j,front+1]=-50;lowerSurface_smoth[j,front+2]=-25.;lowerSurface_smoth[j,front+3]=-12.5
    #lowerSurface_smoth[j,front10::] = gaussian_filter(lowerSurface_smoth[j,front10::],sigma)
    #upperSurface_smoth[j,front::] = gaussian_filter(upperSurface_smoth[j,front::],sigma)
    
#thick_smoth = upperSurface_smoth - lowerSurface_smoth
#lowerSurface_smoth[thick<100.]=0.0
#thick_smoth = upperSurface_smoth - lowerSurface_smoth


# interpolate to coarse grid
xnew=x[::2];ynew=y[::2]
f_thick = interpolate.interp2d(x, y, thick_smoth, kind='cubic')
thick_new=f_thick(xnew, ynew)
x=xnew/1.0e3;y=ynew/1.0e3
jm,im = thick_new.shape

# calve ice thiner than 100 m until x=145 (before the front)
for j in range(jm):
	for i in range(0,146):
		if (thick_new[j,i]<100):
			thick_new[j,i] = 0.0


# compute mass (kg/m^2)
rho_ice = 918.
mass = thick_new * rho_ice
# pressure (kg/(m s^2))
g=9.806
p_ice = mass * g
p_ocean = rho_warm * g * z
#def get_ice_draft(thickness,scale_factor,rho_ice)
#    psurf = scale_factor * thickness

# constrain ocean_thickness
min_thickness = 20.0
thick_new[thick_new<2.] = 0.0
area = np.ones((thick_new.shape))* (xnew[1]-xnew[0]) * (ynew[1]-ynew[0])
area[thick_new==0]=0.0
ice_mask = np.zeros(area.shape);draft = np.zeros(area.shape)
for i in range(im):
    for j in range(jm):
	    if area[j,i] == 0:
	       ice_mask[j,i]=1
               print 'Outside ice shelf!'
            else:
               ind = np.nonzero(z<=-B[j,i])[-1][-1]
	       if (p_ice[j,i] > p_ocean[ind]):
		    ice_mask[j,i]=0
                    draft[j,i] = B[j,i]
		    print 'Grounded!' 
	       else:
		    # find height where ice shelf will float (if not grounded)
		    ice_draft=-np.interp(p_ice[j,i], p_ocean, z)
		    if ice_draft>(B[j,i]+min_thickness):
		       ice_mask[j,i]=0.5
		       draft[j,i] = ice_draft
                       print 'Floating and min_thickness is fine.'
		    else:
		       ice_mask[j,i]=0.25
		       if ((ice_draft-B[j,i])<=min_thickness/4.):
                            thick_new[j,i]=thick_new[j,i] + min_thickness 
                            tmp1=(thick_new[j,i] * rho_ice)*g
		       	    tmp2=-np.interp(tmp1, p_ocean, z)
                            draft[j,i] = tmp2
                       else:
                            thick_new[j,i]=thick_new[j,i] + (min_thickness-(ice_draft-B[j,i]))
                            tmp1=(thick_new[j,i] * rho_ice)*g
		            tmp2=-np.interp(tmp1, p_ocean, z)
		            draft[j,i] = tmp2

ocean_thick = draft - B
ind1,ind2 = np.nonzero((ocean_thick<min_thickness) & (ocean_thick>0.))
for i in range(len(ind1)):
	if(ocean_thick[ind1[i],ind2[i]] < min_thickness/4.):
		thick_new[ind1[i],ind2[i]] = thick_new[ind1[i],ind2[i]] + min_thickness/2.
	else:
		thick_new[ind1[i],ind2[i]] = thick_new[ind1[i],ind2[i]] - 3*min_thickness/4. 
# update pice
mass = thick_new * rho_ice
p_ice = mass * g

# assures connectivity between cells 
# min. thickness of 3 * min_thickness
for i in range(98,120):
   for j in range(jm):
       ice_draft=-np.interp(p_ice[j,i], p_ocean, z)
       if ice_draft>(B[j,i]+4*min_thickness):
          print 'Floating and min_thickness is fine.'
       else:
          thick_new[j,i]=thick_new[j,i] - 4*min_thickness

# manuall changes
thick_new[29,71]=thick_new[29,71] + min_thickness
thick_new[27,71]=thick_new[27,71] - min_thickness
thick_new[28,71]=thick_new[28,71] - min_thickness
for i in range(72,90):
	for j in np.array([7,8,29,30,31,32]):
		thick_new[j,i]=thick_new[j,i] + min_thickness


for j in range(jm):
        for i in range(0,146):
                if (thick_new[j,i]<100):
                        thick_new[j,i] = 0.0

# update mass and pressure
mass = thick_new * rho_ice
p_ice = mass * g

#smooth one mor etime
#sigma = [2,2] # (y,x) the standard deviation of the distribution
#thick_new = gaussian_filter(thick_new,sigma)

# remove ice < 0.1 m
#thick_new[thick_new<0.1] = 0.0

# update area
area = np.ones((thick_new.shape))* (xnew[1]-xnew[0]) * (ynew[1]-ynew[0])
area[thick_new==0]=0.0

#save into netcdf file
# 3D
file3D.variables['area'][:,:] = area[:,:]
file3D.variables['thick'][:,:] = thick_new[:,:]

# 2D (middle of the domain) y = 40 km
for j in range(file2D.variables['area'].shape[0]):
    file2D.variables['area'][j,:] = area[20,:]
    file2D.variables['thick'][j,:] = thick_new[20,:]

file2D.close()
file3D.close()
