import netCDF4
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt

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
plt.show()

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
plt.show()

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

sigma = [6,6] #  the standard deviation of the distribution
#lowerSurface_smoth[:,0:front] = gaussian_filter(lowerSurface[:,0:front],sigma)
lowerSurface_smoth = gaussian_filter(lowerSurface,sigma)
#upperSurface_smoth[:,0:front] = gaussian_filter(upperSurface[:,0:front],sigma)
upperSurface_smoth = gaussian_filter(upperSurface,sigma)

thick_smoth = upperSurface_smoth - lowerSurface_smoth

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
#plt.show()

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

# plot again
plt.figure(figsize=(12,6))
plt.contourf(X,Y,upperSurface_smoth - upperSurface,50)
plt.colorbar()
plt.title('upperSurface (differece)')
plt.xlabel('x [km]')
plt.ylabel('y [km]')

plt.figure(figsize=(12,6))
plt.contourf(X,Y,lowerSurface_smoth - lowerSurface, 50)
plt.colorbar()
plt.title('lowerSurface (difference)')
plt.xlabel('x [km]')
plt.ylabel('y [km]')

plt.figure(figsize=(12,6))
plt.contourf(X,Y,thick_smoth - thick,50)
plt.colorbar()
plt.title('Thickness (differece)')
plt.close('all')

plt.figure(figsize=(12,6))
plt.plot(lowerSurface_smoth[40,:],'b');plt.plot(lowerSurface[40,:],'r')
plt.xlim(300,350)
plt.xlabel('x [km]')
plt.ylabel('z [m]')
plt.title('Lower surface at x = 40 km')

plt.figure(figsize=(12,6))
plt.plot(lowerSurface_smoth[5,:],'b');plt.plot(lowerSurface[5,:],'r')
plt.xlim(300,350)
plt.xlabel('x [km]')
plt.ylabel('z [m]')
plt.title('Lower surface at x = 5 km')

plt.figure(figsize=(12,6))
plt.plot(lowerSurface_smoth[:,100],'b');plt.plot(lowerSurface[:,100],'r')
plt.xlabel('x [km]')
plt.ylabel('z [m]')
plt.title('Lower surface at y = 400 km')


plt.figure(figsize=(12,6))
plt.plot(lowerSurface_smoth[:,200],'b');plt.plot(lowerSurface[:,200],'r')
plt.xlabel('x [km]')
plt.ylabel('z [m]')
plt.title('Lower surface at y = 500 km')
plt.show()

# put into coarse grid
x1=x[::2];y1=y[::2]
thick1=thick_smoth[::2,::2]
height=lowerSurface_smoth[::2,::2]
top=upperSurface_smoth[::2,::2]
#height[thick1<100.]=0.0
#top[thick1<100.]=0.0
thick1 = top - height

# compute area
area = np.ones((thick1.shape))* (x1[1]-x1[0]) * (y1[1]-y1[0])
area[thick1==0.0]=0.0

# compute mass (kg/m^2)
rho_ice = 918.
mass = thick1 * rho_ice
# pressure (kg/(m s^2))
g=9.806
p_ice = mass * g
p_ocean = rho_warm * g * z
#def get_ice_draft(thickness,scale_factor,rho_ice)
#    psurf = scale_factor * thickness

# calculate ice_draft and ocean_thickness
min_thickness = 40.
im,jm = area.shape
ice_draft = np.zeros((im,jm))
ice_draft = np.ma.masked_where(area == 0, ice_draft)
ocean_thickness = np.zeros((im,jm))
ocean_thickness = np.ma.masked_where(area == 0, ocean_thickness)
for i in range(im):
    for j in range(jm):
	    if area[i,j] == 0:
	       ocean_thickness.mask[i,j] = True
	       ice_draft.mask[i,j] = True
            else:
	       if (p_ice[i,j] > p_ocean.max()):
		    ocean_thickness.mask[i,j] = True ; ice_draft[i,j]=B[i,j]
	       else:
		    # find height where ice shelf will float (if not grounded)
		    ice_draft[i,j]=-np.interp(p_ice[i,j], p_ocean, z)
		    if ice_draft[i,j]<=B[i,j]:
		       ice_draft[i,j]=B[i,j]; ocean_thickness.mask[i,j] = True
		    else:
		       ocean_thickness[i,j] = ice_draft[i,j]-B[i,j]

                    if (ocean_thickness[i,j] <= min_thickness):
                            # always force to be grounded
                            thick1[i,j]=thick1[i,j] + 2*min_thickness
                            ocean_thickness[i,j] = min_thickness
			    #if (ocean_thickness[i,j] <= min_thickness*0.5):
                            #   # force to be grounded
			    #   thick1[i,j]=thick1[i,j] + min_thickness*0.5 
			    #   ice_draft[i,j]=B[i,j]; ocean_thickness.mask[i,j] = True
                            #else:
			    #   # force space between ice shelf and ocean of at 
			    #   # least min_thickness	  
			    #   thick1[i,j]=thick1[i,j] - min_thickness*0.5
                            #   mass[i,j] = thick1[i,j] * rho_ice  
                            #   p_ice[i,j] = mass[i,j] * g
			    #   ice_draft[i,j]=-np.interp(p_ice[i,j], p_ocean, z)
			    #   ocean_thickness[i,j] = ice_draft[i,j]-B[i,j]

#save into netcdf file
# 3D
file3D.variables['area'][:,:] = area[:,:].T
file3D.variables['thick'][:,:] = thick1[:,:].T
file3D.variables['height'][:,:] = height[:,:].T

# 2D (middle of the domain)
for i in range(file2D.variables['area'].shape[1]):
    file2D.variables['area'][:,i] = area[20,:]
    file2D.variables['thick'][:,i] = thick1[20,:]
    file2D.variables['height'][:,i] = height[20,:] 

file2D.close()
file3D.close()
