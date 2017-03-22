from remapping import mom_remapping
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import os

j = 20 # section location in y
n = 10 # plot every n

# create dir locally
os.system('mkdir -p PNG')

# static variables
depth = netCDF4.Dataset('ocean_geometry.nc').variables['D'][j,:]
x = netCDF4.Dataset('ocean_geometry.nc').variables['geolon'][j,:]
ave_ssh = netCDF4.Dataset('ISOMIP_IC.nc').variables['ave_ssh'][0,j,:]
z = np.linspace(0,depth.max(),1000)
h = 0.5* (z[0:-1]+z[1::])

# time
tind = np.arange(0,7300,n)
tm = len(tind)
# remapping params
cs = mom_remapping.Remapping_Cs()
# PCM
#cs.remapping_scheme = 0 # PCM
# PPM-H4
cs.remapping_scheme = 2
cs.degree = 2

def remap(cs,h0,e0,h,data0,z):
   km,im = h0.shape
   data1 = np.ma.zeros((len(h),im))
   for i in range(im):
     h1 = np.diff(z)
     if h0[:,i].sum() > 0.01: # ocean
	h1[h<-e0[0,i]] = 0.0 # top
	h1[h>depth[i]] = 0.0 # bottom
	# need to account for SSH and make sure total thicknesses are
	# the same
	dh = h0[:,i].sum() - h1.sum()
        tmp1 = np.nonzero(h1!=0.0)[0]
        if len(tmp1)>0:
	  if dh > 0.0:
	     # correct thickness in the top non-vanished layer
             h1[tmp1[0]] = h1[tmp1[0]] + dh # add
          elif dh < 0.0:
             h1[tmp1[0]] = h1[tmp1[0]] - dh # remove
        else:
           data1[:,i] = np.ma.masked

        # for debugging
        #if h0[:,i].sum() != h1.sum():
	#   print 'WARNING: dh, h0[:,i].sum(), h1.sum()',dh, h0[:,i].sum(), h1.sum()

        # remap
	data1[:,i] = mom_remapping.remapping_core_h(h0[:,i], data0[:,i], h1, cs)
	# mask
	data1[h1==0.0,i] = np.ma.masked;
     else: # land/iceshelf
	data1[:,i] = np.ma.masked

   return data1

# plotting params
temp_levs = np.linspace(-2,1.1,50)
rho_levs = np.linspace(26.94,27.89,25) # ~ every 0.04
X,Z = np.meshgrid(x,h)
# main loop
for t in range(tm):
    time = netCDF4.Dataset('prog.nc').variables['time'][tind[t]]/365. # in years
    print 'Time is:',time
    temp0 = netCDF4.Dataset('prog.nc').variables['temp'][tind[t],:,j,:]
    rho0 = netCDF4.Dataset('prog.nc').variables['rhoinsitu'][tind[t],:,j,:]-1000.
    h0 = netCDF4.Dataset('prog.nc').variables['h'][tind[t],:,j,:]
    e0 = netCDF4.Dataset('prog.nc').variables['e'][tind[t],:,j,:]

    # remap
    temp1 = remap(cs,h0,e0,h,temp0,z)
    rho1 = remap(cs,h0,e0,h,rho0,z)
    # plot
    fig, ax = plt.subplots()
    ct = ax.contourf(X,-Z,temp1,temp_levs, cmap=plt.cm.jet)
    cbar = fig.colorbar(ct, ticks=[-2, -1, 0, 1], orientation = 'horizontal')
    cbar.set_label(r'Temperature [$^o$C]',fontsize=15)
    ax.contour(X,-Z,rho1,rho_levs,colors = 'k',linewidths=0.25)
    ax.set_xlim(x.min(),x.max())
    ax.set_ylim(-720,0)
    ss = str("ax.set_title('Time: %6.4f years')"% (time))
    eval(ss)
    #ax.set_title('Time: '+str(time)+' years')
    ax.set_xlabel('x [km]',fontsize=15)
    ax.set_ylabel('depth [m]',fontsize=15)
    ax.plot(x,ave_ssh,'k-',lw=1)
    ax.plot(x,-depth,'k',lw=1)
    ax.plot(x,e0[0,:],'k',lw=1)
    ss = str("plt.savefig('PNG/temp-%03d.png')"% (t))
    eval(ss)
    plt.close()

