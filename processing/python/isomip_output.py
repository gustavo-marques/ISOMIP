#!/usr/bin/env python

# generate the ISOMIP + diagnostics 
# Gustavo Marques, Aug. 2016

import argparse
from netCDF4 import MFDataset, Dataset
import numpy as np
import matplotlib.pyplot as plt
import warnings
import os

class MyError(Exception):
  """
  Class for error handling
  """
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Compute the isomip+ diagnostics and write them into a netCDF file. This script 
      assumes that files/variables names are consistent with the default ISOMIP/GFDL
      options (see: https://github.com/NOAA-GFDL/MOM6-examples/tree/dev/master/ocean_on       ly/ISOMIP).
      ''',
  epilog='Written by Gustavo Marques, Aug. 2016.')

  parser.add_argument('-type', type=str, default='ocean0',
      help='''The type pf experiment that will computed (ocean0, ocean1 etc). Default is ocean0.''')

  parser.add_argument('-n', type=str, default='Ocean0_COM_MOM6-LAYER',
      help='''The name of the experiment following the ISOMIP+ definition (expt_COM_model). This name is used to save the netCDF file. Default is Ocean0_COM_MOM6-LAYER.''')

  parser.add_argument('--mean4gamma', help='''Compute (and save in a text file) melting over the area where ice base > 300 m and over final six months.''', action="store_true")
 
  optCmdLineArgs = parser.parse_args()
  global name
  name = optCmdLineArgs.n
  driver(optCmdLineArgs)

def driver(args):
   """
   This is where all the action happends, i.e., netcdf is generated/populated and functions for each diagnostic are called.
   """
   if name == 'Ocean0_COM_MOM6-LAYER':
      print("WARNING: exp. name has not been defined, using the default name Ocean0_COM_MOM6-LAYER.")
   
   # load essential variables
   # bedrock
   depth = Dataset('ocean_geometry.nc').variables['D'][:]
   # area under shelf 
   shelf_area = Dataset('MOM_Shelf_IC.nc').variables['shelf_area'][0,:,:]
   # base of STATIC ice shelf, which is ssh(t=0); make it positive
   ice_base = -Dataset('ISOMIP_IC.nc').variables['ave_ssh'][0,:,:]
   ice_base[ice_base<1e-5] = 0.0
   #mask grounded ice and open ocean
   shelf_area = mask_grounded_ice(shelf_area,depth,ice_base) 
   shelf_area = mask_ocean(shelf_area,shelf_area)
   
   x=Dataset('ocean_geometry.nc').variables['lonh'][:]*1.0e3 # in m
   #lonq=Dataset('ocean_geometry.nc').variables['lonq'][:]  
   y=Dataset('ocean_geometry.nc').variables['lath'][:]*1.0e3 # in m
   #latq=Dataset('ocean_geometry.nc').variables['latq'][:]
   #lonqs, latqs = np.meshgrid(lonq,latq)
   x,y = np.meshgrid(x,y)
   
   # create dictionary
   # not for now 
   #data = {'name':name, 'depth':depth, 'area':area, }

   if args.mean4gamma:
      print("Computing mean melt rate to calibrate gamma...")
      melt4gamma(name,shelf_area,ice_base,depth)
   
   # read entire time variable from montly mean
   # it might be useful to add the option of passing
   # just part of the data (e.g., last 6 months)
   time = Dataset('ocean_month.nc').variables['time'][:]
   # create ncfile and zero fields. Each function below will corresponding values
   create_ncfile(name,time,args.type)

   # load additional variables
   ocean_area = Dataset('ocean_geometry.nc').variables['Ah'][:] 
   h = Dataset('ocean_month.nc').variables['h'][:]
   h[h<1.0e-5] = 0.0
   melt = Dataset('ocean_month.nc').variables['melt'][:]
   mass_flux = Dataset('ocean_month.nc').variables['mass_flux'][:]
   # mask open ocean and grounded ice
   ocean_area = mask_grounded_ice(ocean_area,depth,ice_base)
   melt = mask_grounded_ice(melt,depth,ice_base)
   melt = mask_ocean(melt,shelf_area)
   # tracers
   salt = Dataset('ocean_month.nc').variables['salt'][:]   
   temp = Dataset('ocean_month.nc').variables['temp'][:]
   rho = get_rho(salt,temp)

   # compute total volume
   total_volume(ocean_area,h)

   # compute mean melt rate
   mean_melt_rate(melt,shelf_area)
   
   # compute total melt flux
   total_melt_flux(mass_flux)

   # mean temp
   mean_tracer(ocean_area,h,temp,'meanTemperature','C')
   
   # mean salt
   mean_tracer(ocean_area,h,salt,'meanSalinity','PSU')

   # horizontal fields
   
   # bathymetry (negative)
   bathymetry = -mask_grounded_ice(depth,depth,ice_base)
   bathymetry.fill_value=0.0
   saveXY(bathymetry,'bathymetry')

   # meltRate, already masked above
   melt = melt/(3600.*24*365) # in m/s
   saveXY(melt,'meltRate')
   
   # frictionVelocity
   ustar_shelf = Dataset('ocean_month.nc').variables['ustar_shelf'][:]
   # mask open ocean and grounded ice
   ustar_shelf = mask_grounded_ice(ustar_shelf,depth,ice_base)
   ustar_shelf = mask_ocean(ustar_shelf,shelf_area)
   saveXY(ustar_shelf,'frictionVelocity')

   # thermalDriving
   thermal_driving = Dataset('ocean_month.nc').variables['thermal_driving'][:]
   thermal_driving = mask_grounded_ice(thermal_driving,depth,ice_base)
   thermal_driving = mask_ocean(thermal_driving,shelf_area)
   saveXY(thermal_driving,'thermalDriving')

   # halineDriving
   haline_driving = Dataset('ocean_month.nc').variables['haline_driving'][:]
   haline_driving = mask_grounded_ice(haline_driving,depth,ice_base)
   haline_driving = mask_ocean(haline_driving,shelf_area)
   saveXY(haline_driving,'halineDriving')

   # uBoundaryLayer
   u_ml = Dataset('ocean_month.nc').variables['u_ml'][:]
   u_ml = mask_grounded_ice(u_ml,depth,ice_base)
   u_ml = mask_ocean(u_ml,shelf_area)
   saveXY(u_ml,'uBoundaryLayer')

   # vBoundaryLayer
   v_ml = Dataset('ocean_month.nc').variables['v_ml'][:]
   v_ml = mask_grounded_ice(v_ml,depth,ice_base)
   v_ml = mask_ocean(v_ml,shelf_area)
   saveXY(v_ml,'vBoundaryLayer')

   #if (args.type == 'ocean3' or args.type == 'ocean4'):
   # iceDraft
   # will have to works this out when we run these cases
   iceDraft = mask_grounded_ice(ice_base,depth,ice_base)
   iceDraft = mask_ocean(iceDraft,shelf_area)
   saveXY(-iceDraft,'iceDraft')
   
   # data from ocean_month_z
   temp_z = Dataset('ocean_month_z.nc').variables['temp'][:]
   salt_z = Dataset('ocean_month_z.nc').variables['salt'][:]
   
   # data at bottom most cell
   bottomTemperature = get_bottom_data(temp_z)
   bottomSalinity = get_bottom_data(salt_z)
   saveXY(bottomTemperature,'bottomTemperature')
   saveXY(bottomSalinity,'bottomSalinity')

   # XZ y = 40 km (j=20)
   temperatureXZ = temp_z[:,:,20,:]
   salinityXZ = salt_z[:,:,20,:]
   saveXY(temperatureXZ,'temperatureXZ')
   saveXY(salinityXZ,'salinityXZ')

   # YZ x = 520 km (i=100)
   temperatureYZ = temp_z[:,:,:,100]
   salinityYZ = salt_z[:,:,:,100]
   saveXY(temperatureYZ,'temperatureYZ')
   saveXY(salinityYZ,'salinityYZ')

   # barotropic streamfunction
   ubt = Dataset('ocean_month.nc').variables['ubtav'][:]
   vbt = Dataset('ocean_month.nc').variables['vbtav'][:]
   # mask grouded region
   ubt = mask_grounded_ice(ubt,depth,ice_base)
   vbt = mask_grounded_ice(vbt,depth,ice_base)
   psi = get_psi(x,y,ubt,vbt)
   saveXY(psi,'barotropicStreamfunction')

   print('Done!')
   return

def get_psi(x,y,u,v):
    '''
    Loop in time and cCall flowfun to compute the streamfunction psi of a 2D flow.
    '''
    NT,NY,NX = u.shape
    psi = np.zeros(u.shape)
    for t in range(NT):
        psi[t,:,:] = flowfun(x,y,u[t,:],v[t,:])

    return psi 

def flowfun(x,y,u,v,variable='psi'):
	"""
	FLOWFUN  Computes the potential PHI and the streamfunction PSI
	 of a 2-dimensional flow defined by the matrices of velocity
	 components U and V, so that
	       d(PHI)    d(PSI)          d(PHI)    d(PSI)
	  u =  -----  -  ----- ,    v =  -----  +  -----
	        dx        dy              dx        dy
	P = FLOWFUN(x,y,u,v) returns an array P of the same size as u and v,
	which can be the velocity potential (PHI) or the streamfunction (PSI)
	Because these scalar fields are defined up to the integration constant,
	their absolute values are such that PHI[0,0] = PSI[0,0] = 0.
	For a potential (irrotational) flow  PSI = 0, and the Laplacian
	of PSI is equal to the divergence of the velocity field.
	A solenoidal (non-divergent) flow can be described by the
	streamfunction alone, and the Laplacian of the streamfunction
	is equal to the vorticity (curl) of the velocity field.
	The units of the grid coordinates are assumed to be consistent
	with the units of the velocity components, e.g., [m] and [m/s].
	If variable=='psi', the streamfunction (PSI) is returned.
	If variable=='phi', the velocity potential (PHI) is returned.
	Uses function 'cumsimp()' (Simpson rule summation).
	Author: Kirill K. Pankratov, March 7, 1994.
	Source: http://www-pord.ucsd.edu/~matlab/stream.htm
	Translated to Python by Andre Paloczy, January 15, 2015.
        Modified by Gustavo Marques, September, 2016.
	"""
	x,y,u,v = map(np.asanyarray, (x,y,u,v))

	if not x.shape==y.shape==u.shape==v.shape:
		print "Error: Arrays (x, y, u, v) must be of equal shape."
		return

	## Calculating grid spacings.
	dy, _ = np.gradient(y)
	_, dx = np.gradient(x)

	ly, lx = x.shape                         # Shape of the (x,y,u,v) arrays.

	## Now the main computations.
	## Integrate velocity fields to get potential and streamfunction.
	## Use Simpson rule summation (function CUMSIMP).

	## Compute velocity potential PHI (non-rotating part).
	if variable=='phi':
		cx = cumsimp(u[0,:]*dx[0,:])         # Compute x-integration constant
		cy = cumsimp(v[:,0]*dy[:,0])         # Compute y-integration constant
		cx = np.expand_dims(cx, 0)
		cy = np.expand_dims(cy, 1)
		phiy = cumsimp(v*dy) + np.tile(cx, (ly,1))
		phix = cumsimp(u.T*dx.T).T + np.tile(cy, (1,lx))
		phi = (phix + phiy)/2.
		return phi

	## Compute streamfunction PSI (non-divergent part).
	if variable=='psi':
		cx = cumsimp(v[0,:]*dx[0,:])         # Compute x-integration constant
		cy = cumsimp(u[:,0]*dy[:,0])         # Compute y-integration constant
		cx = np.expand_dims(cx, 0)
		cy = np.expand_dims(cy, 1)
		psix = -cumsimp(u*dy) + np.tile(cx, (ly,1))
		psiy = cumsimp(v.T*dx.T).T - np.tile(cy, (1,lx))
		psi = (psix + psiy)/2.
		return psi

def cumsimp(y):
	"""
	F = CUMSIMP(Y)    Simpson-rule column-wise cumulative summation.
	Numerical approximation of a function F(x) such that
	Y(X) = dF/dX.  Each column of the input matrix Y represents
	the value of the integrand  Y(X)  at equally spaced points
	X = 0,1,...size(Y,1).
	The output is a matrix  F of the same size as Y.
	The first row of F is equal to zero and each following row
	is the approximation of the integral of each column of matrix
	Y up to the givem row.
	CUMSIMP assumes continuity of each column of the function Y(X)
	and uses Simpson rule summation.
	Similar to the command F = CUMSUM(Y), exept for zero first
	row and more accurate summation (under the assumption of
	continuous integrand Y(X)).
	Author: Kirill K. Pankratov, March 7, 1994.
	Source: http://www-pord.ucsd.edu/~matlab/stream.htm
	Translated to Python by Andr? Pal?czy, January 15, 2015.
	"""
	y = np.asanyarray(y)

	## 3-point interpolation coefficients to midpoints.
	## Second-order polynomial (parabolic) interpolation coefficients
	## from  Xbasis = [0 1 2]  to  Xint = [.5 1.5]
	c1 = 3/8.
	c2 = 6/8.
	c3 = -1/8.

	if y.ndim==1:
		y = np.expand_dims(y,1)
		f = np.zeros((y.size,1))    # Initialize summation array.
		squeeze_after = True
	elif y.ndim==2:
		f = np.zeros(y.shape)       # Initialize summation array.
		squeeze_after = False
	else:
		print "Error: Input array has more than 2 dimensions."
		return

	if y.size==2:                   # If only 2 elements in columns - simple average.
		f[1,:] = (y[0,:] + y[1,:])/2.
		return f
	else:                           # If more than two elements in columns - Simpson summation.
		## Interpolate values of y to all midpoints.
		f[1:-1,:] = c1*y[:-2,:] + c2*y[1:-1,:] + c3*y[2:,:]
		f[2:,:] = f[2:,:] + c3*y[:-2,:] + c2*y[1:-1,:] + c1*y[2:,:]
		f[1,:] = f[1,:]*2
		f[-1,:] = f[-1,:]*2

		## Simpson (1,4,1) rule.
		f[1:,:] = 2*f[1:,:] + y[:-1,:] + y[1:,:]
		f = np.cumsum(f, axis=0)/6. # Cumulative sum, 6 - denominator from the Simpson rule.

	if squeeze_after:
		f = f.squeeze()

	return f
 
def get_bottom_data(data):
    '''
    Return a 3D array (time,y,x) with data at the bottom most cell of each column. 
    Data is already masked where there is no ocean.
    '''
    NT,NK,NY,NX = data.shape
    bottom_data = np.zeros((NT,NY,NX))
    for j in range(NY):
      for i in range(NX):
          # ice and bottom are static, so we can get indies at t = 0 since they dont change.
          ind = np.nonzero(data[0,:,j,i]!= data.fill_value)[0]
          if ind.any():
             bottom_data[:,j,i] = data[:,ind[-1],j,i]
          else:
             bottom_data[:,j,i] = np.nan

    return np.ma.masked_invalid(bottom_data)

def saveXY(var,varname):
   '''
   Save 2D (x,y) or 3D array (time,x,y) into the netCDF file.
   '''
   #if len(var.shape) == 2:
   ncwrite(name,varname,var) # needs to tanspose it
   #else:
     #var = var.transpose(0, 2, 1)
   #  ncwrite(name,varname,var)
   return

def mean_tracer(area,h,var,varname,units):
   '''
   Compute tracer (temp,salt or any other tracer) averaged over the ocean volume.
   '''
   area = np.resize(area,h.shape)
   tracer = np.zeros(h.shape[0])
   for t in range(len(tracer)):
       vol = (area[t,:] * h[t,:])
       tracer[t] = (var[t,:] * vol).sum() /vol.sum()
       print(varname+' ('+units+') at t='+str(t)+' is:'+str(tracer[t]))

   ncwrite(name,varname,tracer)
   return

def total_melt_flux(mass):
    '''
    Compute the total mass flux of freshwater across the ice-oceaninterface, positive 
    for melting and negative for freezing.
    '''    
    total_mass = np.zeros(mass.shape[0])
    for t in range(len(total_mass)):
        total_mass[t] = mass[t,:].sum()
        print('Total melt flux (kg/s) at t='+str(t)+' is:'+str(total_mass[t]))

    ncwrite(name,'totalMeltFlux',total_mass)
    return

def mean_melt_rate(melt,area):
    '''
    Compute the melt rate, positive for melting negative for freezing, averaged over
    the ice shelf base (area>0).
    '''
    melt = melt/(3600*24.*365) # in m/s
    total_melt = np.zeros(melt.shape[0])
    for t in range(len(total_melt)):
        total_melt[t] = (melt[t,:] * area).sum() /area.sum()
        print('Averaged melt (m/s) at t='+str(t)+' is:'+str(total_melt[t]))

    ncwrite(name,'meanMeltRate',total_melt)
    return

def total_volume(area,h):
    area = area=np.resize(area,h.shape) 
    vol = np.zeros(h.shape[0])
    for t in range(len(vol)):
        vol[t] = (area[t,:] * h[t,:]).sum()
        print('Total vol. (m^3) at t='+str(t)+' is:'+str(vol[t])) 
        
    ncwrite(name,'totalOceanVolume',vol)
    return

def ncwrite(name,var_name,data):
    '''
    Write data (with var_name) into netCDF file (name.nc).
    '''
    file = Dataset(name+'.nc','r+')
    file.variables[var_name][:]=data[:]
    file.close()
    print('*** Variable '+var_name+' was written successfully. \n')
    return
 
def get_rho(S,T):
    '''
    Compute density with a linear EoS and using isomip+ coeffs. 
    '''
    rho0=1027.51; Sref=34.2; Tref=-1.0
    alpha = 3.733e-5; beta = 7.843e-4
    rho = rho0 * (1 - alpha * (T-Tref) + beta * (S-Sref))
    return rho

def melt4gamma(name,area,ice_base,depth):
   '''
   Compute (and save in a text file) melting over the area where ice base > 300 m and over final six months.
   '''

   print('WARNING: have you checked if exp. has reached steady state?')
   # read melt over last 6 months
   melt = Dataset('ocean_month.nc').variables['melt'][-6::,:,:]
   # total area under ice shelf (exclude grounded ice)
   # mask area where ice base <= 300 
   area = np.ma.masked_where(ice_base<=300.,area)
   print('Ice shelf area below 300 m and excluding grouded ice is (m^2): '+str(area.sum()))
   total_melt = np.zeros(melt.shape[0])
   # compute mean melt at each time
   for t in range(len(total_melt)):
       total_melt[t] = ((area * melt[t,:,:]).sum()/area.sum())
       print('Averaged melt at t='+str(t)+' is:'+str(total_melt[t]))
   
   print('Time-averaged <melt> is: '+str(total_melt.mean()))
   return


def mask_ocean(data,area):
   """
   Mask open ocean. Works with 2D or 3D arrays.
   """
   if len(data.shape) == 2: # 2D array
     data = np.ma.masked_where(area==0,data)   

   else: # 3D array
     NZ,NY,NX = data.shape
     area=np.resize(area,(NZ,NY,NX))
     data = np.ma.masked_where(area==0,data)
   
   return  data

def mask_grounded_ice(data,depth,base):
   """
   Mask regions where the ice shelf is grounded (base>=depth). Works with 2D or 3D arrays.
   """
   if len(data.shape) == 2: # 2D array
      data = np.ma.masked_where(base+1.0>=depth, data) # need +1 here      
   else: # 3D array
      NZ,NY,NX = data.shape
      base = np.resize(base,(NZ,NY,NX))
      depth = np.resize(depth,(NZ,NY,NX))
      data = np.ma.masked_where(base+1.0>=depth, data) # need +1 here

   return data

def create_ncfile(exp_name, ocean_time, type): # may add exp_type
   """
   Creates a netcdf file with the fields required for the isomip+ experiments. Different fields are generated based on the type of experiment that is being analyzed (Ocean0, Ocean1 etc).
   """

   # open a new netCDF file for writing.
   ncfile = Dataset(exp_name+'.nc','w',format='NETCDF4')
   # dimensions
   nx = 240 ; ny = 40 ; nz = 144 
   # create dimensions.
   #ncfile.createDimension('nTime', None)
   ncfile.createDimension('nTime', len(ocean_time))
   ncfile.createDimension('nx',nx)
   ncfile.createDimension('ny',ny)
   ncfile.createDimension('nz',nz)

   # create variables, assign units and provide decription
   x = ncfile.createVariable('x',np.dtype('float32').char,('nx',))
   x.units = 'm'; x.description = 'x location of cell centers'
   y = ncfile.createVariable('y',np.dtype('float32').char,('ny',))
   y.units = 'm'; y.description = 'y location of cell centers'
   z = ncfile.createVariable('z',np.dtype('float32').char,('nz',))
   z.units = 'm'; z.description = 'z location of cell centers'
   time = ncfile.createVariable('time',np.dtype('float32').char,('nTime',))
   time.units = 's'; time.description = 'time since start of simulation'
   meanMeltRate = ncfile.createVariable('meanMeltRate',np.dtype('float32').char,('nTime',))
   meanMeltRate.units = 'm/s'; meanMeltRate.description = 'mean melt rate averaged over area of floating ice, positive for melting'
   totalMeltFlux = ncfile.createVariable('totalMeltFlux',np.dtype('float32').char,('nTime',))
   totalMeltFlux.units = 'kg/s'; totalMeltFlux.description = 'total flux of melt water summed over area of floating ice, positive for melting'
   totalOceanVolume = ncfile.createVariable('totalOceanVolume',np.dtype('float32').char,('nTime',))   
   totalOceanVolume.units = 'm^3'; totalOceanVolume.description = 'total volume of ocean'
   meanTemperature = ncfile.createVariable('meanTemperature',np.dtype('float32').char,('nTime',))
   meanTemperature.units = 'deg C'; meanTemperature.description = 'the potential temperature averaged over the ocean volume'
   meanSalinity = ncfile.createVariable('meanSalinity',np.dtype('float32').char,('nTime',))
   meanSalinity.units = 'PSU'; meanSalinity.description = 'the salinity averaged over the ocean volume'
   if (type == 'ocean3' or type == 'ocean4'):
     iceDraft = ncfile.createVariable('iceDraft',np.dtype('float32').char,('nTime','ny','nx'))  
     bathymetry = ncfile.createVariable('bathymetry',np.dtype('float32').char,('nTime','ny','nx'))
   else:
     iceDraft = ncfile.createVariable('iceDraft',np.dtype('float32').char,('ny','nx')) 
     bathymetry = ncfile.createVariable('bathymetry',np.dtype('float32').char,('ny','nx'))  
   iceDraft.units = 'm'; iceDraft.description = 'elevation of the ice-ocean interface'
   bathymetry.units = 'm'; bathymetry.description = 'elevation of the bathymetry'
   meltRate = ncfile.createVariable('meltRate',np.dtype('float32').char,('nTime','ny','nx'))
   meltRate.units = 'm/s'; meltRate.description = 'melt rate, positive for melting'
   frictionVelocity = ncfile.createVariable('frictionVelocity',np.dtype('float32').char,('nTime','ny','nx'))
   frictionVelocity.units = 'm/s'; frictionVelocity.description = 'friction velocity u* used in melt calculations'
   thermalDriving = ncfile.createVariable('thermalDriving',np.dtype('float32').char,('nTime','ny','nx'))
   thermalDriving.units = 'deg C'; thermalDriving.description = 'thermal driving used in the melt calculation'
   halineDriving = ncfile.createVariable('halineDriving',np.dtype('float32').char,('nTime','ny','nx'))
   halineDriving.units = 'PSU';  halineDriving.description = 'haline driving used in the melt calculation'
   uBoundaryLayer = ncfile.createVariable('uBoundaryLayer',np.dtype('float32').char,('nTime','ny','nx'))
   uBoundaryLayer.units = 'm/s'; uBoundaryLayer.description = 'x-velocity in the boundary layer used to compute u*'
   vBoundaryLayer = ncfile.createVariable('vBoundaryLayer',np.dtype('float32').char,('nTime','ny','nx'))
   vBoundaryLayer.units = 'm/s'; vBoundaryLayer.description = 'y-velocity in the boundary layer used to compute u*'
   barotropicStreamfunction = ncfile.createVariable('barotropicStreamfunction',np.dtype('float32').char,('nTime','ny','nx'))
   barotropicStreamfunction.units = 'm^3/s'; barotropicStreamfunction.description = 'barotropic streamfunction'
   overturningStreamfunction = ncfile.createVariable('overturningStreamfunction',np.dtype('float32').char,('nTime','nz','nx'))
   overturningStreamfunction.units = 'm^3/s'; overturningStreamfunction.description = 'overturning (meridional) streamfunction'
   bottomTemperature = ncfile.createVariable('bottomTemperature',np.dtype('float32').char,('nTime','ny','nx'))   
   bottomTemperature.units = 'deg C'; bottomTemperature.description = 'temperature in the bottom grid cell of each ocean column'
   bottomSalinity = ncfile.createVariable('bottomSalinity',np.dtype('float32').char,('nTime','ny','nx'))   
   bottomSalinity.units = 'PSU'; bottomSalinity.description = 'salinity in the bottom grid cell of each ocean column'
   temperatureXZ = ncfile.createVariable('temperatureXZ',np.dtype('float32').char,('nTime','nz','nx'))
   temperatureXZ.units = 'deg C'; temperatureXZ.description = 'temperature slice in x-z plane through the center of the domain (y = 40 km)'
   salinityXZ = ncfile.createVariable('salinityXZ',np.dtype('float32').char,('nTime','nz','nx'))
   salinityXZ.units = 'PSU'; salinityXZ.description = 'salinity slice in x-z plane through the center of the domain (y = 40 km)'
   temperatureYZ = ncfile.createVariable('temperatureYZ',np.dtype('float32').char,('nTime','nz','ny'))
   temperatureYZ.units = 'deg C'; temperatureYZ.description = 'temperature slice in y-z plane through x = 500 km'
   salinityYZ = ncfile.createVariable('salinityYZ',np.dtype('float32').char,('nTime','nz','ny'))
   salinityYZ.units = 'PSU'; salinityYZ.description = 'salinity slice in y-z plane through x = 500 km'

   # write data to coordinate vars.
   x[:] = np.arange(321,800,2)*1.0e3
   y[:] = np.arange(1,80,2)*1.0e3
   z[:] = -np.arange(2.5,720,5) 
   # time since start of simulation
   time[0] = 0.
   time[1::] = ocean_time[0:-1] * 3600 * 24 # in secs
   #time[:] = np.array([0, 2678400, 5097600, 7776000, 1.0368e+07, 1.30464e+07, 1.56384e+07,1.83168e+07, 2.09952e+07, 2.35872e+07, 2.62656e+07, 2.88576e+07])

   # assign zero to variables 
   #meanMeltRate[:] = np.zeros(len(time))
   #totalMeltFlux[:] = np.zeros(len(time))
   # close the file.
   ncfile.close()
   print ('*** SUCCESS creating '+exp_name+'.nc!')
   return

def write_ncfile(exp_name, data): # may add exp_type
   """
   This function writes the fields required for the isomip+ experiments into a pre-exiting netcdf file (exp_name). Different fields (given by the structure data) are added based on the type experiment that is being analyzed (Ocean0, Ocean1 etc).
   """
   print ('*** SUCCESS writing data into '+exp_name+'.nc!')
   return

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
