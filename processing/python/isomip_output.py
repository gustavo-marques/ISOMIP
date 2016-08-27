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

  parser.add_argument('-t','-type', type=str, default='ocean0',
      help='''The type pf experiment that will computed (ocean0, ocean1 etc). Default is ocean0.''')

  parser.add_argument('-n','-name', type=str, default='Ocean0_COM_MOM6-LAYER',
      help='''The name of the experiment following the ISOMIP+ definition (expt_COM_model). Default is Ocean0_COM_MOM6-LAYER.''')

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
   shelf_area = Dataset('INPUT/Ocean1_3D.nc').variables['area'][:]
   # base of STATIC ice shelf, which is ssh(t=0); make it positive
   ice_base = -Dataset('ISOMIP_IC.nc').variables['ave_ssh'][0,:,:]
   ice_base[ice_base<1e-5] = 0.0
   #exclude grounded ice and open ocean
   shelf_area = mask_grounded_ice(shelf_area,depth,ice_base) 
   shelf_area = mask_ocean(shelf_area,shelf_area)
   
   #lonh=Dataset('ocean_geometry.nc').variables['lonh'][:]
   #lonq=Dataset('ocean_geometry.nc').variables['lonq'][:]  
   #lath=Dataset('ocean_geometry.nc').variables['lath'][:]
   #latq=Dataset('ocean_geometry.nc').variables['latq'][:]
   #lonqs, latqs = np.meshgrid(lonq,latq)
   #lons, lats = np.meshgrid(lonh,lath)
   
   # create dictionary
   # not for now 
   #data = {'name':name, 'depth':depth, 'area':area, }

   if args.mean4gamma:
      print("Computing mean melt rate to calibrate gamma...")
      melt4gamma(name,shelf_area,ice_base,depth)
   
   # create ncfile and zero fields. Each function below will corresponding values
   create_ncfile(name)

   # load additional variables
   ocean_area = Dataset('ocean_geometry.nc').variables['Ah'][:] 
   ocean_area = mask_grounded_ice(ocean_area,depth,ice_base)
   h = Dataset('ocean_month.nc').variables['h'][:]
   h[h<1.0e-5] = 0.0
   melt = Dataset('ocean_month.nc').variables['melt'][:]
   # under shelf melt is always > 0. Therefore we can mask melt=0.
   melt = np.ma.masked_where(melt==0,melt)
   mass_flux = Dataset('ocean_month.nc').variables['mass_flux'][:]

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
   print('Done!')
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
   print('Ice shelf area below 300 m and excluding grouded ice is (m^2):',area.sum())
   total_melt = np.zeros(melt.shape[0])
   # compute mean melt at each time
   for t in range(len(total_melt)):
       total_melt[t] = ((area * melt[t,:,:]).sum()/area.sum())
       print('Averaged melt at t='+str(t)+' is:'+str(total_melt[t]))
   
   print('Time-averaged <melt> is: '+str(total_melt.mean()))
   return


def mask_ocean(data,area,D=False):
   """
   Mask open ocean. Works with 2D or 3D arrays (default 2D,   set D=true for 3D).
   """
   if not D:
     data = np.ma.masked_where(area==0,data)   

   else:
     NZ,NY,NX = data.shape
     area=np.resize(area,(NZ,NY,NX))
     data = np.ma.masked_where(area==0,data)
   
   return  data

def mask_grounded_ice(data,depth,base,D=False):
   """
   Mask regions where the ice shelf is grounded (base>=depth). Works with 2D or 3D arrays (default 2D,   set D=true for 3D).
   """
   if not D:
      data = np.ma.masked_where(base+1.0>=depth, data) # need +1 here      
   else:
      NZ,NY,NX = data.shape
      base = np.resize(base,(NZ,NY,NX))
      depth = np.resize(depth,(NZ,NY,NX))
      data = np.ma.masked_where(base+1.0>=depth, data) # need +1 here

   return data

def create_ncfile(exp_name): # may add exp_type
   """
   This function creates a netcdf file with the fields required for the isomip+ experiments. Different fields are generated based on the type experiment that is being analyzed (Ocean0, Ocean1 etc).
   """

   # open a new netCDF file for writing.
   ncfile = Dataset(exp_name+'.nc','w',format='NETCDF4')
   # dimensions
   nx = 240 ; ny = 40 ; nz = 144 
   # create dimensions.
   ncfile.createDimension('nTime', None)
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
   iceDraft = ncfile.createVariable('iceDraft',np.dtype('float32').char,('ny','nx'))   
   iceDraft.units = 'm'; iceDraft.description = 'elevation of the ice-ocean interface'
   bathymetry = ncfile.createVariable('bathymetry',np.dtype('float32').char,('ny','nx'))     
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
   time[:] = np.array([0, 2678400, 5097600, 7776000, 1.0368e+07, 1.30464e+07, 1.56384e+07,                      1.83168e+07, 2.09952e+07, 2.35872e+07, 2.62656e+07, 2.88576e+07])

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
