from netCDF4 import Dataset
from numpy import arange, dtype, array, zeros
"""
These functions are used to create/populate netcdf files.

Gustavo Marques, Aug. 2016.

"""
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
   x = ncfile.createVariable('x',dtype('float32').char,('nx',))
   x.units = 'm'; x.description = 'x location of cell centers'
   y = ncfile.createVariable('y',dtype('float32').char,('ny',))
   y.units = 'm'; y.description = 'y location of cell centers'
   z = ncfile.createVariable('z',dtype('float32').char,('nz',))
   z.units = 'm'; z.description = 'z location of cell centers'
   time = ncfile.createVariable('time',dtype('float32').char,('nTime',))
   time.units = 's'; time.description = 'time since start of simulation'
   meanMeltRate = ncfile.createVariable('meanMeltRate',dtype('float32').char,('nTime',))
   meanMeltRate.units = 'm/s'; meanMeltRate.description = 'mean melt rate averaged over area of floating ice, positive for melting'
   totalMeltFlux = ncfile.createVariable('totalMeltFlux',dtype('float32').char,('nTime',))
   totalMeltFlux.units = 'kg/s'; totalMeltFlux.description = 'total flux of melt water summed over area of floating ice, positive for melting'
   totalOceanVolume = ncfile.createVariable('totalOceanVolume',dtype('float32').char,('nTime',))   
   totalOceanVolume.units = 'm^3'; totalOceanVolume.description = 'total volume of ocean'
   meanTemperature = ncfile.createVariable('meanTemperature',dtype('float32').char,('nTime',))
   meanTemperature.units = 'deg C'; meanTemperature.description = 'the potential temperature averaged over the ocean volume'
   meanSalinity = ncfile.createVariable('meanSalinity',dtype('float32').char,('nTime',))
   meanSalinity.units = 'PSU'; meanSalinity.description = 'the salinity averaged over the ocean volume'
   iceDraft = ncfile.createVariable('iceDraft',dtype('float32').char,('ny','nx'))   
   iceDraft.units = 'm'; iceDraft.description = 'elevation of the ice-ocean interface'
   bathymetry = ncfile.createVariable('bathymetry',dtype('float32').char,('ny','nx'))     
   bathymetry.units = 'm'; bathymetry.description = 'elevation of the bathymetry'
   meltRate = ncfile.createVariable('meltRate',dtype('float32').char,('nTime','ny','nx'))
   meltRate.units = 'm/s'; meltRate.description = 'melt rate, positive for melting'
   frictionVelocity = ncfile.createVariable('frictionVelocity',dtype('float32').char,('nTime','ny','nx'))
   frictionVelocity.units = 'm/s'; frictionVelocity.description = 'friction velocity u* used in melt calculations'
   thermalDriving = ncfile.createVariable('thermalDriving',dtype('float32').char,('nTime','ny','nx'))
   thermalDriving.units = 'deg C'; thermalDriving.description = 'thermal driving used in the melt calculation'
   halineDriving = ncfile.createVariable('halineDriving',dtype('float32').char,('nTime','ny','nx'))
   halineDriving.units = 'PSU';  halineDriving.description = 'haline driving used in the melt calculation'
   uBoundaryLayer = ncfile.createVariable('uBoundaryLayer',dtype('float32').char,('nTime','ny','nx'))
   uBoundaryLayer.units = 'm/s'; uBoundaryLayer.description = 'x-velocity in the boundary layer used to compute u*'
   vBoundaryLayer = ncfile.createVariable('vBoundaryLayer',dtype('float32').char,('nTime','ny','nx'))
   vBoundaryLayer.units = 'm/s'; vBoundaryLayer.description = 'y-velocity in the boundary layer used to compute u*'
   barotropicStreamfunction = ncfile.createVariable('barotropicStreamfunction',dtype('float32').char,('nTime','ny','nx'))
   barotropicStreamfunction.units = 'm^3/s'; barotropicStreamfunction.description = 'barotropic streamfunction'
   overturningStreamfunction = ncfile.createVariable('overturningStreamfunction',dtype('float32').char,('nTime','nz','nx'))
   overturningStreamfunction.units = 'm^3/s'; overturningStreamfunction.description = 'overturning (meridional) streamfunction'
   bottomTemperature = ncfile.createVariable('bottomTemperature',dtype('float32').char,('nTime','ny','nx'))   
   bottomTemperature.units = 'deg C'; bottomTemperature.description = 'temperature in the bottom grid cell of each ocean column'
   bottomSalinity = ncfile.createVariable('bottomSalinity',dtype('float32').char,('nTime','ny','nx'))   
   bottomSalinity.units = 'PSU'; bottomSalinity.description = 'salinity in the bottom grid cell of each ocean column'
   temperatureXZ = ncfile.createVariable('temperatureXZ',dtype('float32').char,('nTime','nz','nx'))
   temperatureXZ.units = 'deg C'; temperatureXZ.description = 'temperature slice in x-z plane through the center of the domain (y = 40 km)'
   salinityXZ = ncfile.createVariable('salinityXZ',dtype('float32').char,('nTime','nz','nx'))
   salinityXZ.units = 'PSU'; salinityXZ.description = 'salinity slice in x-z plane through the center of the domain (y = 40 km)'
   temperatureYZ = ncfile.createVariable('temperatureYZ',dtype('float32').char,('nTime','nz','ny'))
   temperatureYZ.units = 'deg C'; temperatureYZ.description = 'temperature slice in y-z plane through x = 500 km'
   salinityYZ = ncfile.createVariable('salinityYZ',dtype('float32').char,('nTime','nz','ny'))
   salinityYZ.units = 'PSU'; salinityYZ.description = 'salinity slice in y-z plane through x = 500 km'

   # write data to coordinate vars.
   x[:] = arange(321,800,2)*1.0e3
   y[:] = arange(1,80,2)*1.0e3
   z[:] = -arange(2.5,720,5) 
   # time since start of simulation
   time[:] = array([0, 2678400, 5097600, 7776000, 1.0368e+07, 1.30464e+07, 1.56384e+07,                      1.83168e+07, 2.09952e+07, 2.35872e+07, 2.62656e+07, 2.88576e+07])

   # assign zero to variables 
   #meanMeltRate[:] = zeros(len(time))
   #totalMeltFlux[:] = zeros(len(time))
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
