#!/usr/bin/env python

# plot all the terms used to compute melting for a user specified point (i,j)
# Gustavo Marques, Sep. 2016

import argparse
from netCDF4 import MFDataset, Dataset
import numpy as np
import warnings
import os
from datetime import datetime

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Apply a correctoion on the mask under grounded ice.
      ''',
  epilog='Written by Gustavo Marques, Apr. 2019.')

  parser.add_argument('-file_in', type=str, default='Ocean0_COM_MOM6-LAYER.nc',
      help='''Name of input file to be processed. Default is Ocean0_COM_MOM6-LAYER.nc.''')

  parser.add_argument('-file_out', type=str, default='IceOcean1ra_COM_ocean_MOM6_CISM_NCAR_DIVA.nc',
      help='''Name of output file to be written. Default is IceOcean1ra_COM_ocean_MOM6_CISM_NCAR_DIVA.nc.''')

  optCmdLineArgs = parser.parse_args()
  applymask(optCmdLineArgs)

def applymask(arg):
  # load data
  data_in = Dataset(arg.file_in)
  data_out = Dataset(arg.file_out, 'r+')

  # dimensions and time series
  data_out.variables['x'][:] = data_in.variables['x'][160:400]
  data_out.variables['y'][:] = data_in.variables['y'][:]
  data_out.variables['z'][:] = data_in.variables['z'][:]
  data_out.variables['time'][:] = data_in.variables['time'][:]
  #data_out.variables['time'][1:] = data_in.variables['time'][1::] - data_in.variables['time'][1] + 2678400
  data_out.variables['meanMeltRate'][:] = data_in.variables['meanMeltRate'][:]
  data_out.variables['totalMeltFlux'][:] = data_in.variables['totalMeltFlux'][:]
  data_out.variables['totalOceanVolume'][:] = data_in.variables['totalOceanVolume'][:]
  data_out.variables['meanTemperature'][:] = data_in.variables['meanTemperature'][:]
  data_out.variables['meanSalinity'][:] = data_in.variables['meanSalinity'][:]

  # list of variables to be corrected
  variables = ['temperatureYZ','salinityYZ']
  for v in data_in.variables:
    if len(data_in.variables[v][:].shape) > 1:
      print('Processing variable ',v)
      if v in variables:
        tmp = data_in.variables[v][:]
        data_out.variables[v][:] = tmp
      else:
        tmp = data_in.variables[v][:]
        new_tmp = tmp[:,:,160:400]
        data_out.variables[v][:] = new_tmp

    else:
      print('Skipping.. \n')

  # mask  bathymetry
  print('Masking bathymetry...\n')
  bat =  data_in.variables['bathymetry'][:,:,160:400]
  bottomSalinity =  data_in.variables['bottomSalinity'][:,:,160:400]
  tm, jm, im = bat.shape
  for t in range(tm):
    bat[t,:] = np.ma.masked_where(bottomSalinity[t,:].mask == True, bat[t,:])

  data_out.variables['bathymetry'][:,:,:] = bat[:,:,:]

  # : mask zero values
  print('Masking iceDraft...\n')
  iceDraft = data_in.variables['iceDraft'][:,:,160:400]
  for t in range(tm):
    iceDraft[t,:] = np.ma.masked_where(bottomSalinity[t,:].mask == True, iceDraft[t,:])

  data_out.variables['iceDraft'][:,:,:] = iceDraft[:,:,:]
  # change attributes
  data_out.Author = 'Gustavo Marques (gmarques@ucar.edu)'
  data_out.Created = datetime.now().isoformat()

  data_out.sync()
  data_out.close()
  data_in.close()
  print('Done!')

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()

