#!/usr/bin/env python

# plot some diagnostics that I find useful for the isomip setup
# Gustavo Marques, Aug. 2016

try: import argparse
except: raise MyError('This version of python is not new enough. python 2.7 or newer is required.')
try: from netCDF4 import MFDataset, Dataset
except: raise MyError('Unable to import netCDF4 module. Check your PYTHONPATH.\n'
          +'Perhaps try:\n   module load python_netcdf4')
try: import numpy as np
except: raise MyError('Unable to import numpy module. Check your PYTHONPATH.\n'
          +'Perhaps try:\n   module load python_numpy')
try: import matplotlib.pyplot as plt
except: raise MyError('Unable to import matplotlib.pyplot module. Check your PYTHONPATH.\n'
          +'Perhaps try:\n   module load python_matplotlib')
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
      plot useful diagnostics (1- and 2-dimensional data) useful during the 
      development of the the isomip test case.
      ''',
  epilog='Written by Gustavo Marques, Aug. 2016.')

  parser.add_argument('-t','--time', type=int, default=0,
      help='''Time indice (integer). Default is 0.''')

  parser.add_argument('-d','--diagnostic', type=str,default='min_ocean_thickness',
      help='''Specify the desired disgnostic. Valid form are:  min_ocean_thickness;
              total_volume;''')

  optCmdLineArgs = parser.parse_args()
  createPlot(optCmdLineArgs.diagnostic, optCmdLineArgs)


def createPlot(diagnostic, args):
  """
  Generates a plot based on the file/type of diag. specified
  """

  if diagnostic == 'min_ocean_thickness':
     ncfile = open_ncfile('ocean_geometry.nc')
     D = ncfile.variables['D'][:]
     ncfile = open_ncfile('ISOMIP_IC.nc')
     ice_draft = ncfile.variables['ave_ssh'][0,:,:]
     X = ncfile.variables['lonh'][:]; Y = ncfile.variables['lath'][:]
     thick = ice_draft + D
     print 'Min/Max thickness:',thick.min(), thick.max()
     plt.figure()
     plt.contourf(X,Y,thick,np.arange(0.1,20,0.1))
     plt.colorbar()
     plt.xlabel('Long. (km)')
     plt.ylabel('Lat. (km)')
     plt.title('Total (initial) water column thickness (m)')
     plt.show()
 

def open_ncfile(fileName):
  """
  Open netCDF file and return the netcdf object.
  """
    # Open netcdf file
  try: data = MFDataset(fileName, 'r', aggdim='time')
  except:
    try: data = Dataset(fileName, 'r')
    except:
      if os.path.isfile(fileName): raise MyError('There was a problem opening "'+fileName+'".')
      raise MyError('Could not find file "'+fileName+'".')

  return data

def readVariableFromFile(fileName, variableName,dynamic=True, ignoreCoords=False, alternativeNames=None):
  """
  Open netCDF file, find and read the variable meta-information and return both
  the netcdf object and variable object
  """

  if dynamic:
    # split string
    fileName,time = splitFile(fileName)

  # Open netcdf file
  try: data = MFDataset(fileName, 'r', aggdim='time')
  except:
    try: rg = Dataset(fileName, 'r')
    except:
      if os.path.isfile(fileName): raise MyError('There was a problem opening "'+fileName+'".')
      raise MyError('Could not find file "'+fileName+'".')

  if dynamic:
    # If no variable is specified, summarize the file contents and exit
    if not variableName:
       print 'No variable name specified! Specify a varible from the following summary of "'\
          +fileName+'":\n'
       summarizeFile(rg)
       exit(0)

    # Check that the variable is in the file (allowing for case mismatch)
    for v in rg.variables:
      if variableName.lower() == v.lower(): variableName=v ; break
    if not variableName in rg.variables:
      if alternativeNames==None:
        print 'Known variables in file: '+''.join( (str(v)+', ' for v in rg.variables) )
        raise MyError('Did not find "'+variableName+'" in file "'+fileName+'".')
      else:
        for v in alternativeNames:
          if v in rg.variables: variableName=v ; break

    return rg, rg.variables[variableName][time,:]
  else: 
    return rg

# Generate a succinct summary of the netcdf file contents
def summarizeFile(rg):
  dims = rg.dimensions; vars = rg.variables
  print 'Dimensions:'
  for dim in dims:
    oString = ' '+dim+' ['+str(len( dims[dim] ))+']'
    if dim in vars:
      n = len( dims[dim] ); obj = rg.variables[dim]
      if n>5: oString += ' = '+str(obj[0])+'...'+str(obj[n-1])
      else: oString += ' = '+str(obj[:])
      if 'long_name' in obj.ncattrs(): oString += ' "'+obj.long_name+'"'
      if 'units' in obj.ncattrs(): oString += ' ('+obj.units+')'
    print oString
  print; print 'Variables:'
  for var in vars:
    if var in dims: continue # skip listing dimensions as variables
    oString = ' '+var+' [ '; dString = ''
    obj = vars[var]; varDims = obj.dimensions
    for dim in varDims:
      if len(dString)>0: dString += ', '
      dString += dim+'['+str(len( dims[dim] ))+']'
    oString += dString+' ]'
    if 'long_name' in obj.ncattrs(): oString += ' "'+obj.long_name+'"'
    if 'units' in obj.ncattrs(): oString += ' ('+obj.units+')'
    print oString

def splitFile(string):
  """
  Split a string in form of "diagnostic,time_indice" into one string and one integer
  Valid forms are "diagnostic", "diagnostic,1[2,3,4...]" or "diagnostic,-1[-2,-3...]"
  """
  ssplit = string.split(',')
  fname = ssplit[0]
  if not ssplit[1]:
    print('Time indice has not been specified, using the default (1) value... \n')
    time = 1
  else:
    time = int(ssplit[1])
  return fname, time


# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()

