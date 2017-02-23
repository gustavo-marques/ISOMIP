#!/usr/bin/env python
import numpy
from optparse import OptionParser
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from netCDF4 import Dataset

'''
This script can be used to plot ISOMIP+ or MISOMIP1 ocean data
in the standard NetCDF format specified for these MIPs.  The
script produces images for each time index of each field, which
can be used to produce movies if desired.  Image files are skipped if they
exist.  (To regenerate images, the user must first delete the
existing images.)
Usage: plotMISOIP1OceanData.py <in_file.nc> <out_dir>
'''
def plotOcean(fileName, outFolder):
  def plot(varName, label, cmap, scale=None):
    # the file name is the variable followed by the zero-padded time intex
    imageFileName = '%s/%s_%04i.png'%(outFolder, varName, timeIndex)
    if(os.path.exists(imageFileName)):
      # the image exists so we're going to save time and not replot it
      return
  
    # get the variable from the netCDF file
    try:
      var = ncFile.variables[varName]
    except KeyError:
      return
  
    # the axes are 'xy', 'xz', or 'yz'
    axes = '%s%s'%(var.dimensions[2][1],var.dimensions[1][1])
  
    # get the extent based on the axes
    extent = extents[axes]
  
    # aspect ratio
    if(axes == 'xy'):
      # pixels are 1:1
      aspectRatio = None
    else:
      # stretch the axes to fill the plot area
      aspectRatio = 'auto'
  
    field = ncFile.variables[varName][timeIndex,:,:]
    
    missingValue = 9.9692099683868690e36
    if(numpy.ma.amax(field) == missingValue):
      # the array didnt' get masked properly, so do it manually
      field = numpy.ma.masked_array(field, mask=(field == missingValue), dtype=float)
    else:
      # convert from float32 to float64 to avoid overflow/underflow problems
      if hasattr(field, 'mask'):
        field = numpy.ma.masked_array(field, mask=field.mask, dtype=float)
      else:
        field = numpy.array(field, dtype=float)
    

    # scale the variable if a scale is given
    if scale is not None:
      field *= scale
  
    # get the plotting limits from the dictionary of limits we created
    (lower, upper) = limits[varName]
  
    # make a figure
    plt.figure(1, figsize=[16,9], dpi=100, facecolor='w')
  
    # activate the specified colorbar and set the background color
    cmap = plt.get_cmap(cmap)
    cmap.set_bad(backgroundColor)
  
    # clear the figure from the last plot
    plt.clf()
  
    # plot the data as an image
    plt.imshow(field, extent=extent, cmap=cmap, vmin=lower, vmax=upper, 
               aspect=aspectRatio, interpolation='nearest')
  
    plt.colorbar()
    plt.title(label)
  
    if(axes == 'xy'):
      # y axis will be upside down in imshow, which we don't want for xy
      plt.gca().invert_yaxis()
      plt.xlabel('x (km)')
      plt.ylabel('y (km)')
    elif(axes == 'xz'):
      # upside-down y axis is okay
      plt.xlabel('x (km)')
      plt.ylabel('z (m)')
    else:
      # upside-down y axis is okay
      plt.xlabel('y (km)')
      plt.ylabel('z (m)')
    
    if (axes in ['xz','yz']) and (z[0] < z[-1]):
      # flip the y axis
      plt.gca().invert_yaxis()
  
    # save the figure as an image
    plt.tight_layout()
    plt.draw()
    plt.savefig(imageFileName, dpi=100)
    plt.close()
  
  def makeFerretColormap():
    red = numpy.array([[0,0.6],
                       [0.15,1],
                       [0.35,1],
                       [0.65,0],
                       [0.8,0],
                       [1,0.75]])
    
    green = numpy.array([[0,0],
                         [0.1,0],
                         [0.35,1],
                         [1,0]])
    
    
    blue = numpy.array([[0,0],
                       [0.5,0],
                       [0.9,0.9],
                       [1,0.9]])
    
    colorCount = 21
    ferretColorList = numpy.ones((colorCount,4),float)
    ferretColorList[:,0] = numpy.interp(numpy.linspace(0,1,colorCount),red[:,0],red[:,1])
    ferretColorList[:,1] = numpy.interp(numpy.linspace(0,1,colorCount),green[:,0],green[:,1])
    ferretColorList[:,2] = numpy.interp(numpy.linspace(0,1,colorCount),blue[:,0],blue[:,1])
    ferretColorList = ferretColorList[::-1,:]
    
    cmap = colors.LinearSegmentedColormap.from_list('ferret',ferretColorList,N=255)
    return cmap
  
  
  
  try:
    os.makedirs(outFolder)
  except OSError:
    pass

  sPerYr = 365.*24.*60.*60.
  
  # a rainbow colormap based on the default map in the application ferret
  cmap = makeFerretColormap()
  
  # get the filename without the path
  baseName = os.path.basename(fileName)
  # pull out the first part of the file name as the experiment name,
  # following the ISOMIP+ and MISOMIP1 filenaming requirements
  experiment = baseName.split('_')[0]
  
  # make sure this is a valid experiment (i.e. that the file has been
  # named correctly)
  if(experiment not in ['Ocean0','Ocean1','Ocean2','Ocean3','Ocean4',
                        'IceOcean1','IceOcean2']):
    print "Unknown experiment", experiment
    exit(1)
  
  # open the netCDF file with the ISOMIP+ or MISOMIP1 ocean data
  ncFile = Dataset(fileName,'r')
  
  # convert x and y to km
  try:
    x = 1e-3*ncFile.variables['x'][:]
  except KeyError:
    # just use the index (axes will be labeled incorrectly as km)
    nx = len(ncFile.dimensions['nx'])
    x = numpy.arange(nx)
  try:
    y = 1e-3*ncFile.variables['y'][:]
  except KeyError:
    # just use the index (axes will be labeled incorrectly as km)
    ny = len(ncFile.dimensions['ny'])
    y = numpy.arange(ny)

  try:
    # leave z in m
    z = ncFile.variables['z'][:]
  except KeyError:
    # just use the index (axes will be labeled incorrectly as km)
    nz = len(ncFile.dimensions['nz'])
    z = numpy.arange(nz) 

  try:
    times = numpy.array(ncFile.variables['time'][:], dtype=float)/sPerYr
  except KeyError:
    # fill in missing time array with best guess
    nTime = len(ncFile.dimensions['nTime'])
    times = numpy.arange(nTime)/12.0
  print len(times)
  
  
  # the extents of the different plots for use in imshow
  extents = {}
  # the y extent is max then min because the y axis then gets flipped 
  # (imshow is weird that way)
  extents['xy'] = [numpy.amin(x),numpy.amax(x),numpy.amax(y),numpy.amin(y)]
  if(z[0] > z[-1]):
    # we want the z-axis to be flipped so min and max are as one would expect
    extents['xz'] = [numpy.amin(x),numpy.amax(x),numpy.amin(z),numpy.amax(z)]
    extents['yz'] = [numpy.amin(y),numpy.amax(y),numpy.amin(z),numpy.amax(z)]
  else:
    extents['xz'] = [numpy.amin(x),numpy.amax(x),numpy.amax(z),numpy.amin(z)]
    extents['yz'] = [numpy.amin(y),numpy.amax(y),numpy.amax(z),numpy.amin(z)]
    
  
  # set the limits for the colorbars for each field, which are different
  # in some cases between the "colder" experiments (Ocean2 and Ocean4) and
  # "warmer" experiments (all the rest)
  # NOTE: Users should feel free to modify these limits to suit their needs.
  limits = {}
  if(experiment in ['Ocean0', 'Ocean1', 'Ocean3', 'IceOcean1', 'IceOcean2']):
    TLimits = [-2.5,1.1]
    SLimits = [33.6,34.8]
    limits['meltRate'] = [-100.,100.]
    limits['thermalDriving'] = [-1., 1.]
    limits['halineDriving'] = [-10., 10.]
    limits['frictionVelocity'] = [0, 0.02]
    limits['barotropicStreamfunction'] = [-0.3, 0.3]
    limits['overturningStreamfunction'] = [-0.2, 0.2]
  else:
    TLimits = [-2.5,-1.8]
    SLimits = [33.6,34.8]
    limits['meltRate'] = [-5.,5.]
    limits['thermalDriving'] = [-0.2, 0.2]
    limits['halineDriving'] = [-2.0, 2.0]
    limits['frictionVelocity'] = [0, 0.005]
    limits['barotropicStreamfunction'] = [-0.05, 0.05]
    limits['overturningStreamfunction'] = [-0.01,0.01]
  
  vLimits = [-0.5, 0.5]
  limits['uBoundaryLayer'] = vLimits
  limits['vBoundaryLayer'] = vLimits
  limits['bottomTemperature'] = TLimits
  limits['bottomSalinity'] = SLimits
  limits['temperatureXZ'] = TLimits
  limits['salinityXZ'] = SLimits
  limits['temperatureYZ'] = TLimits
  limits['salinityYZ'] = SLimits
  

  
  # light gray for use as an "invalid" background value wherever
  # data has been masked out in the NetCDF file
  backgroundColor = (0.9,0.9,0.9)
  
  for timeIndex in range(len(times)):
    print timeIndex, 'year: ', times[timeIndex]
      
    # make plots for all the fields
    # meltRate needs to be scaled to m/a (as we are used to seeing)
    plot('meltRate', 'melt rate (m/a water equiv.)', cmap, scale=sPerYr)
    plot('thermalDriving', 'thermal driving (C)', cmap)
    plot('halineDriving', 'haline driving (PSU)', cmap)
    plot('frictionVelocity', 'friction velocity (m/s)', cmap)
    plot('bottomTemperature', 'sea-floor temperature (C)', cmap)
    plot('bottomSalinity', 'sea-floor salinity (PSU)', cmap)
    plot('uBoundaryLayer', 'x-velocity (m/s) in the sub-shelf boundary layer', 
         cmap)
    plot('vBoundaryLayer', 'y-velocity (m/s) in the sub-shelf boundary layer', 
         cmap)
    # streamfunctions are scaled to Sv (10^6 m^3/s)
    plot('barotropicStreamfunction', 'barotropic streamfunction (Sv)', cmap, 
         scale=1e-6)
  
    plot('overturningStreamfunction', 'overturning streamfunction (Sv)', cmap, 
         scale=1e-6)
    plot('temperatureXZ', 'temperature (C) sliced through y=40 km', cmap)
    plot('salinityXZ', 'salinity (PSU) sliced through y=40 km', cmap)
    plot('temperatureYZ', 'temperature (C) sliced through x=540 km', cmap)
    plot('salinityYZ', 'salinity (PSU) sliced through x=540 km', cmap)
               
  ncFile.close()

if __name__ == "__main__":
  # we could add some optional command-line argument here but so far none...
  parser = OptionParser()
  options, args = parser.parse_args()

  if(len(args) < 2):
    print "usage: plotMISOMIP1OceanData.py <in_file.nc> <out_dir>"
    exit(1)

  # arguments are the COM file name and the directory for images
  fileName = args[0]
  outFolder = args[1]

  plotOcean(fileName, outFolder)
