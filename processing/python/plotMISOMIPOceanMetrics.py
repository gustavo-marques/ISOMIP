#!/Users/gmarques/anaconda/bin/python

import numpy
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os.path

parser = OptionParser()
# could add optional arguments here, but none so far
options, args = parser.parse_args()

if(len(args) < 1):
  print "Usage: python plotMISOMIPOceanMetrics.py <Ocean1_result.nc> <Ocean2_results.nc> ..."
  exit(1)

files = args

sPerYr = 365.*24.*60.*60.
GTPerKg = 1e-12

fieldNames = ['meanMeltRate', 'totalMeltFlux', 'totalOceanVolume',
              'meanTemperature', 'meanSalinity']
titles = ['mean melt rate (m/a)','total melt flux (GT/a)',
          'total ocean volume (1000 km^3)',
          'mean temperature (C)', 'mean salinity (PSU)']
semiLog = [True, True, False, False, False]
scales = [sPerYr, sPerYr*GTPerKg, 1e-12, 1., 1.]

for fileName in files:
  print fileName
  ncFile = Dataset(fileName,'r')
  baseName = os.path.basename(fileName)
  label = baseName.split('_')[0]
  
  try:
    times = numpy.array(ncFile.variables['time'][:], dtype=float)/sPerYr
  except KeyError:
    # fill in missing time array with best guess
    nTime = len(ncFile.dimensions['nTime'])
    times = numpy.arange(nTime)/12.0

  for fIndex in range(len(fieldNames)):
    fieldName = fieldNames[fIndex]
    field = scales[fIndex]*numpy.array(ncFile.variables[fieldName][:], dtype=float)
    plt.figure(fIndex+1)
    if(semiLog[fIndex]):
      plt.semilogy(times, field, label=label)
    else:
      plt.plot(times, field, label=label)
    
for fIndex in range(len(fieldNames)):
  plt.figure(fIndex+1)
  plt.xlabel('time (a)')
  plt.ylabel(titles[fIndex])
  plt.legend()
  plt.draw()
  plt.savefig('%s.png'%fieldNames[fIndex])
  plt.savefig('%s.pdf'%fieldNames[fIndex])

