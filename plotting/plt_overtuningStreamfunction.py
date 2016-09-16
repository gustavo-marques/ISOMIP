#!/usr/bin/env python

from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
import m6toolbox
from optparse import OptionParser

'''
Plot the overturning streamfunction in z space.
'''

if __name__ == "__main__":
  # we could add some optional command-line argument here but so far none...
  parser = OptionParser()
  options, args = parser.parse_args()

  if(len(args) < 2):
    print "usage: plt_overtuningStreamfunction.py <ncfile.nc> 10"
    exit(1)

  # file name
  fname = args[0]
  # time indice 
  t = int(args[1])

  x=Dataset(fname).variables['nx'][:]
  nx=np.zeros(len(x)+1); nx[0:-1]=x; nx[-1]=x[-1]
  e=Dataset(fname).variables['e'][t,:]
  var=Dataset(fname).variables['overturningStreamfunction'][t,:]

  xCoord, yCoord, zData = m6toolbox.section2quadmesh(nx, e, var, representation='linear')
  plt.figure()
  plt.pcolormesh(xCoord, yCoord, zData)
  plt.colorbar()
  plt.xlabel('x (km)')
  plt.ylabel('Elevation (m)')
  plt.ylim((-720,0))
  #plt.savefig(diags[d]+'.png',bbox_inches='tight')
  plt.show()
