#!/usr/bin/env python

from optparse import OptionParser
from netCDF4 import Dataset

'''
This script reads the salt values from the first file and writes them 
into the second file. 
Usage: make_sponge.py <salt_file.nc> <temp_file.nc>
'''

if __name__ == "__main__":
  # we could add some optional command-line argument here but so far none...
  parser = OptionParser()
  options, args = parser.parse_args()

  if(len(args) < 2):
    print "usage: make_sponge.py <salt_file.nc> <temp_file.nc>"
    exit(1)

  fname1 = args[0]
  fname2 = args[1]
  salt = Dataset(fname1).variables['Salt'][:]
  file = Dataset(fname2,'r+')
  file.variables['Salt'][:]=salt[:]
  file.close()
  print('*** Variable Salt was written successfully. \n')
