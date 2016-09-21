#!/usr/bin/env python

# plot all the terms used to compute melting for a user specified point (i,j)
# Gustavo Marques, Sep. 2016

import argparse
from netCDF4 import MFDataset, Dataset
import numpy as np
import matplotlib.pyplot as plt
import warnings
import os

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Plot all the terms used to compute melting for a user specified point (i,j)
      ''',
  epilog='Written by Gustavo Marques, Sep. 2016.')

  parser.add_argument('-j', type=int, default=0,
      help='''j indice (integer). Default is 0.''')

  parser.add_argument('-i', type=int, default=0,
      help='''i indice (integer). Default is 0.''')

  optCmdLineArgs = parser.parse_args()
  createPlot(optCmdLineArgs)


def createPlot(args):
     """
     Generates a plot based on i,j
     """
     i=args.i; j=args.j
     print('Extracting time series for (i,j) = '+str(i)+', '+str(j))
     # read vars
     time = Dataset('prog.nc').variables['time'][:]
     melt = Dataset('prog.nc').variables['melt'][:,j,i]
     ustar = Dataset('prog.nc').variables['ustar_shelf'][:,j,i]
     #u_ml = ncfile.variables['u_ml'][:,j,i]
     #v_ml = ncfile.variables['v_ml'][:,j,i]
     #h_ml = ncfile.variables['h_ML'][:,j,i]
     exch_vel_t = Dataset('prog.nc').variables['exch_vel_t'][:,j,i]
     #exch_vel_s = ncfile.variables['exch_vel_s'][:,j,i]
     #tfreeze = ncfile.variables['tfreeze'][:,j,i]
     #mass_flux = ncfile.variables['mass_flux'][:,j,i]
     #haline_driving = ncfile.variables['haline_driving'][:,j,i]
     thermal_driving = Dataset('prog.nc').variables['thermal_driving'][:,j,i]
     #h = ncfile.variables['h'][:,:,j,i]
     #SST = ncfile.variables['SST'][:,j,i]
     #SSS = ncfile.variables['SSS'][:,j,i]
     #SSH = ncfile.variables['SSH'][:,j,i]

     # some params
     L = 3.34E5
     rho_i = 918.0
     cw = 3974.0
     rho_sw = 1028.0
     rho_fw = 1000.
     tmp = (rho_sw * cw)/L

     n = 4
     f, ax = plt.subplots(n, sharex=True, figsize=(18, 10))
     # estimated melting
     est_melt = tmp * exch_vel_t * thermal_driving * (86400.0*365.0/rho_fw)
     ustar_melt = tmp * exch_vel_t/ustar * ustar.mean() * thermal_driving * (86400.0*365.0/rho_fw)
     temp_melt = tmp * exch_vel_t * thermal_driving.mean() * (86400.0*365.0/rho_fw)
     print 'recovered melt, min/max',est_melt.min(),est_melt.max()

     ax[0].plot(time, melt,'k', label = 'model',lw=1.5)
     ax[0].plot(time[::4], est_melt[::4],'r*', label = 'estimated')
     ax[0].plot(time, ustar_melt,'b', label = 'melt with <ustar>')
     ax[0].plot(time, temp_melt,'g', label = 'melt with <dT>')
     ax[0].legend(loc = 'best', shadow=True, ncol=4)
     ax[0].set_ylabel('melt [m/y]')
     ax[0].set_title('i = '+str(i)+' and j = '+str(j))

     # gamma_t
     gamma = exch_vel_t/ustar
     gamma = np.around(gamma, decimals=2)
     print 'gamma_t min/max', gamma.min(), gamma.max()

     ax[1].plot(time, gamma, lw=1.5)
     ax[1].set_ylabel('gamma_t')

     # ustar
     ax[2].plot(time, ustar, lw=1.5)
     ax[2].plot(time, ustar.mean()*np.ones(len(time)),'k', lw=1.0)
     ax[2].set_ylabel('u* [m/s]')

     # thermal_driving
     ax[3].plot(time, thermal_driving, lw=1.5)
     ax[3].plot(time, thermal_driving.mean()*np.ones(len(time)),'k', lw=1.0)
     ax[3].set_ylabel(r'$\Delta T$ [C]')

     ax[n-1].set_xlabel('Time (days)')
     ax[n-1].set_xlim(time[0], time[-1])
     plt.show()
 
# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()

