#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 16:55:35 2020

@author: olga
"""


import netCDF4 as nc
import numpy as np

rho_i = 918.0              # density of ice (kg/m^3)
rho_sw = 1028.0            # density of sea water (kg/m^3)
g = 9.8                    # the gravity acceleration (m/s^2)
Bbar = 1.68e8              # ice stiffness parameter ( Pa s^(1/3))
secinyear = 86400.*365.25  # sec per year
n = 3                      # ice rheology exponent
Hgl = 720.                 # maximum ice thickness at the grounding line (m)
Ugl = 300.0/secinyear      # ice velocity at the grounding line (m/s)
Qgl = Ugl*Hgl              # ice flux at the grounding line (m^2/s)
adot=-0.3/secinyear        # accumulation/ablation rate of the ice shelf (negagtive for meelting, m/s)

ny, nx = 4, 240            # 2D ISOMIp geometry
ds = 4000.                 # x spasing (m)
nshelf = 140               # number of ice shelf points
ng = 40                    # number of grounded ice points



delta = 1 - rho_i/rho_sw
Bn=(rho_i*g*delta/8/Bbar)**n
xsh = ds*np.linspace(0,99,num=nshelf)
denom = Qgl**(n+1)*(adot - Bn*Hgl**(n+1)) + Bn*Hgl**(n+1)*(Qgl + adot*xsh)**(n+1) 
Hsh=Hgl*(adot*(Qgl + adot*xsh)**(n+1)/denom )**(1/(n+1))
thick=np.concatenate((Hgl*np.ones(ng),Hsh,np.zeros(nx-nshelf-ng)),axis=None)
h_shelf=np.tile(thick,(ny,1))

area_shelf=np.zeros((ny,nx))
area_shelf[:,:ng+nshelf]=4e6

f_sg = nc.Dataset('Ocean2D_1Dshelf.nc', 'w', format='NETCDF3_CLASSIC')
f_sg.createDimension('ny', ny)
f_sg.createDimension('nx', nx)
area = f_sg.createVariable('area', 'f8', ('ny', 'nx'))
thick = f_sg.createVariable('thick', 'f8', ('ny', 'nx'))
area.units = 'm2'
thick.units = 'm'
f_sg.variables['area'][:] = area_shelf
f_sg.variables['thick'][:] = h_shelf
f_sg.sync()
f_sg.close()