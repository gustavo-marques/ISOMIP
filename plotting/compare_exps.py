from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np

# Compare time series from different runs

paths = ['/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ocean_only/ISOMIP/OCEAN0/layer/COM/RUN','/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ocean_only/ISOMIP/OCEAN0/layer/TPY/g1.0','/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ocean_only/ISOMIP/OCEAN0/layer/TPY/RUN0']

exps = ['COM','TPY1','TPY2']
files = ['Ocean0_COM_MOM6.nc','Ocean0_TPY1_MOM6.nc','Ocean0_TPY2_MOM6.nc']
diags = ['meanMeltRate', 'totalMeltFlux','totalOceanVolume','meanTemperature','meanSalinity']
colors = ['k-','b-','r-']
units = ['m/s','kg/s','m^3','deg C','PSU']
for d in range(len(diags)):
    plt.figure()
    for e in range(len(exps)):
        data = Dataset(paths[e]+'/'+files[e]).variables[diags[d]][:]         
        time = Dataset(paths[e]+'/'+files[e]).variables['time'][:]/(3600*24*365)
        plt.plot(time,data,colors[e],label=files[e],lw=1.5)

    plt.legend(loc='best', shadow=True)
    plt.xlabel('time (years)')
    plt.ylabel(diags[d] + ' ('+units[d]+')')
    plt.savefig(diags[d]+'.png',bbox_inches='tight')
    plt.close()
