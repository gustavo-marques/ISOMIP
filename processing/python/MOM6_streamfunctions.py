import numpy
from netCDF4 import Dataset
import shutil
import argparse

from computeOSF import computeOSF

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--plot", dest="plot", action='store_true',
                    help="Perform plotting at each time step.")
parser.add_argument("-m", "--mom6_file", dest="mom6_file",
                    help="The ROMS input file", required=True)
parser.add_argument("-i", "--input_file", dest="input_file",
                    help="The ISOMIP+ input file", required=True)
parser.add_argument("-o", "--output_file", dest="output_file",
                    help="The ISOMIP+ output file", required=True)
args = parser.parse_args()

# first, make a copy of the input file in which we'll just overwrite the OSF
shutil.copy(args.input_file, args.output_file)

nxOut = 240
nyOut = 40
nzOut = 144
dzOut = 720./nzOut


nzExtra = 6

# the z location of grid-cell corners (or top-bottom inderfaces) on the output
# grid, with extra points at the top to catch anything above 0
zInterfaceOut = -dzOut*numpy.arange(-nzExtra, nzOut+1)

# the z location of grid-cell centers on the output grid
zOut = 0.5*(zInterfaceOut[0:-1] + zInterfaceOut[1:])

ncIn = Dataset(args.mom6_file, 'r')

variables = ncIn.variables

nTime = len(ncIn.dimensions['time'])
nx = len(ncIn.dimensions['xh'])
ny = len(ncIn.dimensions['yh'])
nz = len(ncIn.dimensions['zl'])

assert(nx == nxOut)
assert(ny == nyOut)

# 1D x and y of rho points
x = 1e3*variables['xh'][:]
y = 1e3*variables['yh'][:]

dx = x[1]-x[0]
dy = y[1]-y[0]

# x_u = 1e3*variables['xq'][0:-1]
# y_v = 1e3*variables['yq'][0:-1]

ncOut = Dataset(args.output_file, 'r+')
outVariables = ncOut.variables

for tIndex in range(nTime):
    print "time index {} of {}".format(tIndex, nTime)
    zInterface = variables['e'][tIndex, :, :, :]
    h = variables['h'][tIndex, :, :, :]

    # the last u point is on the boundary and can (should?) be ignored
    u = variables['u'][tIndex, :, :, 0:-1]

    zInterface_u = 0.5*(zInterface[:, :, 0:-1] + zInterface[:, :, 1:])
    h_u = 0.5*(h[:, :, 0:-1] + h[:, :, 1:])

    uMask = numpy.ones((nz, ny, nx-1), bool)

    osfOut = computeOSF(u, uMask, dy=dy*numpy.ones(u.shape),
                        zInterface=zInterface_u,
                        zInterfaceOut=zInterfaceOut, plot=args.plot,
                        xPlot=x, zPlot=zInterface, tIndex=tIndex)

    # skip the extra points we added at the top, since those aren't part of
    # standard ISOMIP+ output
    outVariables['overturningStreamfunction'][tIndex, :, :] = \
        osfOut[nzExtra:, :]

ncIn.close()
ncOut.close()
