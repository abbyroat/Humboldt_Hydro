#!/usr/bin/env python
'''
Script to plot GOLIVE velocity data that was processed by the proc_golive_pt.py
script.  Invoke with -f specifying a .npy file created by the processing script.
Optionally include the -m argument to include MALI subglacial hydrology output.
The MALI output is hardcoded to a path on Cori, but could be adjusted manually
if needed.  In general, it is desirable to include the MALI output data, but the
option to exclude it is provided if those data are not available.
'''

from matplotlib import pyplot as plt
import numpy as np
import sys
import pyproj
import argparse
import netCDF4

# parse command line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.description = __doc__
parser.add_argument("-f", "--file", dest="pdata",
                    help=".npy file containing processed GOLIVE point data",
                    metavar="FILENAME")
parser.add_argument("-m", "--plotmali", dest="plotmali", action='store_true', help="plot mali data")
args = parser.parse_args()

# path to MALI output. Adjust if needed.
MALIdir = "/global/project/projectdirs/piscees/MALI_projects/Humboldt_hydrology"

# Get processed GOLIVE data from .npy file
fname = args.pdata
arr = np.load(fname)

# create a 2x3 set of subplots to plot into
fig, ax = plt.subplots(2, 3, figsize=(15, 7))

# Begin reading and plotting golive point data from npy file
cnt = 0
for i in range(len(arr)):
     #print(i)

     cnt += 1
     v = arr[i,0]
     startDOY = arr[i,1]
     endDOY = arr[i,2]
     midDOY = arr[i,3]
     delt = arr[i,4]
     yr = np.floor(startDOY)

     if delt == 16.0:
         col = 'r'
     elif delt == 32.0:
         col = 'g'
     elif delt == 48.0:
         col = 'c'
     elif delt == 64.0:
         col = '0.6'
     else:
         col = '0.8'

     # only actually plot the years where we have historical RACMO forcing
     if yr>=2014 and yr<=2019:
        sp = int(yr-2014)
        doy = midDOY - np.floor(midDOY)
        a = ax.flatten()[sp] # get subplot for this year
        #a.plot(doy, v, 'o', color=col) # simple point version
        a.plot([startDOY-yr, endDOY-yr], [v, v], '-', color=col) # bar version
print("count={}".format(cnt))

# Finalize plot axes after all velocity plotting is complete
for i in range(6):
   a = ax.flatten()[i]
   a.set_xlabel('Date')
   a.set_ylabel('Velocity (m/yr)')
   a.set_title(f'{int(i+2014)}')
   a.set_xlim((0.15, 0.85))
   a.set_zorder(1); a.patch.set_visible(False) # move velo data to top (easier to see)

# Parse filename UTM and get UTM coords as x0,y0
if '/' in fname:
    fstr = fname.split('/')[-1]
else:
    fstr = fname
x0 = fstr.split('_')[0]
y0 = fstr.split('_')[1].split('.')[0]

# needed helper function for offsetting second righthand axis
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

# Plot MALI output if requested
if args.plotmali:
   # First convert UTM coords to the coord system MALI is using (xm,ym)
   crs_in = 'proj=utm +zone=20 +ellps=WGS84'
   crs_out = '+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=315.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84'
   t = pyproj.Transformer.from_crs(crs_in,crs_out)
   xm, ym = t.transform(x0, y0)

   # loop over the years we have chosen to plot above
   for y in range(2014, 2019+1):
      # Open needed MALI output file and get needed vars
      f = netCDF4.Dataset(f'{MALIdir}/output_{y}.nc')
      xCell = f.variables['xCell'][:]
      yCell = f.variables['yCell'][:]
      ind = np.argmin( (xCell-xm)**2 + (yCell-ym)**2)

      N = f.variables['effectivePressure'][:, ind]

      days = np.arange(0, 365, 1)/365.0
      sp = y-2014
      ax2 = ax.flatten()[sp].twinx() # replicate the velocity axis for this year
      #ax2.plot(days, (N.max()-N)/(N.max()-N).max()*200.0+500.0, '-k') # old version before I added a 2nd axis
      ax2.plot(days, N, '-k')
      ax2.set_ylim((0.0, 2.0e6)) # hard code ylim so we can more easiliy compare plots
      ax2.set_ylabel('N (Pa)')
      ax2.invert_yaxis() # invert y axis so N changes are in same direction as velo changes

      # now plot RACMO runoff on a second right y-axis
      melt = f.variables['externalWaterInput'][:, ind] /1000.0 * 3600.0 * 24.0 # convert to m/d
      ax3 = ax.flatten()[sp].twinx()
      ax3.plot(days, melt, '-b', linewidth=0.5)
      # next lines offset the yaxis from the N one
      ax3.spines["right"].set_position(("axes", 1.2))
      make_patch_spines_invisible(ax3)
      ax3.spines["right"].set_visible(True)
      ax3.tick_params(axis='y', colors='b') # make tick labels same color as plot line
      ax3.set_ylabel('runoff (m/d)')
      ax3.yaxis.label.set_color('b') # make label same color

   # plot location map
   figMap, axMap = plt.subplots(1, 1, figsize=(5, 5))
   # uses final output file, seems ok
   H = f.variables['thickness'][0,:]
   melt = (f.variables['externalWaterInput'][:] /1000.0 * 3600.0 * 24.0).sum(axis=0)
   # plot approx. of ablation zone - places where ice thickness is nonzero and melt is nonzero
   # note we are arbitrarily using the runoff from the final year plotted above - ok for an approximation
   axMap.scatter(xCell, yCell, 3, (H*melt)>0.0)
   axMap.plot(xCell[ind], yCell[ind], 'r*') # plot the point we analyzed with a red star

   axMap.axis('equal')
   # zoom in a bit to where the action is
   axMap.set_xlim((-4.85e5, -2.65e5))
   axMap.set_ylim((-1.191e6, -1.0e6))



fig.suptitle(f'UTM x={x0}, y={y0}') # state the point we used for provenance
fig.tight_layout() # clean up subplots
plt.show()





