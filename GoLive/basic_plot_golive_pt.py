from matplotlib import pyplot as plt
import numpy as np
import sys
import pyproj
import argparse
import netCDF4


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.description = __doc__
parser.add_argument("-f", "--file", dest="pdata",
                    help=".npy file containing processed GOLIVE point data",
                    metavar="FILENAME")
parser.add_argument("-m", "--plotmali", dest="plotmali", action='store_true', help="plot mali data")
args = parser.parse_args()

MALIdir = "/global/project/projectdirs/piscees/MALI_projects/Humboldt_hydrology"


fname = args.pdata
arr = np.load(fname)



fig, ax = plt.subplots(2, 3, figsize=(15, 7))

keepThreshold = 0.0 # fraction of cells that need to be good to process this scene
minSpd = 0.8 # minimum speed in m/d for a location to be included

pth='035'
row='003'
sep='016'
pth='*'; row='*'
sep='*'

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
     if yr>=2014 and yr<=2018:
        sp = int(yr-2014)
        doy = midDOY - np.floor(midDOY)
        #ax.flatten()[sp].plot(doy, v, 'o', color=col)
        ax.flatten()[sp].plot([startDOY-yr, endDOY-yr], [v, v], '-', color=col)
print("count={}".format(cnt))

#f2 = Dataset('Humboldt_2to8km_2014_runoff_RACMO2.3p2.nc','r')
#xCell = f2.variables['xCell'][:]
#yCell = f2.variables['yCell'][:]
#ind = np.argmin(((xCell-xm)**2 + (yCell-ym)**2)**0.5)
#print(ind, xCell[ind], yCell[ind])
#melt = f2.variables['externalWaterInput'][:,:].mean(axis=1)
#print(melt.max())
#plt.plot(np.arange(365), melt/melt.max()*200+100, '-y')

for i in range(5):
   a = ax.flatten()[i]
   a.set_xlabel('Date')
   a.set_ylabel('Velocity (m/yr)')
   a.set_title(f'{int(i+2014)}')
   a.set_xlim((0.15, 0.85))
   a.set_zorder(1)
   a.patch.set_visible(False)

# Parse filename UTM
if '/' in fname:
    fstr = fname.split('/')[-1]
else:
    fstr = fname
x0 = fstr.split('_')[0]
y0 = fstr.split('_')[1].split('.')[0]


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

if args.plotmali:
   crs_in = 'proj=utm +zone=20 +ellps=WGS84'
   crs_out = '+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=315.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84'
   t = pyproj.Transformer.from_crs(crs_in,crs_out)
   xm, ym = t.transform(x0, y0)

   for y in range(2014, 2019):
      f = netCDF4.Dataset(f'{MALIdir}/output_{y}.nc')
      xCell = f.variables['xCell'][:]
      yCell = f.variables['yCell'][:]
      ind = np.argmin( (xCell-xm)**2 + (yCell-ym)**2)

      N = f.variables['effectivePressure'][:, ind]

      days = np.arange(0, 365, 1)/365.0
      sp = y-2014
      ax2 = ax.flatten()[sp].twinx()
      #ax2.plot(days, (N.max()-N)/(N.max()-N).max()*200.0+500.0, '-k')
      ax2.plot(days, N, '-k')
      ax2.set_ylim((0.0, 2.0e6))
      ax2.set_ylabel('N (Pa)')
      ax2.invert_yaxis()

      melt = f.variables['externalWaterInput'][:, ind] /1000.0 * 3600.0 * 24.0
      ax3 = ax.flatten()[sp].twinx()
      ax3.plot(days, melt, '-b', linewidth=0.5)
      ax3.spines["right"].set_position(("axes", 1.2))
      make_patch_spines_invisible(ax3)
      ax3.spines["right"].set_visible(True)
      ax3.tick_params(axis='y', colors='b')
      ax3.set_ylabel('runoff (m/d)')
      ax3.yaxis.label.set_color('b')

   # plot location map
   # uses final output file, seems ok
   a = ax.flatten()[-1]
   H = f.variables['thickness'][0,:]
   melt = (f.variables['externalWaterInput'][:] /1000.0 * 3600.0 * 24.0).sum(axis=0)
   a.scatter(xCell, yCell, 3, (H*melt)>0.0)
   a.plot(xCell[ind], yCell[ind], 'r*')

   a.axis('equal')
   a.set_xlim((-4.85e5, -2.65e5))
   a.set_ylim((-1.191e6, -1.0e6))



fig.suptitle(f'UTM x={x0}, y={y0}')
plt.tight_layout()
plt.show()





