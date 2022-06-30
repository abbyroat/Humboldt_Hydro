from matplotlib import pyplot as plt
import datetime
from netCDF4 import Dataset
import glob
import numpy as np
import numpy.ma as ma


# ITSLIVE in UTM :(
x0=475e3; y0=8862e3 # fast region in peninsula
xp=486e3; y0=8863e3 # u/s of fastest peninsula region
xp=492e3; y0=8861e3 # farther u/s of fastest peninsula region
#xp=481e3; y0=8847e3 # u/s of secondary area


fig, ax = plt.subplots(figsize=(12, 4))

keepThreshold = 0.0 # fraction of cells that need to be good to process this scene
minSpd = 0.8 # minimum speed in m/d for a location to be included

pth='035'
row='003'
sep='016'
pth='*'; row='*'
sep='*'
buf=2

dir='p{}_r{}/'.format(pth, row)
filelist = glob.glob(dir+'L8_{}_{}_{}*.nc'.format(pth, row, sep))

cnt = 0
for infile in sorted(filelist, reverse=True):
  print(infile)
  f = Dataset(infile, 'r')
  x=f.variables['x'][:]
  y=f.variables['y'][:]

  i = np.argmin(np.absolute(x-x0))
  j = np.argmin(np.absolute(y-y0))
  data = f.variables['vv_masked'][j-buf:j+buf+1, i-buf:i+buf+1]
  corr = f.variables['corr'][j-buf:j+buf+1, i-buf:i+buf+1]
  badCorrMask = corr<0.3
  data.mask = np.logical_or(data.mask, badCorrMask)
  #print(goodCorrMask.sum(), buf, (2*buf+1)**2)
  if data.count() > 0 and data.mask.sum() < (2*buf+1)**2 * 0.25:
     cnt += 1
     v = data.mean() * 365.0
     startDOY = float(f.variables['image_pair_times'].start_time_decimal_year)
     endDOY = float(f.variables['image_pair_times'].end_time_decimal_year)
     midDOY = float(f.variables['image_pair_times'].mid_time_decimal_year)

     delt = float(f.variables['image_pair_times'].del_t)
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
     #plt.plot([startDOY, endDOY], [v, v], '-', color=col)
     plt.plot(midDOY, v, 'o', color=col)
print("count={}".format(cnt))

#f2 = Dataset('Humboldt_2to8km_2014_runoff_RACMO2.3p2.nc','r')
#xCell = f2.variables['xCell'][:]
#yCell = f2.variables['yCell'][:]
#ind = np.argmin(((xCell-xm)**2 + (yCell-ym)**2)**0.5)
#print(ind, xCell[ind], yCell[ind])
#melt = f2.variables['externalWaterInput'][:,:].mean(axis=1)
#print(melt.max())
#plt.plot(np.arange(365), melt/melt.max()*200+100, '-y')



plt.xlabel('Date')
plt.ylabel('Velocity (m/yr)')
plt.show()





