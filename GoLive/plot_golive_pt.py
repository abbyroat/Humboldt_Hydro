from matplotlib import pyplot as plt
import numpy as np
import sys

fname = sys.argv[1]
arr = np.load(fname)


fig, ax = plt.subplots(figsize=(12, 4))

keepThreshold = 0.0 # fraction of cells that need to be good to process this scene
minSpd = 0.8 # minimum speed in m/d for a location to be included

pth='035'
row='003'
sep='016'
pth='*'; row='*'
sep='*'

cnt = 0
for i in range(len(arr)):
     print(i)

     cnt += 1
     v = arr[i,0]
     startDOY = arr[i,1]
     endDOY = arr[i,2]
     midDOY = arr[i,3]
     delt = arr[i,4] 
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





