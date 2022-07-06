from matplotlib import pyplot as plt
import datetime
from netCDF4 import Dataset
import glob
import numpy as np
import numpy.ma as ma
import pyproj
import argparse
import csv


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.description = __doc__
parser.add_argument("-d", "--data_dir", dest="data_dir",
                    help="Directory containing GoLive datasets",
                    default='/global/cfs/projectdirs/piscees/MALI_projects/GOLIVE_velocity/',
                    metavar="DIRNAME")
parser.add_argument("-f", "--pointsfile", dest="points_file",
                    help="CSV file containing desired analysis points in lat-lon",
                    default=None,
                    metavar="FILENAME")
parser.add_argument("-p", "--path", dest="pth", help="Landsat path number", default='*')
parser.add_argument("-r", "--row", dest="row", help="Landsat row number", default='*')
parser.add_argument("-s", "--sep", dest="sep", help="Landsat separation number", default='*')

args = parser.parse_args()

data_dir = args.data_dir
pth = str(args.pth)
row = str(args.row)
sep = str(args.sep)
buf = 2
# ITSLIVE in UTM :(
# Below are some reference points
#x0=475e3; y0=8862e3 # fast region in peninsula
#p=486e3; y0=8863e3 # u/s of fastest peninsula region
#xp=492e3; y0=8861e3 # farther u/s of fastest peninsula region
#xp=481e3; y0=8847e3 # u/s of secondary area

# Read in points for processing from CSV file
lats = []
lons = []
with open(args.points_file, 'r', newline='') as points_file:
    reader = csv.reader(points_file, delimiter=',')
    next(reader)  # skip header row of CSV file
    for csv_row in reader:
        lons.append(float(csv_row[0]))
        lats.append(float(csv_row[1]))

keepThreshold = 0.0 # fraction of cells that need to be good to process this scene
minSpd = 0.8 # minimum speed in m/d for a location to be included

# Create object used to transform from lat-lon to UTM
transformer = pyproj.Proj(proj='utm', zone='20', ellps='WGS84', preserve_units=False)

direc='{}p{}_r{}/'.format(data_dir, pth, row)
filelist = glob.glob(direc+'L8_{}_{}_{}*.nc'.format(pth, row, sep))

for lat, lon in zip(lats, lons):
    arr = np.empty((0,5))
    cnt = 0

    # Perform coordinate transformation from lat-lon to UTM
    (x0, y0) = transformer(lon, lat)
    print(f'Processing point {x0}, {y0}')
    for infile in sorted(filelist, reverse=True):
      print('    ', infile)
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
    
         delt = float(f.variables['image_pair_times'].del_t) #diff between start and end day
    
         arr = np.append(arr,[[v,startDOY,endDOY, midDOY, delt],], 0)
    
    print(f"count={cnt}")
    print("data shape =", arr.shape)
    fname = f'{x0:06.0f}_{y0:07.0f}.npy'
    np.save(fname, arr)

