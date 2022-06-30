#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 20:41:02 2021

@author: abbyroat
"""


from __future__ import absolute_import, division, print_function, unicode_literals


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import NearestNDInterpolator 
outputfiles = ['/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2088.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/noChannels/output_2088_noCH.nc']#, '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2091.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2092.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2093.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2094.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2095.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2096.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2097.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2098.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2099.nc']
for file in outputfiles:
    f=Dataset(file,'r') # open a NetCDF file in 'read' mode
    f.set_auto_mask(False)
    
    xCell=f.variables['xCell'][:]  # x coordinate of cell position for all cells
    yCell=f.variables['yCell'][:]  # y coordinate of cell position for all cells
    xEdge = f.variables['xEdge'][:]
    yEdge = f.variables ['yEdge'][:]
    areaCell = f.variables['areaCell'][:]
    dvEdge = f.variables['dvEdge'][:] #nedges
    
    externalWaterInput = f.variables['externalWaterInput'][:]
    daysSinceStart = f.variables['daysSinceStart'][:]
    day = daysSinceStart
    basalMeltInput = f.variables['basalMeltInput'][:]
    channelMelt = f.variables['channelMelt'][:]
    dcEdge = f.variables['dcEdge'][:]
    waterThickness = f.variables['waterThickness'][:]
    effectivePressure = f.variables['effectivePressure'][:]
    thickness = f.variables['thickness'][:]
    waterFlux = f.variables['waterFlux'][:] #m^2/s
    cellMask = f.variables['cellMask'][:]
    
    
    '''
    Masks
    '''
    li_mask_ValueFloating = 4
    floatMask = (cellMask & li_mask_ValueFloating) // li_mask_ValueFloating
    groundedMask = (1 - floatMask) * (thickness > 0.0)
    groundedMask= groundedMask.astype('float')
    groundedMask[groundedMask==0]= np.nan
    #groundedMaskInterpolator = LinearNDInterpolator(np.vstack((xCell, yCell)).T, groundedMask[0,:])
    #groundedMaskOnEdges = groundedMaskInterpolator(np.vstack((xEdge, yEdge)).T)
    
    elevationMask_10km = thickness *1
    elevationMask_10km[elevationMask_10km < 475]='nan' # the 10kms at the front of the glacier
    elevationMask_10km[elevationMask_10km>=475] = 1 
    
    EWI_ab = externalWaterInput *1
    daily_sum_ablation = np.sum(EWI_ab, axis=1)
    max_ablation = np.argmax(daily_sum_ablation) #the index of the day with the largest external water input
    ablation_mask = EWI_ab[max_ablation,:] #index day of max ablation
    ablation_mask[ablation_mask<0.0001] = np.nan
    ablation_mask[ablation_mask>=0.0001] = 1
    
    convert_EWI1 = externalWaterInput*areaCell
    convert_EWI2 = convert_EWI1/1000
    #grounded_EWI_con = convert_EWI2* elevationMask_10km #10km mask
    
    grounded_EWI_con = groundedMask*waterThickness#*convert_EWI2 # grounded Mask
    grounded_EWI_con_noNaN = grounded_EWI_con[:, ~np.isnan(grounded_EWI_con).any(axis=0)]
    gr_convEWI_timeseries = np.sum(grounded_EWI_con_noNaN, axis=1)
    weight_EWI_summer = np.sum(gr_convEWI_timeseries[151:245])
    print(weight_EWI_summer)
    
    
    groundedAreaCell = areaCell * groundedMask[0,:]
    groundedAreaCell = groundedAreaCell[~np.isnan(groundedAreaCell)]
    areaCellmean = np.sum(groundedAreaCell)
    

    waterPressure = ((910.0*9.81*thickness)-effectivePressure)
    #grnd_P_weight = elevationMask_10km*waterPressure*areaCell # 10km mask
    grnd_P_weight = groundedMask*waterPressure*areaCell #grnd mask
    grnd_P_weight_noNaN = grnd_P_weight [:, ~np.isnan(grnd_P_weight).any(axis=0)]
    grnd_P_weight_timeseries = np.sum(grnd_P_weight_noNaN, axis=1)
    weighted_grnd_P_weight_timeseries = grnd_P_weight_timeseries/areaCellmean
    weight_P_summer= np.mean(weighted_grnd_P_weight_timeseries[151:245])
    #print(weight_P_summer)
    
    #water pressure with and without channels 1 year timeseries
    year = np.arange(0, len(gr_convEWI_timeseries))
    plt.plot(year, weighted_grnd_P_weight_timeseries)
    plt.xlabel('day')
    plt.ylabel('water pressure')
    plt.title('water pressure over time chnls vs no chnls (2088)')
    plt.legend(('channels', 'no channels'))
    