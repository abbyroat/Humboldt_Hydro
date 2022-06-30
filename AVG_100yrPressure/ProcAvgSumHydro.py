#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 09:29:27 2022

@author: abbyroat
Extract the averaged summer data from all 100 (99 (06/23/22) because no 2013) years - ploted in ../PlotAvgSumHydro

So everything works now I am duplicating 2012 once I get 2013 that will change.
"""

from __future__ import absolute_import, division, print_function, unicode_literals


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import NearestNDInterpolator 
EP_MEAN_summer= [] #weighted mean effective pressure
EWI_SUM_summer = [] # total ExternalWaterInput for a summer
WP_MEAN_summer= []# mean WaterPressure for a summer
outputfiles = ['/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2000.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2001.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2002.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2003.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2004.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2005.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2006.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2007.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2008.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2009.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2010.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2011.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2012.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2012.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2014.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2015.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2016.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2017.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2018.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2019.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2020.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2021.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2022.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2023.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2024.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2025.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2026.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2027.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2028.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2029.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2030.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2031.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2032.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2033.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2034.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2035.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2036.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2037.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2038.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2039.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2040.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2041.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2042.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2043.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2044.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2045.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2046.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2047.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2048.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2049.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2050.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2051.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2052.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2053.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2054.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2055.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2056.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2057.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2058.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2059.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2060.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2061.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2062.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2063.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2064.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2065.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2066.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2067.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2068.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2069.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2070.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2071.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2072.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2073.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2074.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2075.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2076.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2077.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2078.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2079.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2080.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2081.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2082.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2083.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2084.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2085.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2086.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2087.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2088.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2089.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2090.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2091.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2092.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2093.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2094.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2095.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2096.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2097.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2098.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2099.nc']
#outputfiles = ['/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/noChannels/output_2090.nc']#, '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/noChannels/output_2091.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/noChannels/output_2092.nc']#, '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/noChannels/output_2093.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/noChannels/output_2094.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/noChannels/output_2095.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/noChannels/output_2096.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/noChannels/output_2097.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/noChannels/output_2098.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/noChannels/output_2099.nc']#'/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2090.nc']
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
    
    
    ''' ----------------------------------------------Masks-----------------------------------------------------'''
    
    ''' Mask all grounded ice'''
    li_mask_ValueFloating = 4
    floatMask = (cellMask & li_mask_ValueFloating) // li_mask_ValueFloating
    groundedMask = (1 - floatMask) * (thickness > 0.0)
    groundedMask= groundedMask.astype('float')
    groundedMask[groundedMask==0]= np.nan
    #groundedMaskInterpolator = LinearNDInterpolator(np.vstack((xCell, yCell)).T, groundedMask[0,:])
    #groundedMaskOnEdges = groundedMaskInterpolator(np.vstack((xEdge, yEdge)).T)
    
    '''Mask 10km in front of glacier'''
    elevationMask_10km = thickness *1
    elevationMask_10km[elevationMask_10km < 475]='nan' # the 10kms at the front of the glacier
    elevationMask_10km[elevationMask_10km>=475] = 1 
    
    '''Mask by ablation zone'''
    EWI_ab = externalWaterInput *1
    daily_sum_ablation = np.sum(EWI_ab, axis=1)
    max_ablation = np.argmax(daily_sum_ablation) #the index of the day with the largest external water input
    ablation_mask = EWI_ab[max_ablation,:] #index day of max ablation
    ablation_mask[ablation_mask<0.0001] = np.nan
    ablation_mask[ablation_mask>=0.0001] = 1

    
    '''grounded ice mask'''
    groundedAreaCell = areaCell * groundedMask[0,:]
    groundedAreaCell = groundedAreaCell[~np.isnan(groundedAreaCell)]
    areaCellmean = np.sum(groundedAreaCell)
    
    ablationAreaCell = areaCell * ablation_mask*groundedMask[0,:]
    ablationAreaCell =  ablationAreaCell[~np.isnan( ablationAreaCell)]
    areaCell_ablation =np.sum(ablationAreaCell)
    
    
    '''----------------------------------------MAIN VARIABLES---------------------------------------------------'''
    '''EFFECTIVE PRESSURE - weighted by area of cell - currently calculates MEAN over each summer '''
    grounded_N_tkm = groundedMask * effectivePressure * areaCell
    grounded_N_tkm_noNaN = grounded_N_tkm [:, ~np.isnan(grounded_N_tkm ).any(axis=0)]
    grounded_N_tkm_timeseries = np.sum(grounded_N_tkm_noNaN, axis=1)
    weighted_grounded_N_tkm_timeseries = grounded_N_tkm_timeseries/areaCellmean
    W_N_thkm = np.mean(weighted_grounded_N_tkm_timeseries[151:245]) #mean
    
    '''WATER PRESSURE - MEAN over summer '''
    waterPressure = ((910.0*9.81*thickness)-effectivePressure)
    #grnd_P_weight = elevationMask_10km*waterPressure*areaCell # 10km mask
    grnd_P_weight = groundedMask*waterPressure*areaCell #grnd mask
    grnd_P_weight_noNaN = grnd_P_weight [:, ~np.isnan(grnd_P_weight).any(axis=0)]
    grnd_P_weight_timeseries = np.sum(grnd_P_weight_noNaN, axis=1)
    weighted_grnd_P_weight_timeseries = grnd_P_weight_timeseries/areaCellmean
    weight_P_summer= np.mean(weighted_grnd_P_weight_timeseries[151:245])
    print(weight_P_summer)
    
    ''' EXTERNAL WATER INPUT - calculates the SUM over the summer - NOT MEAN '''
    grounded_EWI_con = groundedMask*((externalWaterInput*areaCell)/1000) # grounded Mask + units to m^3 s^-1
    grounded_EWI_con_noNaN = grounded_EWI_con[:, ~np.isnan(grounded_EWI_con).any(axis=0)]
    gr_convEWI_timeseries = np.sum(grounded_EWI_con_noNaN, axis=1)
    weight_EWI_summer = np.sum(gr_convEWI_timeseries[151:245])
    
    
    '''APPEND ALL IMPORTANT VARIABLES'''
    EP_MEAN_summer.append(W_N_thkm) # mean EffectivePressure for a summer
    WP_MEAN_summer.append(weight_P_summer)# mean WaterPressure for a summer
    EWI_SUM_summer.append(weight_EWI_summer)# total ExternalWaterInput for a summer
    
    
    np.save('EP_MEAN_summer.npy', np.array(EP_MEAN_summer))
    np.save('WP_MEAN_summer.npy', np.array(WP_MEAN_summer))
    np.save('EWI_SUM_summer.npy', np.array(EWI_SUM_summer))
