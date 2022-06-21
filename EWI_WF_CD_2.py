#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:57:36 2022

@author: abbyroat
"""


from __future__ import absolute_import, division, print_function, unicode_literals


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
#from scipy.interpolate import NearestNDInterpolator 
outputfiles = ['/Users/abbyroat/Desktop/Badger_outputs_2021/3km_Darcy/Varying_withRunoff/output_2075.nc']
for file in outputfiles:
    f=Dataset(file,'r') # open a NetCDF file in 'read' mode
    f.set_auto_mask(False)
    
    xCell=f.variables['xCell'][:]  # x coordinate of cell position for all cells
    yCell=f.variables['yCell'][:]  # y coordinate of cell position for all cells
    xEdge = f.variables['xEdge'][:]
    yEdge = f.variables ['yEdge'][:]
    areaCell = f.variables['areaCell'][:]
       
    waterFlux = f.variables['waterFlux'][:] #nEdges and time (m^2/s)
    bedTopography = f.variables['bedTopography'][0,:] #nCells 
    #floatingIce = f.variables['cellMask_floating'][0,:] #ncells
    channelDischarge = f.variables['channelDischarge'][:] # edges and time (m^3/s)
    dvEdge = f.variables['dvEdge'][:] #nedges
    dcEdge = f.variables['dcEdge'][:]
    externalWaterInput = f.variables['externalWaterInput'][:] #kg m^{-2} s^{-1}
    daysSinceStart = f.variables['daysSinceStart'][:]
    day = daysSinceStart
    basalMeltInput = f.variables['basalMeltInput'][:]
    channelMelt = f.variables['channelMelt'][:]
    cEdge = f.variables['dcEdge'][:]
    waterThickness = f.variables['waterThickness'][:]
    effectivePressure = f.variables['effectivePressure'][:]
    thickness = f.variables['thickness'][:]
    waterFlux = f.variables['waterFlux'][:] #m^2/s
    cellMask = f.variables['cellMask'][:]
    edgeMask = f.variables['edgeMask'][:]
    nt= len(f.dimensions['Time']) #number of time levels
    nEdges= len(f.dimensions['nEdges']) #number of edges
    nCells= len(f.dimensions['nCells']) #number of edges
    
    areaEdge = dcEdge*dvEdge
    sumAE=  np.sum(areaEdge)
    
    
    '''
    correct units
    '''
    waterFlux_m3d = np.zeros((nt, nEdges))
    for t in range(nt):
        waterFlux_m3d[t,:]= np.abs(waterFlux[t,:]) *60*60*24* dvEdge #m^2/s to m^3/day
    channelDischarge_m3d = np.abs(channelDischarge)*60*60*24 #* areaEdge# m^3/s to m^3/day
    
    '''
    EdgeMask
    '''
    li_mask_ValueMargin = 8 # This is the last cell with ice. 
    li_mask_ValueDynamicMargin = 16 # This is the last dynamically active cell with ice
    li_mask_ValueGroundingLine= 256
    #marginMask = ((edgeMask & li_mask_ValueGroundingLine) // li_mask_ValueGroundingLine)
    marginMask = np.logical_or( (edgeMask & li_mask_ValueGroundingLine) // li_mask_ValueGroundingLine,  (edgeMask & li_mask_ValueDynamicMargin) // li_mask_ValueDynamicMargin ) # MARGIN MASK W/ DYNAMIC MARGIN
    #marginMask = np.logical_or( (edgeMask & li_mask_ValueGroundingLine) // li_mask_ValueGroundingLine,  (edgeMask & li_mask_ValueMargin) // li_mask_ValueMargin ) # MARGIN MASK W/ REGULAR MARGIN (AKA NOT DYNAMIC)

    '''
    Necessary variables on edges
    '''
    #varying = externalWaterInput
    #externalwaterInput on edges
    '''
    time_EWI = 0
    time_EWI_final = len(externalWaterInput)
    EWIOnEdges = np.zeros(np.shape(waterFlux))
    while time_EWI < time_EWI_final:
        EWIInterpolator = LinearNDInterpolator(np.vstack((xCell, yCell)).T, externalWaterInput[time_EWI,:])
        EWIEdges = EWIInterpolator(np.vstack((xEdge, yEdge)).T)
        EWIOnEdges[time_EWI,:]= EWIEdges
        time_EWI += 1
    '''
        
    
    
    
    #not varying = bedTopo, thickness
    bedTopoInterpolator = LinearNDInterpolator(np.vstack((xCell, yCell)).T, bedTopography)
    bedTopoEdges = bedTopoInterpolator(np.vstack((xEdge, yEdge)).T)
    
    
    thicknessInterpolator = LinearNDInterpolator(np.vstack((xCell, yCell)).T, thickness[0,:])
    thicknessOnEdges = thicknessInterpolator(np.vstack((xEdge, yEdge)).T)
    
    
    '''
    Create Masks
    '''
    #grounded ice mask - mask out all area with no/floating ice
    li_mask_ValueFloating = 4
    floatMask = (cellMask & li_mask_ValueFloating) // li_mask_ValueFloating
    groundedMask = (1 - floatMask) * (thickness > 0.0)
    groundedMask= groundedMask.astype('float')
    groundedMask[groundedMask==0]= 'nan'
    groundedMaskInterpolator = LinearNDInterpolator(np.vstack((xCell, yCell)).T, groundedMask[0,:])
    groundedMaskOnEdges = groundedMaskInterpolator(np.vstack((xEdge, yEdge)).T)
    
    groundedMask_2 = np.logical_not(floatMask) * (thickness > 0.0) # ~ reverses the mask - this is the same as groundedMask
    
    groundedAreaEdge = areaEdge * groundedMaskOnEdges
    groundedAreaEdge = groundedAreaEdge[~np.isnan(groundedAreaEdge)]
    areaEdgeSum = np.sum(groundedAreaEdge)
    
    
    #mask by elevation - make masks for a few elevation bands - likely combined with the grounded ice mask
    
    elevationMask_10km = thicknessOnEdges*1
    elevationMask_10km[elevationMask_10km < 475]='nan' # the 10kms at the front of the glacier
    elevationMask_10km[elevationMask_10km>=475] = 1 
    
    elevationMask_1000m = thicknessOnEdges*1 #this
    elevationMask_1000m[elevationMask_1000m < 1000] = 'nan'
    elevationMask_1000m[elevationMask_1000m>1000]=1
    
    
   
    
    #mask by ablation zone some value of ablation above which is included in the mask - ideas: 1 mask ablation zone for single area on the day of greatest ablation for a single year. 2. mask by area of greatest ablation for each individual year (ie the area masked changes every year 3. Mask by some average of the ablation zone over the entire period - this seems like more trouble than it's worth.
    
    '''
    basic time series and masking variables
    '''
    channelDischarge_m3d_timeseries = np.sum(channelDischarge_m3d, axis=1)
    waterFlux_m3d_timeseries = np.sum(waterFlux_m3d, axis=1)
    cd_time = channelDischarge_m3d_timeseries/sumAE
    
    #grounded mask -THIS IS THE IMPORTANT ONE
    channelDischarge_m3d_grdmsk = channelDischarge_m3d * marginMask
    channelDischarge_m3d_grdmsk_noNaN = channelDischarge_m3d_grdmsk[:, ~np.isnan(channelDischarge_m3d_grdmsk).any(axis=0)]
    channelDischarge_m3d_grdmsk_timeseries = np.sum(channelDischarge_m3d_grdmsk_noNaN, axis=1)
    cummulative_chDis = np.cumsum(channelDischarge_m3d_grdmsk_timeseries)
    
    waterFlux_m3d_grdmsk = waterFlux_m3d * marginMask
    waterFlux_m3d_grdmsk_noNaN = waterFlux_m3d_grdmsk[:, ~np.isnan(waterFlux_m3d_grdmsk).any(axis=0)]
    waterFlux_m3d_grdmsk_timeseries = np.sum(waterFlux_m3d_grdmsk_noNaN, axis=1)
    cummulative_waterFlux = np.cumsum(waterFlux_m3d_grdmsk_timeseries) 
    cummulative_disch = cummulative_waterFlux+cummulative_chDis #cumm disch
    
    #groundedMask weighted
    weightedchannelDischarge_m3d_grdmsk = channelDischarge_m3d * groundedMaskOnEdges * areaEdge
    weightedchannelDischarge_m3d_grdmsk_noNaN = weightedchannelDischarge_m3d_grdmsk[:, ~np.isnan(weightedchannelDischarge_m3d_grdmsk).any(axis=0)]
    weightedchannelDischarge_m3d_grdmsk_timeseries = np.sum(weightedchannelDischarge_m3d_grdmsk_noNaN, axis=1)
    weighted_ChDis = weightedchannelDischarge_m3d_grdmsk_timeseries/areaEdgeSum
    #cummulative_chDis = np.cumsum(weighted_ChDis)
    
    weightedwaterFlux_m3d_grdmsk = waterFlux_m3d * groundedMaskOnEdges * areaEdge
    weightedwaterFlux_m3d_grdmsk_noNaN = weightedwaterFlux_m3d_grdmsk[:, ~np.isnan(weightedwaterFlux_m3d_grdmsk).any(axis=0)]
    weightedwaterFlux_m3d_grdmsk_timeseries = np.sum(weightedwaterFlux_m3d_grdmsk_noNaN, axis=1)
    weighted_watFlux = weightedwaterFlux_m3d_grdmsk_timeseries/areaEdgeSum
    #cummulative_waterFlux = np.cumsum(weighted_watFlux ) 
    
    '''
    Weighted surface runoff- but not weighted lol
    '''
    #this produces the same plot as the code below
    
    weight_EWI = np.zeros((nt, nCells))
    for t in range(nt):
        weight_EWI[t,:]= externalWaterInput[t,:]*areaCell #kg/s
    grounded_EWI = ((weight_EWI/1000)* (60*60*24)* groundedMask_2).sum(axis=1) 
    cum_ewi = np.cumsum(grounded_EWI)
    
    '''
    weight_EWI = externalWaterInput*areaCell
    
    convert_EWI2 = weight_EWI/1000
    grounded_EWI_con = groundedMask*convert_EWI2*60*60*24 #m^3/day groundedMask
    grounded_EWI_con_noNaN = grounded_EWI_con[:, ~np.isnan(grounded_EWI_con).any(axis=0)]
    gr_convEWI_timeseries = np.sum(grounded_EWI_con_noNaN, axis=1)
    cum_ewi = np.cumsum(gr_convEWI_timeseries)
    totalEWI = np.sum(externalWaterInput)
    '''
    '''
    Unweighted surface runoff
    
    weight_EWI = externalWaterInput #kg/ m^2 s^1
    
    convert_EWI2 = weight_EWI/1000 #(divide by density of water  kg/(m^2s)*m)
    grounded_EWI_con = groundedMask*convert_EWI2*60*60*24 #m^3/day groundedMask
    grounded_EWI_con_noNaN = grounded_EWI_con[:, ~np.isnan(grounded_EWI_con).any(axis=0)]
    gr_convEWI_timeseries = np.sum(grounded_EWI_con_noNaN, axis=1)
    cum_ewi = np.cumsum(gr_convEWI_timeseries)
    totalEWI = np.sum(externalWaterInput)
    '''
    '''
    Characterixe extent of channelization
    '''
    
    #Possiblities:
    #Qc/(Qc+Qd)
    chDis_frac_totalDis = channelDischarge_m3d/(channelDischarge_m3d + waterFlux_m3d) #basic equation (over entire year)
    chDis_frac_totalDis_noNaN = chDis_frac_totalDis[:, ~np.isnan(chDis_frac_totalDis).any(axis=0)]
    chDis_frac_totalDis_timeseries = np.mean(chDis_frac_totalDis_noNaN, axis=1)
    
    #with grounded Mask
    chDis_frac_totalDis_grdmsk = channelDischarge_m3d_grdmsk/(channelDischarge_m3d_grdmsk+waterFlux_m3d_grdmsk)
    chDis_frac_totalDis_grdmsk_noNaN = chDis_frac_totalDis_grdmsk[:, ~np.isnan(chDis_frac_totalDis_grdmsk).any(axis=0)]
    chDis_frac_totalDis_grdmsk_timeseries = np.mean(chDis_frac_totalDis_grdmsk_noNaN, axis=1)
    #print(np.mean(chDis_frac_totalDis_grdmsk_timeseries[151:245]))
    
    chDis_frac_frtimeseries = channelDischarge_m3d_grdmsk_timeseries/(channelDischarge_m3d_grdmsk_timeseries+waterFlux_m3d_grdmsk_timeseries)
    
    ChD_f_timeseries_test = channelDischarge_m3d_grdmsk_timeseries/(channelDischarge_m3d_grdmsk_timeseries+waterFlux_m3d_grdmsk_timeseries)
    #
    
    # with 10km mask
    '''
    chDis_frac_totalDis_10kmmsk = channelDischarge_m3d_grdmsk/(channelDischarge_m3d_grdmsk+waterFlux_m3d_grdmsk)
    chDis_frac_totalDis_grdmsk_noNaN = chDis_frac_totalDis_grdmsk[:, ~np.isnan(chDis_frac_totalDis_grdmsk).any(axis=0)]
    chDis_frac_totalDis_grdmsk_timeseries = np.mean(chDis_frac_totalDis_grdmsk_noNaN, axis=1)
    '''
    #Qc/meltFlux For this one I'd have to figure out what melt flux actually is (maybe melt flux is just all the discharge) - that was dumb (I should ask Matt) - this is the external water input  - sum of surface melt over space and time over an entire summer season
    #in the past I have defined the summer season thus: 151:245 -  I do both over the summer like I did for effective pressure
    '''
    #first just with the grounded mask
    #I need to come back to this EWI on edges thing because it's being weird right now and I'm not sure wht
    grounded_EWI = EWIOnEdges * groundedMaskOnEdges
    grounded_EWI_noNaN = grounded_EWI[:, ~np.isnan(grounded_EWI).any(axis=0)]
    EWI_total_timeseries = np.sum(grounded_EWI_noNaN, axis=1)
    meltFlux = np.sum(EWI_total_timeseries[151:245])
    '''
    meltFlux_2010s = [34.7439485098203, 44.16771303951361, 35.716244069454675, 39.51414650843285, 35.51735694526069, 19.721322288221966, 14.486278299406758]
    meltFlux_2070s = [89.10129455262734, 51.9765243412584, 95.71400076303998, 88.95420668023445, 56.563633477100346]#71-75 summer
    mean_ch_o_tot_dis = [0.05114379777909823, 0.03295487412210305, 0.046423899273045396, 0.04605943593625629, 0.03248687564609797] #qc/(qd+qc)
    summer_ChDis_sum = np.sum(channelDischarge_m3d_grdmsk_noNaN[151:245,:])
    #print(summer_ChDis_sum)
    sumChannelDis2070s = [154006707084.71317, 55126111043.852646, 184210820359.5477, 160280252153.45956]#, 67062571750.41909] #71-75
    sumChannelDis2000s = [7338201043.917766, 13949959687.478983, 11835865237.795628, 16598469365.82673, 8458144270.530559, 21277337821.22555, 11823567171.481339, 19963249651.515224, 38275446649.10882, 20719697688.403423] # 2000 - 2009
    sumChannelDis2010s = [22434175093.650608, 27918729133.367847, 18102503331.507816, 25231281305.945732, 21755063819.715992, 9936050297.282625, 7638823113.038562] #2010-2018, no 12/13
    
    
    #ExternalWaterInput - on cells
    summer_N_num_2070s = [952433.1130654175, 1023240.635874282, 931910.1900052043, 964907.3625573837] #71-74
    
    summer_N_num_2010s = [1336002.2551886402, 1265313.8774235856, 1302415.0783568895, 1314226.9087370979, 1329556.9745543324, 1510044.4222078896, 1559800.7537990825] #2010-2018 removed 2012 1178650.843603887,
    
    summer_N_num_2000s = [1448801.0106701024, 1368734.5745848196, 1394439.7121913084, 1343034.1246399514, 1468641.357072073, 1348832.6150821776, 1424013.1857352082, 1352888.273385269, 1239666.7054327382, 1383351.2459341804]# 2000-2009 
    
    

    
    
    #N/meltFlux (where N is ice effective pressure) - which is just effective pressure
    #N/Ndist -this would be comparing a run with channels on to one with channels off (therefore just distributed water pressure)
    '''
    Plots
    '''
    
    oneYear = np.arange(1, len(day)+1)
    '''
    #plt.plot(oneYear, channelDischarge_m3d_grdmsk_timeseries)
    #plt.plot(oneYear, (waterFlux_m3d_grdmsk_timeseries+channelDischarge_m3d_grdmsk_timeseries))
    plt.plot(oneYear, chDis_frac_frtimeseries)
    plt.title('high vs low year channel fraction')
    plt.xlabel('day')
    plt.ylabel('channel discharge fraction')
    plt.legend(('high to avg', 'low to avg'))
    '''
    # IMPORTANT STUFF
    '''
    #fig1, ax1 = plt.subplots()
    oneYear = np.arange(1, len(day)+1)
    
    fig, ax_left=plt.subplots()
    ax_right=ax_left.twinx()
    ax_left.set_xlabel('day')
    ax_left.set_ylabel('discharge (m^3 day^-1)')
    ax_left.set_title(' Channelized and distributed discharge (2000)')
    ax_right.set_ylabel('surface runoff(m^3 day^-1)')
    #plt.plo
    ax_left.plot(oneYear, waterFlux_m3d_grdmsk_timeseries   , label = 'distributed')
    ax_left.plot(oneYear, channelDischarge_m3d_grdmsk_timeseries  , label = 'channelized')
    ax_right.plot(oneYear, gr_convEWI_timeseries, linestyle='dashed', color='green', linewidth='0.8', label = 'surface melt')
    ax_left.legend(loc='upper left')
    ax_right.legend(loc='upper right')
    
    
    fig1, ax1 = plt.subplots()
    ax1.set_xlabel('day')
    ax1.set_ylabel('discharge (m^3 day^-1)')
    ax1.set_title('2075')
    
    ax1.plot(oneYear, cummulative_waterFlux, label='distributed')
    ax1.plot(oneYear, cummulative_chDis, label='channelized')
    ax1.plot(oneYear, cummulative_disch, label='total discharge')
    ax1.plot(oneYear, cum_ewi, linestyle='dashed', label='surfaceMelt')
    
    ax1.legend()
    
    fig.savefig('ax1_figure.png', dpi=300)
    '''
    oneYear = np.arange(1, len(day)+1)
    
    fig, ax_left=plt.subplots()
    ax_right=ax_left.twinx()
    ax_left.set_xlabel('day')
    ax_left.set_ylabel('discharge (m^3 day^-1)')
    ax_left.set_title(' Channelized and distributed discharge (2000)')
    ax_left.set_ylabel('surface runoff(m^3 day^-1)')
    #plt.plo
    ax_left.plot(oneYear, waterFlux_m3d_grdmsk_timeseries   , label = 'distributed')
    ax_left.plot(oneYear, channelDischarge_m3d_grdmsk_timeseries  , label = 'channelized')
    ax_left.plot(oneYear, grounded_EWI, linestyle='dashed', color='green', linewidth='0.8', label = 'surface melt')
    ax_left.legend(loc='upper left')
    ax_right.legend(loc='upper right')
    
    
    fig1, ax1 = plt.subplots()
    ax1.set_xlabel('day')
    ax1.set_ylabel('discharge (m^3 day^-1)')
    ax1.set_title('2000')
    
    ax1.plot(oneYear, cummulative_waterFlux, label='distributed')
    ax1.plot(oneYear, cummulative_chDis, label='channelized')
    ax1.plot(oneYear, cummulative_disch, label='total discharge')
    ax1.plot(oneYear, cum_ewi, linestyle='dashed', label='surfaceMelt')
    
    ax1.legend()
    
    fig.savefig('ax1_figure.png', dpi=300)
    '''
    plt.scatter(xCell, yCell, s=15, c=groundedMask_2[0,:], cmap='Blues')
    margin_waterFlux = margin_waterFlux.astype('float')
    margin_waterFlux[margin_waterFlux == 0] = 'nan' # or use np.nan
    plt.scatter(xEdge,yEdge, s=15,c=margin_waterFlux[300,:], cmap='Reds_r')
    '''
    