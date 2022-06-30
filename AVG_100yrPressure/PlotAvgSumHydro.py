#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 14:14:53 2022

@author: abbyroat
"""

from __future__ import absolute_import, division, print_function, unicode_literals


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import NearestNDInterpolator 
import itertools as it
'''
Import files from ProcAvgSumHydro.py
'''
EP_MEAN_summer = np.load('EP_MEAN_summer.npy')
WP_MEAN_summer = np.load('WP_MEAN_summer.npy')
EWI_SUM_summer = np.load('EWI_SUM_summer.npy')

'''
Divide data into groups by decade
'''

pogi = np.load('WP_MEAN_summer.npy')
mvdp = np.load('EWI_SUM_summer.npy')
EWI_SUM_summer_groups = [EWI_SUM_summer[x:x+10] for x in range(0, len(EWI_SUM_summer), 10)]# EWI 
EP_MEAN_summer_groups = [EP_MEAN_summer[x:x+10] for x in range(0, len(EP_MEAN_summer), 10)]# effective pressure
WP_MEAN_summer_groups = [WP_MEAN_summer[x:x+10] for x in range(0, len(WP_MEAN_summer), 10)]# water pressure


'''
Create lines of best fit
'''
#EFFECTIVE PRESSURE
slopes_EP = []
intercepts_EP = []
for EWI_SUM_summer_group, EP_MEAN_summer_group in zip(EWI_SUM_summer_groups, EP_MEAN_summer_groups):
    slope_EP, intercept_EP = np.polyfit(EWI_SUM_summer_group, EP_MEAN_summer_group, 1)     
    slopes_EP.append(slope_EP)
    intercepts_EP.append(intercept_EP)
# WATER PRESSURE
slopes_WP = []
intercepts_WP = []
for EWI_SUM_summer_group, WP_MEAN_summer_group in zip(EWI_SUM_summer_groups, WP_MEAN_summer_groups):
    slope_WP, intercept_WP = np.polyfit(EWI_SUM_summer_group, WP_MEAN_summer_group, 1)     
    slopes_WP.append(slope_WP)
    intercepts_WP.append(intercept_WP)
        
''' 
plotting
'''
#plot points
# itertools.cycle - specify colors you want to loop over you have to import itertools
fig, ax = plt.subplots()
colors = it.cycle(['mistyrose','coral', 'tomato', 'red', 'violet', 'mediumorchid', 'darkorchid', 'slateblue','blue', 'black'])
for EWI_SUM_summerDecade,EP_MEAN_summerDecade in zip(EWI_SUM_summer_groups,EP_MEAN_summer_groups):
    ax.scatter(EWI_SUM_summerDecade,EP_MEAN_summerDecade, color=next(colors)) 
    #ax.plot(EWI_2000s_w3, EWI_2000s_w3*slope2000_3k + intercept2000_3k, 'mistyrose', label='2000s slope %f'%slope2000_3k)
    ax.set_xlabel("total RACMO summer runoff(m^3 s^-1)")
    ax.set_ylabel('mean summer EFFECTIVE pressure (Pa)')
    ax.set_title('summer surface melt vs mean summer EFFECTIVE pressure')

fig2, ax2 = plt.subplots()
for EWI_SUM_summerDecade,WP_MEAN_summerDecade in zip(EWI_SUM_summer_groups,WP_MEAN_summer_groups):
    ax2.scatter(EWI_SUM_summerDecade,WP_MEAN_summerDecade, color=next(colors)) 
    #ax2.plot(EWI_2000s_w3, EWI_2000s_w3*slope2000_3k + intercept2000_3k, 'mistyrose', label='2000s slope %f'%slope2000_3k)
    ax2.set_xlabel("total RACMO summer runoff(m^3 s^-1)")
    ax2.set_ylabel(' mean summer WATER pressure (Pa)')
    ax2.set_title('summer surface melt vs mean summer WATER pressure')


