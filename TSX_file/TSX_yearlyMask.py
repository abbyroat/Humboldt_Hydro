#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 13:55:18 2021

@author: abbyroat
"""

from __future__ import absolute_import, division, print_function, unicode_literals


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

#outputfiles = ['/Users/abbyroat/Desktop/Badger_outputs_2021/2km/DarcyFlowLaw/SummerRuns/2014_TSX_iceVelocity/TSX_W79.75N_02May14_13May14_11-46-14_vv_v03.0.nc',  '/Users/abbyroat/Desktop/Badger_outputs_2021/2km/DarcyFlowLaw/SummerRuns/2014_TSX_iceVelocity/TSX_W79.75N_26Jun14_07Jul14_11-46-17_vv_v03.0.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/2km/DarcyFlowLaw/SummerRuns/2014_TSX_iceVelocity/TSX_W79.75N_09Aug14_20Aug14_11-46-19_vv_v03.0.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/2km/DarcyFlowLaw/SummerRuns/2014_TSX_iceVelocity/TSX_W79.75N_19Dec14_30Dec14_11-46-19_vv_v03.0.nc']# 2014
#outputfiles = ['/Users/abbyroat/Desktop/Badger_outputs_2021/2km/DarcyFlowLaw/SummerRuns/2017_TSX_iceVelocity/TSX_W79.75N_06May17_17May17_11-46-34_vv_v03.0.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/2km/DarcyFlowLaw/SummerRuns/2017_TSX_iceVelocity/TSX_W79.75N_30Jun17_11Jul17_11-46-36_vv_v03.0.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/2km/DarcyFlowLaw/SummerRuns/2017_TSX_iceVelocity/TSX_W79.75N_04Sep17_15Sep17_11-46-39_vv_v03.0.nc'] # 2017
'''
outputfiles = ['/Users/abbyroat/Desktop/Badger_outputs_2021/2km/DarcyFlowLaw/SummerRuns/2016_TSX_iceVelocity/TSX_W79.75N_16Apr16_27Apr16_11-46-25_vv_v03.0.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/2km/DarcyFlowLaw/SummerRuns/2016_TSX_iceVelocity/TSX_W79.75N_10Jun16_21Jun16_11-46-27_vv_v03.0.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/2km/DarcyFlowLaw/SummerRuns/2016_TSX_iceVelocity/TSX_W79.75N_26Aug16_06Sep16_11-46-32_vv_v03.0.nc', '/Users/abbyroat/Desktop/Badger_outputs_2021/2km/DarcyFlowLaw/SummerRuns/2016_TSX_iceVelocity/TSX_W79.75N_14Dec16_25Dec16_11-46-33_vv_v03.0.nc']
for file_tsx in outputfiles:
    f_tsx=Dataset(file_tsx, 'r')
    f_tsx.set_auto_mask(False)
    x_tsx = f_tsx.variables['x'][:]
    y_tsx = f_tsx.variables['y'][:]
    iceVelocity_tsx = f_tsx.variables['Band1'][:]
    
    iceVelocity_tsx_list = iceVelocity_tsx.flatten()
    iceVelocity_tsx_list[iceVelocity_tsx_list <=20] = 'nan'
    iceVelocity_tsx_list_noNaN = iceVelocity_tsx_list[ ~np.isnan(iceVelocity_tsx_list)]
    mean_iceVelocity_tsx =np.mean(iceVelocity_tsx_list_noNaN)
    
    print(mean_iceVelocity_tsx)
'''

average_tsx = [326.3525, 578.63007, 355.21393, 324.3376]
average_tsx_2017 = [393.44943, 470.0972, 333.6774]
average_tsx_2016 = [372.06464, 423.25034, 308.4964, 351.66403]
    
file = '/Users/abbyroat/Desktop/Badger_outputs_2021/2km/DarcyFlowLaw/SummerRuns/2014_TSX_iceVelocity/TSX_W79.75N_19Dec14_30Dec14_11-46-19_vv_v03.0.nc'
f_mask = Dataset(file, 'r')
f_mask.set_auto_mask(False)

iceVelocity = f_mask.variables['Band1'][:]
np.save('IceVelocity_Dec2014', iceVelocity)


IV_april = np.load('IceVelocity_April2016.npy')
IV_april_list = IV_april*1#.flatten()
IV_april_list[IV_april_list <=20] = 0
IV_april_list[IV_april_list >20] = 1

IV_May = np.load('IceVelocity_May2017.npy')
IV_May_list = IV_May*1#.flatten()
IV_May_list[IV_May_list <=20] = 0
IV_May_list[IV_May_list >20] = 1


IV_june = np.load('IceVelocity_June2017.npy')
IV_june_list = IV_june#1.flatten()
IV_june_list[IV_june_list <=20] = 0
IV_june_list[IV_june_list >20] = 1

IV_august = np.load('IceVelocity_August2016.npy')
IV_august_list = IV_august#.flatten()
IV_august_list[IV_august_list <=20] = 0
IV_august_list[IV_august_list >20] = 1

IV_sept = np.load('IceVelocity_Sept2017.npy')
IV_sept_list = IV_sept#1.flatten()
IV_sept_list[IV_sept_list <=20] = 0
IV_sept_list[IV_sept_list >20] = 1

IV_dec = np.load('IceVelocity_Dec2016.npy')
IV_dec_list = IV_dec#.flatten()
IV_dec_list[IV_dec_list <=20] = 0
IV_dec_list[IV_dec_list >20] = 1

mask_2017 = IV_june_list* IV_May_list*IV_sept_list
np.save('TSXmask_2017', mask_2017)

# when I save all the masks I can just write it as a new variable to one of the .nc files and then I can redo the interpolation.

