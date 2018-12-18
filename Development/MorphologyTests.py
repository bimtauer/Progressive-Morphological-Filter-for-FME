# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 12:04:48 2018

@author: bimta
"""
import numpy as np
from Plotting import mPlotter
from scipy import ndimage
# Understanding erosion, opening and progressive opening

raster = np.loadtxt("../Data/Merged/DEM_0.5_min_filtered_merged - kopia.asc")

raster[raster==-9999] = np.nan

test_area = raster[250:260, 364:374]

# I need a trisurf plotter first and test rasters

# SLope raster
s = np.ones((7,7))
for row in range(len(s)):
    s[:,row] *= (row + 1)
    
# Object raster
o = np.ones((7,7))
o[2:5,2:5] = 4
o[3,5] = 4

#Now I want to compare several operations:
# 1. progressive opening over same surface
# 2. progressive opening with last output
# 3. opening with mxk

def progOpening(raster):
    maxk = 3
    k = 1
    my_dict = {}
    my_dict['Input'] = raster
    while k <= maxk:
        w = k*2 +1
        window = (w,w)
        opened = ndimage.morphology.grey_opening(raster, size = window)
        my_dict[w] = opened
        k += 1          #one step further
    return my_dict

def progOpeningWL(raster):
    maxk = 3
    k = 1
    my_dict = {}
    my_dict['Input'] = raster

    last_surface = np.copy(raster)
    while k <= maxk:
        w = k*2 +1
        window = (w,w)
        this_surface = ndimage.morphology.grey_opening(last_surface, size = window)
        my_dict[w] = this_surface
        last_surface = this_surface
        k += 1          #one step further
    return my_dict

def differences(d1, d2):
    my_dict = {}
    for (k1,v1), (k2,v2) in zip(d1.items(), d2.items()):
        diff = v1 - v2
        my_dict[k1] = diff
    return my_dict

def makeList(raster):
    d1 = progOpening(raster)
    d2 = progOpeningWL(raster)
    d3 = differences(d1, d2)
    l = [d1, d2, d3]
    return l
mPlotter(makeList(test_area))
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            