# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 16:46:15 2018

@author: timbau
"""

import numpy as np
from scipy import ndimage, interpolate
import matplotlib.pyplot as plt
from skimage.restoration import inpaint

raster = np.loadtxt("Data/raster.asc")

c = 0.5


#Exchanges Raster Band Nodata (-9999) with np.nan
def nanFormatter(raster):
    raster = np.where(raster==-9999, np.nan, raster)
    return raster

#Replaces -9999 with 9999 for opening filter
def nanReplacer(raster):
    raster = np.where(raster==-9999, 9999, raster)
    return raster

###############################################################################
# To get rid of holes
def medianFilter(raster):
    median = ndimage.median_filter(raster, size = (3,3))
    output = np.where(output - median < -1, median, output)
    return output

def linearInterpolation(raster):
    mask = np.isnan(raster)
    raster[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), raster[~mask])
    return raster

z1 = normals[:,:,2]
interp = linearInterpolation(z)

x = np.arange(z.shape[0])
y = np.arange(z.shape[1])
f = interpolate.NearestNDInterpolator(x, y, z)
znew = f(x, y)


average_slope = ndimage.uniform_filter(z[2:10,3:10], size = (3,3))

def progressiveMorphologicalFilter(raster, maxk=15, dh0=0.3):
    #Replace -9999 with 9999
    input_raster = nanReplacer(raster)
    lastSurface = np.copy(input_raster)

    mask = np.zeros(input_raster.shape)
    k = 1
    while k <= maxk:
        wk = 2*k+1       #windows size
        window = (wk,wk)  
        
        #Opening
        thisSurface = ndimage.morphology.grey_opening(lastSurface, size=window)
        
        #Increasing maxdh by window size - always the maximum slope to corner point
        dhmax = np.sqrt(k**2 + k**2) * c * dh0
        
        #Only return those below cutoff
        mask = np.where(lastSurface - thisSurface > dhmax, 1, mask)
        
        k += 1          #one step further
    
    output = np.where(mask, np.nan, input_raster)
    output[output==9999]=np.nan
    return output

out = progressiveMorphologicalFilter(raster[575:581,104:110])