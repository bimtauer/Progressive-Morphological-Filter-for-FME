# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 13:52:49 2018

@author: bimta
"""

import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt


raster = np.loadtxt("../Data/Merged/DEM_0.5_min_filtered_merged - kopia.asc")

raster[raster==-9999] = np.nan

input_raster = raster


kernel_radius = 7
c = 0.5

def progressiveMorphologicalfilter(input_raster, slope_threshold, maxk):
        def circularFootprint(cells):
            y,x = np.ogrid[-cells: cells+1, -cells: cells+1]
            footprint = x**2+y**2 <= cells**2
            return footprint

        #nan to 9999
        last_surface = np.where(np.isnan(input_raster), 9999, input_raster)

        #The mask we use to indicate non-ground points
        mask = np.zeros(input_raster.shape)
        k = 1
        while k <= maxk:
            #Generate footprint
            f = circularFootprint(k)
            #Opening
            this_surface = ndimage.morphology.grey_opening(last_surface, footprint=f)
            #Increasing maxdh by footprint radius in cells
            dhmax = k * c * slope_threshold
            #Only mask those below cutoff
            mask = np.where(last_surface - this_surface > dhmax, 1, mask)
            last_surface = this_surface
            k += 1          #one step further

        output = np.where(mask, np.nan, input_raster)
        return output
    
out = progressiveMorphologicalfilter(input_raster, 0.4, 14)

def plotter(raster_list):
    rows = len(raster_list)
    cols = len(raster_list[0])
    fig, ax = plt.subplots(rows, cols ,figsize=(30, 9), sharex = True, sharey = True)

    ind = 0
    for l in raster_list:
        j = 0
        for raster in l:
            #Name = [name for name in globals() if globals()[name] is l][0]
            img = ax[ind,j].imshow(raster, cmap="viridis")
            #ax[ind].set_title('{}'.format(Name))
            #fig.colorbar(img, ax = ax[ind,j], orientation = 'horizontal')
            ax[ind,j].set_axis_off()
            j += 1
        ind += 1
    plt.tight_layout()
    plt.show()
    return

###############################################################################
def progopeningWL(input_raster, maxk):
    last_surface = np.copy(input_raster)
    k = 1
    while k <= maxk:
        #Generate footprint
        w = 2*k + 1
        window = (w,w)
        #Opening
        this_surface = ndimage.morphology.grey_opening(last_surface, size = window)
       
        last_surface = this_surface
        k += 1          #one step further

    return last_surface

def progopening(input_raster, maxk):
    last_surface = np.copy(input_raster)
    k = 1
    while k <= maxk:
        #Generate footprint
        w = 2*k + 1
        window = (w,w)
        #Opening
        this_surface = ndimage.morphology.grey_opening(last_surface, size = window)
       
        k += 1          #one step further

    return this_surface

def opening(input_raster, maxk):
    w = 2*maxk + 1
    window = (w,w)
    opened = ndimage.morphology.grey_opening(input_raster, size = window)
    return opened


progopened_wl = []
progopened = []
opened = []
diffs = []
diffs2 = []
for maxk in range(1,9,1):
    print(maxk)
    propwl = progopeningWL(input_raster, maxk)
    progopened_wl.append(propwl)
    prop = progopening(input_raster, maxk)
    progopened.append(prop)
    op = opening(input_raster, maxk)
    opened.append(op)
    propwl_nan = np.where(np.isnan(propwl), 99, propwl)
    prop_nan = np.where(np.isnan(prop), 99, prop)
    op_nan = np.where(np.isnan(op), 99, op)
    diff = np.where(propwl_nan != prop_nan, 1, 0)
    diffs.append(diff)
    diff2 = np.where(prop_nan != op_nan, 1, 0)
    diffs2.append(diff2)
    
    
    
raster_list = [progopened_wl, progopened, opened, diffs, diffs2]

plotter(raster_list)