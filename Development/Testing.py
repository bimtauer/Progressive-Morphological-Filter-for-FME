# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 10:04:12 2018

@author: bimta
"""

from FME_Scripts.PMF import ProgressiveMorphologicalFilter
import numpy as np

#Old test area aerial and terrestrial combined
#raster = np.loadtxt(r"./Merged/DEM_0.5_min_filtered_merged.asc", skiprows = 6)
#raster[raster == -9999] = np.nan
#raster = raster[560:600, 95:135]    #Problem Trees
#raster = raster[290:360, 360:430]   #Problem House


#New test area inner city terrestrial
raster = np.loadtxt(r"./Data/GEOTIFF2.asc", skiprows = 6)
raster[raster == -9999] = np.nan




#################################################################################

parameters = {'c' : 0.5,                     # The cell size of the input raster
              'kernel_radius' : 8,          # The kernel size of the first filter iteration
              'initial_cutoff' : 0.5,        # The slope threshold of the first filter iteration
              'average_sigma' : 7,           # The gaussian function used to average out local slope
              'dh0' : 0.1,                   # The slope threshold for flat areas
              'hole_cutoff' : -1}          # Threshold for individual holes beneath median in 3x3 kernel

MyFilter = ProgressiveMorphologicalFilter(raster, parameters)

output = MyFilter.filter()
hole_filtered = MyFilter.hole_filtered
interpolated = MyFilter.interpolated
initial = MyFilter.initial_filtered
scaling = MyFilter.scaling_matrix
second = MyFilter.second_filtered

import matplotlib as mpl
import matplotlib.pyplot as plt
def plotter(raster_list):
    nr_plots = len(raster_list)
    fig, ax = plt.subplots(1, nr_plots ,figsize=(30, 9), sharex = True, sharey = True)

    ind = 0
    for raster in raster_list:
        Name = [name for name in globals() if globals()[name] is raster][0]
        img = ax[ind].imshow(raster, cmap='viridis')
        ax[ind].set_title('{}'.format(Name))
        fig.colorbar(img, ax = ax[ind], orientation = 'horizontal')
        ax[ind].set_axis_off()
        ind += 1

    plt.tight_layout()
    plt.show()
    return

plotter([interpolated, hole_filtered, initial, scaling, second, output])

###############################################################################
from scipy import ndimage
def progressiveMorphologicalfilter(maxk, last_surface, slope_threshold):
        openings = []
        diffs = []
        masks = []
        dhmaxs = []
        #The mask we use to indicate non-ground points
        mask = np.ones(last_surface.shape, dtype=bool)
        k = 1
        while k <= maxk:
            #Generate window
            w = k*2 + 1
            window = (w,w)
            #Opening
            this_surface = ndimage.morphology.grey_opening(last_surface, size = window)
            openings.append(this_surface)
            #Increasing dhmax by diagonal windowf extent in cells
            dhmax = np.sqrt(k**2 + k**2) * 0.5 * slope_threshold
            dhmaxs.append(dhmax)
            diff = last_surface - this_surface
            diffs.append(diff)
            mask = np.where(last_surface - this_surface > dhmax, 0, mask)
            masks.append(mask)
            last_surface = this_surface
            k += 1          #one step further
        return [openings, diffs, dhmaxs, masks]



results = progressiveMorphologicalfilter(16, interpolated, scaling)

filtereds = []
for mask in results[3]:
    filtered = np.where(mask, raster, np.nan)
    filtereds.append(filtered)

results[3] = filtereds

def mplotter(lol):
    rows = len(lol)
    cols = len(lol[0])
    fig, ax = plt.subplots(rows, cols ,figsize=(30, 9), sharex = True, sharey = True)
    
    #Consistent scaling
    scaling = np.empty((rows, 2))
    i = 0
    for l in lol:    #Each list of raster
        mins = []
        maxs = []
        for raster in l:     #Each raster
            mn = np.nanmin(raster)
            mx = np.nanmax(raster)
            mins.append(mn)
            maxs.append(mx)
        scaling[i,0] = min(mins)
        scaling[i,1] = max(maxs)
        i += 1
    
    norms = []
    for row in scaling:
        norm = mpl.colors.Normalize(vmin=row[0],vmax=row[1])
        norms.append(norm)

    i = 0
    for l in lol:    #Each list of raster
        j = 0 
        for raster in l:     #Each raster
            img = ax[i,j].imshow(raster, norm = norms[i], cmap='viridis')
            ax[i,j].set_axis_off()
            j += 1
        #fig.colorbar(img, ax = ax[i,:], orientation = 'horizontal')
        i += 1
    
    cols = ['k = {}'.format(col+1) for col in range(cols)]
        
    plt.tight_layout()
    plt.show()

    return


mplotter(results)


