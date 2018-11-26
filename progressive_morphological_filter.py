# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 16:10:57 2018

@author: bimta
"""
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

raster = np.loadtxt("Data/raster.asc")

c = 0.5

def linearInterpolation(raster):
    interpolated = np.where(raster==-9999, np.nan, raster)
    mask = np.isnan(interpolated)
    interpolated[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), raster[~mask])
    return interpolated

#Exchanges Raster Band Nodata (-9999) with np.nan
def nanFormatter(raster):
    raster = np.where(raster==-9999, np.nan, raster)
    return raster

#Replaces -9999 with 9999 for opening filter
def nanReplacer(raster):
    raster = np.where(raster==-9999, 9999, raster)
    return raster

def medianFilter(raster):
    median = ndimage.median_filter(raster, size = (3,3))
    return median

def progressiveMorphologicalFilter(raster, maxk=15, dh0=0.3):
    """
    fig, ax = plt.subplots(1, maxk)
    """    
    
    #Replace -9999 with 9999
    input_raster = nanReplacer(raster)
    processing_raster = np.copy(input_raster)

    
    mask = np.zeros(processing_raster.shape)
    k = 1
    while k <= maxk:
        wk = 2*k+1       #windows size
        window = (wk,wk)  
        
        #Opening
        opened = ndimage.morphology.grey_opening(processing_raster, size=window)
        
        #Increasing maxdh by window size - always the maximum slope to corner point
        dhmax = np.sqrt(k**2 + k**2) * c * dh0
        
        #Only return those below cutoff
        mask = np.where(processing_raster - opened > dhmax, 1, mask)
        """
        ax[k-1].imshow(mask, cmap="viridis")
        ax[k-1].set_title('At k = {}'.format(k))
        ax[k-1].set_axis_off()
        """
        k += 1          #one step further
    
    """
    plt.tight_layout()
    plt.show()
    """
    output = np.where(mask, np.nan, input_raster)
    output[output==9999]=np.nan
    
    
    #Filter holes
    median = medianFilter(output)
    output = np.where(output - median < -1, median, output)
    
    return output

#out = progressiveMorphologicalFilter(raster)

###############################################################################
from Slope_Erosion_Experiment import maxSlopeFilter
import itertools
# Test Areas:
Trees = raster[500:600,100:200]
Tunnels = raster[315:415, 370:470]
Holes = raster[200:300, 50:150]

def comparisonPlot(raster1, raster2, raster3, raster4):
    
    fig, ax = plt.subplots(2, 4 ,figsize=(30, 9), sharex = 'col', sharey = 'col')
    
    rasters = [raster1, raster2, raster3, raster4]
    index = 0
    for raster in rasters:
        Name = [name for name in globals() if globals()[name] is raster][0]
        morphout = progressiveMorphologicalFilter(raster)
        #erosionout = maxSlopeFilter(raster, (9,9), 0.3)
        #ys = itertools.cycle((morphout,erosionout))
        
        img1 = ax[0,index].imshow(nanFormatter(raster), cmap="viridis")
        ax[0,index].set_title('{} Original'.format(Name))
        ax[0,index].set_axis_off()
        
        img2 = ax[1,index].imshow(morphout, cmap="viridis")
        ax[1,index].set_title('{} Morph Filtered'.format(Name))
        fig.colorbar(img2, ax = ax[1,index])
        ax[1,index].set_axis_off()

        
        
        plt.tight_layout()
        plt.show()
        
        index += 1
    return ys



comparisonPlot(raster, Trees, Tunnels, Holes)

"""
#ys = itertools.cycle((morphout,morphout2))

fig, ax = plt.subplots(1, 4 ,figsize=(30, 9), sharex = 'col', sharey = 'col')
def onclick(event):
    for index in range(4):
        ax[index].clear()
        img2 = ax[index].imshow(next(ys), cmap="viridis")
        #ax[1,index].set_title('{} Filtered'.format(Name))
        #fig.colorbar(img2, ax = ax[index])
        #ax[index].set_axis_off()
    fig.canvas.draw()

cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.tight_layout()
plt.show()
"""


"""
###############################################################################
# Evaluation
import pandas as pd
def Evaluator():
    results_list = []
    for k in range(1, 31):
        for d in range(2, 7):
            d = d/10
            out = progressiveMorphologicalFilter(raster, maxk = k, dh0 = d)
            minz = np.nanmin(out)
            maxz = np.nanmax(out)
            nan = np.isnan(out).sum()
            mylist = [k, d, minz, maxz, nan]
            results_list.append(mylist)
    
    return results_list

evaluation = Evaluator()

df = pd.DataFrame(evaluation)
df.columns = ['maxk', 'dh0', "min", "max", "nan"]

df.pivot(index='maxk', columns='dh0', values='max')

import seaborn as sns; sns.set()
fig, ax = plt.subplots(2)
ax = sns.heatmap(df.pivot(index='maxk', columns='dh0', values='max'))
#ax[1] = sns.heatmap(df.pivot(index='maxk', columns='dh0', values='nan'))
plt.show()  
"""
