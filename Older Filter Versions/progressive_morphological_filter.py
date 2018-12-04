# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 16:10:57 2018

@author: bimta
"""
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from skimage.restoration import inpaint

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

# Inpainting for interpolation - takes very long
def imageInpainter(raster):
    mask = np.where(raster == -9999, 1, 0)
    raster_result = inpaint.inpaint_biharmonic(raster, mask)
    return raster_result

# To get rid of holes
def medianFilter(raster):
    median = ndimage.median_filter(raster, size = (3,3))
    return median

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
        #lastSurface = thisSurface
        
        
        
        k += 1          #one step further
    
    output = np.where(mask, np.nan, input_raster)
    output[output==9999]=np.nan
    """
    #Filter holes
    median = medianFilter(output)
    output = np.where(output - median < -1, median, output)
    """
    return output

out = progressiveMorphologicalFilter(raster[575:581,104:110])

###############################################################################
# Visualizing results for test areas
    
# Test Areas:
def testAreas(raster):
    Trees = raster[500:600,100:200]
    Tunnels = raster[315:415, 370:470]
    Holes = raster[200:300, 50:150]
    return raster, Trees, Tunnels, Holes

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
        
        index += 1
    plt.tight_layout()
    plt.show()
    return 

showRaster, Trees, Tunnels, Holes = testAreas(raster)

comparisonPlot(showRaster, Trees, Tunnels, Holes)

###############################################################################
#Stepwise

def stepwisePlot(raster, maxk):
    
    fig, ax = plt.subplots(3, maxk ,figsize=(30, 9), sharex = 'row', sharey = 'row')
    
    index = 0
    for k in range(1, maxk+1):
        
        wk = 2*k+1       #windows size
        window = (wk,wk) 
        dhmax = np.sqrt(k**2 + k**2) * c * dh0
        
        out = progressiveMorphologicalFilter(raster, k)
        op = ndimage.morphology.grey_opening(raster, size= window)
        
        
        img = ax[0,index].imshow(nanFormatter(raster), cmap="viridis")
        for i in range(raster.shape[0]):
            for j in range(raster.shape[1]):
                text = ax[0,index].text(j, i, format(raster[i, j], '.2f'),
                               ha="center", va="center", color="w")
        fig.colorbar(img, ax = ax[0,index])
        ax[0,index].set_axis_off()
        
        
        img1 = ax[1,index].imshow(nanFormatter(out), cmap="viridis")
        for i in range(raster.shape[0]):
            for j in range(raster.shape[1]):
                text = ax[1,index].text(j, i, format(raster[i,j]-op[i,j], '.2f'),
                               ha="center", va="center", color="w")
        ax[1,index].set_title('For k: {}'.format(k))
        fig.colorbar(img1, ax = ax[1,index])
        ax[1,index].set_axis_off()
        

        img2 = ax[2,index].imshow(nanFormatter(op), cmap="viridis")
        for i in range(raster.shape[0]):
            for j in range(raster.shape[1]):
                text = ax[2,index].text(j, i, format(op[i, j], '.2f'),
                               ha="center", va="center", color="w")
        ax[2,index].set_title('Cutoff: {}'.format(dhmax))
        fig.colorbar(img1, ax = ax[2,index])
        ax[2,index].set_axis_off()
        index += 1
        
    plt.tight_layout()
    plt.show()
    return 

stepwisePlot(raster[575:581,104:110], 3)




###############################################################################
# Interactive Plot
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


import seaborn as sns; sns.set()
fig, (ax,ax2) = plt.subplots(ncols=2)
#fig.subplots_adjust(wspace=0.01)
sns.heatmap(df.pivot(index='maxk', columns='dh0', values='max'), cmap="rocket", ax=ax, cbar=False)
fig.colorbar(ax.collections[0], ax=ax,location="right", use_gridspec=False, pad=0.2)
sns.heatmap(df.pivot(index='maxk', columns='dh0', values='nan'), cmap="icefire", ax=ax2, cbar=False)
fig.colorbar(ax2.collections[0], ax=ax2,location="right", use_gridspec=False, pad=0.2)
ax2.yaxis.tick_right()
ax2.tick_params(rotation=0)
plt.show()
"""