# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 22:07:10 2018

@author: bimta
"""
import numpy as np
import sys
sys.path.insert(0,r'C:\Users\bimta\Documents\Arbeit\Malmö_Stad\Git\dem-filter')
from FME_Scripts.PMF import ProgressiveMorphologicalFilter as Filter
import matplotlib.pyplot as plt


raster = np.loadtxt("../Data/Merged/DEM_0.5_min_filtered_merged - kopia.asc")
raster[raster==-9999] = np.nan
raster = raster[200:500, 0:300]
"""
raster = np.ones((7,7))
i = 1
for row in raster: 
    row += i/10
    i += 1
raster[3,3] = 1.6
"""
parameters = {'c' : 0.6,                     # The cell size of the input raster
              'kernel_radius' : 10,           # The kernel size of the opening
              'initial_cutoff' : 0.4,        # The slope threshold of the first filtr iteration
              'average_sigma' : 8,           # The gaussian function used to average out local slope
              'dh0' : 0.1,                   # The slope threshold for flat areas
              'hole_cutoff' : -0.2}          # Threshold for individual holes beneath median in 3x3 kernel

MyFilter = Filter(raster, parameters)

import time
start = time.process_time()
final = MyFilter.filter()
end = time.process_time()        
print("Done after {0:.2f} seconds.".format(end - start))

initial = MyFilter.initial_filtered
scaling = MyFilter.scaling_matrix
final = MyFilter.final_filtered


def plotter(raster_list):
    nr_plots = len(raster_list)
    fig, ax = plt.subplots(1, nr_plots ,figsize=(30, 9), sharex = True, sharey = True)

    ind = 0
    for raster in raster_list:
        Name = [name for name in globals() if globals()[name] is raster][0]
        img = ax[ind].imshow(raster, cmap="viridis")
        ax[ind].set_title('{}'.format(Name))
        fig.colorbar(img, ax = ax[ind], orientation = 'horizontal')
        ax[ind].set_axis_off()
        ind += 1

    plt.tight_layout()
    plt.show()
    return

from scipy import ndimage
opened = ndimage.morphology.grey_opening(raster, size = (30,30))

plotter([raster, initial, scaling, final])

"""
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
def TINinterpolation(raster):
    def rasterToList(raster):
        leni = raster.shape[0]
        lenj = raster.shape[1]
        i = np.arange(0, leni)
        j = np.arange(0, leni)
        jj, ii = np.meshgrid(i, j)
    
        x = jj.reshape(leni * lenj,)
        y = ii.reshape(leni * lenj,)
        z = raster.reshape(leni * lenj,)
    
        mask = np.isnan(z)
    
        xknown = x[~mask]
        yknown = y[~mask]
        zknown = z[~mask]
    
        coordinates = np.stack((yknown,xknown)).T
        tointerpolate = np.stack((y,x)).T
        return coordinates, zknown, tointerpolate, leni, lenj
    
    coordinates, values, tointerpolate, leni, lenj = rasterToList(raster)
    
    tri = Delaunay(coordinates)
    Interpolator = LinearNDInterpolator(tri, values)
    interpolated = Interpolator(tointerpolate).reshape(leni,lenj)[1:-1,1:-1]
    return interpolated


interpolated = TINinterpolation(new_out)

plotter([new_out, interpolated])
"""