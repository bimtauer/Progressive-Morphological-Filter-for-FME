# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 08:35:21 2018

@author: bimta
"""

import time
import numpy as np
import matplotlib.pyplot as plt

#raster = np.loadtxt("Data/raster.asc")

#TODO: Implement kdTree distance computation in place of loop

def iwd(x,y,v,grid,power):
    start = time.process_time()
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            distance = np.sqrt((x-i)**2+(y-j)**2)
            if (distance**power).min()==0: 
                grid[i,j] = v[(distance**power).argmin()]
            else:
                total = np.sum(1/(distance**power))
                grid[i,j] = np.sum(v/(distance**power)/total)
    end = time.process_time()       
    print("Done after {0:.2f} seconds.".format(end - start))
    return grid

def findReplace(where, find, replace):
    if np.isnan(find):
        output = np.where(np.isnan(where), replace, where)
    else:
        output = np.where(where == find, replace, where)
    return output

#Own data
#test = raster[0:300, 0:300]
#test = findReplace(raster, -9999, np.nan)

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
    
    z = z[~mask]
    x = x[~mask]
    y = y[~mask]
    
    grid = np.zeros((leni, lenj))
    return x, y, z, grid

x, y, z, grid = rasterToList(output[0:200, 0:200])

grid = iwd(x,y,z,grid,2)

np.count_nonzero(~np.isnan(test))
#Out 396317
np.count_nonzero(np.isnan(test))
#Out 93683
###############################################################################





"""
img1 = plt.imshow(grid.T,origin='lower',interpolation='nearest',cmap='jet')
#plt.scatter(x,y,c='black',s=120)
plt.xlim(0,grid.shape[0])
plt.ylim(0,grid.shape[1])
plt.colorbar(img1)
plt.grid()
plt.show()
"""