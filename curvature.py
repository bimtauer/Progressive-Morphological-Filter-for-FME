# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 13:08:55 2018

@author: timbau
"""
#For Reference
"""
with open("GEOTIFF.asc") as data:
    content = data.readlines()
"""    

"""
FME Structure start:

import fmeobjects

class FeatureProcessor(object):
    def __init_(self):
	pass

    def input(self, feature):
	self.featureList.append(feature)
	pass

    def close(self):
	pass

"""

import numpy as np
import time
import pandas as pd


#Manually deleted first 6 text lines
raster = np.loadtxt("raster.asc")

#Grid size
a = 0.5

def normalVectorEstimation(raster):
    #Creating 3d matrix to store resulting vectors
    out = np.zeros((raster.shape[0], raster.shape[1], 3))
    
    #Padding Raster
    padded = np.pad(raster, 1, "median")
    
    rows = padded.shape[0]
    cols = padded.shape[1]
    
    print("Beginning normal vector estimation...")
    start = time.clock()
    for i in range(1, rows-1):
        for j in range(1, cols-1):
            #Surrounding points
            top = padded[i-1,j]
            bot = padded[i+1,j]
            left = padded[i,j-1]
            right = padded[i,j+1]
            
            # calculating surface normal based on 4 surrounding points:
            normal = np.array([(left - right), (top - bot), (2*a)]) 
            mag = np.sqrt(normal.dot(normal))
            unit_normal = normal / mag
            
            out[i-1, j-1] = unit_normal
    
    end = time.clock()        
    print("Done after {0:.2f} seconds.".format(end - start))
    return out


normals = normalVectorEstimation(raster)

# I also want a reliable ascii out writer

"""from astropy.io import ascii
ascii.write(slopes, "test.asc")

"""

""" In resulting document remove col lines and add this before points:
    
ncols 700
nrows 700
xllcorner 116549.99999999983
yllcorner 6164149.999999888
cellsize 0.5000000000021828
nodata_value -9999
"""

# Need an x y coordinate grid

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import seaborn as sns

def sectionVisualizer(yy, xx, smoothing=1): # ii, jj
    #expects tuple like (0,3), (0,3) as input
    
    #Replace -9999 with NaN
    raster[raster==-9999]=np.nan
    
    #Setting up points (and vectors) collections possible length
    points = np.zeros(( int((xx[1]-xx[0])/smoothing)*int((yy[1]-yy[0])/smoothing) , 6))
    
    #Counting and looping
    pointID = 0
    for i in range(yy[0], yy[1], smoothing):
        for j in range(xx[0], xx[1], smoothing):
            points[pointID, 0] = j * a                  #X Point coordinate
            points[pointID, 1] = i * a                  #Y Point coordinate
            points[pointID, 2] = raster[i,j]            #Z Point coordinate
            points[pointID, 3] = normals[i, j, 0]       #X component of normal
            points[pointID, 4] = normals[i, j, 1]       #Y component of normal
            points[pointID, 5] = normals[i, j, 2]       #Z component of normal
            
            pointID += 1

    #Plotting    
    fig = plt.figure()
    
    ax = fig.gca(projection='3d') 
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_zlabel('Elevation (m)')
    
    #Scaling of axes
    #ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([2, 2, 1, 1]))
    #ax.set_zticks([0,10,20])
    
    #Scaling for colors
    mynorm = plt.Normalize(vmin=0, vmax=10)
    
    #Plot surface
    ax.plot_trisurf(points[:, 0], points[:, 1], points[:, 2],  cmap=plt.cm.viridis, norm=mynorm, linewidth=0.2, antialiased=True)
    #ax.bar3d(points[:, 1], points[:, 0], points[:, 2], 0, a, a,  shade=True)
    #surf=ax.plot_trisurf(points[:, 1], points[:, 0], points[:, 2],  cmap=plt.cm.viridis, linewidth=0.2, antialiased=True)
    #fig.colorbar(surf)
    
    #Plot vectors
    av = fig.gca(projection='3d')
    av.quiver(points[:, 0], points[:, 1], points[:, 2], points[:, 3], points[:, 4], points[:, 5], length= 0.2, color= "red")
    plt.show()    
    return

sectionVisualizer((680,700),(0,20))
 
        
    
    
""" Not necessary
xx, yy = np.meshgrid(x, y)         
r3coordinates = np.dstack((xx,yy,raster))
"""

            
