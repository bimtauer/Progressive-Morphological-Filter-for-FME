# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 20:57:18 2018

@author: bimta
"""

import matplotlib.pyplot as plt
import numpy as np

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


def mPlotter(raster_lod):
    rows = len(raster_lod)
    cols = len(raster_lod[0])
    fig, ax = plt.subplots(rows, cols,figsize=(30, 9), sharex = True, sharey = True)
    i = 0
    for d in raster_lod:
        j = 0
        for key, raster in d.items():
            name = key
            ax[i,j].imshow(raster, cmap="winter")
            ax[i,j].set_title('{}'.format(name))
            ax[i,j].set_axis_off()
            for (j2,i2),label in np.ndenumerate(raster):
                ax[i,j].text(i2,j2, "{:.1f}".format(label),ha='center',va='center')
            j += 1
        #ax[i,:].annotate('test',size='large', ha='right', va='center')
        i += 1
    plt.tight_layout()
    plt.show()
    return


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