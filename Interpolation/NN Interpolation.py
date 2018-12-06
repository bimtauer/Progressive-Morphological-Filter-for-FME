# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 10:20:46 2018

@author: bimta
"""

from scipy.interpolate import NearestNDInterpolator
import numpy as np
import matplotlib.pyplot as plt

raster = np.loadtxt("../Outputs/MorphFilteredtestIn.asc")

#Own data
test = np.where(raster == -9999, np.nan, raster)

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

    #xunknown = x[mask]
    #yunknown = y[mask]

    tointerpolate = np.stack((y,x)).T

    return coordinates, zknown, tointerpolate, leni, lenj

coordinates, values, tointerpolate, leni, lenj = rasterToList(test)

""" Usage of NN:
x : (Npoints, Ndims) ndarray of floats

    Data point coordinates.
y : (Npoints,) ndarray of float or complex

    Data values.
"""

Interpolator = NearestNDInterpolator(coordinates, values)
interpolation = Interpolator(tointerpolate).reshape(leni,lenj)

interpolated = interpolation - test
interpolated = findReplace(interpolated, np.nan, 1)

def plotter(raster_list):
    nr_plots = len(raster_list)
    fig, ax = plt.subplots(1, nr_plots ,figsize=(30, 9), sharex = True, sharey = True)

    ind = 0
    for raster in raster_list:
        Name = [name for name in globals() if globals()[name] is raster][0]
        img = ax[ind].imshow(raster, cmap="viridis")
        ax[ind].set_title('{}'.format(Name))
        """
        for i in range(raster.shape[0]):
            for j in range(raster.shape[1]):
                text = ax[i].text(j, i, format(raster[i, j], '.1f'), ha="center", va="center", color="w")
        """

        fig.colorbar(img, ax = ax[ind], orientation = 'horizontal')
        ax[ind].set_axis_off()

        ind += 1

    plt.tight_layout()
    plt.show()
    return

plotter([test,interpolation, interpolated])
"""
###############################################################################
#Example
lat = np.linspace(2, 6, 10)
lon = np.linspace(5, 9, 14)

latM, lonM = np.meshgrid(lat, lon)  # M is for Matrix

dataM = np.sin(latM)*np.cos(lonM)  # example of data, Matrix form


points = np.array((latM.flatten(), lonM.flatten())).T
print( points.shape )
# >>> (140, 2)

f_nearest = NearestNDInterpolator(points, dataM.flatten())
f_nearest(5, 5)
"""
