# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 09:58:10 2018

@author: Tim Bauer (tim.bauer@gmx.us)
"""

import numpy as np
from scipy.interpolate import NearestNDInterpolator

#Estimates the z component of the surface normal for each raster cell using
#the left, right, top, and bottom cell as approximators. Input is supposed to
#be free of nans
def zComponentEstimation(raster):
    raster = np.pad(raster, 1, "reflect")
    left = raster[:,:-2]
    right = raster[:,2:]
    top = raster[:-2,:]
    bot = raster[2:,:]
    x = (left - right)[1:-1,:]
    y = (top - bot)[:,1:-1]
    magnitudes = np.sqrt(x**2 + y**2 + 1)
    z = 1/magnitudes
    return z

#Takes a raster containing nan as input and interpolates all holes with NN
def nnInterpolation(raster):
    def rasterToList(raster):
        leni = raster.shape[0]
        lenj = raster.shape[1]

        #Retrieve indexes as meshgrid of coordinates
        i = np.arange(0, leni)
        j = np.arange(0, leni)
        jj, ii = np.meshgrid(i, j)

        #Stack array rows
        x = np.hstack(jj)
        y = np.hstack(ii)
        z = np.hstack(raster)
        mask = np.isnan(z)

        #Mask nan
        xknown = x[~mask]
        yknown = y[~mask]
        zknown = z[~mask]

        coordinates = np.stack((yknown,xknown)).T
        tointerpolate = np.stack((y,x)).T

        return coordinates, zknown, tointerpolate, leni, lenj

    coordinates, values, tointerpolate, leni, lenj = rasterToList(raster)
    Interpolator = NearestNDInterpolator(coordinates, values)
    interpolation = Interpolator(tointerpolate).reshape(leni,lenj)
    return interpolation

""" Example usage:
interpolated_raster = nnInterpolation(raster_with_nan)
z = zComponentEstimation(interpolated_raster)
"""
