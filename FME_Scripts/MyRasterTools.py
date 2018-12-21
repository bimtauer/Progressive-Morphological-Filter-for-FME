# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 09:58:10 2018

@author: Tim Bauer (tim.bauer@gmx.us)
"""

import numpy as np
from scipy.interpolate import NearestNDInterpolator
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator

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

#Estimates rise/run for a 8 cell kernel following method from:
# http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-slope-works.htm
def slopeEstimation(raster, c):
    raster = np.pad(raster, 1, "edge")
    x1 = raster[ :-2,:]
    x2 = raster[1:-1,:]
    x3 = raster[2:  ,:]
    xs = (x1 + 2*x2 + x3)
    dzx = xs[:,2:] - xs[:,:-2]
    dzdx = dzx/8*c
    y1 = raster[:, :-2]
    y2 = raster[:,1:-1]
    y3 = raster[:,2:  ]
    ys = (y1 + 2*y2 + y3)
    dzy = ys[2:,:] - ys[:-2,:]
    dzdy = dzy/8*c
    rise_run = np.sqrt((dzdx**2) + (dzdy**2))
    return rise_run

#Takes a raster containing nan as input and interpolates all holes with NN
def nnInterpolation(raster):
    def rasterToList(raster):
        leni = raster.shape[0]
        lenj = raster.shape[1]

        #Retrieve indexes as meshgrid of coordinates
        i = np.arange(0, leni)
        j = np.arange(0, lenj)
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
def tinInterpolation(raster):
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
    interpolated = Interpolator(tointerpolate).reshape(leni,lenj)
    return interpolated
