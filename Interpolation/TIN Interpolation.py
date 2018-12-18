# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 09:25:44 2018

@author: timbau
"""

import numpy as np
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator

"""
You can give the Delaunay triangulation to scipy.interpolate.LinearNDInterpolator 
together with the set of Z-values, and it should do the job for you.
"""
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


interpolated = TINinterpolation(raster)

plotter([raster, interpolated])

import plotly.plotly as py
import plotly.graph_objs as go

# Read data from a csv

data = [
    go.Surface(
        z=test
    )
]
layout = go.Layout(
    scene = dict(
                    xaxis = dict(
                        nticks=7, range = [0,700],),
                    yaxis = dict(
                        nticks=7, range = [0,700],),
                    zaxis = dict(
                        nticks=3, range = [-5,50],),)
)
fig = go.Figure(data=data, layout=layout)
py.iplot(fig, filename='Test Area')