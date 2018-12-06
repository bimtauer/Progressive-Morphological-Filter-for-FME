# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 09:25:44 2018

@author: timbau
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay


raster = np.loadtxt("../Outputs/MorphFilteredtestIn.asc")

#Own data
test = np.where(raster == -9999, np.nan, raster)[200:210, 200:210]

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

points = np.column_stack((coordinates, values))
tri = Delaunay(points)
indices = tri.simplices
vertices = points[indices]
"""
You can give the Delaunay triangulation to scipy.interpolate.LinearNDInterpolator 
together with the set of Z-values, and it should do the job for you.
"""

"""
Also zou could estimate normals from the vertices/simplices
"""



from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri



fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1, projection='3d')
plot_trisurf(points[:,0], points[:,1], points[:,2], cmap=plt.cm.Spectral)
plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot_trisurf(points[:,0], points[:,1], points[:,2], cmap=plt.cm.Spectral)

plt.show()


import plotly.plotly as py
import plotly.graph_objs as go

import pandas as pd

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