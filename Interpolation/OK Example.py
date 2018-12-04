# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 11:31:28 2018

@author: bimta
"""

from __future__ import division 
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as LA

raster = np.loadtxt("Data/raster.asc")


def OK(x,y,v,variogram,grid):
    cov_angulos = np.zeros((x.shape[0],x.shape[0]))
    cov_distancias = np.zeros((x.shape[0],x.shape[0]))
    K = np.zeros((x.shape[0]+1,x.shape[0]+1))
    for i in range(x.shape[0]-1):
        cov_angulos[i,i:]=np.arctan2((y[i:]-y[i]),(x[i:]-x[i]))
        cov_distancias[i,i:]=np.sqrt((x[i:]-x[i])**2+(y[i:]-y[i])**2)
    for i in range(x.shape[0]):
        for j in range(x.shape[0]):
            if cov_distancias[i,j]!=0:
                amp=np.sqrt((variogram[1]*np.cos(cov_angulos[i,j]))**2+(variogram[0]*np.sin(cov_angulos[i,j]))**2)
                K[i,j]=v[:].var()*(1-np.e**(-3*cov_distancias[i,j]/amp))
    K = K + K.T
    K[-1,:] = 1
    K[:,-1] = 1
    K[-1,-1] = 0

    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
             distancias = np.sqrt((i-x[:])**2+(j-y[:])**2)
             angulos = np.arctan2(i-y[:],j-x[:])
             amplitudes = np.sqrt((variogram[1]*np.cos(angulos[:]))**2+(variogram[0]*np.sin(angulos[:]))**2)
             M = np.ones((x.shape[0]+1,1))
             M[:-1,0] = v[:].var()*(1-np.e**(-3*distancias[:]/amplitudes[:]))
             W = LA.solve(K,M)
             grid[i,j] = np.sum(W[:-1,0]*(v[:]))
    return grid

np.random.seed(123433789) # GIVING A SEED NUMBER FOR THE EXPERIENCE TO BE REPRODUCIBLE
grid = np.zeros((10,10),dtype='float32') # float32 gives us a lot precision
x,y = np.random.randint(0,100,10),np.random.randint(0,100,10) # CREATE POINT SET.
v = np.random.randint(0,10,10) # THIS IS MY VARIABLE

test = raster[0:10, 0:10]
test = findReplace(test, -9999, np.nan)


i = np.arange(0, test.shape[0])
j = np.arange(0, test.shape[1])
jj, ii = np.meshgrid(i, j)

x = jj.reshape(test.shape[0]*test.shape[1],)
y = ii.reshape(test.shape[0]*test.shape[1],)
z = test.reshape(test.shape[0]*test.shape[1],)

mask = np.isnan(z)

z = z[~mask]
x = x[~mask]
y = y[~mask]




grid = OK(x,y,z,(50,30),grid)
plt.imshow(grid.T,origin='lower',interpolation='nearest',cmap='jet')
plt.scatter(x,y,c=z,s=120)
plt.xlim(0,grid.shape[0])
plt.ylim(0,grid.shape[1])
plt.grid()
plt.show()