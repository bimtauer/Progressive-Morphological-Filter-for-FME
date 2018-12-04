# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 14:46:57 2018

@author: bimta
"""

from pykrige.uk import UniversalKriging
import numpy as np
import matplotlib.pyplot as plt


raster = np.loadtxt("Data/raster.asc")


def findReplace(where, find, replace):
    if np.isnan(find):
        output = np.where(np.isnan(where), replace, where)
    else:
        output = np.where(where == find, replace, where)
    return output




def rasterToList(raster):
    raster = findReplace(raster, -9999, np.nan)
    leni = raster.shape[0]
    lenj = raster.shape[1]
    
    i = np.arange(0, leni)
    j = np.arange(0, lenj)
    jj, ii = np.meshgrid(i, j)
    
    origx = jj.reshape(leni*lenj,)
    origy = ii.reshape(leni*lenj,)
    origz = raster.reshape(leni*lenj,)
    
    mask = np.isnan(z)
    
    z = origz[~mask]
    x = origx[~mask]
    y = origy[~mask]
    out = np.array([x,y,z]).T.astype('float16')
    i = i[mask].astype('float32')
    j = j[mask].astype('float32')
    return out, j, i

data, gridx, gridy = rasterToList(raster[0:100,0:100])

# CSV out
#np.savetxt("points.csv", data, delimiter=",")

"""
data = np.array([[1.0, 2.0, -10.00],
                 [1.5, 1.5, 30.00],
                 [2.0, 2.0, 15.00]])

gridx = np.arange(0.0, 5.5, 0.5)
gridy = np.arange(0.0, 5.5, 0.5)
"""


# Create the ordinary kriging object. Required inputs are the X-coordinates of
# the data points, the Y-coordinates of the data points, and the Z-values of the
# data points. Variogram is handled as in the ordinary kriging case.
# drift_terms is a list of the drift terms to include; currently supported terms
# are 'regional_linear', 'point_log', and 'external_Z'. Refer to
# UniversalKriging.__doc__ for more information.
UK = UniversalKriging(data[:, 0], data[:, 1], data[:, 2], variogram_model='linear', drift_terms=['regional_linear'])
# Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
# grid of points, on a masked rectangular grid of points, or with arbitrary points
z, ss = UK.execute('grid', gridx, gridy)


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

test = raster[0:100,0:100]
test = findReplace(test, -9999, np.nan)
plotter([test,z])

np.savetxt("kriging_interpolated.csv", z, delimiter=",")
np.savetxt("krieging_ss.csv", ss, delimiter=",")





from pykrige.ok import OrdinaryKriging

# Create the ordinary kriging object. Required inputs are the X-coordinates of
# the data points, the Y-coordinates of the data points, and the Z-values of the
# data points. If no variogram model is specified, defaults to a linear variogram
# model. If no variogram model parameters are specified, then the code automatically
# calculates the parameters by fitting the variogram model to the binned
# experimental semivariogram. The verbose kwarg controls code talk-back, and
# the enable_plotting kwarg controls the display of the semivariogram.
OK = OrdinaryKriging(data[:, 0], data[:, 1], data[:, 2], variogram_model='linear', verbose=False , enable_plotting=False)
# Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
# grid of points, on a masked rectangular grid of points, or with arbitrary points.
# (See OrdinaryKriging.__doc__ for more information.)
z, ss = OK.execute('grid', gridx, gridy)
# Writes the kriged grid to an ASCII grid file.
kt.write_asc_grid(gridx, gridy, z, filename="output.asc")
