# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 10:04:12 2018

@author: bimta
"""

from FME_Scripts.PMF import ProgressiveMorphologicalFilter
import numpy as np

raster = np.loadtxt(r"./Merged/DEM_0.5_min_filtered_merged.asc", skiprows = 6)
raster[raster == -9999] = np.nan

parameters = {'c' : 0.5,                     # The cell size of the input raster
              'kernel_radius' : 15,          # The kernel size of the first filter iteration
              'initial_cutoff' : 0.4,        # The slope threshold of the first filter iteration
              'average_sigma' : 7,           # The gaussian function used to average out local slope
              'dh0' : 0.1,                   # The slope threshold for flat areas
              'hole_cutoff' : -0.2}          # Threshold for individual holes beneath median in 3x3 kernel

MyFilter = ProgressiveMorphologicalFilter(raster, parameters)

output = MyFilter.filter()