# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 14:51:33 2018

@author: Tim Bauer (tim.bauer@gmx.us)

Description:

This is a morphological filter to extract ground level points from unfiltered
elevation rasters. It can be used to generate a DTM from an unclassified LIDAR
pointcloud, assuming that this pointcloud has already been turned into an evenly
spaced minimum elevation raster (picking minimum z for each cell) which is used
as input.

The algorithm is build on the methodology described in Pingle, Clarke, McBride
(2013). The main filtering proceedure is comparison of the input raster to a
morphological opening of that input raster with a stepwise increasing window.
It discards non-ground points if the difference between original point and
opening exceeds a slope threshold. The algorithm is thus similar to the one
described by Vosselmann (2000) which utilizes morphological erosion. The
advantage of using opening instead of erosion is the improved preservation of
distinct ground features such as tunnels or bridges.

The algorithm proceeds in two steps: First an initial filtering with a
conservative slope threshold (parameter: initial_cutoff) removes the most
prominent artefacts. The resulting provisional surface is interpolated with a
nearest neighbor algorithm. On this interpolated surface the z component of the
normal vectors is estimated for each cell and then averaged out
(parameter: average_sigma) in order to receive an average slope for each area.
This average slope is transformed into a scaling matrix which serves as slope
cutoff in the second iteration. This way it is ensured that the subsequent
round of filtering uses a maximal slope cutoff in evenly sloped areas, while
performing a more moderate filtering along slopes. The parameter scaling_factor
times 0.1 is then equal to the maximum slope tolerated in steep areas.

The resulting raster is finally filtered for single pixel holes which can be
assumed to be measurement errors.


Sources:

Pingle, Clarke, McBride (2013). "An Improved Simple Morphological Filter for the
    Terrain Classification of Airborne LIDAR Data". ISPRS Journal of Photogrammetry
    and Remote Sensing, 77, 21-30.

Vosselman, G. (2000). Slope based filtering of laser altimetry data. International
    Archives of Photogrammetry and Remote Sensing, 33(B3/2; PART 3), 935-942.
"""


import numpy as np
from scipy import ndimage
from .MyRasterTools import zComponentEstimation, nnInterpolation

################################################################################
# The Filter Class

# Parameters
parameters = {'c' : 0.5,                     # The cell size of the input raster
              'initial_size' : 2,            # The kernel size of the first filter iteration
              'initial_cutoff' : 0.4,        # The slope threshold of the first filtr iteration
              'average_sigma' : 7,           # The gaussian function used to average out local slope
              'dh0' : 0.1,                   # The slope threshold for flat areas
              'dhmax' : 1.0,                 # The slope threshold at maximum slope
              'final_size' : 15,             # The kernel size of the second filter iteration
              'hole_cutoff' : -0.2}          # Threshold for individual holes beneath median in 3x3 kernel
              #'nan_value' : #TODO: Get from FME with FMEBand.getNodataValue()

class ProgressiveMorphologicalFilter():
    def __init__(self, input_raster, parameters):
        for key, value in parameters.items():
            setattr(self, key, value)

        self.input_raster = input_raster
        self.scaling_factor = self.dhmax/self.dh0

    # To get rid of holes
    def medianFilter(self, input_raster):
        filter_raster = np.where(np.isnan(input_raster), 9999, input_raster)
        median = ndimage.median_filter(filter_raster, size = (3,3))
        output = np.where(filter_raster - median < self.hole_cutoff, np.nan, input_raster)
        return output

    # The main filter algorithm
    def progressiveMorphologicalfilter(self, raster, slope_threshold, maxk):
        #Replace -9999 with 9999
        input_raster = np.where(raster == -9999, 9999, raster)
        lastSurface = np.copy(input_raster)

        mask = np.zeros(input_raster.shape) #The mask we use to indicate non-ground points
        k = 1
        while k <= maxk:
            wk = 2*k+1       #windows size
            window = (wk,wk)

            #Opening
            thisSurface = ndimage.morphology.grey_opening(lastSurface, size=window)

            #Increasing maxdh by window size - always the maximum slope to corner point
            dhmax = np.sqrt(k**2 + k**2) * self.c * slope_threshold

            #Only return those below cutoff
            mask = np.where(lastSurface - thisSurface > dhmax, 1, mask)

            #Continue with last output
            #last_surface = this_surface                #Deactivated because results show better performance without

            k += 1          #one step further

        output = np.where(mask, np.nan, input_raster)
        #Reintroduce nan
        output[output==9999]=np.nan
        return output

    def scalingMatrix(self, raster):
        interpolated_raster = nnInterpolation(raster)
        z = zComponentEstimation(interpolated_raster)
        # average normals out with wide gaussian filter
        average_slope = ndimage.gaussian_filter(z, self.average_sigma)
        normalized = np.copy(average_slope)
        #Scaling average slope to 0 - 1
        normalized = normalized - normalized.min()
        normalized *= 1/normalized.max()
        normalized = 1 - normalized
        #Decrease minor slopes
        normalized = normalized ** 2
        #Introduce scaling factor
        normalized *= self.scaling_factor - 1
        normalized += 1
        normalized *= self.dh0
        return normalized

    def filter(self):
        self.initial_filtered = self.progressiveMorphologicalfilter(self.input_raster, self.initial_cutoff, self.initial_size)
        self.scaling_matrix = self.scalingMatrix(self.initial_filtered)
        self.final_filtered = self.progressiveMorphologicalfilter(self.input_raster, self.scaling_matrix, self.final_size)
        self.hole_filtered = self.medianFilter(self.final_filtered)
        return self.hole_filtered
