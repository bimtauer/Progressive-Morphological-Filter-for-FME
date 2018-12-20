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

The raster is also initally filtered for single pixel holes which can
be assumed to be measurement errors.


Sources:

Pingle, Clarke, McBride (2013). "An Improved Simple Morphological Filter for the
    Terrain Classification of Airborne LIDAR Data". ISPRS Journal of Photogrammetry
    and Remote Sensing, 77, 21-30.

Vosselman, G. (2000). Slope based filtering of laser altimetry data. International
    Archives of Photogrammetry and Remote Sensing, 33(B3/2; PART 3), 935-942.
"""


import numpy as np
from scipy import ndimage
from .MyRasterTools import slopeEstimation, nnInterpolation, tinInterpolation

################################################################################
"""
# Parameters (Example)
parameters = {'c' : 0.5,                     # The cell size of the input raster
              'kernel_radius' : 15,          # The kernel size of the first filter iteration
              'initial_cutoff' : 0.4,        # The slope threshold of the first filter iteration
              'average_sigma' : 7,           # The gaussian function used to average out local slope
              'dh0' : 0.1,                   # The slope threshold for flat areas
              'hole_cutoff' : -0.2}          # Threshold for individual holes beneath median in 3x3 kernel
"""

# The Filter Class
class ProgressiveMorphologicalFilter():
    def __init__(self, input_raster, parameters):
        for key, value in parameters.items():
            setattr(self, key, value)
        self.input_raster = input_raster

    # To get rid of holes
    def medianFilter(self, input_raster):
        #The mask we use to indicate holes
        mask = np.ones(input_raster.shape, dtype=bool)
        median = ndimage.median_filter(input_raster, size = (3,3))
        mask = np.where(input_raster - median < self.hole_cutoff, 0, mask)
        return mask

    # The main filter algorithm
    def progressiveMorphologicalfilter(self, last_surface, slope_threshold):

        #The mask we use to indicate non-ground points
        mask = np.ones(last_surface.shape, dtype=bool)
        final = int(self.kernel_radius/self.c)
        k = 1
        while k <= final:
            #Generate window
            w = k*2 + 1
            window = (w,w)
            #Opening
            this_surface = ndimage.morphology.grey_opening(last_surface, size = window)
            #Increasing dhmax by diagonal window extent in cells
            dhmax = np.sqrt(k**2 + k**2) * self.c * slope_threshold
            #Only mask cells below cutoff
            mask = np.where(last_surface - this_surface > dhmax, 0, mask)
            last_surface = this_surface
            k += 1          #one step further
        return mask

    def scalingMatrix(self, input_raster):
        z = slopeEstimation(input_raster, self.c)
        # assume nans are even surfaces
        z = np.where(np.isnan(z), 0, z)
        # average slope out with wide gaussian filter
        average_slope = ndimage.gaussian_filter(z, self.average_sigma)
        average_slope += self.dh0
        return average_slope

    def filter(self):
        print("Filtering holes...")
        hole_mask = self.medianFilter(self.input_raster)
        self.hole_filtered = np.where(hole_mask, self.input_raster, np.nan)
        print("Beginning initial filtering...")
        initial_mask = self.progressiveMorphologicalfilter(self.hole_filtered, self.initial_cutoff)
        self.initial_filtered = np.where(initial_mask, self.input_raster, np.nan)
        print("Computing average slope...")
        self.scaling_matrix = self.scalingMatrix(self.initial_filtered)
        print("Final filtering with adaptive threshold...")
        second_mask = self.progressiveMorphologicalfilter(self.hole_filtered, self.scaling_matrix)
        final_mask = hole_mask & initial_mask & second_mask
        self.final_filtered = np.where(final_mask, self.input_raster, np.nan)
        return self.final_filtered
