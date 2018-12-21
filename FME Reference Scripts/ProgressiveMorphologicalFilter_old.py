# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 14:51:33 2018

@author: Tim Bauer

Description:

This is a morphological filter to extract ground level points from unfiltered
elevation rasters. It can be used to generate a DTM from an unclassified LIDAR
pointcloud, assuming that this pointcloud has already been turned into an
evenly spaced minimum elevation raster (picking minimum z for each cell).

The algorithm is build on the methodology described in Pingle, Clarke, McBride
(2013). The main filtering proceedure is comparison of the input raster to a
morphological opening of that input raster with a stepwise increasing window,
discarding points as non-ground if the difference between original point and
opening exceeds a slope threshold. The algorithm is thus similar to the one
described by Vosselmann (2000) which utilizes morphological erosion. The
advantage of using opening instead of erosion is the improved preservation of
distinct ground features such as tunnels or bridges.

The algorithm proceeds in two steps: First an initial filtering with a
conservative slope threshold (parameter: initial_cutoff) removes the most
prominent artefacts. For the resulting provisional surface surface normal
vectors are estimated and averaged with a gaussian kernel
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
import time #TODO: Remove
from scipy import ndimage
import matplotlib.pyplot as plt #TODO: Remove

raster = np.loadtxt("../Data/raster.asc") #TODO: Remove

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

    #A find and replace function for arrays that can handle np.nan
    def findReplace(self, where, find, replace):
        if np.isnan(find):
            output = np.where(np.isnan(where), replace, where)
        else:
            output = np.where(where == find, replace, where)
        return output

    # To get rid of holes
    def medianFilter(self, input_raster):
        filter_raster = np.where(np.isnan(raster), 9999, raster)
        median = ndimage.median_filter(filter_raster, size = (3,3))
        output = np.where(filter_raster - median < self.hole_cutoff, np.nan, input_raster)
        return output

    # The main filter algorithm
    def progressiveMorphologicalfilter(self, raster, slope_threshold, maxk):
        #Replace -9999 with 9999
        input_raster = self.findReplace(raster, -9999, 9999)
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

    def normalVectorEstimation(self, raster):
        #Creating 3d matrix to store resulting vectors
        out = np.zeros((raster.shape[0], raster.shape[1]))

        #TODO: Include interpolation here or make better method

        #If input contains nan, replace with -9999 so that nans don't mess up median padding
        raster = np.where(np.isnan(raster), -9999, raster)
        #Padding Raster
        padded = np.pad(raster, 1, "median")
        #Reintroduce nan to not compute normals on the border to Nan
        padded = np.where(padded==-9999, np.nan, padded)

        rows = padded.shape[0]
        cols = padded.shape[1]

        print("Beginning normal vector estimation...")
        start = time.process_time()
        for i in range(1, rows-1):                  #Replace with vectorized solution in future to make faster
            for j in range(1, cols-1):
                #Surrounding points
                top = padded[i-1,j]
                bot = padded[i+1,j]
                left = padded[i,j-1]
                right = padded[i,j+1]

                # calculating surface normal based on 4 surrounding points:
                normal = np.array([(left - right), (top - bot), (2*self.c)])
                mag = np.sqrt(normal.dot(normal))
                unit_normal = normal / mag

                out[i-1, j-1] = unit_normal[2] #only z component of normal

        end = time.process_time()
        print("Done after {0:.2f} seconds.".format(end - start))
        return out

    def scalingMatrix(self, raster):
        # Assume Nan is even surface:
        z = self.normalVectorEstimation(raster)
        z = np.where(np.isnan(z), 1, z)
        # average normals out with wide gaussian filter
        average_slope = ndimage.gaussian_filter(z, self.average_sigma)
        normalized = np.copy(average_slope)
        normalized = normalized - normalized.min()
        normalized *= 1/normalized.max()
        normalized = 1 - normalized
        normalized = normalized ** 2   #Decrease minor slopes
        normalized *= self.scaling_factor - 1
        normalized += 1
        normalized *= self.dh0
        return normalized

    def filter(self):
        self.initial_filtered = self.progressiveMorphologicalfilter(self.input_raster, self.initial_cutoff, self.initial_size)
        self.scaling_matrix = self.scalingMatrix(self.initial_filtered)
        self.final_filtered = self.progressiveMorphologicalfilter(self.input_raster, self.scaling_matrix, self.final_size)
        hole_filtered = self.medianFilter(self.final_filtered)
        return hole_filtered


###############################################################################
def plotter(raster_list):
    nr_plots = len(raster_list)
    fig, ax = plt.subplots(1, nr_plots ,figsize=(30, 9), sharex = True, sharey = True)

    ind = 0
    for raster in raster_list:
        Name = [name for name in globals() if globals()[name] is raster][0]
        img = ax[ind].imshow(raster, cmap="viridis")
        ax[ind].set_title('{}'.format(Name))
        fig.colorbar(img, ax = ax[ind], orientation = 'horizontal')
        ax[ind].set_axis_off()
        ind += 1

    plt.tight_layout()
    plt.show()
    return
###############################################################################
"""
MyFilter = ProgressiveMorphologicalFilter(raster, parameters)
output = MyFilter.filter()
initial_filtered = MyFilter.initial_filtered
scaling_matrix = MyFilter.scaling_matrix
final_filtered = MyFilter.final_filtered

# Plot all
plotter([initial_filtered, scaling_matrix, final_filtered, output])


test = np.where(output < 0, output, 1)
plotter([output, test])



from writers import asciiOut
asciiOut(final, "MorphFiltered_ins{}_inc{}_avg{}_dh0{}_sca{}_hol{}".format(initial_size,
         initial_cutoff, average_sigma, dh0, scaling_factor, hole_cutoff))

"""
