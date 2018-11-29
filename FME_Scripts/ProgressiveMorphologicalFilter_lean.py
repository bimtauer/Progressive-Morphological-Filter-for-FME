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
from scipy import ndimage
from curvature import normalVectorEstimation



raster = np.loadtxt("Data/raster.asc")

# Parameters
parameters = {'c': 0.5,
              'initial_size': 2,
              'initial_cutoff': 0.4,
              'average_sigma': 7,
              'dh0': 0.1,
              'scaling_factor': 11,
              'hole_cutoff': -0.5}


class ProgressiveMorphologicalFilter():
    def __init__(self, input_raster, parameters):
        self.input_raster = input_raster
        self.c = parameters['c']
        self.initial_size = parameters['initial_size']
        self.initial_curoff = parameters['initial_cutoff']
        self.average_sigma = parameters['average_sigma']
        self.dh0 = parameters['dh0']
        self.scaling_factor = parameters['scaling_factor']
        self.hole_cutoff = parameters['hole_cutoff']
        
#Exchanges Raster Band Nodata (-9999) with np.nan
    def nanFormatter(raster):
        raster = np.where(raster==-9999, np.nan, raster)
        return raster
    #The Inverse
    def nanReFormatter(raster):
        raster = np.where(np.isnan(raster), -9999, raster)
        return raster
    
    #Replaces -9999 with 9999 for opening filter
    def nanReplacer(raster):
        raster = np.where(raster==-9999, 9999, raster)
        return raster
    
    # To get rid of holes
    def medianFilter(raster, cutoff):
        median = ndimage.median_filter(raster, size = (3,3))
        output = np.where(raster - median < cutoff, np.nan, raster)
        return output
    
    # The main filter algorithm
    def progressiveMorphologicalfilter(self, raster, slope_threshold, maxk=15):
        #Replace -9999 with 9999
        input_raster = nanReplacer(raster)
        lastSurface = np.copy(input_raster)
    
        mask = np.zeros(input_raster.shape)
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
            
            k += 1          #one step further
        
        output = np.where(mask, np.nan, input_raster)
        #output = medianFilter(output)
        output[output==9999]=np.nan
        return output
    
    #Slope estimation from surface unit normals z component transformed into radian angle:
    def zExtraction(prov_filtered):
        # Replace Nan with -9999 again temporary solution- figure interpolation at some point
        prov_filtered = np.where(np.isnan(prov_out), -9999, prov_filtered)
        # 2. Calculate normals for this surface
        z = normalVectorEstimation(prov_filtered)[:,:,2]
        return z
    
    def slopeEstimation(z, s):    
        # Assume Nan is even surface:
        z = np.where(np.isnan(z), 1, z)
        
        # average normals out with wide gaussian filter
        average_slope = ndimage.gaussian_filter(z, sigma=s)
        return average_slope
    
   #TODO: remove, integrate into first one
    def secondprogressiveMorphologicalFilter(raster, scaling, maxk=15):
        #Replace -9999 with 9999
        input_raster = nanReplacer(raster)
        lastSurface = np.copy(input_raster)
    
        mask = np.zeros(input_raster.shape)
        k = 1
        while k <= maxk:
            wk = 2*k+1       #windows size
            window = (wk,wk)  
            
            #Opening
            thisSurface = ndimage.morphology.grey_opening(lastSurface, size=window)
            
            #Scaling the maximum slope tolerance by average area slope
            
            dhmax = np.sqrt(k**2 + k**2) * c * scaling
            
            #Only return those below cutoff
            mask = np.where(lastSurface - thisSurface > dhmax, 1, mask)
            
            k += 1          #one step further
        
        output = np.where(mask, np.nan, input_raster)
        output[output==9999]=np.nan
        return output

###############################################################################
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
###############################################################################


def scalingMatrix(average_slope, dh0, scaling):
    normalized = np.copy(average_slope)
    normalized = normalized - normalized.min()
    normalized *= 1/normalized.max() 
    normalized = 1 - normalized
    normalized = normalized ** 2   #Decrease minor slopes
    normalized *= scaling - 1
    normalized += 1
    normalized *= dh0
    return normalized



initial_raster = nanFormatter(raster)

prov_out = progressiveMorphologicalFilter(raster, initial_size, initial_cutoff)

z = zExtraction(prov_out)

average_slope = slopeEstimation(z, average_sigma)

scaling = scalingMatrix(average_slope, dh0, scaling_factor)

almost = secondprogressiveMorphologicalFilter(raster, scaling)

final = medianFilter(almost, hole_cutoff)


# Plot all
plotter([initial_raster, prov_out, scaling, almost, final])


from writers import asciiOut
asciiOut(final, "MorphFiltered_ins{}_inc{}_avg{}_dh0{}_sca{}_hol{}".format(initial_size, 
         initial_cutoff, average_sigma, dh0, scaling_factor, hole_cutoff))


