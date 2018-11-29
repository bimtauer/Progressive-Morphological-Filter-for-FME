# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 16:46:15 2018

@author: timbau
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from curvature import normalVectorEstimation



raster = np.loadtxt("Data/raster.asc")
c = 0.5

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

#Still unsure about this one
def linearInterpolation(raster):
    mask = np.isnan(raster)
    raster[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), raster[~mask])
    return raster

def progressiveMorphologicalFilter(raster, maxk=15, dh0=0.3):
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
        dhmax = np.sqrt(k**2 + k**2) * c * dh0 
        
        #Only return those below cutoff
        mask = np.where(lastSurface - thisSurface > dhmax, 1, mask)
        
        k += 1          #one step further
    
    output = np.where(mask, np.nan, input_raster)
    #output = medianFilter(output)
    output[output==9999]=np.nan
    return output


""" From Pingle, Clarke, McBride (2013), "An Improved Simple Morphological Filter
for the Terrain Classification of Airborne LIDAR Data":
    
"While many authors use a
single value for the elevation threshold, we suggest that a second
parameter be used to increase the threshold on steep slopes, transforming
the threshold to a slope-dependent value. The total permissible
distance is then equal to a fixed elevation threshold plus
the scaling value multiplied by the slope of the DEM at each LIDAR
point."
"""

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

"""
fig = plt.figure(figsize=(30, 30))
fig.add_subplot(1,1,1)
plt.imshow(z)
fig.add_subplot(2,2,2)
plt.imshow(average_slope)
plt.show()
"""
# 3. call prg morph with slope average as scaling
"""
maxk = 15
k = 1
dh0=0.3
while k <= maxk:
    dhmax = np.sqrt(k**2 + k**2) * c * dh0 
    k+= 1
    print(dhmax)
"""
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

#from writers import asciiOut

#asciiOut(out, "newMorphFiltered")

###############################################################################
# Interpolation Experiments
"""
#Inpainting Nan in original raster
def inpainter(input_raster):
    mask = np.where(np.isnan(input_raster), 1, 0)
    start = time.process_time()
    inpainted = inpaint.inpaint_biharmonic(scaled, mask)
    end = time.process_time()        
    print("Done after {0:.2f} seconds.".format(end - start))
    return inpainted

raster_inpainted = inpainter(raster)


# Doesn't remove all Nan
def CVinpainter(input_raster):
    mask = np.where(np.isnan(input_raster), 1, 0)
    start = time.process_time()
    dst = cv2.inpaint(input_raster.astype("float32"),mask.astype("uint8"), 3, cv2.INPAINT_NS)
    end = time.process_time()        
    print("Done after {0:.2f} seconds.".format(end - start))
    return dst

rad_CVinpainted = CVinpainter(rad)

#Interpolation
def Interpolation(a):
    x, y = np.indices(a.shape)
    interp = np.array(a)
    interp[np.isnan(interp)] = griddata((x[~np.isnan(a)], y[~np.isnan(a)]), # points we know
                                          a[~np.isnan(a)],                    # values we know
                                          (x[np.isnan(a)], y[np.isnan(a)]))   # points to interpolate
    return interp

test = Interpolation(rad)


#Now I want a function that calculates average rad over a rolling window
average_slope = ndimage.uniform_filter(rad_CVinpainted, size = (31,31))

weights = np.array([[0.5,0.7,0.5],[0.7,1,0.7],[0.5,0.7,0.5]])
average_slope = ndimage.convolve(rad_CVinpainted, weights)

fig = plt.figure(figsize=(30, 30))
fig.add_subplot(1,1,1)
plt.imshow(rad)
fig.add_subplot(2,2,2)
plt.imshow(average_slope)
plt.show()

interp = linearInterpolation(z)

x = np.arange(z.shape[0])
y = np.arange(z.shape[1])
f = interpolate.NearestNDInterpolator(x, y, z)
znew = f(x, y)


average_slope = ndimage.uniform_filter(z[2:10,3:10], size = (3,3))

"""
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

# Parameters
#parameters = 
initial_size = 2
initial_cutoff = 0.4
average_sigma = 7
dh0 = 0.1
scaling_factor = 11
hole_cutoff = -0.5

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


