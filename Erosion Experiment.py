# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 10:41:20 2018

@author: bimta
"""

import numpy as np
from scipy import ndimage

#Manually deleted first 6 text lines
raster = np.loadtxt("raster.asc")

#Slope tolerance delta H max - tolerated increase in H for distance between two points
dhmax = 0.3

#Grid size
a = 0.5

#Kernel Size   
kernel = (3,3) #only uneven numbers!

def maxSlopeErosion(input, kernel, dhmax):
    def setupStructureElement(kernel, dhmax):
        #Throw error if kernel set wrong
        if kernel[0]%2 == 0 or kernel[1]%2 == 0:
            raise ValueError("The kernel shape cannot contain even numbers!")
        
        structure_element = np.zeros(kernel) #same size
        for i in range(kernel[0]):
            for j in range(kernel[1]):
                centerx = (kernel[1]-1)/2
                centery = (kernel[0]-1)/2
                dx = a * (j - centerx)
                dy = a * (i - centery)
                value = -dhmax * np.sqrt(dx**2+dy**2) #following Vosselman(2000)
                structure_element[i,j] = value
                
        return structure_element

    s = setupStructureElement(kernel, dhmax)
    
    output = ndimage.morphology.grey_erosion(input, structure = s)
    return output

def maxSlopeFilter(input, kernel, dhmax):
    #Replace -9999 with NaN
    input[input==-9999]=np.nan
    #Call Erosion
    eroded = maxSlopeErosion(input, kernel, dhmax)
    #Only return those below cutoff
    output = np.where(input <= eroded, input, np.nan)
    return output

filtered = maxSlopeFilter(raster, kernel, dhmax)
np.nanmin(filtered)
        