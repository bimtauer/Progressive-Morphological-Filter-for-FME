# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 10:41:20 2018

@author: bimta
"""

import numpy as np
from scipy import ndimage
from writers import asciiOut 

#Manually deleted first 6 text lines
raster = np.loadtxt("raster.asc")
raster2 = np.loadtxt("MinMaxComparisonResult.asc")

###############################################################################
# Specify main parameters:

#Slope tolerance delta H max - tolerated increase in H for distance between two points
dhmax = 0.3
#Grid size
a = 0.5
#Kernel Size   
kernel = (3,3) #only uneven numbers!

###############################################################################

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

def maxSlopeFilter(input_raster, kernel, dhmax):
    #Replace -9999 with 9999 to avoid padding issue
    input_raster[input_raster==-9999]=9999

    #Call Erosion
    eroded = maxSlopeErosion(input_raster, kernel, dhmax)
    

    #Only return those below cutoff
    output = np.where(input_raster <= eroded, input_raster, np.nan)
    #Reintroduce nan
    output[output==9999]=np.nan
    return output

###############################################################################

# Output Raster
filtered = maxSlopeFilter(raster, kernel, dhmax)

# Out to Ascii  
asciiOut(filtered, "Slope_Erosion_filtered_{}kernel_{}dhmax".format(str(kernel[0])+str(kernel[1]), dhmax))


###############################################################################
#Examples:
"""
s = np.array([[2,1,1],
              [2,1,2],
              [2,1,2]])

t = np.zeros((7,7), dtype=np.int)
t[1:6, 1:6] = 3
t[4,4] = 2; a[2,3] = 1

test_out = ndimage.morphology.grey_erosion(t, structure = s)

###############################################################################
a = np.zeros((7,7), dtype=np.int)
a[1:6, 1:6] = 3
a[4,4] = 2; a[2,3] = 1
a


test_out = ndimage.morphology.grey_erosion(a, size = (3,3))

t = raster[0:5,0:5]
t = np.pad(t, 1, "constant")
test_out2 = ndimage.morphology.grey_erosion(t, size = (3,3))
"""