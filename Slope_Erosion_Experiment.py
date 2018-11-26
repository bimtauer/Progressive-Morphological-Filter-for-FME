
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

raster = np.loadtxt("Data/raster.asc")
 
dhmax = 0.3
        
#Grid size
a = 0.5

#Kernel Size   
kernel = (9,9)

#Exchanges Raster Band Nodata (-9999) with np.nan
def nanFormatter(raster):
    raster = np.where(raster==-9999, np.nan, raster)
    return raster

#Replaces -9999 with 9999 for opening filter
def nanReplacer(raster):
    raster = np.where(raster==-9999, 9999, raster)
    return raster


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

def maxSlopeFilter(input_raster, kernel = (9,9), dhmax = 0.3):
    #Replace -9999 with NaN
    input_raster = nanReplacer(input_raster)

    #Call Erosion
    eroded = maxSlopeErosion(input_raster, kernel, dhmax)

    #Only return those below cutoff
    output = np.where(input_raster <= eroded, input_raster, np.nan)
    output[output==9999]=np.nan
    return output
    

###############################################################################        
    
Trees = raster[500:600,100:200]
Tunnels = raster[310:410, 370:470]
Holes = raster[200:300, 50:150]

def comparisonPlot(raster1, raster2, raster3):
    
    fig, ax = plt.subplots(2, 3, figsize=(50, 15))
    
    rasters = [raster1, raster2, raster3]
    index = 0
    for raster in rasters:
        Name = [name for name in globals() if globals()[name] is raster][0]
        out = maxSlopeFilter(raster, kernel, dhmax)
        img1 = ax[0,index].imshow(nanFormatter(raster), cmap="viridis")
        ax[0,index].set_title('{} Original'.format(Name))
        ax[0,index].set_axis_off()
        img2 = ax[1,index].imshow(out, cmap="viridis")
        ax[1,index].set_title('{} Filtered'.format(Name))
        fig.colorbar(img2, ax = ax[1,index])
        ax[1,index].set_axis_off()

        index += 1
        
            
    # Plotting   
    plt.tight_layout()
    plt.show()
    return

comparisonPlot(Trees, Tunnels, Holes)      