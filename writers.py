# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 13:35:38 2018

@author: bimta
"""
import numpy as np

#Write out ascii grids in proper format for our Ribban test area
def asciiOut(input_raster, name):
    ncols = input_raster.shape[1]
    nrows = input_raster.shape[0]
    xllcorner = 116549.99999999983
    yllcorner = 6164149.999999888
    cellsize = 0.5000000000021828
    nodata_value = -9999
    
    #Replace np.nan with raster nodata value
    input_raster[np.isnan(input_raster)]= nodata_value
    
    with open("./Outputs/{}.asc".format(name),"w") as TheFile:
        TheFile.write("ncols {}\n".format(ncols))
        TheFile.write("nrows {}\n".format(nrows))
        TheFile.write("xllcorner {}\n".format(xllcorner))
        TheFile.write("yllcorner {}\n".format(yllcorner))
        TheFile.write("cellsize {}\n".format(cellsize))
        TheFile.write("nodata_value {}\n".format(nodata_value))
          
        for i in range(0,nrows):
            for j in range(0,ncols):
                TheFile.write(str(input_raster[i,j]))
                TheFile.write(" ")
            TheFile.write("\n")
    return
