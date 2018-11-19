"""
Works only if input raster is REAL64
"""

import fmeobjects
import numpy as np

def processFeature(feature):
    raster = feature.getGeometry()
    rasterProperties = raster.getProperties()
    numRows, numCols = rasterProperties.getNumRows(), rasterProperties.getNumCols()

    band = raster.getBand(0)
    bandProperties = band.getProperties()
    tile = fmeobjects.FMEReal64Tile(numRows, numCols)

    values = []
    for data in band.getTile(0, 0, tile).getData():
        values += data
    tiles_matrix = np.array([values])
    tiles_matrix = tiles_matrix.reshape(numRows, numCols) #From this point on i have my raster
    print(tiles_matrix)
    return
