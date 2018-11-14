# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 20:48:03 2018

@author: bimta
"""

# PythonCreator Script Example: Create a Feature containing a Raster
# The raster will have a single band with UINT8 interpretation.
import fmeobjects
 
# Define a concrete class derived from FMEBandTilePopulator class.
# An instance of this class will be used to create a tile
# and populate it to a band (FMEBand instance).
class MyTilePopulator(fmeobjects.FMEBandTilePopulator):
    def __init__(self, dataArray):
        self.dataArray = dataArray
        
    # Implement 'clone' method.
    # It will be called multiple times while creating a new band.
    def clone(self):
        return MyTilePopulator(self.dataArray)
        
    # Implement 'getTile' method.
    # You can create a new tile containing desired contents here.
    # It's not essential to use the parameters: startRow, startCol, tile.
    def getTile(self, startRow, startCol, tile):
        numRows, numCols = len(self.dataArray), len(self.dataArray[0])
        newTile = fmeobjects.FMEReal64Tile(numRows, numCols)
        newTile.setData(self.dataArray)      
        return newTile
        
    # The following two methods won't be called while creating a new band.
    # It seems not to be essential to implement these methods in this case,
    # although the API doc says "This method must be implemented
    # in the FMEBandTilePopulator subclass".
    def setDeleteSourceOnDestroy(self, deleteFlag):
        pass
    def setOutputSize(self, rows, cols):
        return (rows, cols)
        
class FeatureCreator(object):
    def __init__(self):
        pass
    
    def input(self, feature):
        pass
    
    def close(self):
        # Contents of a tile for a band to be created.
        # A list of row data, each element is a list of column values.
        dataArray = [
            [  0, 128,   0, 128,   0, 128,   0],
            [128,   0, 128,   0, 128,   0, 128],
            [  0, 128,   0, 128,   0, 128,   0],
            [128,   0, 128,   0, 128,   0, 128],
            [  0, 128,   0, 128,   0, 128,   0],
        ]
        
        # Properties of a raster to be created.
        numRows, numCols = len(dataArray), len(dataArray[0]) # resolution
        xSpacing, ySpacing = 10.0, 10.0 # cell spacing in ground units
        xCellOrigin, yCellOrigin = 0.5, 0.5 # cell origin coordinates
        xOrigin, yOrigin = 0.0, numRows * ySpacing # left-top coordinates
        xRotation, yRotation = 0.0, 0.0 # rotation angle in degrees
        
        # Create a new raster.
        rasterProperties = fmeobjects.FMERasterProperties(numRows, numCols,
            xSpacing, ySpacing, xCellOrigin, yCellOrigin, xOrigin, yOrigin,
            xRotation, yRotation)
        raster = fmeobjects.FMERaster(rasterProperties)
        
        # Create a new band and append it to the raster.
        # It's optional to specify Nodata value when creating a band.
        bandTilePopulator = MyTilePopulator(dataArray)
        bandName = 'My Real64 Band' # can be set to empty.
        
        bandProperties = fmeobjects.FMEBandProperties(bandName,
            fmeobjects.FME_INTERPRETATION_REAL64,
            fmeobjects.FME_TILE_TYPE_FIXED,
            numRows, numCols)
        
        
        #nodataValue = fmeobjects.FMEReal64Tile(1, 1)
        #nodataValue.setData([[0]])
        
        band = fmeobjects.FMEBand(bandTilePopulator, rasterProperties, bandProperties) # nodataValue
        raster.appendBand(band)
        
        # Create and output a feature containing the raster created above.
        feature = fmeobjects.FMEFeature()
        feature.setGeometry(raster)
        print("Message: All fine until here")
        self.pyoutput(feature)