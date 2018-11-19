"""
This is a raster reader/writer template to be used in FME
"""

import fmeobjects
import numpy as np
# Raster reading and writing functions based on code from Takashi: https://knowledge.safe.com/questions/38000/python-fme-objects-api-for-raster-manipulation.html

class RasterReader():
    def __init__(self):
        pass

    def tileAndInterpretation(self, interpretation, numRows, numCols):
        if interpretation == fmeobjects.FME_INTERPRETATION_REAL64:
            return (fmeobjects.FMEReal64Tile(numRows, numCols), 'REAL64')
        elif interpretation == fmeobjects.FME_INTERPRETATION_REAL32:
            return (fmeobjects.FMEReal32Tile(numRows, numCols), 'REAL32')
        ######################################################
        # Add other interpretations here if necessary.
        ######################################################
        else:
            return (None, 'Unsupported')

    # Returns a dictionary containing the statistics of specified band.
    def retrieveTiles(self, band, numRows, numCols):
        # Get the band properties.
        bandProperties = band.getProperties()
        interpretation = bandProperties.getInterpretation()

        # Create a tile that will be used to get cell values of the band.
        # The interpretation, number of rows, and number of columns
        # have to be consistent with the band properties.
        tile, interpret = self.tileAndInterpretation(interpretation,
            numRows, numCols)

        if tile != None:
            # Get all the cell values except Nodata as a list.
            #------> Check if that is actually a problem!!! Doesn't seem to be
            values = []
            for data in band.getTile(0, 0, tile).getData():
                values += data

        return values

    def read(self, feature):
        raster = feature.getGeometry()
        if isinstance(raster, fmeobjects.FMERaster):
            rasterProperties = raster.getProperties()

            #Reading the tiles
            tiles = self.retrieveTiles(raster.getBand(0), rasterProperties.getNumRows(), rasterProperties.getNumCols())

            #Formatting output
            tiles_matrix = np.array([tiles])
            tiles_matrix = tiles_matrix.reshape(rasterProperties.getNumRows(), rasterProperties.getNumCols()) #From this point on i have my raster

        return tiles_matrix

class FeatureProcessor(object):
    def __init__(self):
        # Specify main parameters:
        # TBC
        pass

    def input(self, feature):

        #1. Pass feature to reader and receive extracted tiles
        MyRasterReader = RasterReader()
        raster = MyRasterReader.read(feature)
        print(raster)
        #2. Call some operation on these tiles and receive output

        #self.filtered = self.maxSlopeFilter(tiles_matrix, self.kernel, self.dhmax)
        self.pyoutput(feature)
        return

    def close(self):
        pass
