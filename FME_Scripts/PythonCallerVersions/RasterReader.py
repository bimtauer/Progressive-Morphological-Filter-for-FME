import fmeobjects
import numpy as np
# Raster reading and writing functions based on code from Takashi: https://knowledge.safe.com/questions/38000/python-fme-objects-api-for-raster-manipulation.html

class RasterReader:
    def __init__(self, feature):
        raster = feature.getGeometry()
        #Only proceed if input is raster
        if isinstance(raster, fmeobjects.FMERaster):
            rasterProperties = raster.getProperties()
            self.numRows, self.numCols = rasterProperties.getNumRows(), rasterProperties.getNumCols()
            self.band = raster.getBand(0)

            bandProperties = self.band.getProperties()
            interpretation = bandProperties.getInterpretation() #Not necessary
            self.tile, interpret = self.tileAndInterpretation(interpretation, self.numRows, self.numCols)
            print("Tile is:", self.tile)
            print("Raster interpretion is {}".format(interpret))
        else:
            raise TypeError("Input is not a raster.")
        return

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

    def read(self, feature):
        if self.tile != None:
            values = []
            for data in self.band.getTile(0, 0, self.tile).getData():
                values += data
        else:
            raise TypeError("Data interpretation not yet added to tileAndInterpretation.")
        tiles_matrix = np.array([values])
        tiles_matrix = tiles_matrix.reshape(self.numRows, self.numCols) #From this point on i have my raster
        return tiles_matrix

class FeatureProcessor(object):
    def __init__(self):
        # Specify main parameters:
        # TBC
        pass

    def input(self, feature):

        #1. Pass feature to reader and receive extracted tiles
        MyRasterReader = RasterReader(feature)
        raster = MyRasterReader.read(feature)
        print(raster)
        #2. Call some operation on these tiles and receive output

        #self.filtered = self.maxSlopeFilter(tiles_matrix, self.kernel, self.dhmax)
        self.pyoutput(feature)
        return

    def close(self):
        pass
