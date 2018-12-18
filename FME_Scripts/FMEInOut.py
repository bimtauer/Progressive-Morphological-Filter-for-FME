import fmeobjects
import numpy as np
# Raster reading and writing functions based on code from Takashi: https://knowledge.safe.com/questions/38000/python-fme-objects-api-for-raster-manipulation.html

class RasterReader:
    def __init__(self, feature):
        raster = feature.getGeometry()
        #Only proceed if input is raster
        if isinstance(raster, fmeobjects.FMERaster):
            self.rasterProperties = raster.getProperties()
            self.numRows, self.numCols = self.rasterProperties.getNumRows(), self.rasterProperties.getNumCols()
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
        return tiles_matrix, self.rasterProperties

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
        newTile = fmeobjects.FMEReal64Tile(numRows, numCols)                      # <------------
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

class RasterWriter():
    def __init__(self):
        # Properties of a raster to be created.
        pass

    def write(self, dataArray, rasterProperties):
        # Taking over our input raster's properties
        raster = fmeobjects.FMERaster(rasterProperties)

        #Clean the output array
        dataArray[np.isnan(dataArray)]=-9999.0
        dataArray = dataArray.tolist()

        bandTilePopulator = MyTilePopulator(dataArray)                 # <------------ changed class name

        bandName = "Modified"
        bandProperties = fmeobjects.FMEBandProperties(bandName,
            fmeobjects.FME_INTERPRETATION_REAL64,                         # <----------- changed from UInt8
            fmeobjects.FME_TILE_TYPE_FIXED,
            rasterProperties.getNumRows(), rasterProperties.getNumCols())
        nodataValue = fmeobjects.FMEReal64Tile(1, 1)                      # <----------- changed from UInt8
        nodataValue.setData([[-9999.0]])
        band = fmeobjects.FMEBand(bandTilePopulator, rasterProperties, bandProperties, nodataValue)
        raster.appendBand(band)
        return raster
