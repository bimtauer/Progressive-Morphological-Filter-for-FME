import fmeobjects
import numpy as np



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


    def write(self, dataArray):
        numRows, numCols = len(dataArray), len(dataArray[0]) # resolution
        xSpacing, ySpacing = 10.0, 10.0 # cell spacing in ground units
        xCellOrigin, yCellOrigin = 0.5, 0.5 # cell origin coordinates
        xOrigin, yOrigin = 0.0, numRows * ySpacing # left-top coordinates
        xRotation, yRotation = 0.0, 0.0 # rotation angle in degrees

        # Create a new raster.
        self.rasterProperties = fmeobjects.FMERasterProperties(numRows, numCols,
            xSpacing, ySpacing, xCellOrigin, yCellOrigin, xOrigin, yOrigin,
            xRotation, yRotation)
        # Contents of a tile for a band to be created.
        # A list of row data, each element is a list of column values.
        # Taking over our input raster's properties
        raster = fmeobjects.FMERaster(self.rasterProperties)

        #Clean the output array
        dataArray[np.isnan(dataArray)]=-9999.0
        dataArray = dataArray.tolist()

        bandTilePopulator = MyTilePopulator(dataArray)                 # <------------ changed class name

        bandName = "Slope Erosion Filtered"
        bandProperties = fmeobjects.FMEBandProperties(bandName,
            fmeobjects.FME_INTERPRETATION_REAL64,                         # <----------- changed from UInt8
            fmeobjects.FME_TILE_TYPE_FIXED,
            self.rasterProperties.getNumRows(), self.rasterProperties.getNumCols())
        nodataValue = fmeobjects.FMEReal64Tile(1, 1)                      # <----------- changed from UInt8
        nodataValue.setData([[-9999.0]])
        band = fmeobjects.FMEBand(bandTilePopulator,
            self.rasterProperties, bandProperties, nodataValue)
        raster.appendBand(band)
        return raster


class FeatureCreator(object):
    def __init__(self):
        # Specify main parameters:
        #TBC
        pass

    def input(self, feature):
        pass


    def close(self):

        dataArray = np.array([[ 0.0, 4.5, 6.2],
                    [ 1.2, 3.4, 6.9]])


        print(dataArray)
        #3. Package that output as a raster again
        MyWriter = RasterWriter()
        outputRaster = MyWriter.write(dataArray)


        # Create and output a feature containing the raster created above.
        feature = fmeobjects.FMEFeature()
        feature.setGeometry(outputRaster)
        print("Saving raster output.")
        self.pyoutput(feature)
