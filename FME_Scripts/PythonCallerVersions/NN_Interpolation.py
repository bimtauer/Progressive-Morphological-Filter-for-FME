import numpy as np
import fmeobjects
from scipy.interpolate import NearestNDInterpolator


################################################################################
# Reading/Writing for FME raster

class RasterReader:
    def __init__(self, feature):
        raster = feature.getGeometry()
        #Only proceed if input is raster
        if isinstance(raster, fmeobjects.FMERaster):
            self.rasterProperties = raster.getProperties()
            self.numRows, self.numCols = self.rasterProperties.getNumRows(), self.rasterProperties.getNumCols()
            self.band = raster.getBand(0)
            self.bandProperties = self.band.getProperties()
            interpretation = self.bandProperties.getInterpretation()
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
        tiles_matrix = tiles_matrix.reshape(self.numRows, self.numCols) #From this point on i have my numpy array
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

class RasterWriter:
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
                                                    rasterProperties.getNumRows(),
                                                    rasterProperties.getNumCols())
        nodataValue = fmeobjects.FMEReal64Tile(1, 1)                      # <----------- changed from UInt8
        nodataValue.setData([[-9999.0]])
        band = fmeobjects.FMEBand(bandTilePopulator, rasterProperties, bandProperties, nodataValue)
        raster.appendBand(band)
        return raster

################################################################################
# FME Main class that calls all the other ones

class FeatureProcessor(object):
    def __init__(self):
        # Specify main parameters:
        # TBC
        pass

    def input(self, feature):

        #1. Pass feature to reader and receive extracted tiles
        MyRasterReader = RasterReader(feature)
        data, self.rasterProperties = MyRasterReader.read(feature)
        #2. Turn input raster into lists of the known points coordinates and values
        data = np.where(data == -9999, np.nan, data)
        def rasterToList(raster):
            leni = raster.shape[0]
            lenj = raster.shape[1]
            i = np.arange(0, leni)
            j = np.arange(0, leni)
            jj, ii = np.meshgrid(i, j)

            x = jj.reshape(leni * lenj,)
            y = ii.reshape(leni * lenj,)
            z = raster.reshape(leni * lenj,)

            mask = np.isnan(z)

            xknown = x[~mask]
            yknown = y[~mask]
            zknown = z[~mask]

            coordinates = np.stack((yknown,xknown)).T
            tointerpolate = np.stack((y,x)).T

            return coordinates, zknown, tointerpolate, leni, lenj

        coordinates, values, tointerpolate, leni, lenj = rasterToList(data)

        Interpolator = NearestNDInterpolator(coordinates, values)
        self.output = Interpolator(tointerpolate).reshape(leni,lenj)
        print(self.output)
        return

    def close(self):
        MyWriter = RasterWriter()
        outputRaster = MyWriter.write(self.output, self.rasterProperties)

        # Create and output a feature containing the raster created above.
        feature = fmeobjects.FMEFeature()
        feature.setGeometry(outputRaster)
        print("Saving raster output.")
        self.pyoutput(feature)
        pass
