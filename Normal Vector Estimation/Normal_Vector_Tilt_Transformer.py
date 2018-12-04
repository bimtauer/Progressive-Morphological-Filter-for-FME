
"""
This code is to be pasted into an FME PythonCaller which receives a raster as
input. Currently supported band interpretations include Real64 and Real32. Other
formats can be added.

The code first unpacks the fme raster object and loads it into a numpy array.
On this array it estimates surface unit normal vectors based on the nearest 4 points
to every point in the raster. The z values of those normals are stored as a measure
of local curvature, to be used for edge detection.

The output is again an FME raster object
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
            values = []
            for data in band.getTile(0, 0, tile).getData():
                values += data

        return values

###############################################################################

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

###############################################################################

class NormalEstimator(object):
    def __init__(self):
        # Specify main parameters:

        #Grid size
        self.a = 0.5

    # Returns a tuple (tile object, interpretation name).
    # This method returns (None, 'Unsupported')
    # when the specified interpretation was not supported,


    ###########################################################################
    #My code
    def normalVectorEstimation(self, raster):
    #Creating 3d matrix to store resulting vectors
        out = np.zeros((raster.shape[0], raster.shape[1],3))

        #Padding Raster
        padded = np.pad(raster, 1, "median")

        rows = padded.shape[0]
        cols = padded.shape[1]

        print("Beginning unit surface normal vector estimation...")
        #start = time.clock()
        for i in range(1, rows-1):
            for j in range(1, cols-1):
                #Surrounding points
                top = padded[i-1,j]
                bot = padded[i+1,j]
                left = padded[i,j-1]
                right = padded[i,j+1]

                # calculating surface normal based on 4 surrounding points:
                normal = np.array([(left - right), (top - bot), (2 * self.a)])
                mag = np.sqrt(normal.dot(normal))
                unit_normal = normal / mag

                out[i-1, j-1] = unit_normal  #A 3d matrix containing x,y,z of the normals

        #end = time.clock()
        #print("Done after {0:.2f} seconds.".format(end - start))
        return out

    ###########################################################################

    def input(self, feature):
        raster = feature.getGeometry()
        if isinstance(raster, fmeobjects.FMERaster):
            self.rasterProperties = raster.getProperties()

            MyRasterReader = RasterReader()
            #for band 0
            tiles = MyRasterReader.retrieveTiles(raster.getBand(0), self.rasterProperties.getNumRows(), self.rasterProperties.getNumCols())
            tiles_matrix = np.array([tiles])
            tiles_matrix = tiles_matrix.reshape(self.rasterProperties.getNumRows(), self.rasterProperties.getNumCols()) #From this point on i have my raster

            # Output Raster
            self.normals = self.normalVectorEstimation(tiles_matrix)
            print("Unit surface normal vector estimation done. Beginning output.")

        return

    ###########################################################################
    def close(self):
        # Contents of a tile for a band to be created.
        # A list of row data, each element is a list of column values.
        print("Converting to raster.")
        # Taking over our input raster's properties
        raster = fmeobjects.FMERaster(self.rasterProperties)

        #Clean the output array
        normals = self.normals.astype(float)
        normals[np.isnan(normals)]=-9999.0

        #For each coordinate one band
        for i in range(normals.shape[2]):
            coordinateArray = normals[:,:,i].tolist()

            # Create a new band and append it to the raster.
            # It's optional to specify Nodata value when creating a band.

            bandTilePopulator = MyTilePopulator(coordinateArray)                 # <------------ changed class name

            if i == 0:
                bandName = "X band"
            elif i == 1:
                bandName = "Y band"
            elif i == 2:
                bandName = "Z band"
            else:
                raise ValueError("Too many values to unpack")

            bandProperties = fmeobjects.FMEBandProperties(bandName,
                fmeobjects.FME_INTERPRETATION_REAL64,                         # <----------- changed from UInt8
                fmeobjects.FME_TILE_TYPE_FIXED,
                self.rasterProperties.getNumRows(), self.rasterProperties.getNumCols())
            nodataValue = fmeobjects.FMEReal64Tile(1, 1)                      # <----------- changed from UInt8
            nodataValue.setData([[-9999.0]])
            band = fmeobjects.FMEBand(bandTilePopulator,
                self.rasterProperties, bandProperties, nodataValue)
            raster.appendBand(band)

        # Create and output a feature containing the raster created above.
        feature = fmeobjects.FMEFeature()
        feature.setGeometry(raster)
        print("Saving raster output.")
        self.pyoutput(feature)
