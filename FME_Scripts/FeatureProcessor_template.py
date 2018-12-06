import fme
import fmeobjects
import sys
sys.path.insert(0,r'\\intra.malmo.se\dfs\gemensamt\Projekt\lidar-data (josbir)\FME\Tim\Git\dem-filter')
from FME_Scripts.InOut import RasterReader, RasterWriter

class FeatureProcessor(object):
    def __init__(self):
        pass
    def input(self,feature):
        MyReader = RasterReader(feature)
        data, self.rasterProperties = MyReader.read(feature)

        #2. Call some operation on data and assign to output
        print("Beginning raster operation.")
        self.output = data
        print("Completed operation.")
        return
    def close(self):
        # Package that output as a raster again
        MyWriter = RasterWriter()
        outputRaster = MyWriter.write(self.output, self.rasterProperties)

        # Create and output a feature containing the raster created above.
        feature = fmeobjects.FMEFeature()
        feature.setGeometry(outputRaster)
        print("Saving raster output.")
        self.pyoutput(feature)
        return
