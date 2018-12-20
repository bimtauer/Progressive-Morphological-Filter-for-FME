import fme
import fmeobjects
import numpy as np
import sys
# The path where the python scripts are currently stored:
sys.path.insert(0,r'\\intra.malmo.se\dfs\gemensamt\Projekt\lidar-data (josbir)\FME\DEM Generation from Lidar (Project with Tim, Oct.18 - Dec. 18)\Git\dem-filter')
from FME_Scripts.FMEInOut import RasterReader, RasterWriter
from FME_Scripts.PMF import ProgressiveMorphologicalFilter

class FeatureProcessor(object):
    def __init__(self):
        pass
    def input(self,feature):
        # 1. Retrieve raster tiles as numpy array
        MyReader = RasterReader(feature)
        data, self.rasterProperties = MyReader.read(feature)
        print("Beginning filtering:")
        # 2. Create an instance of ProgressiveMorphologicalFilter with the respective parameters
        parameters = {'c' : 0.5,                     # The cell size of the input raster
                      'kernel_radius' : 15,          # The kernel size of the first filter iteration in meters
                      'initial_cutoff' : 0.5,        # The slope threshold of the first filter iteration
                      'average_sigma' : 7,           # The gaussian function used to average out local slope
                      'dh0' : 0.1,                   # The slope threshold for flat areas
                      'hole_cutoff' : -0.2}          # Threshold for individual holes beneath median in 3x3 kernel
        MyFilter = ProgressiveMorphologicalFilter(data, parameters)
        filtered = MyFilter.filter()
        filtered[np.isnan(filtered)] = -9999.0    # The FME nodata value of the raster  #TODO: make dynamic
        self.output = filtered
        print("Completed filtering operation.")
        return
    def close(self):
        # 3. Package that output as a raster again
        MyWriter = RasterWriter()
        outputRaster = MyWriter.write(self.output, self.rasterProperties)

        # 4. Create and output a feature containing the raster created above.
        feature = fmeobjects.FMEFeature()
        feature.setGeometry(outputRaster)
        print("Saving raster output.")
        self.pyoutput(feature)
        return
