# Progressive Morphological Filter

This is a morphological filter to extract ground level points from unfiltered
elevation rasters. It can be used to generate a DTM from a minimum elevation raster.

The algorithm is build on the methodology described in [Pingle, Clarke, McBride
(2013)](https://www.researchgate.net/publication/258333806_An_Improved_Simple_Morphological_Filter_for_the_Terrain_Classification_of_Airborne_LIDAR_Data). The main filtering proceedure is comparison of the input raster to a
morphological opening of that input raster with a stepwise increasing window.
It discards non-ground points if the difference between original point and
opening exceeds a slope threshold.
