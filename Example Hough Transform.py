import numpy as np

from skimage.transform import (hough_line, hough_line_peaks,
                               probabilistic_hough_line)
from skimage.feature import canny
from skimage import data
import scipy

import matplotlib.pyplot as plt
from matplotlib import cm

z = normals[:,:,2]

tilt = np.arccos(z)

selection = tilt[200:400,250:450]

reduced = np.clip(selection-0.3, 0, selection.max())

scaled = reduced * 255.0/reduced.max()

image = scaled

image_closed = scipy.ndimage.morphology.grey_closing(image, size=(3,3))

image_opened = scipy.ndimage.morphology.grey_opening(image, size=(3,3))

fig, axes = plt.subplots(1, 3, figsize=(30, 5), sharex=True, sharey=True)
ax = axes.ravel()

ax[0].imshow(image, cmap=cm.gray)
ax[0].set_title('Input image')

ax[1].imshow(image_closed, cmap=cm.gray)
ax[1].set_title('Closed')

ax[2].imshow(image_opened, cmap=cm.gray)
ax[2].set_title('Opened')

for a in ax:
    a.set_axis_off()

plt.tight_layout()
plt.show()

###############################################################################
#Trying with raster as input image
image = raster[200:400,250:450]
image = image  * 255.0/image.max()

################################################################################
#
# Line finding using the Probabilistic Hough Transform
edges = canny(image, 2, 1, 25)
lines = probabilistic_hough_line(image, threshold=10, line_length=5,
                                 line_gap=3)

# Generating figure 2
fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)
ax = axes.ravel()

ax[0].imshow(image, cmap=cm.gray)
ax[0].set_title('Input image')

ax[1].imshow(edges, cmap=cm.gray)
ax[1].set_title('Canny edges')

ax[2].imshow(edges * 0)
for line in lines:
    p0, p1 = line
    ax[2].plot((p0[0], p1[0]), (p0[1], p1[1]))
ax[2].set_xlim((0, image.shape[1]))
ax[2].set_ylim((image.shape[0], 0))
ax[2].set_title('Probabilistic Hough')

for a in ax:
    a.set_axis_off()

plt.tight_layout()
plt.show()