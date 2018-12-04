# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 20:57:18 2018

@author: bimta
"""

import matplotlib.pyplot as plt
import numpy as np
import itertools

from progressive_morphological_filter import progressiveMorphologicalFilter
morphout = progressiveMorphologicalFilter(Trees)

morphout2 = Holes

fig = plt.figure()
ax = fig.add_subplot(111)
y1 = np.random.rand(10)
y2 = np.random.rand(10)
ys = itertools.cycle((morphout,morphout2))


def onclick(event):
    ax.clear()
    ax.imshow(next(ys))
    fig.canvas.draw()

cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()