## imports ##
import ehtim as eh
import HoughTransform
import Pontifex
import Hyperion
import numpy as np
import matplotlib.pyplot as plt
import random
import scipy.ndimage

from Tkinter import *
import tkFileDialog
from tkFileDialog import askopenfilename

## ask for a file ##
imageFilePath = askopenfilename()

## open the file in ehtim ##
im = eh.image.load_image(imageFilePath)

## display image ##
im.display()
raw_input()


## calculate the gradient of the image ##
testIm = im
imvecArray = testIm.imvec.reshape(testIm.xdim, testIm.ydim)
gradientImvecArray = np.gradient(imvecArray)
grad = (gradientImvecArray[0]**2 + gradientImvecArray[1]**2)**0.5
locs = np.where(grad > np.max(grad)/2.)

## display the gradient ##
plt.imshow(grad > np.max(grad)/2., interpolation='gaussian', cmap='gray')
plt.show()

print (len(list(locs[0])), len(list(locs[1])))

# ## alternative? calculate center with scipy ##
# (cx, cy) = Pontifex.get_non_colinear_orthocenter(im, npoints=len(locs)/10, points_tuple=(list(locs[0]), list(locs[1])), return_r=False)
# print "CENTER:", cx, cy

# rslocs = [[], []]
# numPoints = len(locs[0])
# for i in range(numPoints):
#     rIdx = random.choice(range(len(locs[1])))
#     rslocs[0].append(locs[0][rIdx])
#     rslocs[1].append(locs[1][rIdx])

# resolution=100

# print float(np.max(locs[0]) - np.min(locs[0] ))

# ## define the circle ##
# def paramshadow(x,y, r):
#         return (((x-cx)**2) + ((y-cy)**2)) - r**2

# # HoughTransformObject = HoughTransform.HoughTransform((rslocs[0], rslocs[1]), 
# #             [
# #                 ('h', resolution, float(np.min(rslocs[0])), float(np.max(rslocs[0]))), 
# #                 ('k', resolution, float(np.min(rslocs[1])), float(np.max(rslocs[1]))), 
# #                 ('r', resolution, 0.0, float(np.max(locs[0]) - np.min(locs[0] )))
# #             ], 
# #             paramshadow
# #         )
# HoughTransformObject = HoughTransform.HoughTransform((rslocs[0], rslocs[1]), 
#             [ 
#                 ('r', resolution, 0.0, float(np.max(locs[0]) - np.min(locs[0] )))
#             ], 
#             paramshadow
#         )

# res = HoughTransformObject.get_estimation(threaded='single', title=r'\textbf{Est. of center (h,k) and radius (r) in pixels}', show=True)

# # cx = res[0][0]
# # cy = res[1][0]
# rx = (testIm.psize/eh.RADPERUAS)*res[0][0]
# ry = (testIm.psize/eh.RADPERUAS)*res[0][0]

# if testIm.psize > 1e-6: 
#     print "probably in radians"
#     rx *= eh.RADPERUAS*(1/0.0174532925199)
#     ry *= eh.RADPERUAS*(1/0.0174532925199)


(cx, cy, rx, ry) = Hyperion.get_inner_circle(testIm, 0, 0, thresh=0.0)
print "Estimated radius:", rx

## plot the resulting circle back on the original image ##
plt.imshow(imvecArray, cmap='afmhot', interpolation='gaussian')
fig = plt.gcf()
ax = fig.gca()
circle1 = plt.Circle((cy, cx), rx, color='r', fill=False)
ax.add_artist(circle1)
plt.show()
