# -*- coding: utf-8 -*-

import sys
import glob
import math
import copy
import random
import Hyperion
import numpy as np 
import ehtim as eh
import matplotlib.pyplot as plt 
from scipy.signal import argrelextrema
from matplotlib.patches import Ellipse
from scipy.interpolate import UnivariateSpline 

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 16})
plt.rcParams['axes.linewidth'] = 2 #set the value globally
plt.rcParams["font.weight"] = "bold"

# Global variables ##########################
VERBOSE_FLAG = False
valve = None
pi = np.pi




# Classes ###################################

class AnalysisObject(object):

    def __init__(self, image, intervals=75, scale='log'):
        self.image = image
        self.intervals = intervals
        self.scale = scale

        # things that we will eventually want to fill
        self.radius_dictionary = {
            'fluxcenter':0,
            'circumcenter':0,
            'combined':0,
            'longdist':0,
            'hough':0,
        }
        self.error_dictionary = {
            'systematic':0,
            'corruption':0,
            'asymmetry':0
        }
        self.points = {
            'xcoord':[],
            'ycoord':[]
        }
        self.centers = {
            'fluxcenter':(0,0),
            'circumcenter':(0,0),
            'combined':(0,0)
        }
        self.parameters = {
            'intervals':self.intervals,
            'npoints':self.intervals,
            'iterations':0
        }
        self.intensityProfiles = {
            0:Hyperion.get_horizontal_peaks(image, return_center_row=True),
            1:Hyperion.get_horizontal_peaks(image, angle= np.pi/2., return_center_row=True),
            2:Hyperion.get_horizontal_peaks(image, angle= np.pi, return_center_row=True),
            3:Hyperion.get_horizontal_peaks(image, angle= 3*np.pi/2., return_center_row=True)
        }



    def pix_to_uas(self, pix):
        if self.image.psize > 1e-6: 
            return pix*(self.image.psize)*(1/0.0174532925199)
        
        return pix*(self.image.psize/eh.RADPERUAS)

    def uas_to_pix(self, uas):
        return uas/(self.image.psize/eh.RADPERUAS)


    def get_radius_hough_transform(self):
        radius = Hyperion.get_inner_circle(self.image, 0, 0)[2]
        print "ESTIMATED RADIUS", radius
        self.radius_dictionary['hough'] = radius
        return radius

    def get_radius_FLUX_METHOD(self):
        shadow_estimate, error = get_shadow_size_from_image(self.image, 
                                                            interval=self.parameters['intervals'], 
                                                            display=False)
        self.radius_dictionary['fluxcenter'] = self.pix_to_uas(shadow_estimate)
        return self.pix_to_uas(shadow_estimate)

    def get_radius_CIRCUMCENTER_METHOD(self):
        shadow_estimate = get_non_colinear_orthocenter(self.image, npoints=self.parameters['npoints'], return_r=True)
        self.radius_dictionary['circumcenter'] = self.pix_to_uas(shadow_estimate)
        return self.pix_to_uas(shadow_estimate)

    def iterative_peak_finding(self):
        x, y = iterative_shadow_estimate(self.image, iterations=self.parameters['iterations'], display=False, interval=self.parameters['intervals'], scale=self.scale)
        # x = np.asarray(x)
        # y = np.asarray(y)
        # print len(x)
        # x = x[np.where(y > 60)]
        # y = y[np.where(y > 60)]
        # print len(x)
        self.points['xcoord'] = x
        self.points['ycoord'] = y
        shadow_estimate = get_non_colinear_orthocenter(self.image, npoints=self.parameters['npoints'], return_r=True, points_tuple=(self.points['xcoord'], self.points['ycoord']))
        self.radius_dictionary['combined'] = self.pix_to_uas(shadow_estimate)
        self.get_center_CIRCUMCENTER_METHOD()
        return self.pix_to_uas(shadow_estimate)

    def get_center_CIRCUMCENTER_METHOD(self):
        (cx, cy) = get_non_colinear_orthocenter(self.image, npoints=self.parameters['npoints'], points_tuple=(self.points['xcoord'], self.points['ycoord']))
        self.centers['circumcenter'] = (cx, cy)

    def get_radius_LONGDIST_METHOD(self):
        x, y = iterative_shadow_estimate(self.image, iterations=self.parameters['iterations'], display=False, interval=self.parameters['intervals'], scale=self.scale)
        self.points['xcoord'] = x
        self.points['ycoord'] = y
        diameters = []
        for idx, x in enumerate(self.points['xcoord']):
            max_dist = 0
            start_point = (x, self.points['ycoord'][idx])
            for e_idx, e_x in enumerate(self.points['xcoord']):
                end_point = (e_x, self.points['ycoord'][e_idx])
                dist = get_distance_between_two_points(start_point[0], start_point[1], end_point[0], end_point[1])
                if dist > max_dist:
                    max_dist = dist
                else:
                    continue
            diameters.append(max_dist)
        self.radius_dictionary['longdist'] = np.median(diameters)/2.

    def account_for_offset(self, key):
        print "centers: ", self.centers[key][0], self.centers[key][1]
        self.points['xcoord'] = np.array(self.points['xcoord']) - self.centers[key][1]
        self.points['ycoord'] = np.array(self.points['ycoord']) - self.centers[key][0]

    def centerImageGaussian(self):
        global valve
        fov = 150.*eh.RADPERUAS
        npix = 128
        base_image = eh.image.Image(np.zeros((npix,npix)), fov/npix, eh.RA_DEFAULT, eh.DEC_DEFAULT, rf=230e9)

        print "making a gaussian, just for you"
        base_image = base_image.add_gauss(1.0, [40.*eh.RADPERUAS, 40.*eh.RADPERUAS, 0.*180./np.pi, 0.*eh.RADPERUAS, 0.*eh.RADPERUAS])
        im_array = [self.image]
        im_array = base_image.align_images(im_array)[0]
        self.image = im_array[0]


    ### function attributes specific to the model comparison pipeline ###
    def get_asymmetry(self):
        x, y = iterative_shadow_estimate(self.image, iterations=self.parameters['iterations'], display=False, interval=self.parameters['intervals'], scale=self.scale)
        self.points['xcoord'] = x
        self.points['ycoord'] = y
        diameters = []
        for idx, x in enumerate(self.points['xcoord']):
            max_dist = 0
            start_point = (x, self.points['ycoord'][idx])
            for e_idx, e_x in enumerate(self.points['xcoord']):
                end_point = (e_x, self.points['ycoord'][e_idx])
                dist = get_distance_between_two_points(start_point[0], start_point[1], end_point[0], end_point[1])
                if dist > max_dist:
                    max_dist = dist
                else:
                    continue
            diameters.append(max_dist)

        return np.std(diameters)

    def get_mean_width(self):
        widths = []
        for theta in np.linspace(0, 2*np.pi, self.intervals):
            (c, e) = Hyperion.get_horizontal_peaks(self.image, angle=theta, return_error=True)
            widths.append(e)
        print widths
        print np.median(widths)
        return np.median(widths)


    def get_flux_inside_radius(self):
        r = self.uas_to_pix(self.radius_dictionary['combined'])
        center = self.centers['combined']
        h = self.image.xdim
        w = self.image.ydim
        mask = create_circular_mask(h, w, center=center, radius=r)

        masked_img = self.image.imvec.reshape((h, w)).copy()
        masked_img[~mask] = 0

        return np.sum(masked_img) / np.sum(self.image.imvec)

    def get_flux_outside_radius(self):
        r = self.uas_to_pix(self.radius_dictionary['combined'])
        center = self.centers['combined']
        h = self.image.xdim
        w = self.image.ydim
        mask = create_circular_mask(h, w, center=center, radius=r)

        masked_img = self.image.imvec.reshape((h, w)).copy()
        masked_img[mask] = 0

        return np.sum(masked_img) / np.sum(self.image.imvec)

    ### end function attributes specific to the model comparison pipeline ###


# Functions #################################
def print_verbose(message):
    """
        Prints a message if the verbose flag is set. To 
        set VERBOSE_FLAG to true, use SET_VERBOSE()
        Args:
            message: message to print if verbose flag is on
    """
    global VERBOSE_FLAG
    if VERBOSE_FLAG:
        print message

def SET_VERBOSE(flag=False):
    """
        Function to either pass a specific flag to VERBOSE_FLAG
        or to just switch the state
    """
    global VERBOSE_FLAG
    if flag:
        VERBOSE_FLAG = flag
        return 1
    if not flag:
        VERBOSE_FLAG = not VERBOSE_FLAG

def get_center_from_matrix(vec, thresh):
    """ 
        given a vector and a threshold, return the center of the matrix
        Args:
            vec:    one dimensional vector of intensities
            thresh: percentage of flux to mask
        Returns:
            center: tuple containing xy coordinate of center
    """
    vec = vec.reshape((int(np.sqrt(vec.shape[0])), int(np.sqrt(vec.shape[0]))))
    (x,y) = np.nonzero(vec > np.max(vec)*float(thresh))
    center = (int(np.mean(x)), int(np.mean(y)))
    print_verbose(center)
    return center

def get_linearity(points):
    # points is  a list of lists, where each sublist is a point [x,y]
    coords = map(list, zip(*points))
    p = np.polyfit(map(float, coords[0]), map(float,coords[1]), 1)
    chi_squared = np.sum((np.polyval(p, coords[0]) - coords[1]) ** 2)
    return chi_squared



def get_non_colinear_orthocenter(im, npoints=10, points_tuple=False, return_r=False):
    global_xs = []
    global_ys = []
    global_error = []

    (cx, cy, rx, ry) = Hyperion.get_outer_circle(im)
    (inner_cx, inner_cy, radiusx, radiusy) = Hyperion.get_inner_circle(im, rx, ry)

    # plt.imshow(im.imvec.reshape(im.xdim, im.ydim))
    # plt.show()
    # raw_input()

    interv = 75
    for i in range(0,int(interv)+3):
        angle = 2*i*pi/float(int(interv))
        # print "ANGLE: ", angle/pi
        (xs, ys), error = get_single_point_max(im, angle=angle, get_center=(cx,cy), flux_thresh=0.1, ring_thresh=0.0, return_error=True, bounds=False)
        if xs == False:
            return False, False
        global_error.append(error)
        for x in xs:
            global_xs.append(x)
        for y in ys:
            global_ys.append(y)

    if points_tuple:
        global_xs = points_tuple[0]
        global_ys = points_tuple[1]

    lftrng = 0.05
    rtrng = 0.05
    points_to_delete = []
    for i, x in enumerate(global_xs):
        x = x
        y = global_ys[i]
        if x >= int((im.xdim/2.)-lftrng*im.xdim) and x < int((im.xdim/2.)+lftrng*im.xdim) and y > int((im.ydim/2.)-rtrng*im.ydim) and y <= int((im.ydim/2.)+rtrng*im.ydim):
            points_to_delete.append(i)

    print_verbose("DELETING {0} BAD POINTS".format(len(points_to_delete)))
    for i, point in enumerate(points_to_delete):
        del global_xs[point-i]
        del global_ys[point-i]

    cx_sum, cy_sum, r = [], [], []
    points_permutations, chi_squares = [], []

    for k in range(1,npoints+1):
        ## select three points at random
        choices = range(0, len(global_xs)-1)
        indices = []

        for i in range(3):
            idx = random.choice(choices)
            indices.append(idx)
            choices.remove(idx)

        ## get the three points
        A = [global_xs[indices[0]], global_ys[indices[0]]]
        B = [global_xs[indices[1]], global_ys[indices[1]]]
        C = [global_xs[indices[2]], global_ys[indices[2]]]

        points_permutations.append([A, B, C])
        chi_squares.append(get_linearity([A, B, C]))

    max_indices = np.argsort(chi_squares)[5:]
    # print chi_squares
    # print np.where(chi_squares == np.max(chi_squares))
    # print max_indices
    # raw_input()
    k = 0
    while 1:
        k+=1
        try:
            for idx in max_indices:
                A = points_permutations[idx][0]
                B = points_permutations[idx][1]
                C = points_permutations[idx][2]

                D = 2*( A[0]*(B[1]-C[1]) + B[0]*(C[1] - A[1]) + C[0]*(A[1] - B[1]))
                # print D

                cx = int((1./D)*((A[0]**2 + A[1]**2)*(B[1] - C[1]) +  (B[0]**2 +B[1]**2)*(C[1] - A[1]) +  (C[0]**2 + C[1]**2)*(A[1] - B[1])))

                cy = int((1./D)*((A[0]**2 + A[1]**2)*(C[0] - B[0]) +  (B[0]**2 +B[1]**2)*(A[0] - C[0]) +  (C[0]**2 + C[1]**2)*(B[0] - A[0])))

                cx_sum.append(cx)
                cy_sum.append(cy)
                r.append(np.linalg.norm(np.asarray(A)-np.asarray([float(cx), float(cy)])))
        except ValueError:
            print "VALUEERROR"
            if k >= 10: break
            continue
        except OverflowError:
            print "VALUEERROR"
            if k >= 10: break
            continue
        break

        # print np.linalg.norm(np.asarray(A-)cannp.asarray([float(cx_sum/k), float(cy_sum/k)]))

    cx = float(np.median(cx_sum))
    cy = float(np.median(cy_sum))

    # for idx in max_indices:
    #     A = points_permutations[idx][0]
    #     r.append(np.linalg.norm(np.asarray(A)-np.asarray([float(cx), float(cy)])))

    # print_verbose( "CENTER:"float(np.median(cx_sum)), float(np.median(cy_sum)))
    # print "CENTER:", float(cx/npoints), float(cy/npoints)
    if return_r:
        return np.median(r)

    return float(np.median(cx)), float(np.median(cy))






def get_orthocenter(im, npoints=10, points_tuple=False, return_r=False):
    global_xs = []
    global_ys = []
    global_error = []

    (cx, cy, rx, ry) = Hyperion.get_outer_circle(im)
    (inner_cx, inner_cy, radiusx, radiusy) = Hyperion.get_inner_circle(im, rx, ry)

    # plt.imshow(im.imvec.reshape(im.xdim, im.ydim))
    # plt.show()
    # raw_input()

    interv = 75
    for i in range(0,int(interv)+3):
        angle = 2*i*pi/float(int(interv))
        # print "ANGLE: ", angle/pi
        (xs, ys), error = get_single_point_max(im, angle=angle, get_center=(cx,cy), flux_thresh=0.1, ring_thresh=0.0, return_error=True, bounds=False)
        if xs == False:
            return False, False
        global_error.append(error)
        for x in xs:
            global_xs.append(x)
        for y in ys:
            global_ys.append(y)

    if points_tuple:
        global_xs = points_tuple[0]
        global_ys = points_tuple[1]

    cx_sum, cy_sum, r = [], [], []

    for k in range(1,npoints+1):
        ## select three points at random
        choices = range(0, len(global_xs)-1)
        indices = []

        for i in range(3):
            idx = random.choice(choices)
            indices.append(idx)
            choices.remove(idx)

        ## get the three points
        A = [global_xs[indices[0]], global_ys[indices[0]]]
        B = [global_xs[indices[1]], global_ys[indices[1]]]
        C = [global_xs[indices[2]], global_ys[indices[2]]]

        D = 2*( A[0]*(B[1]-C[1]) + B[0]*(C[1] - A[1]) + C[0]*(A[1] - B[1]))
        # print D

        cx = int((1./D)*((A[0]**2 + A[1]**2)*(B[1] - C[1]) +  (B[0]**2 +B[1]**2)*(C[1] - A[1]) +  (C[0]**2 + C[1]**2)*(A[1] - B[1])))

        cy = int((1./D)*((A[0]**2 + A[1]**2)*(C[0] - B[0]) +  (B[0]**2 +B[1]**2)*(A[0] - C[0]) +  (C[0]**2 + C[1]**2)*(B[0] - A[0])))

        cx_sum.append(cx)
        cy_sum.append(cy)
        r.append(np.linalg.norm(np.asarray(A)-np.asarray([float(cx), float(cy)])))
        # print np.linalg.norm(np.asarray(A)-np.asarray([float(cx_sum/k), float(cy_sum/k)]))

    # print "CENTER:", float(cx/npoints), float(cy/npoints)
    if return_r:
        return np.median(r)

    return float(np.median(cx)), float(np.median(cy))

def get_single_point_max(image, thresh=0.0, get_center=True, angle=0, flux_thresh = 0.0, ring_thresh=0.0, return_error=False, bounds=False, scale='log'):
    """ 
        get a single horizontal peak flux from 
        two sides of a centerpoint

        !!! IF YOU WANT get_center TO BE ANYTHING BUT TRUE
            YOU MUST SET IT EQUAL TO AN (X,Y) TUPLE!
        !!!
    """
    ''' sanitize input '''
    thresh = float(thresh)
    print_verbose("Threshold set to {0}".format(thresh))
    im = image
    im = im.rotate(angle)

    ''' get center of image if get_center is flagged '''
    if get_center ==  True:    
        (CENTER_X, CENTER_Y) = get_orthocenter(im)
        # (CENTER_X, CENTER_Y) = get_center_from_matrix(im.imvec, thresh)
        print_verbose("Got center as: x={0}, y={1}".format(CENTER_X, CENTER_Y))
    else:
        (CENTER_X, CENTER_Y) = get_center


    # lftrng = 0.0005
    # rtrng = 0.15
    # for x in range(int((im.xdim/2.)-lftrng*im.xdim), int((im.ydim/2.)+rtrng*im.xdim)):
    #     for y in range(int((im.xdim/2.)-rtrng*im.xdim), int((im.ydim/2.)+rtrng*im.xdim)):
    #         im.imvec.reshape((im.xdim, im.ydim))[x][y] = np.min(im.imvec)

    vec = im.imvec

    if scale=='linear':
        vec = vec
    elif scale=='log':
        vec = np.log(vec)

    def add_center(a):
        return a#+CENTER_X
    if bounds:
        bounds = map(add_center, bounds)
    # reshape the image vector based on the squareroot of the length
    vec = vec.reshape((int(np.sqrt(vec.shape[0])), int(np.sqrt(vec.shape[0]))))
    print_verbose("Vec shape now: {0}".format(vec.shape))
    # print CENTER_X, CENTER_Y, "CENTER"
    center_row = vec[CENTER_X,:]
    center_row = center_row[:CENTER_Y]

    # plt.plot(center_row)
    # plt.title("Cross sectional intensity profile")
    # plt.xlabel('Position')
    # plt.ylabel('Intensity')
    # plt.show()

    # center_row = smooth(center_row, 10)
    splice = vec[CENTER_X,:]
    mx = np.median(np.where(center_row == np.max(center_row)))

    if scale == 'log':
        mx = argrelextrema(center_row, np.greater)[0]
        xx = 0
        for x in mx:
            if center_row[x] > center_row[xx]:
                xx = x
        mx = xx
    try:
        if bounds:
            splice = center_row[len(center_row)-int(bounds[1]):len(center_row)-int(bounds[0])]
            mx = np.where(splice == np.max(splice))[0] + len(center_row)-int(bounds[1])
    except ValueError:
        print "TOO MANY ITERATIONS--try reducing the number of iterations by 1."
        return (False, False), False

    # print "MX", mx
    # print vec.size
    # error = FWHM(range(len(vec[CENTER_X,:])), vec[CENTER_X,:])
    k = ((len(np.where(vec > flux_thresh*np.max(vec))[0])/float(vec.size)))
    # error = FWHM(range(len(vec[CENTER_X,:])), np.array(vec[CENTER_X,:]))*k
    error = 0
    # print "ERROR", error
    return convert_coordinates(mx, CENTER_Y, angle, CENTER_X, CENTER_Y), error

def convert_coordinates(x1, y1, theta, cx, cy):
    hx1 = x1-cx

    # x1 = x1-(x1+x2)/2.
    hypotenuse = math.sqrt(x1**2 + y1**2)
    x1prime = hx1*math.cos(theta)
    y1prime = hx1*math.sin(theta) 
    # print x1prime+cx, y1prime +cy
    return ([x1prime+cx], [y1prime+cy])

def FWHM(X,Y):
    half_max = max(Y) / 2.
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - Y[])
    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    # plt.plot(X,d) #if you are interested
    # plt.show()
    #find the left and right most indexes
    left_idx = np.where(d > 0)[0][0]
    try:
        right_idx = np.where(d < 0)[-1][0]
    except IndexError:
        right_idx = len(Y) - 1
    # print left_idx
    # print right_idx
    return (X[right_idx] - X[left_idx])/2. #return the difference (full width)

def get_rho(n, r):
    return 2.*r*np.sin(np.pi/float(n+2))

def corruption_error(xs, ys, vec, estr):
    rho = get_rho(len(xs), estr)
    # print "RHO", rho
    # print xs, ys
    # xs = [x[0] for x in xs]
    # ys = [y[0] for y in ys]
    # return round(round((((len(np.where(vec > 0.1*np.max(vec))[0])/float(vec.size)))*((np.mean(np.absolute(np.diff((xs))))+np.mean(np.absolute(np.diff((ys)))))/2)),2)/rho, 2)

    k = 1#(len(np.where(vec > 0.0*np.max(vec))[0])/float(vec.size))/rho
    return round(k*np.mean(np.sort([get_distance_between_two_points(xs[i], ys[i],xs[i+1], ys[i+1]) for i in range(len(xs)-1)])[:])/rho, 2)


def get_distance_between_two_points(x1,y1,x2,y2):
    dist = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return dist

def get_shadow_size(image, points, center):
    radii = []
    for i, x in enumerate(points[0]):
        radii.append(get_distance_between_two_points(x, points[1][i], points[0][i-int((len(points[0]))/2.)], points[1][i-int((len(points[0]))/2.)])/2.)
        # radii.append(get_distance_between_two_points(x, points[1][i], center[0], center[1]))

    return (np.mean(radii), np.std(radii))

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def get_shadow_size_from_image(im, interval=75, display=True, bounds=False, return_points=False, scale='log'):

    global_xs = []
    global_ys = []
    global_error = []

    (cx, cy, rx, ry) = Hyperion.get_outer_circle(im)
    (inner_cx, inner_cy, radiusx, radiusy) = Hyperion.get_inner_circle(im, rx, ry)

    # (cx, cy) = get_orthocenter(im)

    print cx, cy
    interv = interval
    for i in range(0,int(interv)+3):
        angle = 2*i*pi/float(int(interv))
        # print "ANGLE: ", angle/pi
        (xs, ys), error = get_single_point_max(im, angle=angle, get_center=(cx,cy), flux_thresh=0.1, ring_thresh=0.0, return_error=True, bounds=bounds, scale=scale)
        if xs == False:
            if return_points: return False, False, False, False
            return False, False
        global_error.append(error)
        for x in xs:
            global_xs.append(x)
        for y in ys:
            global_ys.append(y)

    new_gxs = [np.mean([global_xs[i], global_xs[i+2]]) for i in range(len(global_xs)-3)]
    new_gys = [np.mean([global_ys[i], global_ys[i+2]]) for i in range(len(global_xs)-3)]

    (shadow_estimate, asymmetry) = get_shadow_size(im, (global_xs, global_ys), (cx, cy))
    # print asymmetry/im.xdim

    corr_error = corruption_error(global_xs, global_ys, im.imvec, shadow_estimate/(im.psize/eh.RADPERUAS))
    # print "CORR_ERROR", corr_error


    t = np.linspace(0, 2*pi, 100)
    sys_error = round((100*np.mean(global_error)/(im.xdim)),2)

    if display:
        plt.imshow(im.imvec.reshape(im.xdim, im.ydim), cmap='plasma', interpolation='gaussian')
        plt.colorbar()
        plt.scatter(
                    new_gxs,
                    new_gys, 
                    s=40, 
                    marker='+',
                    c='royalblue', 
                    label='Photon ring FEX \n $\sigma_s= ${0}\% \n $\sigma_c= ${1} $\\rho$'.format(
                                sys_error, 
                                corr_error
                            )
                )
        plt.grid(color='gray',linestyle='--')
        plt.title(r'\textbf{Feature extraction, simulated reconstruction of Sgr A*}')
        plt.xlabel(r'\textbf{X position (pixels)}')
        plt.ylabel(r'\textbf{y position (pixels)}')
        plt.legend()
        plt.show()

    # build dictionary of errors:
    error_dict = {
        'se':sys_error,
        'ce':corr_error,
        'as':asymmetry
    }

    if return_points:
        return shadow_estimate, error_dict, global_xs, global_ys

    return shadow_estimate, error_dict


def iterative_shadow_estimate(im, iterations=1, display=True, interval=75, scale='log'):

    ## bootstrap and get the initial mean value. we will take the first set of bounds from the asymmetry
    std_coefficient = 5
    xs = []
    ys = []
    shadows = []
    corr_errors = []
    shadow_estimate, error, x, y = get_shadow_size_from_image(im, display=False, return_points=True, interval=interval, scale=scale)
    shadows.append(shadow_estimate)
    corr_errors.append(error['ce'])
    xs.append(x)
    ys.append(y)
    for it in range(iterations):
        # std_coefficient = 5*(5./iterations)-5*it/iterations
        std_coefficient = 1.5*iterations-it
        # print "STD", std_coefficient
        # print "BOUDNS", int(shadow_estimate-std_coefficient*error['as']), int(shadow_estimate+std_coefficient*error['as'])
        shadow_estimate, error, x,y = get_shadow_size_from_image(im, bounds=[shadow_estimate-std_coefficient*error['as'], shadow_estimate+std_coefficient*error['as']], display=display, return_points=True, interval=interval, scale=scale)
        if shadow_estimate:
            xs.append(x)
            ys.append(y)
            shadows.append(shadow_estimate)
            corr_errors.append(error['ce'])
            continue
        break


    min_corr = np.where(corr_errors == np.min(corr_errors))[-1][-1]
    lftrng = 0.05
    rtrng = 0.05

    # print int((im.xdim/2.)-lftrng*im.xdim), int((im.xdim/2.)+lftrng*im.xdim)
    # print int((im.ydim/2.)-rtrng*im.ydim), int((im.ydim/2.)+rtrng*im.ydim)
    # for x in range(int((im.xdim/2.)-lftrng*im.xdim), int((im.ydim/2.)+rtrng*im.xdim)):
    #     for y in range(int((im.ydim/2.)-rtrng*im.ydim), int((im.ydim/2.)+rtrng*im.ydim)):
    #         if float(x) in map(np.floor, xs[min_corr]):
    #             if float(y) in map(np.floor, ys[min_corr]):
    #                 whereinx = np.where(map(np.floor, xs[min_corr]) == x)
    #                 whereiny = np.where(map(np.floor, ys[min_corr]) == y)
    #                 if(len(set(whereinx[0]).intersection(whereiny[0]))) > 0:
    #                     common_element = list(set(whereinx).intersection(whereiny))
    #                     print common_element
    #                     del xs[min_corr][common_element]
    #                     del ys[min_corr][common_element]
    points_to_delete = []
    for i, x in enumerate(xs[min_corr]):
        x = x
        y = ys[min_corr][i]
        if x >= int((im.xdim/2.)-lftrng*im.xdim) and x < int((im.xdim/2.)+lftrng*im.xdim) and y > int((im.ydim/2.)-rtrng*im.ydim) and y <= int((im.ydim/2.)+rtrng*im.ydim):
            points_to_delete.append(i)

    print "DELETING {0} BAD POINTS".format(len(points_to_delete))
    for i, point in enumerate(points_to_delete):
        del xs[min_corr][point-i]
        del ys[min_corr][point-i]

    # plt.plot(xs[min_corr])
    # plt.plot(ys[min_corr])
    # plt.show()
    # print xs, ys
    return xs[min_corr], ys[min_corr]

    print 'MEAN ANGULAR ESTIMATE:', shadows[np.where(corr_errors == np.min(corr_errors))[-1][-1]]*(im.psize/eh.RADPERUAS)
    print 'CIRCUMCENTER ESTIMATE:', get_orthocenter(im, points_tuple=(x,y), return_r=True)*(im.psize/eh.RADPERUAS)

    # print np.asarray(shadows)*(im.psize/eh.RADPERUAS)
    print corr_errors
    ## return final shadow size
    return shadows[np.where(corr_errors == np.min(corr_errors))[-1][-1]]



### functions specific to the model comparison pipeline ###

def create_circular_mask(h, w, center=None, radius=None):

    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask