# -*- coding: utf-8 -*-

import sys
import glob
import math
import copy
import random
import thread
import Hyperion
import matplotlib
import numpy as np
import ehtim as eh
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.interpolate import spline  
from scipy.signal import argrelextrema
from matplotlib.patches import Ellipse
from scipy.interpolate import UnivariateSpline 

import itertools
import seaborn as sns
plt.rc('text', usetex=True)
plt.rcParams['axes.linewidth'] = 2 #set the value globally
plt.rc('font', family='serif')
plt.rcParams["font.weight"] = "bold"
plt.rcParams.update({'font.size': 16})
cycol = itertools.cycle(sns.color_palette())
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

# Global variables ##########################
VERBOSE_FLAG = False
valve = None
pi = np.pi


np.random.seed(2)

# Classes ###################################
class HoughTransform(object):
    ''' 
        Performs the Hough Transform on a series of points, a given formula, and an established parameter space. Has a wide variety of knobs and levers, and its most powerful feature is parallelization. Designed to work in concert with the PontiFEX library.

        Attributes:
                  pset: a TUPLE or LIST containing the X and Y coordinates of the points to test
            parameters: LIST of TUPLES containing: [(1, 2, 3, 4), (1, 2, 3, 4) ... ] where:
                            1: STRING with the name of the parameter (ex: theta)
                            2: INTEGER with the resolution of the parameter axis in 1 (ex: 100)
                            3: FLOAT with the minimum value of the parameter range (ex: 0)
                            4: FLOAT with the maximum value of the parameter range (ex: 2*pi)
               formula: a FUNCTION that takes parameter inputs and returns ZERO. the first two inputs for the formula should be X and Y 
    ''' 

    def __init__(self, pset, parameters, formula):
        ''' 
            Initialize the attributes for the HoughTransform object
        '''
        # input attributes
        self.pset = pset
        self.parameters = parameters
        self.formula = formula

        # define empty accumulator -- to initialize it, run self.empty_accumulator(self.parameters)
        self.accumulator = None
        self.accumulator_flat = None

        # thresholds
        self.threshold = None
        self.dist_threshold = None

        # resulting distributions -- key will be the parameter, value will be a list of results
        self.distribution_dictionary = {}


    def empty_accumulator(self, parameters):
        '''
            creates an empty accumulator pased on the parameters provided in the initial attribution
        '''

        # create the tuples representing the length of each axis
        resolution = tuple([param_set[1] for param_set in parameters])
        self.accumulator = np.zeros(resolution)

        # define empty lists for each distribution
        for param_set in parameters:
            self.distribution_dictionary[param_set[0]] = []

    def single_thread_parameter_search(self, point, parameters):
        ''' 
            perform a parameter search with an arbitrary number of parameters over arbitrary ranges at arbitrary resolution
        '''

        # generate the parameter ranges
        ranges = [np.linspace(param_set[2], param_set[3], param_set[1]) for param_set in parameters]

        # generate all possible combinations of the parameters
        permutations = list(itertools.product(*ranges))

        # define threshold
        self.get_threshold(parameters)

        # temporarily flatten accumulator
        self.accumulator_flat = self.accumulator.flatten()

        for idx, p in enumerate(permutations):
            if abs(self.formula(point[0], point[1], *p)) <= self.threshold:
                self.accumulator_flat[idx] += 1
            else:
                continue


        # self.iterate_permutations(point, permutations)

        self.accumulator = self.accumulator_flat.reshape(self.accumulator.shape)



    def iterate_permutations(self, point, permutations, const=0):
        for idx, p in enumerate(permutations):
            if abs(self.formula(point[0], point[1], *p)) <= self.threshold:
                self.accumulator_flat[idx+int(const)] += 1
            else:
                continue


    def get_threshold(self, parameters):
        ''' 
            get the threshold based on a parameter list
        '''
        diffs = [abs(param_set[3]-param_set[2])/float(param_set[1]) for param_set in parameters]
        self.threshold = np.sqrt(np.sum(np.array(diffs)**2))

    def get_estimation(self, threaded='single', progress=True, title='', show=True, get2dist=False):
        '''
            run thruogh the parameter space and update the accumulator for all points
        '''
        self.empty_accumulator(self.parameters)

        if progress:
            toolbar_width = len(self.pset[0])
            sys.stdout.write("[%s]" % (" " * toolbar_width))
            sys.stdout.flush()
            sys.stdout.write("\b" * (toolbar_width+1))
        for idx, x in enumerate(self.pset[0]):
            if threaded=='single':
                self.single_thread_parameter_search((x, self.pset[1][idx]), self.parameters)
            elif threaded=='multi':
                self.multi_thread_parameter_search((x, self.pset[1][idx]), self.parameters)
            if progress:
                sys.stdout.write("-")
                sys.stdout.flush()
        sys.stdout.write("\n")

        self.accumulator = self.accumulator_flat.reshape(self.accumulator.shape)

        # find where accumulator peaks
        if self.dist_threshold == None:
            self.dist_threshold = np.max(self.accumulator)
        
        print "MAX", np.max(self.accumulator)
        locs = np.where(self.accumulator >= np.max(self.accumulator)/2.)

        # plt.hist2d(locs[0], locs[1], bins=5)
        # plt.imshow(self.accumulator)
        # plt.colorbar()
        # plt.show()
        # raw_input()


        returns = []
        # generate plots
        for idx, param_set in enumerate(self.parameters):

            fig = matplotlib.figure.Figure()
            ax = matplotlib.axes.Axes(fig, (0,0,0,0))
            binwidth = (param_set[3]-param_set[2])/100.
            # print np.linspace(param_set[2], param_set[3], param_set[1])[locs[idx]]
            n, bins, patches = ax.hist(np.linspace(param_set[2], param_set[3], param_set[1])[locs[idx]], bins=100, density=True)#np.arange(param_set[2], param_set[3] + binwidth, binwidth))
            del ax, fig
            (mu, sigma) = norm.fit(np.linspace(param_set[2], param_set[3], param_set[1])[locs[idx]])

            y = norm.pdf( bins, mu, sigma)
            # print bins, mu, sigma
            plt.plot(bins, y, linewidth=2, label=str(param_set[0]), c=cycol.next())
            plt.axvline(x=mu, c='black', linestyle='--', alpha=0.3)
            # plt.title(param_set[0])
            # plt.xlim(param_set[2], param_set[3])
            returns.append((mu, sigma))
        plt.title(title)
        plt.legend()

        if get2dist:
            sampling_matrix = []
            for x, ii in enumerate(self.accumulator):
                for y, kk in enumerate(ii):
                    for jj in range(int(kk)):
                        sampling_matrix.append((x, y))
            return sampling_matrix

        if not show:
            plt.clf()
            return returns
        plt.show()


        return returns

    def multi_thread_parameter_search(self, point, parameters, threads=10):
        ''' 
            perform a parameter search with an arbitrary number of parameters over arbitrary ranges at arbitrary resolution, using multithreading to split up the work
        '''

        # generate the parameter ranges
        ranges = [np.linspace(param_set[2], param_set[3], param_set[1]) for param_set in parameters]

        # generate all possible combinations of the parameters
        permutations = list(itertools.product(*ranges))

        # define threshold
        self.get_threshold(parameters)

        # temporarily flatten accumulator
        self.accumulator_flat = self.accumulator.flatten()

        chunks = list(split(permutations, threads))

        for idx, chunk in enumerate(chunks):
            thread.start_new_thread( self.iterate_permutations, (point, chunk, (len(chunk[0])*idx)-1) )
            # self.iterate_permutations(point, chunk, const=threads*idx)

# Functions #################################

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))


