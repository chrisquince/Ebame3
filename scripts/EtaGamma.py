#!/usr/bin/env python

import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
import scipy as sp
import scipy.misc as spm
import math
import argparse
import cPickle
import logging

from operator import mul, div, eq, ne, add, ge, le, itemgetter
from itertools import izip
from itertools import compress
from numpy import array, log, exp
from scipy.special import gammaln
from scipy.optimize import minimize_scalar
from numpy.random import RandomState
from scipy.stats import chi2
from collections import defaultdict


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("eta_file", help="gene assignments")

    parser.add_argument("gamma_file", help="relative frequency")

    parser.add_argument("tau_file", help="relative frequency")

    args = parser.parse_args()

    eta = p.read_csv(args.eta_file, header=0, index_col=0)

    gamma = p.read_csv(args.gamma_file, header=0, index_col=0)

    tau = p.read_csv(args.tau_file, header=0, index_col=0)
    #import ipdb; ipdb.set_trace()
    tau_names = list(tau.columns.values)

    tau_names1 = tau_names[1:]

    tau_idxs = tau_names1[::4]
    tau_ints = [int(x) for x in tau_idxs]
    tau_ints_array = np.array(tau_ints)
    tau_ints_array /= 4

#,Position,4,5,6,7,12,13,14,15,20,21,22,23

    eta_matrix = eta.as_matrix()
    eta_matrix = eta_matrix[:,tau_ints_array]    
    gamma_matrix = gamma.as_matrix()
    gamma_matrix = gamma_matrix[:,tau_ints_array]
    #import ipdb; ipdb.set_trace()
    
    G = eta_matrix.shape[1]
    Z = eta_matrix.shape[0]
    assert G == gamma_matrix.shape[1]

    for g in range(G):
    
        for h in range(g):

            etaDist = np.sum(eta_matrix[:,g] != eta_matrix[:,h])
            
            etaDist = etaDist/float(Z)
            
            gammaDist = np.sum(np.absolute(gamma_matrix[:,g] - gamma_matrix[:,h]))
            
            gammaDist /= np.sum(gamma_matrix[:,g]) + np.sum(gamma_matrix[:,h])
            
            print str(g) + "," + str(h) + "," + str(etaDist) + "," + str(gammaDist)

if __name__ == "__main__":
    main(sys.argv[1:])
