"""
Code from Mark Lindeman starting to play with a weighted
Stouffer combination instead of Fisher's method.
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:03:54 2020

@author: Mark Lindeman
"""

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm, chi2
import math


def wtd_stouffer(pvalues, weights=np.array([0.5, 0.5])):
    assert len(pvalues) == len(weights)
    pvalues = [min(0.9999, p) for p in pvalues]
    zp = np.dot(norm.ppf(pvalues), weights) / math.sqrt(sum(weights**2))
    return(norm.cdf(zp))


# This optimizes weights when we can reliably predict, or hedge, what the p values 
# will/may be. Of course the weights must be set before the p values are known! 

# Note that this family of functions presently specializes on two p-values,
# while wtd_stouffer generalizes. Convenient because minimize optimizes on
# a single parameter, but not essential.

def opt_stouffer_p(pvalue1, pvalue2):
    opt = lambda x: wtd_stouffer(np.array([pvalue1, pvalue2]), 
                                 np.array([x, 1 - x]))
    return minimize(opt, 0.5, bounds=[(0, 1)]).fun  # doesn't give wt


def opt_stouffer_wt(pvalue1, pvalue2):
    opt = lambda x: wtd_stouffer(np.array([pvalue1, pvalue2]), 
                                 np.array([x, 1 - x]))
    return minimize(opt, 0.5, bounds=[(0, 1)]).x 


def opt_stouffer(pvalue1, pvalue2):
    opt = lambda x: wtd_stouffer(np.array([pvalue1, pvalue2]), 
                                 np.array([x, 1 - x]))
    result = minimize(opt, 0.5, bounds=[(0, 1)])
    return {"p": result.fun[0], "wts": [result.x[0], 1 - result.x[0]]}


def fisher_pvalue(pvalues):
    """
    Shorthand for fisher_combined_pvalue
    """
    if np.any(np.array(pvalues)==0):
        return 0
    obs = -2*np.sum(np.log(pvalues))
    return 1-chi2.cdf(obs, df=2*len(pvalues))


def print_ps(p1, p2):
    print("input p values:", p1, p2)
    print("Fisher:", fisher_pvalue([p1, p2]))
    print("Stouffer", wtd_stouffer([p1, p2]))
    res = opt_stouffer(p1, p2)
    print("optimal weighted Stouffer", res["p"])
    print("weights:", res["wts"])

"""
print_ps(0.1, 0.1)
print_ps(0.14, 0.16)
print_ps(0.05, 0.3)
print_ps(0.01, 0.5)
print_ps(0.02, 0.7)
print_ps(0.02, 1)
"""
