"""
Functions for a bayesian probability ratio test comparison audit.
actually maybe not bayesian, why bother with bayesian?
ah, past oliver, glad you asked
if you just take alt hypothesis to be as announed then your alt dist has
probabilty zero for anything less than all matches. that is silly. you
want to allow for matches

Oliver Broadrick 2020
"""

import math
import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt

def generate_comparison_distributions(k, n, N, alpha, Vw, Vl, null_margin=0):
    """
    Generates the probability distributions over number of matches for both the null and alternative hypothesis of a comparison audit.
    This is proof of concept so I've kept things simple in the following ways:
        Only mismatches considered are 2-vote overstatements (This could be extended to multiple variables for 2-vote and 1-vote overstatements and understatements)
        The alternative half of the prior is all one point (This could be exteneded first to uniform distribution but also to any efficiency-maximizing distribution)

    Parameters:
        k : number of ballots that match their CVR or are understatements
            all others are 2-vote overstatements
        N : total ballots
        alpha : risk limit
        Vw : reported votes for winner
        Vl : reported votes for loser
        null_margin : the margin in votes under the null
                    Optional: defaults to 0 (used for stratified application)

    Returns:
        float : pvalue?
    """
    # reported margin from passed tallies
    reported_margin = Vw - Vl
    
    # b is 2-vote overstatements under the null
    #   reported_margin - 2b <= null_margin
    b_null = math.floor((reported_margin - null_margin) / 2)
    print("b_null:",b_null)
    # maximum number of matches under the null
    #   x = N - 2b
    x_null = N - 2 * b_null
    print("x_null:",x_null)

    # risk limiting prior (0's, .5 at x_null, then uniform)
    prior = np.concatenate([np.zeros(x_null), np.array([.5]), np.full(N - x_null, .5 / (N - x_null))])

    # generate distributions
    dist_range = range(0, n + 1)

    # null dist is only affected by one point on prior
    null_dist = binom.pmf(dist_range, n, x_null / N)

    # alt dist affected by uniform prior
    alt_dist = np.empty(n + 1)
    for k in dist_range:
        for x in range(x_null + 1, N + 1):
            alt_dist[k] = binom.pmf(k, n, x / N) * prior[x]
    # normalize
    alt_dist = alt_dist / sum(alt_dist)


    """
    # A WAY TO DO THIS WITHOUT A UNIFORM PRIOR(JUST ONE POINT)
    # fraction of ballots which we assume to be mismatches under null
    # why have a uniform prior when this is less work!
    # 0 <= fractional_mismatches <= (N - x_null) / N
    fractional_mismatches = .0001 
    x_alt = N - math.ceil(N * fractional_mismatches)

    dist_range = range(0, n + 1)
    alt_dist = binom.pmf(dist_range, n, x_alt / N)
    null_dist = binom.pmf(dist_range, n, x_null / N)
    """

    # plot for viewing pleasure
    plt.scatter(dist_range, alt_dist, label='alt', marker='o', color='b')
    plt.scatter(dist_range, null_dist, label='null', marker='x', color='r')
    plt.show()

    alt_tail = sum(alt_dist[k:])
    null_tail = sum(null_dist[k:])

    pvalue = null_tail / alt_tail

    print("pvalue:",pvalue)
    return pvalue
 

###############   TEST THINGS WOOHOO   #########################

k = 100
n = 100
N = 1000
alpha = .1
fractional_margin = .01
Vw = math.ceil(N * (1 + fractional_margin) / 2)
Vl = N - Vw
null_margin = 0

generate_comparison_distributions(k, n, N, alpha, Vw, Vl, null_margin)

