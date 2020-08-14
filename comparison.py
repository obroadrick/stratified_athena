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
from matplotlib import cm
import time

def generate_comparison_dists(n, alpha, Vw, Vl, null_margin=0, plot=False):
    """
    Generates the first round probability distributions over number 
    of matches for both the null and alternative hypothesis of 
    a comparison audit.

    For simplicity (this is just proof of concept) I'm considering 
    all overstatement errors 2-vote overstatements

    Parameters:
        n : sample size
        N : total ballots
        alpha : risk limit
        Vw : reported votes for winner
        Vl : reported votes for loser
        null_margin : the margin in votes under the null
                    Optional: defaults to 0 (used for stratified application)
        plot : Optional: if true, will plot the distributions (default false)

    Returns: dict with
        dist_range : values of k to which the distributions correspond
        alt_dist : probability distribution over matches under the alternative hypothesis
        null_dist : probability distribution over matches under the null hypothesis
    """
    # relevant ballots
    N = Vw + Vl

    # reported margin from passed tallies
    reported_margin = Vw - Vl
    
    # b is 2-vote overstatements under the null
    #   reported_margin - 2b <= null_margin
    b_null = math.floor((reported_margin - null_margin) / 2)
    #print("b_null:",b_null)

    # maximum number of matches under the null
    #   x = N - 2b
    x_null = N - 2 * b_null
    #print("x_null:",x_null)

    # risk limiting prior (0's, .5 at x_null, then uniform)
    prior = np.concatenate([np.zeros(x_null), np.array([.5]), np.full(N - x_null, .5 / (N - x_null))])

    """
    # plot for viewing pleasure
    plt.scatter(range(N+1), prior, label='prior (comparison)', marker='o', color='b')
    plt.show()
    """

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

    if (plot):
        # plot for viewing pleasure
        plt.scatter(dist_range, alt_dist, label='alt', marker='o', color='b')
        plt.scatter(dist_range, null_dist, label='null', marker='x', color='r')
        plt.show()

    """
    # CODE FOR COMPUTING PVALUES FOR THIS COMPARISON AUDIT (tail ratio) (FOR TESTING... IT WORKED)
    alt_tail = sum(alt_dist[k:])
    null_tail = sum(null_dist[k:])

    pvalue = null_tail / alt_tail

    print("pvalue:",pvalue)
    return pvalue
    """

    return {
        'dist_range': dist_range,
        'alt_dist': alt_dist,
        'null_dist': null_dist
    }
 
def generate_polling_dists(n, alpha, Vw, Vl, null_margin=0, plot=False):
    """
    Generates first round distributions over winner votes for 
    ballot polling audit.

    Parameters:
        k : number of votes for the winner in the sample (rest for loser)
        N : total ballots
        alpha : risk limit
        Vw : reported votes for winner
        Vl : reported votes for loser
        null_margin : the margin in votes assumed under the null
                Optional: default to 0 (for stratified application)
        plot : Optional: if true, will plot the distributions (default false)

    Returns: dict with
        dist_range : values of k to which the distributions correspond
        alt_dist : probability distribution over winner votes under the alternative hypothesis
        null_dist : probability distribution over winner votes under the null hypothesis
    """
    # relevant ballots
    N = Vw + Vl

    # winner votes under the null
    x = (N + null_margin) / 2

    # don't bother with bayesian prior other than .5 and .5 at x and Vw
    # could easily extend to accept other priors, same as in comparison stratum
    dist_range = range(0, n + 1)
    alt_dist = binom.pmf(dist_range, n, Vw / N)
    null_dist = binom.pmf(dist_range, n, x / N)

    """
    # CODE FOR COMPUTING PVALUE OF THIS POLLING AUDIT (tail ratio)
    alt_tail = sum(alt_dist[k:])
    null_tail = sum(null_dist[k:])

    pvalue = null_tail / alt_tail

    return pvalue
    """

    # plot if option set true
    if (plot):
        # plot for viewing pleasure
        plt.scatter(dist_range, alt_dist, label='alt', marker='o', color='b')
        plt.scatter(dist_range, null_dist, label='null', marker='x', color='r')
        plt.show()

    return {
        'dist_range': dist_range,
        'alt_dist': alt_dist,
        'null_dist': null_dist
    }

def generate_joint_dist(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, null_margin1=0, plot=False):
    """
    Generates the joint probability distribution for the first round
    of a 2-strata (comparison and polling) audit.

    Parameters:
        N_w1 : winner votes in comparison stratum
        N_l1 : loser votes in comparison stratum
        N_w2 : winner votes in polling stratum
        N_lw : loser votes in polling stratum
        n1 : comparison sample size (first round)
        n2 : polling sample size (first round)
        alpha : risk limit
        null_margin1 : the margin in votes assumed under the null in the comparison stratum
                        (defaults to zero)
        plot : Optional: if true, will plot the distributions (default false)

    Returns:
        alt_joint : joint probability distribution under the alternative hypothesis
        null_joint : joint probability distribution under the null hypothesis
    """

    # generate comparison dists
    comparison_dists = generate_comparison_dists(n1, alpha, N_w1, N_l1, null_margin1)
    dist_range_comparison = comparison_dists['dist_range']
    alt_dist_comparison = comparison_dists['alt_dist']
    null_dist_comparison = comparison_dists['null_dist']

    # generate polling dists
    polling_dists = generate_polling_dists(n2, alpha, N_w2, N_l2, -1 * null_margin1)
    dist_range_polling = polling_dists['dist_range']
    alt_dist_polling = polling_dists['alt_dist']
    null_dist_polling = polling_dists['null_dist']

    # compute joint dists
    alt_joint = np.empty((len(dist_range_comparison), len(dist_range_polling)))
    null_joint = np.empty((len(dist_range_comparison), len(dist_range_polling)))
    for k_c in dist_range_comparison:
        for k_p in dist_range_polling:
            alt_joint[k_c][k_p] = alt_dist_comparison[k_c] * alt_dist_polling[k_p]
            null_joint[k_c][k_p] = null_dist_comparison[k_c] * null_dist_polling[k_p]

    if (plot):
        # plot for viewing pleasure
        fig = plt.figure()
        axes = fig.gca(projection='3d')

        X, Y = np.meshgrid(dist_range_polling, dist_range_comparison)

        surf1 = axes.plot_surface(X, Y, alt_joint, label='joint dist under alt', cmap=cm.coolwarm, linewidth=0, antialiased=False)
        # matplotlib error fixed with these two lines
        surf1._facecolors2d=surf1._facecolors3d
        surf1._edgecolors2d=surf1._edgecolors3d

        surf2 =axes.plot_surface(X, Y, null_joint, label='joint dist under null', cmap=cm.magma, linewidth=0, antialiased=False)
        # matplotlib error fixed with these two lines
        surf2._facecolors2d=surf2._facecolors3d
        surf2._edgecolors2d=surf2._edgecolors3d

        axes.legend(fontsize=20)

        axes.set_xlabel('Polling Stratum: Winner Ballots', fontsize=20, labelpad=24)
        axes.set_ylabel('Comparison Stratum: Matching Ballots', fontsize=20, labelpad=24)
        axes.set_zlabel('Probability', fontsize=20, labelpad=24)

        plt.setp(axes.get_xticklabels(), fontsize=18)
        plt.setp(axes.get_yticklabels(), fontsize=18)
        plt.setp(axes.get_zticklabels(), fontsize=18)

        plt.show()

    return {
        'alt_joint': alt_joint,
        'null_joint': null_joint
    }

def compute_pvalue_for_dist(k_c, k_p, alt_joint, null_joint):
    """
    Compute the pvalue for the passed joint distribution.

    Parameters:
        k_c: matches in comparison sample
        k_p: winner votes in polling sample
        alt_joint: joint
        alt_joint: the joint probability distribution under the alternative hypothesis
        null_joint: the joint probability distribution under the null hypothesis

    Returns:
        float: pvalue (ratio of corners)
    """
    # sum the corners of the joint dists
    alt_corner = 0
    for k in range(k_c, len(alt_joint)):
        alt_corner += sum(alt_joint[k][k_p:])
    null_corner = 0
    for k in range(k_c, len(null_joint)):
        null_corner += sum(null_joint[k][k_p:])

    # pvalue is ratio of the corners
    return null_corner / alt_corner

 

def find_kmin_pairs(alpha, alt_joint, null_joint):
    """
    Finds all valid kmin pairs for the passed joint distributions 
    and risk limit.

    Parameters:
        alpha: risk limit
        alt_joint: the joint probability distribution under the alternative hypothesis
        null_joint: the joint probability distribution under the null hypothesis

    Returns: dict with
        k_mins_c: values of kmin for comparison stratum
        k_mins_p: corresponding values of kmin for polling stratum
    """

    k_mins_c = []
    k_mins_p = []
    for k_c in range(len(alt_joint)):
        for k_p in range(len(alt_joint[0])):

            # test if this pair of k's meet the stopping condition (you know, the stopping condition that I'm making up)
            pvalue = compute_pvalue_for_dist(k_c, k_p, alt_joint, null_joint)

            if not math.isnan(pvalue) and pvalue <= alpha:
                # add these k mins to the lists
                k_mins_c.append(k_c)
                k_mins_p.append(k_p)

                """
                # print for viewing pleasure
                print("k_c:",k_c," k_p:",k_p)
                print("pvalue:",pvalue)
                """

                # break after first k_p min is found for this value of k_c
                break

    return {
        'k_mins_c': k_mins_c,
        'k_mins_p': k_mins_p
    }

def plot_joint_dist_with_kmins(alt_joint, null_joint, k_mins_c, k_mins_p):
    """
    Plots the joint probability distribution as well as the kmin line.

    Parameters:
        alt_joint : joint probability distribution under the alternative hypothesis
        null_joint : joint probability distribution under the null hypothesis
        k_mins_c: values of kmin for comparison stratum
        k_mins_p: corresponding values of kmin for polling stratum
    """
    fig = plt.figure()
    axes = fig.gca(projection='3d')

    dist_range_polling = range(len(alt_joint[0] + 1))
    dist_range_comparison = range(len(alt_joint + 1))
    X, Y = np.meshgrid(dist_range_polling, dist_range_comparison)

    surf1 = axes.plot_surface(X, Y, alt_joint, label='joint dist under alt', cmap=cm.coolwarm, linewidth=0, antialiased=False)
    # matplotlib error fixed with these two lines
    surf1._facecolors2d=surf1._facecolors3d
    surf1._edgecolors2d=surf1._edgecolors3d

    surf2 =axes.plot_surface(X, Y, null_joint, label='joint dist under null', cmap=cm.magma, linewidth=0, antialiased=False)
    # matplotlib error fixed with these two lines
    surf2._facecolors2d=surf2._facecolors3d
    surf2._edgecolors2d=surf2._edgecolors3d

    for k_min_c, k_min_p in zip(k_mins_c, k_mins_p):
        z = [0,.06]
        axes.plot([k_min_c,k_min_c], [k_min_p,k_min_p], z, linewidth=3, color='g')

    axes.legend(fontsize=20)

    axes.set_xlabel('Polling Stratum: Winner Ballots', fontsize=20, labelpad=24)
    axes.set_ylabel('Comparison Stratum: Matching Ballots', fontsize=20, labelpad=24)
    axes.set_zlabel('Probability', fontsize=20, labelpad=24)

    plt.setp(axes.get_xticklabels(), fontsize=18)
    plt.setp(axes.get_yticklabels(), fontsize=18)
    plt.setp(axes.get_zticklabels(), fontsize=18)

    plt.show()

def maximize_joint_pvalue(k_c, k_p, N_w1, N_l1, N_w2, N_l2, n1, n2, alpha):
    """
    Maximizes the joint pvalue for the given sample over all
    possible allocations of error.

    Parameters:
        k_c: matches in comparison sample
        k_p: winner votes in polling sample
        N_w1 : winner votes in comparison stratum
        N_l1 : loser votes in comparison stratum
        N_w2 : winner votes in polling stratum
        N_lw : loser votes in polling stratum
        n1 : comparison sample size (first round)
        n2 : polling sample size (first round)
        alpha : risk limit

    Return: dict with
        pvalue: maximum pvalue found
    """



 
###############   TEST THINGS WOOHOO   #########################

"""
k = 100
n = 100
N = 1000
alpha = .1
fractional_margin = .01
Vw = math.ceil(N * (1 + fractional_margin) / 2)
Vl = N - Vw
null_margin = 0

generate_comparison_dists(k, n, N, alpha, Vw, Vl, null_margin)
"""

N_w1 = 600
N_l1 = 400
N_w2 = 600
N_l2 = 400
n1 = 100
n2 = 150
alpha = .1
null_margin1 = 0
plot = False

joint_dist_results = generate_joint_dist(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, null_margin1, plot)
alt_joint = joint_dist_results['alt_joint']
null_joint = joint_dist_results['null_joint']

kmin_results = find_kmin_pairs(alpha, joint_dist_results['alt_joint'], joint_dist_results['null_joint'])
k_mins_c = kmin_results['k_mins_c']
k_mins_p = kmin_results['k_mins_p']

print(k_mins_c)
print(k_mins_p)

plot_joint_dist_with_kmins(alt_joint, null_joint, k_mins_c, k_mins_p)









