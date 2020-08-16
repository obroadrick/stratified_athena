"""
Functions for a new method of stratified audit with 2 strata (ballot
comparison and ballot polling) that finds a sharper p-value by 
analytically finding the overall pvalue, rather than combining
independent p-values.

Oliver Broadrick 2020
"""

import math
import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt
from matplotlib import cm
import time

def generate_comparison_dists(n, Vw, Vl, null_margin=0, plot=False):
    """
    Generates the first round probability distributions over number 
    of matches for both the null and alternative hypothesis of 
    a comparison audit.

    For simplicity (this is just proof of concept) I'm considering 
    all overstatement errors 2-vote overstatements

    Parameters:
        n : sample size
        N : total ballots
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
    if b_null < 0:
        # when b_null is less than 0, understatement errors are required to make this null_margin possible
        # so for this case of only considering 2-vote overstatements, let's just say no errors
        b_null = 0
    #print("b_null:",b_null)

    # maximum number of matches under the null
    #   x = N - b
    x_null = N - b_null
    #print("x_null:",x_null)

    # risk limiting prior (0's, .5 at x_null, then uniform)
    if x_null == N:
        prior = np.concatenate([np.zeros(x_null), np.array([1])])
    else:
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
            alt_dist[k] += binom.pmf(k, n, x / N) * prior[x]
    # normalize
    alt_dist = alt_dist / sum(alt_dist)

    """
    # plot for viewing pleasure
    plt.scatter(dist_range, alt_dist, label='alt dist (comparison)', marker='o', color='b')
    plt.show()
    """

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
 
def generate_polling_dists(n, Vw, Vl, null_margin=0, plot=False):
    """
    Generates first round distributions over winner votes for 
    ballot polling audit.

    Parameters:
        k : number of votes for the winner in the sample (rest for loser)
        N : total ballots
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

def generate_joint_dist(N_w1, N_l1, N_w2, N_l2, n1, n2, null_margin1=0, plot=False):
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
        null_margin1 : the margin in votes assumed under the null in the comparison stratum
                        (defaults to zero)
        plot : Optional: if true, will plot the distributions (default false)

    Returns:
        alt_joint : joint probability distribution under the alternative hypothesis
        null_joint : joint probability distribution under the null hypothesis
    """

    # generate comparison dists
    comparison_dists = generate_comparison_dists(n1, N_w1, N_l1, null_margin1)
    dist_range_comparison = comparison_dists['dist_range']
    alt_dist_comparison = comparison_dists['alt_dist']
    null_dist_comparison = comparison_dists['null_dist']

    # generate polling dists
    polling_dists = generate_polling_dists(n2, N_w2, N_l2, -1 * null_margin1)
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

        surf2 =axes.plot_surface(X, Y, null_joint, label='joint dist under null', cmap=cm.seismic, linewidth=0, antialiased=False)
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

def compute_winner_vote_bounds(N_w1, N_l1, N_w2, N_l2, k_c=0, k_p=0):
    """
    Computes upper and lower bounds for the number of winner votes in 
    the comparison stratum, x1. 

    Parameters:
        N_w1 : winner votes in comparison stratum
        N_l1 : loser votes in comparison stratum
        N_w2 : winner votes in polling stratum
        N_lw : loser votes in polling stratum
        k_c : matches in comparison sample
            Optional: defaults to zero
        k_p : winner votes in polling sample
            Optional: defaults to zero

    Return: dict with
        x1_l : lower bound on x1
        x1_u : upper bound on x1
    """
    
    N1 = N_w1 + N_l1
    N2 = N_w2 + N_l2
    N = N1 + N2

    winner_votes_null = math.floor(N / 2)

    x1_l = max(0, winner_votes_null - N2, k_c)
    x1_u = min(N1 - k_p, winner_votes_null)

    return {
        'x1_l' : x1_l,
        'x1_u' : x1_u
    }

def maximize_joint_pvalue(k_c, k_p, N_w1, N_l1, N_w2, N_l2, n1, n2, lower_bound=None, upper_bound=None):
    """
    Maximizes the joint pvalue for the given sample by searching
    over all possible allocations of winner votes under the null
    hypothesis.

    Parameters:
        k_c : matches in comparison sample
        k_p : winner votes in polling sample
        N_w1 : reported winner votes in comparison stratum
        N_l1 : reported loser votes in comparison stratum
        N_w2 : reported winner votes in polling stratum
        N_lw : reported loser votes in polling stratum
        n1 : comparison sample size (first round)
        n2 : polling sample size (first round)
        x1_l : lower bound for winner votes in comparison stratum under the null
        x1_u : upper bound for winner votes in comparison stratum under the null

    Return: dict with
        pvalue : maximum pvalue found
    """
    
    if lower_bound is None or upper_bound is None:
        # compute bounds for search
        bounds = compute_winner_vote_bounds(N_w1, N_l1, N_w2, N_l2, k_c, k_p)
        lower_bound = bounds['x1_l']
        upper_bound = bounds['x1_u']

    """
    print("lower_bound:",lower_bound)
    print("upper_bound:",upper_bound)
    """

    # define a step size and generate lists for testing
    divisions = 50 # could tweak this for effiency
    step_size = math.ceil((upper_bound - lower_bound) / divisions)
    test_xs = np.arange(lower_bound, upper_bound, step_size)
    pvalues = np.empty_like(test_xs, dtype=float)

    for i in range(len(test_xs)):
        # compute null margin based on winner votes
        x1 = test_xs[i]
        null_margin1 = x1 - (N_w1 + N_l1 - x1)

        # generate joint distributions
        dists = generate_joint_dist(N_w1, N_l1, N_w2, N_l2, n1, n2, null_margin1)
        alt_joint = dists['alt_joint']
        null_joint = dists['null_joint']

        # compute pvalue
        pvalue = compute_pvalue_for_dist(k_c, k_p, alt_joint, null_joint)
        pvalues[i] = pvalue
   
        """
        # print for viewing pleasure
        print("N_w1:",N_w1)
        print("N_l1:",N_l1)
        print("N_w2:",N_w2)
        print("N_l2:",N_l2)
        print("x1:",x1)
        print("null_margin1:",null_margin1)
        print("pvalue:",pvalue)
        print("i+1:",i+1,"of",len(test_xs))
        """


    # get the maximum pvalue found
    max_index = np.argmax(pvalues)
    max_pvalue = pvalues[max_index]
    max_x = test_xs[max_index]
    """
    print("max_index:",max_index)
    print("max_pvalue:",max_pvalue)
    print("max_x:",max_x)
    """

    # if step_size has reached 1, search is over
    if step_size == 1:
        return {
            'pvalue' : max_pvalue,
            'x1' : max_x,
            'search_iterations' : 1
        }
    else:
        # set bounds for refined search
        x1_l = max_x - 2 * step_size
        x1_u = max_x + 2 * step_size

        # perform refined search
        refined = maximize_joint_pvalue(k_c, k_p, N_w1, N_l1, N_w2, N_l2, n1, n2, x1_l, x1_u)

        # increase iterations by 1 and return results
        refined['search_iterations'] += 1
        return refined
   
###############   TEST THINGS WOOHOO   #########################

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
k_c = 90
k_p = 125

joint_dist_results = generate_joint_dist(N_w1, N_l1, N_w2, N_l2, n1, n2, null_margin1, plot)
alt_joint = joint_dist_results['alt_joint']
null_joint = joint_dist_results['null_joint']

kmin_results = find_kmin_pairs(alpha, joint_dist_results['alt_joint'], joint_dist_results['null_joint'])
k_mins_c = kmin_results['k_mins_c']
k_mins_p = kmin_results['k_mins_p']

print(k_mins_c)
print(k_mins_p)

plot_joint_dist_with_kmins(alt_joint, null_joint, k_mins_c, k_mins_p)
"""

"""
N_w1 = 60
N_l1 = 40
N_w2 = 60
N_l2 = 40
n1 = 10
n2 = 15
alpha = .1
k_c = 9
k_p = 12


results = maximize_joint_pvalue(k_c, k_p, N_w1, N_l1, N_w2, N_l2, n1, n2, lower_bound=None, upper_bound=None)
print(results)
"""


"""
LETS TRY TO REPLICATE A TRIAL AND BEAT IT!
      "percent_polling": 0.05,
      "N_relevant": 104000,
      "N_w": 53000,
      "N_l": 51000,
      "N_2": 5200,
      "N_1": 98800,
      "N_w1": 50350,
      "N_l1": 48450,
      "N_w2": 2650,
      "N_l2": 2550,
      "minerva_round_size": 84,
      "minerva_combined_pvalue": 0.08636737929003713,
      "minerva_comparison_pvalue": 0.2410385629811352,
      "minerva_polling_pvalue": 0.07063014160664928,
      "minerva_alloc_lambda": 0.19456450000000156,
      "r2bravo_round_size": 208,
      "r2bravo_combined_pvalue": 0.08866865456278672,
      "r2bravo_comparison_pvalue": 0.01884495212609485,
      "r2bravo_polling_pvalue": 0.933506058896115,
      "r2bravo_alloc_lambda": 0.5421737000000015
#AND ANOTHER
      "percent_polling": 0.1,
      "N_relevant": 104000,
      "N_w": 53000,
      "N_l": 51000,
      "N_2": 10400,
      "N_1": 93600,
      "N_w1": 47700,
      "N_l1": 45900,
      "N_w2": 5300,
      "N_l2": 5100,
      "minerva_round_size": 314,
      "minerva_combined_pvalue": 0.09476826805247196,
      "minerva_comparison_pvalue": 0.5156528394976762,
      "minerva_polling_pvalue": 0.03707432824218867,
      "minerva_alloc_lambda": 0.0858475999999846,
      "r2bravo_round_size": 740,
      "r2bravo_combined_pvalue": 0.09216942784585869,
      "r2bravo_comparison_pvalue": 0.02481337562856968,
      "r2bravo_polling_pvalue": 0.7440985264790503,
      "r2bravo_alloc_lambda": 0.4781420999999825
 
"""
N_w1 = 47700
N_l1 = 45900
N_w2 = 5300
N_l2 = 5100
n1 = 750
n2 = 314
alpha = .1
k_c = n1
k_p = 171

start = time.time()

results = maximize_joint_pvalue(k_c, k_p, N_w1, N_l1, N_w2, N_l2, n1, n2, lower_bound=None, upper_bound=None)


print(results)


print("time:",(time.time()-start)/60,"minutes")













