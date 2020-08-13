"""
This script plots an example of a round size that achieves a 90% probability of stopping
in a stratified audit using Minerva for the pollign stratum.
Includes the pvalues and probabilities for all possible values of k and includes
a marker of the value of kmax.
Now also includes both Fisher and Stouffer combining function numbers.

Oliver Broadrick 2020
"""

from simulations import minerva_pvalue_direct_count
from round_sizes import compute_dist_over_pvalues
from scipy.stats import binom
import math
import matplotlib.pyplot as plt
import numpy as np
from fisher_stouffer import wtd_stouffer

# these numbers correspond to a 53000/51000 contest with 
#   5% of the relevant ballots in the polling stratum
N_w1 = 50350
N_l1 = 48450
N_1 = N_w1 + N_l1
N_w2 = 2650
N_l2 = 2550
N_2 = N_w2 + N_l2
n1 = 750
n2 = 84
alpha = .1
stopping_probability = .9
kmax = math.floor(binom.ppf(1 - stopping_probability, n2, N_w2 / N_2))

############STOUFFER##################
# throw stouffer function at it 
weight1 = .5
weights1 = np.array([weight1, 1 - weight1])
weight2 = .4
weights2 = np.array([weight2, 1 - weight2])
# define stouffer functions for these weights
def stouffer1(pvalues):
    np.array(pvalues)
    return wtd_stouffer(pvalues, weights1)
def stouffer2(pvalues):
    np.array(pvalues)
    return wtd_stouffer(pvalues, weights2)
############STOUFFER##################

# compute distribution over pvalues for stouffer and fisher
fisher_results = compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None)
stouffer_results1 = compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None, combine_func=stouffer1)
#stouffer_results2 = compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None, combine_func=stouffer2)
"""
Contents of results:
    "possible_winner_votes":possible_winner_votes,
    "dist_over_winner_votes":dist_over_winner_votes,
    "pvalues":pvalues
"""

# put everything in a pretty plot
fig = plt.figure(figsize=(20,10))
axes = fig.add_subplot(111)
axes.scatter(fisher_results['possible_winner_votes'], fisher_results['pvalues'],linestyle="solid",color="blue", label="pvalue (fishers)")
axes.scatter(fisher_results['possible_winner_votes'], fisher_results['alloc_lambda'],color="c", label="lambda (fishers)", marker="x")
axes.scatter(fisher_results['possible_winner_votes'], fisher_results['dist_over_winner_votes'],linestyle="solid",color="red", label="probability of k under alt")
axes.scatter(stouffer_results1['possible_winner_votes'], stouffer_results1['pvalues'],linestyle="solid",color="g", label="pvalue (stouffer "+str(weights1)+")")
axes.scatter(stouffer_results1['possible_winner_votes'], stouffer_results1['alloc_lambda'],linestyle="solid",color="y", label="lambda (stouffer "+str(weights1)+")", marker="x")
#axes.scatter(stouffer_results2['possible_winner_votes'], stouffer_results2['pvalues'],linestyle="solid",color="c", label="pvalue (stouffer "+str(weights2)+")")
axes.plot([kmax,kmax], [0, 1], linestyle="dashed", label="kmax="+str(kmax)+" (.9 quantile to right)")
axes.plot([0,n2],[alpha,alpha], linestyle="dashed", label="alpha")
plt.legend(loc='upper right', fontsize=20)
plt.setp(axes.get_xticklabels(), fontsize=18)
plt.setp(axes.get_yticklabels(), fontsize=18)
plt.show()

