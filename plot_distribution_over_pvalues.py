"""
This script plots a probability distribution over pvalues for given round size
in the ballot polling stratum of a two strata audit.
These numbers are an example of a round size that achieves a 90% probability of stopping
in an audit using Minerva for the pollign stratum.

Oliver Broadrick 2020
"""

from simulations import minerva_pvalue_direct_count
from round_sizes import compute_dist_over_pvalues
from scipy.stats import binom
import math
import matplotlib.pyplot as plt
from fisher_stouffer import wtd_stouffer
import numpy as np

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
weight = .5
weights = np.array([weight, 1 - weight])
print("\nweights:",weights)

# define a stouffer function that uses the weights defined above
def stouffer(pvalues):
    np.array(pvalues)
    return wtd_stouffer(pvalues, weights)
############STOUFFER##################

# compute distribution over pvalues 
results = compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None, combine_func=stouffer)
"""
Contents of results:
    "possible_winner_votes":possible_winner_votes,
    "dist_over_winner_votes":dist_over_winner_votes,
    "pvalues":pvalues
"""
# prob of pvalue less than alpha
prob_sum = 0
for pvalue, prob in zip(results['pvalues'], results['dist_over_winner_votes']):
    if (pvalue <= alpha):
        prob_sum += prob
print(prob_sum)

# put everything in a pretty plot
fig = plt.figure(figsize=(20,10))
axes = fig.add_subplot(111)
axes.scatter(results['pvalues'], results['dist_over_winner_votes'],linestyle="solid",color="blue")
#axes.plot([kmax,kmax], [0, 1], linestyle="dashed", label="kmax (.9 quantile to right)")
axes.plot([alpha,alpha],[0,.1], linestyle="dashed", label="alpha (Pr[pvalue<=alpha]="+str(round(prob_sum, 4))+")")
axes.set_xlabel('pvalues (combined)', fontsize=20)
axes.set_ylabel('probability (under alt)', fontsize=20)
plt.legend(loc='upper right', fontsize=20)
plt.setp(axes.get_xticklabels(), fontsize=18)
plt.setp(axes.get_yticklabels(), fontsize=18)
plt.show()

