"""
This script plots an example of a round size that achieves a 90% probability of stopping
in a stratified audit using Minerva for the pollign stratum.
Includes the pvalues and probabilities for all possible values of k and includes
a marker of the value of kmax.

Oliver Broadrick 2020
"""

from simulations import minerva_pvalue_direct_count
from round_sizes import compute_dist_over_pvalues
from scipy.stats import binom
import math
import matplotlib.pyplot as plt

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

# compute distribution over pvalues 
results = compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None)
"""
Contents of results:
    "possible_winner_votes":possible_winner_votes,
    "dist_over_winner_votes":dist_over_winner_votes,
    "pvalues":pvalues
"""

# put everything in a pretty plot
fig = plt.figure(figsize=(20,10))
axes = fig.add_subplot(111)
axes.scatter(results['possible_winner_votes'], results['pvalues'],linestyle="solid",color="blue", label="pvalue produced by this k")
axes.scatter(results['possible_winner_votes'], results['dist_over_winner_votes'],linestyle="solid",color="red", label="probability of k under alt")
axes.plot([kmax,kmax], [0, 1], linestyle="dashed", label="kmax (.9 quantile to right)")
axes.plot([0,n2],[alpha,alpha], linestyle="dashed", label="alpha")
plt.legend(loc='upper right', fontsize=20)
plt.setp(axes.get_xticklabels(), fontsize=18)
plt.setp(axes.get_yticklabels(), fontsize=18)
plt.show()

