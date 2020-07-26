import numpy as np
from simulation_functions import simulate_fisher_combined_audits
from fishers_combination import calculate_lambda_range
import math
from round_sizes import compute

alpha = 0.1

# overall numbers (same as suite example 1)
N_w = 45500 + 7500
N_l = 49500 + 1500
N_relevant = N_w + N_l
print ("N_relevant: "+str(N_relevant))
N_w_fraction = N_w / N_relevant

# division into strata 
# (testing multiple stratum proportions, all with same N_w_fraction fraction of winner votes)
percent_polling = .1
N_2 = math.ceil(percent_polling * N_relevant) # arbitrary whether this is ceil or floor (should just look up python round to nearest whole number)
print ("N_2: "+str(N_2))
N_1 = N_relevant - N_2
print ("N_1: "+str(N_1))
N_w1 = math.ceil(N_w_fraction * N_1) # arbitrary whether this is ceil or floor (should just look up python round to nearest whole number)
print ("N_w1: "+str(N_w1))
N_l1 = N_1 - N_w1
print ("N_l1: "+str(N_l1))
N_w2 = N_w - N_w1
print ("N_w2: "+str(N_w2))
N_l2 = N_2 - N_w2
assert (N_l2 + N_l1 == N_l) # sanity check (can remove after first successful run)
assert (N_w2 + N_w1 == N_w) # sanity check (can remove after first successful run)
assert (N_1 + N_2 == N_relevant) # sanity check (can remove after first successful run)
print ("N_l2: "+str(N_l2))
print ("N_w: "+str(N_w))
print ("N_w: "+str(N_w1+N_w2))
print ("N_l: "+str(N_l))
print ("N_l: "+str(N_l1+N_l2))
margin = N_w1 + N_w2 - N_l1 - N_l2

np.random.seed(18124328)

n1 = 750 # same for all tests, same as in suite example
n2 = 620

reps = 3000

results = compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None)

         "possible_winner_votes":possible_winner_votes
        "dist_over_winner_votes":dist_over_winner_votes
        "pvalues":pvalues
 
fig = plt.figure(figsize=(20,10))
fig.suptitle('round size vs polling stratum size (as percentages)', fontsize=20)
ax = fig.add_subplot(111)#numrows, numcols, num of subplot being referenced
ax.scatter(percent_pollings, minerva_round_sizes, color='b', marker='o', label='Minerva')
ax.scatter(percent_pollings, r2bravo_round_sizes, color='r', marker='x', label='R2 Bravo')
ax.set_xlabel('polling stratum as percent of relevant ballots', fontsize=20)
ax.set_ylabel('round size', fontsize=20)
plt.legend(loc='upper left', fontsize=20)
plt.setp(ax.get_xticklabels(), fontsize=18)
plt.setp(ax.get_yticklabels(), fontsize=18)
plt.show()



