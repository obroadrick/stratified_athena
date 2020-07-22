import numpy as np
from simulation_functions import simulate_fisher_combined_audits, compute_dist_over_pvalues
from fishers_combination import calculate_lambda_range
import math
import matplotlib.pyplot as plt

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
n2 = 310

results = compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None)

possible_winner_votes = results["possible_winner_votes"]
dist_over_winner_votes = results["dist_over_winner_votes"]
pvalues = results["pvalues"]

index = -1
# find the index of the first pvalue less than the risk limit
for i,pvalue in zip(range(0, n2 + 1), pvalues):
    if (pvalue <= alpha):
        index = i
        break
prob_stop = sum(dist_over_winner_votes[index:])

print("probability of getting less than .1 (alpha) pvalue " + str(prob_stop))

plt.plot(possible_winner_votes,pvalues)
plt.show()
plt.plot(pvalues,dist_over_winner_votes)
plt.show()


