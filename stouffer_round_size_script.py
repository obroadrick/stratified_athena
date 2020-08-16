"""
Rather than using Fisher's method for pvalue combination as in SUITE, 
this script uses Stouffer's method. Exploring the method and the weights
used.
This script generates computes and stores in a json file the polling stratum
round sizes that achieve a desired probability of stopping under the alternative
hypothesis that the election is truly as announced. It produces both
Minerva and R2 Bravo round sizes for various polling stratum sizes (as a 
proportion of relevant ballots).

Oliver Broadrick 2020
"""

from round_sizes import find_sample_size_for_stopping_prob_efficiently, find_sample_size_for_stopping_prob_efficiently_r2bravo
import math
import matplotlib.pyplot as plt
import json
from fisher_stouffer import wtd_stouffer
import numpy as np

# track all data in a struct, output to a json file
data = {}
data['audits'] = []

# file name for data 
data_file_name = 'data/data_stouffer_1.txt'

# risk limit (same as suite example 1)
alpha = 0.1

# overall relevant ballot tallies
fractional_margin = .02
N_relevant = 104000
N_w = round(N_relevant * (1 + fractional_margin) / 2)
N_l = 104000 - N_w
N_w_fraction = N_w / N_relevant

# divide the ballots into strata
# all strata will have the same margin as the overall contest
# define the proportion of relevant ballots in the polling stratum
proportion_polling = .5
N_2 = round(proportion_polling * N_relevant)
N_1 = N_relevant - N_2
N_w1 = round(N_w_fraction * N_1)
N_l1 = N_1 - N_w1
N_w2 = N_w - N_w1
N_l2 = N_2 - N_w2
margin = N_w1 + N_w2 - N_l1 - N_l2

# sanity check
assert (N_l2 + N_l1 == N_l)
assert (N_w2 + N_w1 == N_w)
assert (N_1 + N_2 == N_relevant)

# print for viewing pleasure
print("\nfractional_margin: "+str(fractional_margin))
print("proportion_polling: "+str(proportion_polling))
print ("N_relevant: "+str(N_relevant))
print ("N_w: "+str(N_w))
print ("N_l: "+str(N_l))
print ("N_1: "+str(N_1))
print ("N_w1: "+str(N_w1))
print ("N_l1: "+str(N_l1))
print ("N_2: "+str(N_2))
print ("N_w2: "+str(N_w2))
print ("N_l2: "+str(N_l2))

# comparison stratum round size (fixed)
n1 = 750

# desired probability of stopping under the alternative hypothesis that 
#   the contest truly is as announced
stopping_probability = .9

# define right bounds for round size search (defaults to stratum size)
minerva_right = None
r2bravo_right = None

weight = .05
while weight < 1:
    weights = np.array([weight, 1 - weight])
    print("\nweights:",weights)

    # define a stouffer function that uses the weights defined above
    def stouffer(pvalues):
        return wtd_stouffer(pvalues, weights)

    # obtain and print the minerva round size along with pvalues and lambda
    minerva_results = find_sample_size_for_stopping_prob_efficiently(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None, right=minerva_right, combine_func=stouffer, stouffers=True)
    print ("minerva_round_size: "+str(minerva_results['round_size']))
    print("combined_pvalue: "+str(minerva_results['combined_pvalue']))
    print("comparison pvalue: "+str(minerva_results['comparison_pvalue']))
    print("polling pvalue: "+str(minerva_results['polling_pvalue']))
    print("alloc_lambda: "+str(minerva_results['alloc_lambda']))

    # obtain and print the r2bravo round size along with pvalues and lambda
    r2bravo_results = find_sample_size_for_stopping_prob_efficiently_r2bravo(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None, right=r2bravo_right, combine_func=wtd_stouffer, stouffers=True)
    print ("r2bravo_round_size: "+str(r2bravo_results['round_size']))
    print("combined_pvalue: "+str(r2bravo_results['combined_pvalue']))
    print("comparison pvalue: "+str(r2bravo_results['comparison_pvalue']))
    print("polling pvalue: "+str(r2bravo_results['polling_pvalue']))
    print("alloc_lambda: "+str(r2bravo_results['alloc_lambda']))

    # add this data to the data structure
    data['audits'].append({
        'weight1':weight,
        'weight2':1 - weight,
        'proportion_polling':proportion_polling,
        'N_relevant':N_relevant,
        'N_w':N_w,
        'N_l':N_l,
        'N_2':N_2,
        'N_1':N_1,
        'N_w1':N_w1,
        'N_l1':N_l1,
        'N_w2':N_w2,
        'N_l2':N_l2,
        'minerva_round_size':minerva_results['round_size'],
        'minerva_combined_pvalue':minerva_results['combined_pvalue'],
        'minerva_comparison_pvalue':minerva_results['comparison_pvalue'],
        'minerva_polling_pvalue':minerva_results['polling_pvalue'],
        'minerva_alloc_lambda':minerva_results['alloc_lambda'],
        'r2bravo_round_size':r2bravo_results['round_size'],
        'r2bravo_combined_pvalue':r2bravo_results['combined_pvalue'],
        'r2bravo_comparison_pvalue':r2bravo_results['comparison_pvalue'],
        'r2bravo_polling_pvalue':r2bravo_results['polling_pvalue'],
        'r2bravo_alloc_lambda':r2bravo_results['alloc_lambda']
    })

    # update the file each loop (for convenience of checking progress)
    with open(data_file_name, 'w') as outfile:
        json.dump(data, outfile, indent=2)

    # increase to next weight
    weight += .05



