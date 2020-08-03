"""
This script is used to generate the polling round size for a 50%-50% 2-strata
audit for various comparisons stratum sample sizes.

Oliver Broadrick 2020
"""

from round_sizes import find_sample_size_for_stopping_prob_efficiently, find_sample_size_for_stopping_prob_efficiently_r2bravo
import math
import matplotlib.pyplot as plt
import json

# track all data in one structure, output to a json file
data = {}
data['audits'] = []

# loop through the various comparison stratum sample sizes
for n1 in range(100, 1001, 10):
    # risk limit
    alpha = 0.1

    # overall relevant ballot tallies
    N_w = 53000
    N_l = 104000 - N_w
    N_relevant = N_w + N_l
    N_w_fraction = N_w / N_relevant

    # divide the ballots into strata
    # all strata will have the same margin as the overall contest
    percent_polling = .5
    N_2 = round(percent_polling * N_relevant)
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
    print("\npercent_polling: "+str(percent_polling))
    print ("N_relevant: "+str(N_relevant))
    print ("N_w: "+str(N_w))
    print ("N_l: "+str(N_l))
    print ("N_1: "+str(N_1))
    print ("N_w1: "+str(N_w1))
    print ("N_l1: "+str(N_l1))
    print ("N_2: "+str(N_2))
    print ("N_w2: "+str(N_w2))
    print ("N_l2: "+str(N_l2))
    print ("n1: "+str(n1))
 
    # desired probability of stopping under the alternative hypothesis that 
    #   the contest truly is as announced
    stopping_probability = .9

    # obtain and print the minerva round size along with pvalues and lambda
    minerva_results = find_sample_size_for_stopping_prob_efficiently(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None, right = 10*N_2)
    print ("minerva_round_size: "+str(minerva_results['round_size']))
    print("combined_pvalue: "+str(minerva_results['combined_pvalue']))
    print("comparison pvalue: "+str(minerva_results['comparison_pvalue']))
    print("polling pvalue: "+str(minerva_results['polling_pvalue']))
    print("alloc_lambda: "+str(minerva_results['alloc_lambda']))

    # obtain and print the r2bravo round size along with pvalues and lambda
    r2bravo_results = find_sample_size_for_stopping_prob_efficiently_r2bravo(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None, right = 10*N_2)
    print ("r2bravo_round_size: "+str(r2bravo_results['round_size']))
    print("combined_pvalue: "+str(r2bravo_results['combined_pvalue']))
    print("comparison pvalue: "+str(r2bravo_results['comparison_pvalue']))
    print("polling pvalue: "+str(r2bravo_results['polling_pvalue']))
    print("alloc_lambda: "+str(r2bravo_results['alloc_lambda']))

    # add this data to the data structure
    data['audits'].append({
        'percent_polling':percent_polling,
        'N_relevant':N_relevant,
        'N_w':N_w,
        'N_l':N_l,
        'N_2':N_2,
        'N_1':N_1,
        'N_w1':N_w1,
        'N_l1':N_l1,
        'N_w2':N_w2,
        'N_l2':N_l2,
        'n1':n1,
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
    with open('data/data_stratified_vary_comparisons_TEST_RIGHTS.txt', 'w') as outfile:
        json.dump(data, outfile, indent=2)



