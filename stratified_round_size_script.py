import numpy as np
from simulation_functions import simulate_fisher_combined_audits
from round_sizes import compute_dist_over_pvalues, find_sample_size_for_stopping_prob_efficiently, find_sample_size_for_stopping_prob_efficiently_r2bravo
from fishers_combination import calculate_lambda_range
import math
import matplotlib.pyplot as plt
import time
import json

data = {}
data['audits'] = []

for percent_polling in [.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95]:
    print("\npercent_polling: "+str(percent_polling))
    alpha = 0.1

    # overall numbers (same as suite example 1)
    N_w = 45500 + 7500
    N_l = 49500 + 1500
    N_relevant = N_w + N_l
    print ("N_relevant: "+str(N_relevant))
    N_w_fraction = N_w / N_relevant

    # division into strata 
    # (testing multiple stratum proportions, all with same N_w_fraction fraction of winner votes)
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

    n1 = 1500 # same for all tests, same as in suite example

    stopping_probability = .9

    minerva_results = find_sample_size_for_stopping_prob_efficiently(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None)
    print ("minerva_round_size: "+str(minerva_results['round_size']))
    print("combined_pvalue: "+str(minerva_results['combined_pvalue']))
    print("comparison pvalue: "+str(minerva_results['comparison_pvalue']))
    print("polling pvalue: "+str(minerva_results['polling_pvalue']))
    print("alloc_lambda: "+str(minerva_results['alloc_lambda']))

    r2bravo_results = find_sample_size_for_stopping_prob_efficiently_r2bravo(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None)
    print ("r2bravo_round_size: "+str(r2bravo_results['round_size']))
    print("combined_pvalue: "+str(r2bravo_results['combined_pvalue']))
    print("comparison pvalue: "+str(r2bravo_results['comparison_pvalue']))
    print("polling pvalue: "+str(r2bravo_results['polling_pvalue']))
    print("alloc_lambda: "+str(r2bravo_results['alloc_lambda']))


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

    #update the file each time (hopefully will write over?
    with open('data/data_double_comparisons.txt', 'w') as outfile:
        json.dump(data, outfile, indent=2)



