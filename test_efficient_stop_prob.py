import numpy as np
from simulation_functions import simulate_fisher_combined_audits
from round_sizes import compute_dist_over_pvalues, find_sample_size_for_stopping_prob, find_sample_size_for_stopping_prob_efficiently, find_sample_size_for_stopping_prob_efficiently_r2bravo
from fishers_combination import calculate_lambda_range
import math
import matplotlib.pyplot as plt
import time
import json

data = {}
data['audits'] = []

for percent_polling in [.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65]:
    print("percent_polling: "+str(percent_polling))
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

    n1 = 750 # same for all tests, same as in suite example

    stopping_probability = .9

    minerva_round_size = find_sample_size_for_stopping_prob_efficiently(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None)
    print ("minerva_round_size: "+str(minerva_round_size))

    r2bravo_round_size = find_sample_size_for_stopping_prob_efficiently_r2bravo(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None)
    print ("r2bravo_round_size: "+str(r2bravo_round_size))
 
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
        'minerva_round_size':minerva_round_size,
        'r2bravo_round_size':r2bravo_round_size
    })

    #update the file each time (hopefully will write over?
    with open('data_with_r2bravo_backup.txt', 'w') as outfile:
        json.dump(data, outfile, indent=2)



