import numpy as np
from simulation_functions import simulate_fisher_combined_audits
from round_sizes import compute_dist_over_pvalues, find_sample_size_for_stopping_prob_efficiently, find_sample_size_for_stopping_prob_r2bravo
from fishers_combination import calculate_lambda_range
import math
import matplotlib.pyplot as plt
import time
import json

data = {}
data['audits'] = []

alpha = 0.1

# overall numbers (same as suite example 1)
N_w = 45500 + 7500
N_l = 49500 + 1500
N_relevant = N_w + N_l

n1 = 750 # same for all tests, same as in suite example

stopping_probability = .9

round_size = find_sample_size_for_stopping_prob_r2bravo(stopping_probability, N_w, N_l, alpha, underlying=None)
print ("round_size: "+str(round_size))

data['audits'].append({
    'N_relevant':N_relevant,
    'N_w':N_w,
    'N_l':N_l,
    'round_size':round_size
})

#update the file each time (hopefully will write over?
with open('data_r2bravo.txt', 'w') as outfile:
    json.dump(data, outfile, indent=2)



