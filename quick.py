import numpy as np
from simulation_functions import simulate_fisher_combined_audits
from fishers_combination import calculate_lambda_range
import math

alpha = 0.1

# overall numbers (same as suite example 1)
N_w1 = 45500
N_w2 = 7500
N_l1 = 49500
N_l2 = 1500
N_1 = N_w1 + N_l1
N_2 = N_w2 + N_l2
margin = N_w1 + N_w2 - N_l1 - N_l2

np.random.seed(41782714)

n1 = 750 # same for all tests, same as in suite example
n2 = 200

reps = 100

results = simulate_fisher_combined_audits(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, reps=reps, verbose=False, feasible_lambda_range=calculate_lambda_range(N_w1, N_l1, N_1, N_w2, N_l2, N_2), underlying=None)
 
print(results)



