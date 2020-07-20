import numpy as np
from simulation_functions import estimate_round_size_for_stopping_prob
from fishers_combination import calculate_lambda_range
import math

alpha = 0.1

N_w1 = 45500
N_l1 = 49500
N_w2 = 7500
N_l2 = 1500

np.random.seed(49768284)

prob = .9

results = estimate_round_size_for_stopping_prob(prob, N_w1, N_l1, N_w2, N_l2, alpha)

print(results)

