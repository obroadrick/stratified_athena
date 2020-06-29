import matplotlib
import matplotlib.pyplot as plt

import math
import numpy as np
import numpy.random
import scipy as sp
import scipy.stats

from contest import ContestType
from contest import Contest

from ballot_comparison import findNmin_ballot_comparison_rates
from hypergeometric import trihypergeometric_optim, simulate_ballot_polling_power
from fishers_combination import simulate_fisher_combined_audit_minerva, simulate_fisher_combined_audit, calculate_lambda_range

alpha = 0.1

N_w1 = 45500
N_l1 = 49500
N_w2 = 7500
N_l2 = 1500
N1 = 100000
N2 = 10000
margin = (N_w1 + N_w2 - N_l1 - N_l2)

np.random.seed(20280584)

n1 = 500
n2 = 180

minerva_ish_result = simulate_fisher_combined_audit_minerva(N_w1, N_l1, N1, N_w2, N_l2, N2, n1, n2, alpha, reps=1000, feasible_lambda_range=calculate_lambda_range(N_w1, N_l1, N1, N_w2, N_l2, N2))

arlo_ish_result = simulate_fisher_combined_audit(N_w1, N_l1, N1, N_w2, N_l2, N2, n1, n2, alpha, reps=1000, feasible_lambda_range=calculate_lambda_range(N_w1, N_l1, N1, N_w2, N_l2, N2))

print("percent of simulations that stopped in arlo: "+str(arlo_ish_result * 100)+"%")
print("percent of simulations that stopped in minerva: "+str(minerva_ish_result * 100)+"%")
print("percent improvement: "+str(100 - arlo_ish_result / minerva_ish_result * 100)+"%")

candidates = ["Alice","Bob"]
winners = ["Alice"]
losers = ["Bob"]
stratum_sizes = [N1,N2]









