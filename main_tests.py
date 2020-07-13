import matplotlib
import matplotlib.pyplot as plt

import math
import numpy as np
import numpy.random
import scipy as sp
import scipy.stats
from scipy.stats import binom

from contest import ContestType
from contest import Contest

from ballot_comparison import findNmin_ballot_comparison_rates
from hypergeometric import trihypergeometric_optim, simulate_ballot_polling_power
from fishers_combination import calculate_lambda_range
from simulation_functions import simulate_fisher_combined_audits

alpha = 0.1

N_w1 = 45500
N_l1 = 49500
N_w2 = 7500
N_l2 = 1500
N1 = N_w1 + N_w2
N2 = N_w2 + N_w2
margin = (N_w1 + N_w2 - N_l1 - N_l2)

np.random.seed(4976828)

n1 = 500
n2 = 200

reps = 1000

results = simulate_fisher_combined_audits(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, reps=reps, feasible_lambda_range=calculate_lambda_range(N_w1, N_l1, N1, N_w2, N_l2, N2), underlying=math.floor((N_w2 + N_l2) / 2))

print(results)



