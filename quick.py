import math
from scipy.stats import binom

N_w1 = 50350
N_l1 = 48450
N_w2 = 2650
N_l2 = 2550
n1 = 750
n2 = 84
alpha = .1

N_w1 = 47700
N_l1 = 45900
N_w2 = 5300
N_l2 = 5100
n1 = 750
n2 = 314
alpha = .1
k_c = n1


q = .9
n = n2
p = N_w2 / (N_w2 + N_l2)
print(binom.ppf(q, n, p))



