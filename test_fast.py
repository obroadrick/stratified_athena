"""
Quickly confirm that a binomial (and hypergeometric)
distribution wtih = .99 means no errors allowed

Oliver Broadrick 2020
"""

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binom



n = 100
p = .99
dist = []
ks = range(0,101)
for k in ks:
    dist.append(binom.pmf(k, n, p))

plt.scatter(ks, dist)
plt.show()


