"""
Short script for testing Stouffer weighted combination.

Oliver Broadrick 2020
"""

import math
import matplotlib.pyplot as plt
import numpy as np
from fisher_stouffer import wtd_stouffer

# throw stouffer function at it 
weight = .5
weights = np.array([weight, 1 - weight])

# define stouffer functions for these weights
def stouffer(pvalues):
    np.array(pvalues)
    return wtd_stouffer(pvalues, weights)

# combine a bunch
pvals = np.linspace(0,1,100)
stouffers = []
for pval in pvals:
    stouffers.append(stouffer([pval,pval]))


# put everything in a plot for viewing pleasure
fig = plt.figure(figsize=(20,10))
axes = fig.add_subplot(111)
axes.scatter(pvals,stouffers,label="stouffers combined of two of same pval from x axis")
plt.legend(loc='upper right', fontsize=20)
plt.setp(axes.get_xticklabels(), fontsize=18)
plt.setp(axes.get_yticklabels(), fontsize=18)
plt.show()

