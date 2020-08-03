"""
This script will plot stopping probability
for various round sizes by accessing the 
stored data.

Oliver Broadrick 2020
"""

import json
import matplotlib.pyplot as plt

# open desired stratified data json file
with open('data/data_stopping_prob_for_round_sizes.txt') as json_file:
    data = json.load(json_file)
    #pprint.pprint(data)

# acquire desired data from struct
n2s = data['n2s']
stopping_probs = data['stop_probs']
min_n2 = data['min_n2']
max_n2 = data['max_n2']

# plot it as well
fig = plt.figure(figsize=(20,10))
fig.suptitle('Stopping Probability for Round Sizes', fontsize=20)
ax = fig.add_subplot(111)#numrows, numcols, num of subplot being referenced
ax.plot(n2s, stopping_probs, linestyle='solid', color='b', marker='o', label='Minerva')
ax.plot([min_n2, max_n2], [.9, .9], linestyle='dashed', color='b', label='.9 probability of stopping')
ax.set_xlabel('Polling Stratum Round Sizes', fontsize=20)
ax.set_ylabel('Probability of Stopping', fontsize=20)
plt.legend(loc='upper left', fontsize=20)
plt.setp(ax.get_xticklabels(), fontsize=18)
plt.setp(ax.get_yticklabels(), fontsize=18)
plt.show()


