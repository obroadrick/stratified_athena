"""
This script is used to generate a plot of polling stratum round size 
over various comparison stratum sample sizes. This should give a sense
for the number of comparisons necessary to achieve the desired efficiency,
after which point more comparisons are unneccessary for sufficiently large
polling stratum sizes.

Oliver Broadrick 2020
"""

import json
import pprint
import matplotlib.pyplot as plt

# open desired stratified data json file
with open('data/data_stratified_vary_comparisons.txt') as json_file:
    data = json.load(json_file)
    #pprint.pprint(data)

# go through the stratified file and get desired data
n1s = []
minerva_round_sizes = []
r2bravo_round_sizes = []
for audit in data['audits']:
    # margin and percent polling same for all...
    overall_margin = round((audit['N_w1'] + audit['N_w2'] - audit['N_l1'] - audit['N_l2']) / audit['N_relevant'], 4)
    percent_polling = audit['percent_polling']

    # get and add to lists the number of comparisons and the polling stratum
    #   round sizes for both of Minerva and R2 Bravo
    n1 = audit['n1']
    minerva_round_size = audit['minerva_round_size']
    r2bravo_round_size = audit['r2bravo_round_size']

    n1s.append(n1)
    minerva_round_sizes.append(minerva_round_size)
    r2bravo_round_sizes.append(r2bravo_round_size)

# make a pretty plot of the data
fig = plt.figure(figsize=(20,10))
fig.suptitle('Round Sizes for Varying Comparison Sample Sizes (for margin: '+str(overall_margin)+', and percent polling: '+str(percent_polling)+')', fontsize=20)
ax = fig.add_subplot(111)#numrows, numcols, num of subplot being referenced
ax.scatter(n1s, minerva_round_sizes, color='b', marker='o', label='Minerva')
ax.scatter(n1s, r2bravo_round_sizes, color='r', marker='x', label='R2 Bravo')
ax.set_xlabel('Comparison Stratum Sample Size (number of comparisons)', fontsize=20)
ax.set_ylabel('Polling Stratum First Round Size (90% stopping probability)', fontsize=20)
plt.legend(loc='upper right', fontsize=20)
plt.setp(ax.get_xticklabels(), fontsize=18)
plt.setp(ax.get_yticklabels(), fontsize=18)
plt.show()


