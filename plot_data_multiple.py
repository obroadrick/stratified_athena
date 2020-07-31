"""
This script is used to plot the first round sizes for stratified audits using both
Minerva and R2 Bravo in the polling stratum for different polling stratum sizes.
It plots data not just from one test, but multiple. (For the sake of comparison.)
It also includes the first round sizes for contest-wide audits, both Minerva
and R2 Bravo. 

Oliver Broadrick 2020
"""

import json
import pprint
import matplotlib.pyplot as plt

# open desired stratified data json file for baseline
with open('data/data_stratified.txt') as json_file:
    data = json.load(json_file)
    #pprint.pprint(data)

# go through the stratified file and get desired data
N_2s = []
minerva_round_sizes = []
r2bravo_round_sizes = []
percent_pollings = []
for audit in data['audits']:
    overall_margin = round((audit['N_w1'] + audit['N_w2'] - audit['N_l1'] - audit['N_l2']) / audit['N_relevant'], 4)
    percent_polling = audit['percent_polling']

    N_2 = audit['N_2']
    minerva_round_size = audit['minerva_round_size']
    r2bravo_round_size = audit['r2bravo_round_size']

    N_2s.append(N_2)
    percent_pollings.append(percent_polling)
    minerva_round_sizes.append(minerva_round_size)
    r2bravo_round_sizes.append(r2bravo_round_size)

# open desired stratified data json file for comparison
with open('data/data_stratified_smaller_margin_TEST.txt') as json_file:
    data_2 = json.load(json_file)
    #pprint.pprint(data)

# go through the stratified file and get desired data
N_2s_2 = []
minerva_round_sizes_2 = []
r2bravo_round_sizes_2  = []
percent_pollings_2 = []
for audit in data_2['audits']:
    overall_margin_larger = round((audit['N_w1'] + audit['N_w2'] - audit['N_l1'] - audit['N_l2']) / audit['N_relevant'], 4)
    percent_polling = audit['percent_polling']

    N_2 = audit['N_2']
    minerva_round_size = audit['minerva_round_size']
    r2bravo_round_size = audit['r2bravo_round_size']

    N_2s_2.append(N_2)
    percent_pollings_2.append(percent_polling)
    minerva_round_sizes_2.append(minerva_round_size)
    r2bravo_round_sizes_2.append(r2bravo_round_size)

# open the contest-wide data
with open('data/data_contest_wide.txt') as json_file:
    contest_wide_data = json.load(json_file)
    #pprint.pprint(data)

# obtain desired contest-wide data
minerva_round_size_no_stratification = contest_wide_data['audits'][0]['minerva_round_size']
r2bravo_round_size_no_stratification = contest_wide_data['audits'][0]['r2bravo_round_size']

# open the contest-wide data for larger margin
with open('data/data_contest_wide_smaller_margin.txt') as json_file:
    contest_wide_data_larger_margin = json.load(json_file)
    #pprint.pprint(data)

# obtain desired contest-wide data for larger margin
minerva_round_size_no_stratification_larger_margin = contest_wide_data_larger_margin['audits'][0]['minerva_round_size']
r2bravo_round_size_no_stratification_larger_margin = contest_wide_data_larger_margin['audits'][0]['r2bravo_round_size']

# make a pretty plot of the data
fig = plt.figure(figsize=(20,10))
fig.suptitle('Round Sizes for Varying Polling Stratum Sizes (for margins: '+str(overall_margin)+', '+str(overall_margin_larger)+')', fontsize=20)
ax = fig.add_subplot(111)#numrows, numcols, num of subplot being referenced
ax.scatter(percent_pollings, minerva_round_sizes, color='b', marker='o', label='Minerva (.0192 margin)')
ax.scatter(percent_pollings, r2bravo_round_sizes, color='r', marker='x', label='R2 Bravo (.0192 margin)')
ax.scatter(percent_pollings_2, minerva_round_sizes_2, color='g', marker='o', label='Minerva (.0096 margin)')
ax.scatter(percent_pollings_2, r2bravo_round_sizes_2, color='y', marker='x', label='R2 Bravo (.0096 margin)')
ax.set_xlabel('Polling Stratum Size (as percent of relevant ballots)', fontsize=20)
ax.set_ylabel('First Round Size (90% stopping probability)', fontsize=20)
ax.plot([0,1],[minerva_round_size_no_stratification,minerva_round_size_no_stratification], label='Contest-Wide Minerva Audit (.0192 margin)', linestyle='dashed', color='b')
ax.plot([0,1],[r2bravo_round_size_no_stratification,r2bravo_round_size_no_stratification], label='Contest-Wide R2 Bravo Audit (.0192 margin)', linestyle='dashed', color='r')
ax.plot([0,1],[minerva_round_size_no_stratification_larger_margin,minerva_round_size_no_stratification_larger_margin], label='Contest-Wide Minerva Audit (.0096 margin)', linestyle='dashed', color='g')
ax.plot([0,1],[r2bravo_round_size_no_stratification_larger_margin,r2bravo_round_size_no_stratification_larger_margin], label='Contest-Wide R2 Bravo Audit (.0096 margin)', linestyle='dashed', color='y')
plt.legend(loc='upper left', fontsize=20)
plt.setp(ax.get_xticklabels(), fontsize=18)
plt.setp(ax.get_yticklabels(), fontsize=18)
plt.show()


