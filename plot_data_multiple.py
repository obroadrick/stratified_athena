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

# begin plot
fig = plt.figure(figsize=(20,10))
fig.suptitle('Round Sizes for Varying Polling Stratum Sizes (for margins: 2% to 6%)', fontsize=20)
ax = fig.add_subplot(111)#numrows, numcols, num of subplot being referenced

# array of data files to be plotted
file_name_start = 'data/data_stratified_'
file_name_end = '_percent_margin.txt'
file_name_start_contest_wide = 'data/data_contest_wide_'
file_name_end_contest_wide = '_percent_margin.txt'

# loop through margins
for margin, color in zip(range(2,7), ['r','b','g','c','m']):
    # construct file names
    file_name = file_name_start + str(margin) + file_name_end
    file_name_contest_wide = file_name_start_contest_wide + str(margin) + file_name_end_contest_wide

    # open desired stratified data json file for baseline
    with open(file_name) as json_file:
        data = json.load(json_file)

    # go through the stratified file and get desired data
    N_2s = []
    minerva_round_sizes = []
    r2bravo_round_sizes = []
    percent_pollings = []
    for audit in data['audits']:

        # gather data
        overall_margin = round((audit['N_w1'] + audit['N_w2'] - audit['N_l1'] - audit['N_l2']) / audit['N_relevant'], 4)
        percent_polling = audit['percent_polling']
        N_2 = audit['N_2']
        minerva_round_size = audit['minerva_round_size']
        r2bravo_round_size = audit['r2bravo_round_size']
        N_2s.append(N_2)
        percent_pollings.append(percent_polling)
        minerva_round_sizes.append(minerva_round_size)
        r2bravo_round_sizes.append(r2bravo_round_size)
        
    # plot data
    #ax.scatter(percent_pollings, minerva_round_sizes, color=color, marker='o', label='Stratified Minerva ('+str(margin)+'% margin)')
    ax.scatter(percent_pollings, r2bravo_round_sizes, color=color, marker='x', label='Stratified R2 Bravo ('+str(margin)+'% margin)')

    # open the contest-wide data
    with open(file_name_contest_wide) as json_file:
        contest_wide_data = json.load(json_file)

    # obtain desired contest-wide data
    minerva_round_size_no_stratification = contest_wide_data['audits'][0]['minerva_round_size']
    r2bravo_round_size_no_stratification = contest_wide_data['audits'][0]['r2bravo_round_size']

    # plot contest-wide lines
    #ax.plot([0,1],[minerva_round_size_no_stratification,minerva_round_size_no_stratification], label='Contest-Wide Minerva ('+str(margin)+'% margin)', linestyle='dashed',color=color)
    ax.plot([0,1],[r2bravo_round_size_no_stratification,r2bravo_round_size_no_stratification], label='Contest-Wide R2 Bravo ('+str(margin)+'% margin)', linestyle='dashed',color=color)

# add labels and adjust font sizes, then show
ax.set_xlabel('Polling Stratum Size (as proportion of relevant ballots)', fontsize=20)
ax.set_ylabel('First Round Size (90% stopping probability)', fontsize=20)
plt.legend(loc='upper left', fontsize=20)
plt.setp(ax.get_xticklabels(), fontsize=18)
plt.setp(ax.get_yticklabels(), fontsize=18)
plt.show()


