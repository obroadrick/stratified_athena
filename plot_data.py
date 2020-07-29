import json
import pprint
import matplotlib.pyplot as plt

with open('data_with_r2bravo_3.txt') as json_file:
    data = json.load(json_file)
    #pprint.pprint(data)

with open('data_minerva.txt') as json_file:
    minerva_data = json.load(json_file)
    #pprint.pprint(data)

with open('data_r2bravo.txt') as json_file:
    r2bravo_data = json.load(json_file)
    #pprint.pprint(data)

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

minerva_round_size_no_stratification = minerva_data['audits'][0]['round_size']
r2bravo_round_size_no_stratification = r2bravo_data['audits'][0]['round_size']

fig = plt.figure(figsize=(20,10))
fig.suptitle('Round Sizes for Varying Polling Stratum Sizes (for margin: '+str(overall_margin)+')', fontsize=20)
ax = fig.add_subplot(111)#numrows, numcols, num of subplot being referenced
ax.scatter(percent_pollings, minerva_round_sizes, color='b', marker='o', label='Minerva')
ax.scatter(percent_pollings, r2bravo_round_sizes, color='r', marker='x', label='R2 Bravo')
ax.set_xlabel('Polling Stratum Size (as percent of relevant ballots)', fontsize=20)
ax.set_ylabel('First Round Size (90% stopping probability)', fontsize=20)
ax.plot([0,1],[minerva_round_size_no_stratification,minerva_round_size_no_stratification], label='Contest-Wide Minerva Audit', linestyle='dashed', color='b')
ax.plot([0,1],[r2bravo_round_size_no_stratification,r2bravo_round_size_no_stratification], label='Contest-Wide R2 Bravo Audit', linestyle='dashed', color='r')
plt.legend(loc='upper left', fontsize=20)
plt.setp(ax.get_xticklabels(), fontsize=18)
plt.setp(ax.get_yticklabels(), fontsize=18)
plt.show()


