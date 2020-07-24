import json
import pprint
import matplotlib.pyplot as plt

with open('data_with_r2bravo.txt') as json_file:
    data = json.load(json_file)
    #pprint.pprint(data)


N_2s = []
minerva_round_sizes = []
r2bravo_round_sizes = []
percent_pollings = []
for audit in data['audits']:
    percent_polling = audit['percent_polling']

    N_2 = audit['N_2']
    minerva_round_size = audit['minerva_round_size']
    r2bravo_round_size = audit['r2bravo_round_size']

    N_2s.append(N_2)
    percent_pollings.append(percent_polling)
    minerva_round_sizes.append(minerva_round_size)
    r2bravo_round_sizes.append(r2bravo_round_size)

fig = plt.figure(figsize=(20,10))
fig.suptitle('round size vs polling stratum size (as percentages)', fontsize=20)
ax = fig.add_subplot(111)#numrows, numcols, num of subplot being referenced
ax.scatter(percent_pollings, minerva_round_sizes, color='b', marker='o', label='Minerva')
ax.scatter(percent_pollings, r2bravo_round_sizes, color='r', marker='x', label='R2 Bravo')
ax.set_xlabel('polling stratum as percent of relevant ballots', fontsize=20)
ax.set_ylabel('round size', fontsize=20)
plt.legend(loc='upper left', fontsize=20)
plt.setp(ax.get_xticklabels(), fontsize=18)
plt.setp(ax.get_yticklabels(), fontsize=18)
plt.show()


