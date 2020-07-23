import json
import pprint
import matplotlib.pyplot as plt

with open('data_eff.txt') as json_file:
    data = json.load(json_file)
    pprint.pprint(data)


N_2s = []
round_sizes = []
round_size_as_percents = []
percent_pollings = []
for audit in data['audits']:
    percent_polling = audit['percent_polling']

    N_2 = audit['N_2']
    round_size = audit['round_size']
    round_size_as_percent = round_size / N_2

    N_2s.append(N_2)
    round_sizes.append(round_size)
    percent_pollings.append(percent_polling)
    round_size_as_percents.append(round_size_as_percent)

plt.plot(N_2s, round_sizes, 'bo')
plt.show()

fig = plt.figure()
plt.plot(percent_pollings, round_size_as_percents, 'bo')
fig.suptitle('round size vs polling stratum size (as percentages)', fontsize=20)
plt.xlabel('polling stratum as percent of relevant ballots', fontsize=20)
plt.ylabel('round size as percent of stratum size', fontsize=20)
plt.show()
    


