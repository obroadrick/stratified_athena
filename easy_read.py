import json
import pprint
import matplotlib.pyplot as plt

with open('data_with_r2bravo_3.txt') as json_file:
    data = json.load(json_file)
    #pprint.pprint(data)

with open('data_minerva.txt') as json_file:
    minerva_data = json.load(json_file)
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

for i,j,k in zip(minerva_round_sizes, r2bravo_round_sizes, percent_pollings):
    print(str(k)+": "+str(i)+" "+str(j))

