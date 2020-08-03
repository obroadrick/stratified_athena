"""
This script computes and outputs to json file stopping
probabilities for various polling stratum first round
sizes in a 2-strata audit with fixed comparison sample
size.

Oliver Broadrick 2020
"""

from round_sizes import compute_stopping_probability
import json

# risk limit (same as suite example 1)
alpha = 0.1

# overall relevant ballot tallies (same as suite example 1)
fractional_margin = .02
N_relevant = 104000
N_w = round(N_relevant * (1 + fractional_margin) / 2)
N_l = 104000 - N_w
N_w_fraction = N_w / N_relevant

# divide the ballots into strata
# all strata will have the same margin as the overall contest
percent_polling = .05
N_2 = round(percent_polling * N_relevant)
N_1 = N_relevant - N_2
N_w1 = round(N_w_fraction * N_1)
N_l1 = N_1 - N_w1
N_w2 = N_w - N_w1
N_l2 = N_2 - N_w2
margin = N_w1 + N_w2 - N_l1 - N_l2

# sanity check
assert (N_l2 + N_l1 == N_l)
assert (N_w2 + N_w1 == N_w)
assert (N_1 + N_2 == N_relevant)

# print for viewing pleasure
print("\nfractional_margin: "+str(fractional_margin))
print("percent_polling: "+str(percent_polling))
print ("N_relevant: "+str(N_relevant))
print ("N_w: "+str(N_w))
print ("N_l: "+str(N_l))
print ("N_1: "+str(N_1))
print ("N_w1: "+str(N_w1))
print ("N_l1: "+str(N_l1))
print ("N_2: "+str(N_2))
print ("N_w2: "+str(N_w2))
print ("N_l2: "+str(N_l2))

# fixed comparison stratum sample size
n1 = 750

# create desired n2 range
min_n2 = 0
max_n2 = 120
n2s = range(min_n2, max_n2 + 1)
stopping_probs = []

# create a structure for storing results
data = {
    'percent_polling':percent_polling,
    'N_relevant':N_relevant,
    'N_w':N_w,
    'N_l':N_l,
    'N_2':N_2,
    'N_1':N_1,
    'N_w1':N_w1,
    'N_l1':N_l1,
    'N_w2':N_w2,
    'N_l2':N_l2,
    'min_n2': min_n2,
    'max_n2': max_n2,
    'n2s':[],
    'stop_probs':[]
}

# loop through n2s and compute stopping probabilities
for n2 in n2s:
    # compute stopping prob for this n2 
    prob_stop = compute_stopping_probability(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None)
    stopping_probs.append(prob_stop)

    # print for viewing pleasure
    print('n2: '+str(n2)+'  prob_stop: '+str(prob_stop))

    # append to list in data struct
    data['n2s'].append(n2)
    data['stop_probs'].append(prob_stop)

# output the data struct to a file
with open('data/data_stopping_prob_for_round_sizes.txt', 'w') as outfile:
    json.dump(data, outfile, indent=2)

