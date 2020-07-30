"""
This script is used to find the number of ballots required for the first round of
both a Minerva and R2 Bravo audit of the contest. These numbers are used for
comparison to the stratified numbers.

Oliver Broadrick 2020
"""

from round_sizes import find_sample_size_for_stopping_prob_minerva, find_sample_size_for_stopping_prob_r2bravo
import json

data = {}
data['audits'] = []

alpha = 0.1

# overall numbers (same as suite example 1)
N_w = 64000
N_l = 40000
N_relevant = N_w + N_l

n1 = 750

stopping_probability = .9

minerva_round_size = find_sample_size_for_stopping_prob_minerva(stopping_probability, N_w, N_l, alpha, underlying=None)
print ("minerva_round_size: "+str(minerva_round_size))
r2bravo_round_size = find_sample_size_for_stopping_prob_r2bravo(stopping_probability, N_w, N_l, alpha, underlying=None)
print ("r2bravo_round_size: "+str(r2bravo_round_size))

data['audits'].append({
    'N_relevant':N_relevant,
    'N_w':N_w,
    'N_l':N_l,
    'minerva_round_size':minerva_round_size,
    'r2bravo_round_size':r2bravo_round_size
})

#update the file each time (hopefully will write over?
with open('data/data_contest_wide_larger_margin.txt', 'w') as outfile:
    json.dump(data, outfile, indent=2)

