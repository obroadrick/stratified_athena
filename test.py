# quick script to test something
# not sticking around long term
"""
Need to see how many times greater each round size is compared to previous round
size in the usual data ive been getting (for defining right bounds for 
efficiency/ not having too low a bound
"""


import json
import pprint
import matplotlib.pyplot as plt

with open('data/data_stratified_smaller_margin.txt') as json_file:
    data = json.load(json_file)

minerva_round_size = 1
r2bravo_round_size = 1
minerva_round_size_prev = 1
r2bravo_round_size_prev = 1

for audit in data['audits']:
    minerva_round_size = audit['minerva_round_size']
    r2bravo_round_size = audit['r2bravo_round_size']

    print(round(minerva_round_size / minerva_round_size_prev,4))
    print(round(r2bravo_round_size / r2bravo_round_size_prev,4))

    minerva_round_size_prev = minerva_round_size
    r2bravo_round_size_prev = r2bravo_round_size

