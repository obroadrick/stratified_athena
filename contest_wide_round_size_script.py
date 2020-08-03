"""
This script is used to find the number of ballots required for the first round of
both a Minerva and R2 Bravo audit of the contest. These numbers are used for
comparison to the stratified numbers.

Oliver Broadrick 2020
"""

from round_sizes import find_sample_size_for_stopping_prob_minerva, find_sample_size_for_stopping_prob_r2bravo
import json

alpha = 0.1

# loop through the various margins for which I would like results (run overnight)
for margin in [.01, .02, .03, .04, .05, .06, .07, .08]:
    # create new data struct for each margin
    data = {}
    data['audits'] = []

    # file for data named based on margin
    data_file_name = 'data/data_contest_wide_'+str(round(margin * 100))+'_percent_margin.txt'

    # tallies (based on margin)
    N_relevant = 104000
    N_w = round(N_relevant * (1 + margin) / 2)
    N_l = 104000 - N_w
    N_w_fraction = N_w / N_relevant
 
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
    with open(data_file_name, 'w') as outfile:
        json.dump(data, outfile, indent=2)

