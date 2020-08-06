"""
For several states of interest (tight margins), this script generates computes 
and stores in a json file the polling stratum round sizes that achieve a .9 
probability of stopping under the alternative hypothesis that the election is 
truly as announced. It produces both Minerva and R2 Bravo round sizes for various 
polling stratum sizes (as a percentage of relevant ballots).

Oliver Broadrick 2020
"""

"""
States of interest:
    Rhode Island
    Wisconsin
    Michigan
    Pennsylvania
    Ohio
    North Carolina
"""

from round_sizes import find_sample_size_for_stopping_prob_efficiently, find_sample_size_for_stopping_prob_efficiently_r2bravo
import math
import matplotlib.pyplot as plt
import json
import pprint

# make an array of the states of interest
states_of_interest =   [
        'RhodeIsland',
        'Wisconsin',
        'Michigan',
        'Pennsylvania',
        'Ohio',
        'NorthCarolina'
    ]

# open desired stratified data json file
with open('data/2016_one_round_all.json') as json_file:
    data_2016 = json.load(json_file)
    #pprint.pprint(data)

# loop through the states
for state in data_2016:
    if state in states_of_interest:

        # print this states data
        print("\n", state, ":")
        pprint.pprint(data_2016[state])

        # track all data in a struct for each state, output to a json file
        data = data_2016[state]
        data.update({'audits': []})

        # file for data named based on margin
        data_file_name = 'data/round_sizes_'+state+'.txt'

        # risk limit of 10%
        alpha = 0.1

        # overall relevant ballot tallies from state data
        results = data_2016[state]['contests']['presidential']['results']
        N_w = max(results)
        N_l = min(results)
        N_relevant = N_w + N_l
        N_w_fraction = N_w / N_relevant
        fractional_margin = (N_w - N_l) / N_relevant
        print("N_w: "+str(N_w))
        print("N_l: "+str(N_l))
        print("N_relevant: "+str(N_relevant))
        print("fractional_margin(mine): "+str(fractional_margin))

        # loop through the various percent sizes of the polling stratum
        for percent_polling in [.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95]:

            # divide the ballots into strata
            # all strata will have the same margin as the overall contest
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
            print("percent_polling: "+str(percent_polling))
            print ("N_1: "+str(N_1)+"  N_w1: "+str(N_w1)+"  N_l1: "+str(N_l1))
            print ("N_2: "+str(N_2)+"  N_w2: "+str(N_w2)+"  N_l2: "+str(N_l2))
         
            # comparison stratum round size (fixed)
            n1 = 750

            # desired probability of stopping under the alternative hypothesis that 
            #   the contest truly is as announced
            stopping_probability = .9

            # define right bounds for round size search
            minerva_right = None
            r2bravo_right = None
            
            # obtain and print the minerva round size along with pvalues and lambda
            minerva_results = find_sample_size_for_stopping_prob_efficiently(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None, right=minerva_right)
            print ("minerva_round_size: "+str(minerva_results['round_size']))
            print("combined_pvalue: "+str(minerva_results['combined_pvalue']))
            print("comparison pvalue: "+str(minerva_results['comparison_pvalue']))
            print("polling pvalue: "+str(minerva_results['polling_pvalue']))
            print("alloc_lambda: "+str(minerva_results['alloc_lambda']))

            # obtain and print the r2bravo round size along with pvalues and lambda
            r2bravo_results = find_sample_size_for_stopping_prob_efficiently_r2bravo(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None, right=r2bravo_right)
            print ("r2bravo_round_size: "+str(r2bravo_results['round_size']))
            print("combined_pvalue: "+str(r2bravo_results['combined_pvalue']))
            print("comparison pvalue: "+str(r2bravo_results['comparison_pvalue']))
            print("polling pvalue: "+str(r2bravo_results['polling_pvalue']))
            print("alloc_lambda: "+str(r2bravo_results['alloc_lambda']))

            # add this data to the data structure
            data['audits'].append({
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
                'minerva_round_size':minerva_results['round_size'],
                'minerva_combined_pvalue':minerva_results['combined_pvalue'],
                'minerva_comparison_pvalue':minerva_results['comparison_pvalue'],
                'minerva_polling_pvalue':minerva_results['polling_pvalue'],
                'minerva_alloc_lambda':minerva_results['alloc_lambda'],
                'r2bravo_round_size':r2bravo_results['round_size'],
                'r2bravo_combined_pvalue':r2bravo_results['combined_pvalue'],
                'r2bravo_comparison_pvalue':r2bravo_results['comparison_pvalue'],
                'r2bravo_polling_pvalue':r2bravo_results['polling_pvalue'],
                'r2bravo_alloc_lambda':r2bravo_results['alloc_lambda']
            })

            # update the file each loop (for convenience of checking progress)
            with open(data_file_name, 'w') as outfile:
                json.dump(data, outfile, indent=2)



