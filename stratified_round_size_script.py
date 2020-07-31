from round_sizes import find_sample_size_for_stopping_prob_efficiently, find_sample_size_for_stopping_prob_efficiently_r2bravo
import math
import matplotlib.pyplot as plt
import json

# track all data in one structure, output to a json file
data = {}
data['audits'] = []

# loop through the various percent sizes of the polling stratum
for percent_polling in [.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95]:
    # define a right bounds for the round size searchs based on the previous size found
    if (percent_polling == .05):
        minerva_right = None
        r2bravo_right = None
    elif (percent_polling < .25):
        # times 4 is sufficient here (I tested it)
        minerva_prev = data['audits'][len(data['audits']) - 1]['minerva_round_size']
        r2bravo_prev = data['audits'][len(data['audits']) - 1]['r2bravo_round_size']
        minerva_right = 4 * minerva_prev
        print(minerva_right)
        r2bravo_right = 4 * r2bravo_prev
        print(r2bravo_right)
    else:
        # times 2 is sufficient here (I tested it)
        minerva_prev = data['audits'][len(data['audits']) - 1]['minerva_round_size']
        r2bravo_prev = data['audits'][len(data['audits']) - 1]['r2bravo_round_size']
        minerva_right = 2 * minerva_prev
        print(minerva_right)
        r2bravo_right = 2 * r2bravo_prev
        print(r2bravo_right)

    # risk limit (same as suite example 1)
    alpha = 0.1

    # overall relevant ballot tallies (same as suite example 1)
    N_w = 52500
    N_l = 104000 - N_w
    N_relevant = N_w + N_l
    N_w_fraction = N_w / N_relevant

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
    print("\npercent_polling: "+str(percent_polling))
    print ("N_relevant: "+str(N_relevant))
    print ("N_w: "+str(N_w))
    print ("N_l: "+str(N_l))
    print ("N_1: "+str(N_1))
    print ("N_w1: "+str(N_w1))
    print ("N_l1: "+str(N_l1))
    print ("N_2: "+str(N_2))
    print ("N_w2: "+str(N_w2))
    print ("N_l2: "+str(N_l2))
 
    # comparison stratum round size (fixed)
    n1 = 1500

    # desired probability of stopping under the alternative hypothesis that 
    #   the contest truly is as announced
    stopping_probability = .9

    
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
    with open('data/data_stratified_smaller_margin_TEST.txt', 'w') as outfile:
        json.dump(data, outfile, indent=2)



