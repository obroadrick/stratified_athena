import time
import numpy as np
import scipy as sp
import scipy.stats
import scipy.optimize
from ballot_comparison import ballot_comparison_pvalue
from sprt import ballot_polling_sprt
import matplotlib.pyplot as plt
import numpy.testing
from contest import ContestType
from contest import Contest
from minerva_s import Minerva_S
from fishers_combination import create_modulus, maximize_fisher_combined_pvalue, calculate_lambda_range
from scipy.stats import binom
import math
import matplotlib.pyplot as plt
from simulations import minerva_pvalue_direct_count, r2bravo_pvalue_direct_count

def compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None):
    """
    Hopefully compute a distribution over possible pvalues
    """

    N_1 = N_w1 + N_l1
    N_2 = N_w2 + N_l2
    margin = N_w1 + N_w2 - N_l1 - N_l2
    
    feasible_lambda_range=calculate_lambda_range(N_w1, N_l1, N_1, N_w2, N_l2, N_2)

    possible_winner_votes = range(0, n2 + 1)
    dist_over_winner_votes = binom.pmf(possible_winner_votes, n2, N_w2 / N_2)

    pvalues = []

    for k, pr_k in zip(possible_winner_votes, dist_over_winner_votes):

        cvr_pvalue = lambda alloc: ballot_comparison_pvalue(n=n1, gamma=1.03905, \
                                       o1=0, u1=0, o2=0, u2=0,
                                   reported_margin=margin, N=N_1,
                                   null_lambda=alloc)

        mod = create_modulus(n1, n2, k, n2 - k, N_1, margin, 1.03905)

        nocvr_pvalue = lambda alloc: \
            minerva_pvalue_direct_count(winner_votes=k, n=n2, popsize=N_2, alpha=alpha, \
                                Vw=N_w2, Vl=N_l2, \
                                null_margin=(N_w2-N_l2) - alloc*margin)

        pvalue = maximize_fisher_combined_pvalue(N_w1, N_l1, \
                               N_1, N_w2, N_l2, N_2, \
                               pvalue_funs=[cvr_pvalue, nocvr_pvalue], \
                               modulus=mod, alpha=alpha, \
                               feasible_lambda_range=feasible_lambda_range)['max_pvalue']

        pvalues.append(pvalue)
        #print("for k="+str(k)+"   pval="+str(pvalue))

    return {
        "possible_winner_votes":possible_winner_votes,
        "dist_over_winner_votes":dist_over_winner_votes,
        "pvalues":pvalues
    }
        

def compute_stopping_probability(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None):
    """
    Compute the stopping probability for the given sample sizes.
    Computes the full probability distribution over pvalues.
    """

    results = compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None)

    possible_winner_votes = results["possible_winner_votes"]
    dist_over_winner_votes = results["dist_over_winner_votes"]
    pvalues = results["pvalues"]

    # find the index of the first pvalue that passes the stopping condition
    for i,pvalue in zip(range(0, n2 + 1), pvalues):
        if (pvalue <= alpha):
            index = i
            break
    prob_stop = sum(dist_over_winner_votes[index:])

    return prob_stop

def find_sample_size_for_stopping_prob_efficiently(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None):
    """
    This function will also compute minimum round size for the 
    passed stopping probability, but it will do so much more 
    efficiently. At each point in the search only one pvalue 
    will be computed. Should have done it this way to begin with.

    Uses Minerva for the ballot polling stratum.
    """

    N_1 = N_w1 + N_l1
    N_2 = N_w2 + N_l2
    margin = N_w1 + N_w2 - N_l1 - N_l2

    feasible_lambda_range = calculate_lambda_range(N_w1, N_l1, N_1, N_w2, N_l2, N_2)

    left = 1
    right = round(N_2 / 2)
     
    while(1):
        n2 = math.ceil((left + right) / 2 )

        # compute the 1 - stopping_probability quantile of the alt dist
        # kmax where pr[k >= kmax | alt] = stopping_probability
        # floor because we need to ensure at least a stopping_probability prob of stopping
        kmax = math.floor(binom.ppf(1 - stopping_probability, n2, N_w2 / N_2))

        # compute pvalue for this kmax
        cvr_pvalue = lambda alloc: ballot_comparison_pvalue(n=n1, gamma=1.03905, \
                                       o1=0, u1=0, o2=0, u2=0,
                                   reported_margin=margin, N=N_1,
                                   null_lambda=alloc)

        mod = create_modulus(n1, n2, kmax, n2 - kmax, N_1, margin, 1.03905)

        nocvr_pvalue = lambda alloc: \
            minerva_pvalue_direct_count(winner_votes=kmax, n=n2, popsize=N_2, alpha=alpha, \
                                Vw=N_w2, Vl=N_l2, \
                                null_margin=(N_w2-N_l2) - alloc*margin)

        combination_results = maximize_fisher_combined_pvalue(N_w1, N_l1, \
                               N_1, N_w2, N_l2, N_2, \
                               pvalue_funs=[cvr_pvalue, nocvr_pvalue], \
                               modulus=mod, alpha=alpha, \
                               feasible_lambda_range=feasible_lambda_range)
        pvalue = combination_results['max_pvalue']
        pvalue_comparison = combination_results['pvalue1']
        pvalue_polling = combination_results['pvalue2']
        alloc_lambda = combination_results['allocation lambda']
       
        # update binary search bounds
        if (pvalue > alpha):
            left = n2
        elif (pvalue <= alpha):
            right = n2
 
        # when and right converge, right is the minimum round size that achieves stopping_probability
        if (left == right - 1 and n2 == right):
            if (right == round(N_2 / 3)):
                print("more than one third of the stratum votes are required")
            print(combination_results['refined'])
            return  {
                        "round_size":right,
                        "combined_pvalue":pvalue,
                        "comparison_pvalue":pvalue_comparison,
                        "polling_pvalue":pvalue_polling,
                        "alloc_lambda":alloc_lambda
                    }


def find_sample_size_for_stopping_prob_minerva(stopping_probability, N_w, N_l, alpha, underlying=None):
    """
    Finds the first round size that achieves the passed stopping_probability
    for a Minerva audit (with no stratification). 
    """

    N = N_w + N_l 

    left = 1
    right = round(N / 2)
     
    while(1):
        n = math.ceil((left + right) / 2)

        # compute the 1 - stopping_probability quantile of the alt dist
        # kmax where pr[k >= kmax | alt] = stopping_probability
        # floor because we need to ensure at least a stopping_probability prob of stopping
        kmax = math.floor(binom.ppf(1 - stopping_probability, n, N_w / N))

        # compute pvalue for this kmax
        pvalue = minerva_pvalue_direct_count(winner_votes=kmax, n=n, popsize=N, alpha=alpha, Vw=N_w, Vl=N_l, null_margin=0)

        # update binary search bounds
        if (pvalue > alpha):
            left = n
        elif (pvalue <= alpha):
            right = n
 
        # when and right converge, right is the minimum round size that achieves stopping_probability
        if (left == right - 1):
            if (right == round(N / 2)):
                print("more than one third of the stratum votes are required")
            return right

def find_sample_size_for_stopping_prob_r2bravo(stopping_probability, N_w, N_l, alpha, underlying=None):
    """
    Finds the first round size that achieves the passed stopping_probability
    for an R2 Bravo audit (with no stratification). 
    """

    N = N_w + N_l 

    left = 1
    right = round(N / 2)
     
    while(1):
        n = math.ceil((left + right) / 2)

        # compute the 1 - stopping_probability quantile of the alt dist
        # kmax where pr[k >= kmax | alt] = stopping_probability
        # floor because we need to ensure at least a stopping_probability prob of stopping
        kmax = math.floor(binom.ppf(1 - stopping_probability, n, N_w / N))

        # compute pvalue for this kmax
        pvalue = r2bravo_pvalue_direct_count(winner_votes=kmax, n=n, popsize=N, alpha=alpha, Vw=N_w, Vl=N_l, null_margin=0)

        # update binary search bounds
        if (pvalue > alpha):
            left = n
        elif (pvalue <= alpha):
            right = n
 
        # when and right converge, right is the minimum round size that achieves stopping_probability
        if (left == right - 1):
            if (right == round(N / 2)):
                print("more than one third of the stratum votes are required")
            return right


def find_sample_size_for_stopping_prob_efficiently_r2bravo(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None):
    """
    This function will also compute minimum round size for the 
    passed stopping probability, but it will do so much more 
    efficiently. At each point in the search only one pvalue 
    will be computed. Should have done it this way to begin with.
    """

    N_1 = N_w1 + N_l1
    N_2 = N_w2 + N_l2
    margin = N_w1 + N_w2 - N_l1 - N_l2

    feasible_lambda_range=calculate_lambda_range(N_w1, N_l1, N_1, N_w2, N_l2, N_2)

    left = 1
    right = round(N_2 / 2)
     
    while(1):
        n2 = math.ceil((left + right) / 2)

        # compute the 1 - stopping_probability quantile of the alt dist
        # kmax where pr[k >= kmax | alt] = stopping_probability
        # floor because we need to ensure at least a stopping_probability prob of stopping
        kmax = math.floor(binom.ppf(1 - stopping_probability, n2, N_w2 / N_2))

        # compute pvalue for this kmax
        cvr_pvalue = lambda alloc: ballot_comparison_pvalue(n=n1, gamma=1.03905, \
                                       o1=0, u1=0, o2=0, u2=0,
                                   reported_margin=margin, N=N_1,
                                   null_lambda=alloc)

        mod = create_modulus(n1, n2, kmax, n2 - kmax, N_1, margin, 1.03905)

        nocvr_pvalue = lambda alloc: \
            r2bravo_pvalue_direct_count(winner_votes=kmax, n=n2, popsize=N_2, alpha=alpha, \
                                Vw=N_w2, Vl=N_l2, \
                                null_margin=(N_w2-N_l2) - alloc*margin)

        combination_results = maximize_fisher_combined_pvalue(N_w1, N_l1, \
                               N_1, N_w2, N_l2, N_2, \
                               pvalue_funs=[cvr_pvalue, nocvr_pvalue], \
                               modulus=mod, alpha=alpha, \
                               feasible_lambda_range=feasible_lambda_range)

        pvalue = combination_results['max_pvalue']
        pvalue_comparison = combination_results['pvalue1']
        pvalue_polling = combination_results['pvalue2']
        alloc_lambda = combination_results['allocation lambda']

        
        # update binary search bounds
        if (pvalue > alpha):
            left = n2
        elif (pvalue <= alpha):
            right = n2
 
        # when and right converge, right is the minimum round size that achieves stopping_probability
        if (left == right - 1 and n2 == right):
            if (right == round(N_2 / 3)):
                print("more than one third of the stratum votes are required")
            return  {
                        "round_size":right,
                        "combined_pvalue":pvalue,
                        "comparison_pvalue":pvalue_comparison,
                        "polling_pvalue":pvalue_polling,
                        "alloc_lambda":alloc_lambda
                    }



