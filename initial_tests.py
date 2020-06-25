'''
Oliver Broadrick - Testing stratified Athena class audits.
'''
from audit_s import Audit_S
from contest import Contest
from contest import ContestType
from minerva_s import Minerva_S
import scipy
import numpy as np
import scipy.optimize as optimize
      
##############################################
### utility functions for testing purposes ###
##############################################

def fisher_combined_pvalue(pvalues):
    """Find the p-value for Fisher's combined test statistic

    Parameters:
        pvalues: Array of p-values to combine

    Returns:
        float: p-value for Fisher's combined test statistic
    """
    if np.any(np.array(pvalues)==0):
        return 0
    obs = -2*np.sum(np.log(pvalues))
    return 1-scipy.stats.chi2.cdf(obs, df=2*len(pvalues))

def find_stratum_pvalue(num_strata, alpha, tol: float = .001):
    """ Finds the required stratum-level p-value to 
    achieve the passed overall_pvalue in S strata.

    Parameters:
        tol: tolerance for overall_pvalue (default .001)

    Returns:
        float: the stratum-level p-value 
    """
    left = 0
    right = 1

    while (1):
        mid = (left + right) / 2
        left_combined = fisher_combined_pvalue(np.full(num_strata, left))
        right_combined = fisher_combined_pvalue(np.full(num_strata, right))
        mid_combined = fisher_combined_pvalue(np.full(num_strata, mid))
        if (mid_combined > alpha):
            right = mid
        elif (mid_combined < (alpha - tol)):
            left = mid
        else:
            return mid

def create_audits_for_lambda(strata, overall_contest, lambda_, stratum_alpha):
    """Create minerva audits with null margins according to the passed values for lambda.

    Returns:
        []: list of the audits for the strata.
    """
    stratum_audits = []
    winner = 'A'
    loser = 'B'
    overall_margin = overall_contest.tally[winner] - overall_contest.tally[loser]
    for stratum_contest, stratum_lambda in zip (strata, lambda_):
        reported_stratum_margin = stratum_contest.tally[winner] - stratum_contest.tally[loser]
        overstatement_quota = stratum_lambda * overall_margin
        # null margin is the margin in votes under the null hypothesis
        null_margin = reported_stratum_margin - overstatement_quota
        stratum_audit = Minerva_S(stratum_alpha, 1.0, stratum_contest, null_margin)
        stratum_audits.append(stratum_audit)
    return stratum_audits
 
def predict_value_of_lambda(strata,overall_contest):
    """Predict values of lambda by, yanno, totally guessing.

    Returns:
        []: values of lambda.
    """
    lambda_ = []
    i = 1
    winner = 'A'
    loser = 'B'
    overall_margin = overall_contest.tally[winner] - overall_contest.tally[loser]
    for stratum in strata:
        reported_stratum_margin = stratum.tally[winner] - stratum.tally[loser]
        stratum_lambda = reported_stratum_margin / overall_margin
        #print("stratum "+str(i)+" lambda value: "+str(stratum_lambda))
        lambda_.append(stratum_lambda)
        i += 1

    #print("sum of lambda values: "+str(sum(lambda_)))
    return lambda_

def overall_pvalue_for_lambda_and_sample(strata, overall_contest, sample, lambda_, alpha):
    """Computes the overall combined pvalue for the past strata, lambda, and sample.
    
    Parameters:
        sample: a list of sample dicts with sample size and votes for winner in each sample
    """
    stratum_alpha = find_stratum_pvalue(len(strata), alpha)
    strata_audits = create_audits_for_lambda(strata, overall_contest, lambda_, stratum_alpha)
    pvalues = []
    for stratum_sample, stratum_audit in zip(sample, strata_audits):
        stratum_audit.rounds.append(stratum_sample['sample_size'])
        stratum_audit.current_dist_null()
        stratum_audit.current_dist_reported()
        #stratum_audit.compute_min_winner_ballots([stratum_sample['sample_size']])
        #print("stratum kmin: " + str(stratum_audit.min_winner_ballots[0]))
        #print("winner votes: " + str(stratum_sample['winner_votes']))
        stratum_pvalue = stratum_audit.stopping_condition(stratum_sample['winner_votes'])['pvalue']
        pvalues.append(stratum_pvalue)
    print("\nresults for lambda="+str(lambda_))
    print("stratum pvalues: "+str(pvalues))
    combined = fisher_combined_pvalue(pvalues)
    print("combined pvalue: "+str(combined)+"\n")
    return combined

def overall_stopping_condition(strata, overall_contest, sample, alpha):
    """Maximizes pvalue over all possible values of lambda,
    and uses the maximum pvalue to test the overall stopping condition.

    Stopping condition: combined_pvalue < alpha
    """
    # this function can be used to mazimize the pvalue over all values of lambda
    # it will take as param all but one value of lambda 
    def overall_pvalue_for_lambda(lambda_):
        # last value of lambda is 1 - (sum of the others) so that they all sum to 1
        lambda_ = np.append(lambda_, 1 - sum(lambda_))
        return -1 * overall_pvalue_for_lambda_and_sample(strata, overall_contest, sample, lambda_, alpha)

    if (len(strata) > 1):

        # need to place bounds on the values of lambda
        

        initial_guess = np.full(len(strata) - 1, 1 / len(strata))

        result = optimize.minimize(overall_pvalue_for_lambda, initial_guess)
        lambda_ = np.append(result.x,1-sum(result.x))
        overall_pvalue = overall_pvalue_for_lambda_and_sample(strata, overall_contest, sample, lambda_, alpha)
        print("\nwinning lambda: "+str(lambda_))
        print("winning lambda combined pvalue: "+str(overall_pvalue)+"\n")
    else:
        lambda_ = [1]
        print("\nwinning lambda: "+str(lambda_))
        overall_pvalue = overall_pvalue_for_lambda_and_sample(strata, overall_contest, sample, lambda_, alpha)
        print("winning lambda combined pvalue: "+str(overall_pvalue)+"\n")
    return overall_pvalue < alpha







    
##############################################
#######      some actual testing!      #######
##############################################
"""
stratum1 = Contest(45500+49500, {
    'A': 45500,
    'B': 49500
}, 1, ['A'], ContestType.PLURALITY)

stratum2 = Contest(7500+1500, {
    'A': 7500,
    'B': 1500
}, 1, ['A'], ContestType.PLURALITY)

overall_contest = Contest(45500+49500+7500+1500, {
    'A': 45500+7500,
    'B': 49500+1500
}, 1, ['A'], ContestType.PLURALITY)
"""

stratum1 = Contest(100000, {
    'A': 60000,
    'B': 40000
}, 1, ['A'], ContestType.PLURALITY)

stratum2 = Contest(50000, {
    'A': 32000,
    'B': 18000
}, 1, ['A'], ContestType.PLURALITY)

overall_contest = Contest(100000+50000, {
    'A': 32000+60000,
    'B': 18000+40000
}, 1, ['A'], ContestType.PLURALITY)


strata = [stratum1, stratum2]

alpha = .05

sample = [  
            {'sample_size': 600, 'winner_votes': 360}, 
            {'sample_size': 300, 'winner_votes': 280}  
         ]

sample_one_strata = [  
                        {'sample_size': 600+300, 'winner_votes': 360+280}  
                    ]

# minerva for overall (1 stratum essentially)
overall_minerva = Minerva_S(alpha, 1.0, overall_contest)
overall_minerva.rounds.append(sample[0]['sample_size']+sample[1]['sample_size'])
overall_minerva.current_dist_null()
overall_minerva.current_dist_reported()
print("if just one big athena audit, pvalue: "+str(overall_minerva.stopping_condition(sample[0]['winner_votes']+sample[1]['winner_votes'])['pvalue']))

overall_stopping_condition(strata, overall_contest, sample, alpha)
overall_stopping_condition([overall_contest], overall_contest, sample_one_strata, alpha)


