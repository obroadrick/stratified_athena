'''
Oliver Broadrick - Testing stratified Athena class audits.
'''
from audit_s import Audit_S
from contest import Contest
from contest import ContestType
from minerva_s import Minerva_S
import scipy
import numpy as np

class Stratified():
    def __init__(self, strata: [], overall_contest: Contest, alpha: float):
        self.strata = strata
        self.overall_contest = overall_contest
        self.alpha = alpha
        self.num_strata = len(strata)
        self.stratum_pvalue = self.find_stratum_pvalue()
       
    def fisher_combined_pvalue(self, pvalues):
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

    def find_stratum_pvalue(self, tol: float = .001):
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
            left_combined = self.fisher_combined_pvalue(np.full(self.num_strata, left))
            right_combined = self.fisher_combined_pvalue(np.full(self.num_strata, right))
            mid_combined = self.fisher_combined_pvalue(np.full(self.num_strata, mid))
            if (mid_combined > self.alpha):
                right = mid
            elif (mid_combined < (self.alpha - tol)):
                left = mid
            else:
                return mid
     
    def create_audits_for_lambda(self, lambda_: []):
        """Creates a list of audits for some value of lambda.
        
        Parameters:
            lambda_: values of lambda for each stratum
        """
        stratum_audits = []
        winner = 'A'
        loser = 'B'
        overall_margin = self.overall_contest.tally[winner] - self.overall_contest.tally[loser]
        print("overall margin: "+str(overall_margin))
        for stratum_contest, stratum_lambda in zip(strata,lambda_):
            null_margin = stratum_lambda * overall_margin
            stratum_audit = Minerva_S(self.stratum_pvalue, 1.0, stratum_contest, null_margin)
            stratum_audits.append(stratum_audit)
        return stratum_audits

    def predict_value_of_lambda(self):
        """Predict values of lambda by, yanno, totally guessing.

        Returns:
            []: values of lambda.
        """
        lambda_ = []
        i = 1
        for stratum in self.strata:
            stratum_lambda = stratum.contest_ballots / self.overall_contest.contest_ballots
            print("stratum "+str(i)+" lambda value: "+str(stratum_lambda))
            lambda_.append(stratum_lambda)
            i += 1

        print("sum of lambda values: "+str(sum(lambda_)))
        return lambda_

# TESTING

stratum1 = Contest(45500+49500, {
    'A': 49500,
    'B': 45500
}, 1, ['A'], ContestType.PLURALITY)

stratum2 = Contest(9000, {
    'A': 7500,
    'B': 1500
}, 1, ['A'], ContestType.PLURALITY)

overall_contest = Contest(45500+49500+9000, {
    'A': 45500+7500,
    'B': 49500+1500
}, 1, ['A'], ContestType.PLURALITY)
strata = []
strata.append(stratum1)
strata.append(stratum2)

alpha = .05

stratified = Stratified(strata, overall_contest, alpha)

strata_audits = stratified.create_audits_for_lambda(stratified.predict_value_of_lambda())

# check out what's going on by comparing kmins of stratified to overall minerva

overall_minerva = Minerva_S(alpha, 1.0, overall_contest) # null margin defaults to 0 in (tie)

total_drawn = 0
total_min = 0
i = 1
for stratum_audit in strata_audits:
    print ("----------------")
    print ("Stratum "+str(i))
    print ("winner votes: "+str(stratum_audit.contest.tally['A'])+"  loser votes: " + str(stratum_audit.contest.tally['B']))
    print ("required pvalue: "+str(stratum_audit.alpha))
    print ("null margin: "+str(stratum_audit.null_margin))
    stratum_drawn = 300
    stratum_audit.compute_min_winner_ballots([stratum_drawn])
    stratum_kmin = stratum_audit.min_winner_ballots[0]
    total_drawn += stratum_drawn
    total_min += stratum_kmin
    i+=1
    print ("drawn: "+str(stratum_drawn))
    print ("kmin: "+str(stratum_kmin))
    print ("----------------")
print("total drawn: "+str(total_drawn))
print("total min: "+str(total_min))

print ("\n\n----------------")
print("overall:")
draw = 600
overall_minerva.compute_min_winner_ballots([600])
print("drawn:"+str(600))
print("kmin:"+str(overall_minerva.min_winner_ballots[0]))
"""
ok so this initial test gives a lower kmin overall to the stratified 
audit which is not what i was exactly thinking would happen, but no
big deal because i'm guessing that my lambda prediction is just
not as good as i thought it might be
"""
"""
so no the plan is to come up with a more robust way of testing this stuff
potentially using the code from corla to assist
"""
"""
yanno what now that its been a couple minutes i think this result makes
sense and is about what i should have expected.
gonna also print pvalues now
"""

