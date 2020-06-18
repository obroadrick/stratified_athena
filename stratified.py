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
        #overall_margin = self.overall_contest.tally[          #TODO
        for stratum_contest, stratum_lambda in zip(strata,lambda_):
            null_margin = stratum_lambda * overall_margin
            stratum_audit = Minerva_S(self.stratum_pvalue, 1.0, stratum_contest, null_margin)
            stratum_audits.append(stratum_audit)
        return stratum_audits

# TESTING
stratum1 = Contest(100000, {
    'A': 60000,
    'B': 40000
}, 1, ['A'], ContestType.PLURALITY)

stratum2 = Contest(50000, {
    'A': 28000,
    'B': 22000
}, 1, ['A'], ContestType.PLURALITY)

overall = Contest(150000, {
    'A': 28000+60000,
    'B': 22000+40000
}, 1, ['A'], ContestType.PLURALITY)

strata = []
strata.append(stratum1)
strata.append(stratum2)
stratified = Stratified(strata, overall, .05)

print(stratified.stratum_pvalue)

