'''
Oliver Broadrick
Testing stratified Athena class audits.
Will use contest to represent an individual stratum.
'''

from audit_s import Audit_S
import election
from contest import Contest
from contest import ContestType
from minerva_s import Minerva_S

stratum1 = Contest(100000, {
    'A': 60000,
    'B': 40000
}, 1, ['A'], ContestType.PLURALITY)

stratum2 = Contest(50000, {
    'A': 28000,
    'B': 22000
}, 1, ['A'], ContestType.PLURALITY)

# check that null dists can account for an overstatement quota
audit1 = Minerva_S(.05, 1.0, stratum1,0)
audit1.compute_min_winner_ballots([100,200])
print(audit1.min_winner_ballots)

audit2 = Minerva_S(.05, 1.0, stratum1,10000)
audit2.compute_min_winner_ballots([100,200])
print(audit2.min_winner_ballots)

audit3 = Minerva_S(.05, 1.0, stratum2,-10000)
audit3.compute_min_winner_ballots([100,200])
print(audit3.min_winner_ballots)

overall = Contest(150000, {
    'A': 28000+60000,
    'B': 22000+40000
}, 1, ['A'], ContestType.PLURALITY)

strata = []
strata.append(stratum1)
strata.append(stratum2)

