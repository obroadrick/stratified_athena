from audit_s import Audit_S
from contest import Contest
from contest import ContestType
from minerva_s import Minerva_S
import scipy
import numpy as np
 
alpha = .1

contest = Contest(100000, {
    'A': 60000,
    'B': 40000
}, 1, ['A'], ContestType.PLURALITY)

for null_margin in [-3000,3000]:
    print("null_margin: "+str(null_margin))

    audit = Minerva_S(alpha, 1.0, contest, null_margin)

    n = 500

    for k in range(240,280):
        audit.rounds.append(n)

        audit.current_dist_reported()
        audit.current_dist_null()

        result = audit.stopping_condition(k)

        print("k: "+str(k)+", result: "+str(result))

