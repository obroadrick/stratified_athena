import numpy as np
from simulation_functions import simulate_fisher_combined_audits
from fishers_combination import calculate_lambda_range
import math

alpha = 0.1

# overall numbers (same as suite example 1)
N_w = 45500 + 7500
N_l = 49500 + 1500
N_relevant = N_w + N_l

# division into strata 
percent_polling = .1
N_w_fraction = .8
N_2 = round(percent_polling * N_relevant)
N_1 = N_relevant - N_2
N_w2 = round(N_w_fraction * N_2)
N_l2 = N_2 - N_w2
N_w1 = N_w - N_w2
N_l1 = N_1 - N_w1
margin = N_w1 + N_w2 - N_l1 - N_l2

# print the numbers
print ("N_relevant: "+str(N_relevant))
print ("N_1: "+str(N_1))
print ("N_w1: "+str(N_w1))
print ("N_l1: "+str(N_l1))
print ("N_2: "+str(N_2))
print ("N_w2: "+str(N_w2))
print ("N_l2: "+str(N_l2))
print ("N_w: "+str(N_w))
print ("N_l: "+str(N_l))
assert (N_l2 + N_l1 == N_l) # sanity check (can remove after first successful run)
assert (N_w2 + N_w1 == N_w) # sanity check (can remove after first successful run)
assert (N_1 + N_2 == N_relevant) # sanity check (can remove after first successful run)

np.random.seed(14123272)

n1 = 750 # same for all tests, same as in suite example
n2 = 295

reps = 1000

results = simulate_fisher_combined_audits(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, reps=reps, verbose=False, feasible_lambda_range=calculate_lambda_range(N_w1, N_l1, N_1, N_w2, N_l2, N_2), underlying=None)

print(results)


"""
left = 1
right = N_2 / 8
while(1):
    # current value to test
    n2 = math.ceil((left + right) / 2)

    # run simulation for this value
    results = simulate_fisher_combined_audits(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, reps=reps, verbose=False, feasible_lambda_range=calculate_lambda_range(N_w1, N_l1, N_1, N_w2, N_l2, N_2), underlying=None)
    stop = results['minerva']

    # print results for fun
    print("n2: "+str(n2)+"  stop: "+str(stop)+"  reps: "+str(reps))

    # update binary search bounds
    if (stop < .9):
        left = n2
    elif (stop < .9 + .01):
        # require more reps to be confident before stopping
        if (reps == 200):
            print("stopping at "+str(n2))
            break
    else:
        right = n2

    # make sure the bounds make sense
    if (left >= right):
        print("stopping at left>=right with left:"+str(n2))
        break

    # increase reps when getting closer
    if (stop > .9 - .1 and stop < .9 + .1):
        reps = 150
    if (stop > .9 - .05 and stop < .9 + .05):
        reps = 200
"""




