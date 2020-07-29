# stratified_athena

This code is messy and ugly. Slowly working towards things making more sense, but for now I more focused on getting things working and producing meaningful results, less focused on code ready to be integrated into r2b2....

A lot of this code is pulled straight from [the SUITE work](https://github.com/pbstark/CORLA18/tree/master/code) and a few files are pulled from [the minerva branch of r2b2](https://github.com/gwexploratoryaudits/r2b2/tree/minerva).

All current work is on [the s2 branch of this repo](https://github.com/obroadrick/stratified_athena/tree/s2), s2 meaning 2 strata (CVR and no-CVR).

Please ask me any questions that come up.

# files
### audit_s.py
Audit class from r2b2, modified to handle null hypotheses of any margin. 

###ballot_comparison.py
Ballot comparison file from SUITE.

###contest.py
Contest file from r2b2.

###data
Directory with json files containing results from various tests, simulations, and computations of round sizes.

###dist_pvalues.py
Incomplete script for creating a plot of the probability distribution over possible pvalues for a stratified audit.

###easy_read.py
Simple script for printing data from the json data files.

###election.py
Election file from r2b2.

###fishers_combination.py
Fishers combination functions from SUITE.

###hypergeometric.py
Hypergeometric functions from SUITE. (Not used?)

###minerva_round_size_script.py
Script that computes the round size for a contest-wide Minerva audit with no stratification.

###minerva_s.py
Minerva file from r2b2 modified to handle a null hypothesis with nonzero margin.

###old
Some files that are no longer used from when I first started exploring this topic, including some functions (incomplete) for more than two strata.

###process_data.py
Script that produces a plot of round sizes from the json data files.

###pvalue_dist_script.py
Incomplete script that creates a plot of the probability distribution over pvalues for a stratified audit.

###__pycache__
This is a python thing.

###r2bravo_round_size_script.py
Script that computes the round size for a contest-wide R2 Bravo audit with no stratification.

###README.md
This file.

###round_sizes.py
Functions for computing a first round size which accomplishes a given probability of stopping under the alternative hypothesis, that the election is truly as announced.

###simulation_functions.py
Functions for simulating stratified audits. Also includes my own functions for computing the pvalue for one-round Minerva audits and R2 Bravo audits.

###simulation_script.py
Script for running simulations of stratified audits.

###sprt.py
Functions for R2 Bravo from SUITE.

###stratified_round_size_script.py
Script to compute the round size which accomplishes a given probability of stopping for various one-round stratified audits.

