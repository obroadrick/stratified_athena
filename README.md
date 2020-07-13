# stratified_athena

This code is messy and ugly. Slowly working towards things making more sense, but for now I more focused on getting things working and producing meaningful results, less focused on code ready to be integrated into r2b2....

A lot of this code is pulled straight from [the SUITE work](https://github.com/pbstark/CORLA18/tree/master/code) and a few files are pulled from [the minerva branch of r2b2](https://github.com/gwexploratoryaudits/r2b2/tree/minerva).

All current work is on [the s2 branch of this repo](https://github.com/obroadrick/stratified_athena/tree/s2), s2 meaning 2 strata (CVR and no-CVR).

Please ask me any questions that come up.

# files
## simulation_functions.py
This file has functions that can be used to run simulations of audits. This is a main focus right now.
## minerva_s.py and audit_s.py 
These files are adapted from the minerva branch of r2b2 (Grant's code) to handle null-margins, a concept unique to stratified audits.
## main_tests.py
This is the script that runs simulations. I change this frequently for whatever I'm currently testing.
