"""
Small script for getting comfortable with the 2016 presidential election json data.

Oliver Broadrick 2020
"""

"""
States mentioned today:
    Rhode Island
    Wisconsin       # who knew that it's wisonsIn not wisonsOn?? not me...
    Michigan
    Pennsylvania
    Ohio
    North Carolina
"""

import json
import pprint

# make an array of the state names that i care about
states_i_care_about =   [
        'RhodeIsland',
        'Wisconsin',
        'Michigan',
        'Pennsylvania',
        'Ohio',
        'NorthCarolina'
    ]

# open desired stratified data json file
with open('data/2016_one_round_all.json') as json_file:
    data = json.load(json_file)
    #pprint.pprint(data)

# print only the states I currently care about
for state in data:
    if state in states_i_care_about:
        print("\n", state, ":")
        pprint.pprint(data[state])

