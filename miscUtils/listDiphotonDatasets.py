#!/usr/bin/env python3

from __future__ import print_function, division

import tmGeneralUtils, tmDASUtils

query_string = "/*DoubleEMEnriched*/*/MINIAODSIM"
list_of_datasets_matching_query = tmDASUtils.get_list_of_datasets(query_string)
event_counts = {}
for dataset in list_of_datasets_matching_query:
    print("Getting nEvents for dataset: {d}".format(d=dataset))
    event_counts[dataset] = tmDASUtils.get_number_of_events_in_dataset(dataset)

tmGeneralUtils.prettyPrintDictionary(inputDict=event_counts, keyPrintOrder=sorted(event_counts.keys(), key=(lambda dataset: event_counts[dataset])))
