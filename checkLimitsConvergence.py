#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, sys, commonFunctions

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--inputROOTFile', required=True, help='Name of input ROOT file containing observed and expected limits.',type=str)
inputArguments = inputArgumentsParser.parse_args()

limitsConverge = True
try:
    expectedUpperLimit, expectedUpperLimitOneSigmaDown, expectedUpperLimitOneSigmaUp, observedUpperLimit = commonFunctions.get_expected_and_observed_limits_from_combine_output(combineOutputFilePath=inputArguments.inputROOTFile)
    if (abs((observedUpperLimit/expectedUpperLimit)-1.0) >= 0.2): # Observed limit deviates from expected limit by more than 20%
        limitsConverge = False
except ValueError:
    limitsConverge = False

if (limitsConverge): print("true")
else: print("false")
