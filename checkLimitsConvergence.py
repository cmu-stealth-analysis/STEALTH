#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, sys, commonFunctions

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--inputROOTFile', required=True, help='Name of input ROOT file containing observed and expected limits.',type=str)
inputArguments = inputArgumentsParser.parse_args()

ACCEPTABLE_LOWER_RATIO = 0.2
ACCEPTABLE_UPPER_RATIO = 5.0
UPPER_LIMIT_CHECK_RELAXATION_THRESHOLD = 0.1

limitsConverge = True
try:
    expectedUpperLimit, expectedUpperLimitOneSigmaDown, expectedUpperLimitOneSigmaUp, observedUpperLimit = commonFunctions.get_expected_and_observed_limits_from_combine_output(combineOutputFilePath=inputArguments.inputROOTFile)
    ratiosToCheck = [observedUpperLimit/expectedUpperLimit, expectedUpperLimitOneSigmaDown/expectedUpperLimit, expectedUpperLimitOneSigmaUp/expectedUpperLimit]

    # for low expected limit values, the combine tool sometimes seems to not converge... doesn't seem to matter if we're away from the expected limit contours
    # the 50% and "one sigma up" expected limits still seem to be reasonable, so we can use them to check for convergence
    if (expectedUpperLimit < UPPER_LIMIT_CHECK_RELAXATION_THRESHOLD): ratiosToCheck = [expectedUpperLimitOneSigmaUp/expectedUpperLimit]

    for ratio in ratiosToCheck:
        if ((ratio < ACCEPTABLE_LOWER_RATIO) or (ratio > ACCEPTABLE_UPPER_RATIO)):
            limitsConverge = False
            break
except ValueError:
    limitsConverge = False

if (limitsConverge): print("true")
else: print("false")
