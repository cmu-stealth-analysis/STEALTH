#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, sys, commonFunctions

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--inputROOTFile', required=True, help='Name of input ROOT file containing observed and expected limits.',type=str)
inputArguments = inputArgumentsParser.parse_args()

isInExpectedFormat = True
try:
    limits = commonFunctions.read_limits_from_combine_output(inputArguments.inputROOTFile)
except ValueError:
    isInExpectedFormat = False

if (isInExpectedFormat): print("true")
else: print("false")
