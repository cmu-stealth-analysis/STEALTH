#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, sys, commonFunctions

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--inputROOTFile', required=True, help='Name of input ROOT file containing observed and expected limits.',type=str)
inputArguments = inputArgumentsParser.parse_args()

try:
    bestFitValue = commonFunctions.get_best_fit_from_MultiDim_output(multiDimOutputFilePath=inputArguments.inputROOTFile)
except ValueError:
    sys.exit("Error: multidim fit output not in expected format.")

print("{bFV:.6f}".format(bFV=bestFitValue))
