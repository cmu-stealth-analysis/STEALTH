#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, argparse, sys, commonFunctions

inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination from input MC in real data.')
inputArgumentsParser.add_argument('--inputROOTFile', required=True, help='Name of input ROOT file containing observed and expected limits.', type=str)
inputArgumentsParser.add_argument('--rMaxPassed', required=True, help='rMax passed while calling combine.', type=float)
inputArguments = inputArgumentsParser.parse_args()

TOLERANCE=0.01

bestFitValue = -1.0
convergentBestFitFound = True
try:
    bestFitValue = commonFunctions.get_best_fit_from_MultiDim_output(multiDimOutputFilePath=inputArguments.inputROOTFile)
    if not(bestFitValue < TOLERANCE):
        if ((abs((2.0*bestFitValue/(inputArguments.rMaxPassed)) - 1.0)) < TOLERANCE):
            convergentBestFitFound = False
except ValueError:
    convergentBestFitFound = False

if (convergentBestFitFound): print("true")
else: print("false")
