#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, os, sys, argparse
inputArgumentsParser = argparse.ArgumentParser(description='Run Higgs combine tool on already generated datacards.')
inputArgumentsParser.add_argument('--dataCardsDirectory', default="analysis/dataCards", help='Path to directory containing already generated datacards.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="analysis/combineToolOutputs", help='Path to directory containing already generated datacards.',type=str)
inputArgumentsParser.add_argument('--MCTemplate', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArguments = inputArgumentsParser.parse_args()

generatedMCTemplate = ROOT.TFile(inputArguments.MCTemplate)
h_MCTemplate = generatedMCTemplate.Get("h_susyMasses_template")
for gluinoMassBin in range(1, 1+h_MCTemplate.GetXaxis().GetNbins()):
    for neutralinoMassBin in range(1, 1+h_MCTemplate.GetYaxis().GetNbins()):
        if not(int(0.5 + h_MCTemplate.GetBinContent(gluinoMassBin, neutralinoMassBin)) == 1): continue
        gluinoMass = int(0.5 + h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin))
        neutralinoMass = int(0.5 + h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin))

        os.system("./runCombineToolHelper.sh {dataCardsDirectory} {gluinoMass} {neutralinoMass} {outputDirectory}".format(dataCardsDirectory=inputArguments.dataCardsDirectory, gluinoMass=gluinoMass, neutralinoMass=neutralinoMass, outputDirectory=inputArguments.outputDirectory))
generatedMCTemplate.Close()
