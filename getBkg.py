#!/usr/bin/python

from __future__ import print_function, division

import os, sys, ROOT, argparse

import tmROOTUtils

inputArgumentsParser = argparse.ArgumentParser(description='Run STEALTH selection.')
inputArgumentsParser.add_argument('--nSTBins', default=5, help='Number of sT bins.',type=int)
inputArgumentsParser.add_argument('--sTMin', default=900., help='Min value of sT to plot.',type=float)
inputArgumentsParser.add_argument('--sTMax', default=3000., help='Max value of sT.',type=float)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max number of jets.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=2, help='Number of jets w.r.t. which to normalize the sT distributions for other jets.',type=int)
inputArgumentsParser.add_argument('--inputFilesSuffix', required=True, help='Prefix for input files.',type=str)
inputArgumentsParser.add_argument('--outputFilesSuffix', default="", help='Prefix for input files.',type=str)
inputArguments = inputArgumentsParser.parse_args()

n_jets_min = 2
n_jets_max = inputArguments.nJetsMax
inputFilesSuffix = inputArguments.inputFilesSuffix
outputFilesSuffix = inputArguments.outputFilesSuffix
if not outputFilesSuffix:
    outputFilesSuffix = inputArguments.inputFilesSuffix

# Global stores for minimum and maximum y-values, to use while plotting.
yMin = 0.
yMax = 0.

sTBackgroundFunctionNames = ["inversePowerLaw", "logModulatedInversePowerLaw", "inverseExponential"]
sTBackgroundFunctionDefinitions = {
    "inversePowerLaw": "[0]/TMath::Power(x/13000.,[1])",
    "logModulatedInversePowerLaw": "[0]/TMath::Power(x/13000.,[1]*TMath::Log(x))",
    "inverseExponential": "[0]/TMath::Exp([1]*x/13000.)"
}
sTBackgroundLineColors = {
    "inversePowerLaw": ROOT.kBlue,
    "logModulatedInversePowerLaw": ROOT.kBlack,
    "inverseExponential": ROOT.kRed
}

sTBackgroundLabels = {
    "inversePowerLaw": "1/S_{T}^{p_{0}}",
    "logModulatedInversePowerLaw": "1/S_{T}^{p_{1}lnS_{t}}",
    "inverseExponential": "1/e^{p_{2}S_{T}}"
}

legendParameters = {
    "x1": 0.72,
    "y1": 0.65,
    "x2": 0.86,
    "y2": 0.85
}

if (set(sTBackgroundFunctionNames) != set(sTBackgroundFunctionDefinitions.keys())):
    sys.exit("Unknown functions defined for sT background fit.")

sTBackgroundFunctions = {}
for sTBackgroundFunctionName in sTBackgroundFunctionNames:
    sTBackgroundFunctions[sTBackgroundFunctionName] = ROOT.TF1(sTBackgroundFunctionName, sTBackgroundFunctionDefinitions[sTBackgroundFunctionName], inputArguments.sTMin, inputArguments.sTMax)

def setYRange(listOfInputHistograms):
    yMax = 1.05*tmROOTUtils.getMaxValueFromListOfHistograms(listOfInputHistograms)
    yMin = 0.95*tmROOTUtils.getMinValueFromListOfHistograms(listOfInputHistograms)

def plotHistogramWithFit(sTHistogram, normScale, nJets):
    canvas = ROOT.TCanvas("c_{nJets}Jets".format(nJets=nJets), "c_{nJets}Jets".format(nJets=nJets), 1024, 768)
    canvas.SetBorderSize(0);
    canvas.SetFrameBorderMode(0)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetLogy()
    canvas.cd()

    sTHistogram.SetTitle("Double-fake-#gamma, {nJets} jets".format(nJets=nJets))
    sTHistogram.GetYaxis().SetRangeUser(yMin,yMax)
    sTHistogram.GetXaxis().SetTitle("S_{T} [GeV]")
    sTHistogram.GetXaxis().SetTitleOffset(1.1)
    sTHistogram.GetXaxis().SetRangeUser(inputArguments.sTMin,inputArguments.sTMax)
    sTHistogram.SetMarkerStyle(20)
    sTHistogram.SetMarkerSize(.9)
    sTHistogram.Draw("E")
    legend = ROOT.TLegend(legendParameters["x1"], legendParameters["y1"], legendParameters["x2"], legendParameters["y2"])
    legend.AddEntry(sTHistogram, "Background", "LP")

    sTBackgroundFunctionClones = {}
    for sTBackgroundFunctionName in sTBackgroundFunctionNames:
        sTBackgroundFunction = sTBackgroundFunctions[sTBackgroundFunctionName]
        sTBackgroundFunctionClone = sTBackgroundFunction.Clone()
        sTBackgroundFunctionClone.SetParameter(0,float(sTBackgroundFunctionClone.GetParameter(0))/normScale)
        sTBackgroundFunctionClone.SetLineWidth(2)
        sTBackgroundFunctionClone.SetLineColor(sTBackgroundLineColors[sTBackgroundFunctionName])
        sTBackgroundFunctionClone.Draw("SAME")
        canvas.Update()
        legend.AddEntry(sTBackgroundFunctionClone, sTBackgroundLabels[sTBackgroundFunctionName], "LP")
        sTBackgroundFunctionClones[sTBackgroundFunctionName] = sTBackgroundFunctionClone
        
    legend.SetBorderSize(0);
    legend.Draw("SAME")
    canvas.Update()

    canvas.Print("analysis/backgroundFit_%s_%djet.png"%(outputFilesSuffix,nJets))

## MAIN ##
def main():

    sTHistograms = {}
    inputFile = ROOT.TFile("analysis/hSTs_{inputFilesSuffix}.root".format(inputFilesSuffix=inputFilesSuffix),"READ")
    for nJets in range(n_jets_min,n_jets_max+1):
        hist_name = 'h_st_' + str(nJets) + 'Jets'
        sTHistograms[hist_name] = ROOT.gDirectory.Get(hist_name)

    setYRange(sTHistograms.values())

    # Estimate backgrounds analytically
    # i = 0
    hist_name_norm = 'h_st_' + str(inputArguments.nJetsNorm) + 'Jets'
    sTHistogramNorm = sTHistograms[hist_name_norm]
    fitResultsOutputFile = open('analysis/fit_choice_params_%s.dat'%(outputFilesSuffix),'w+')
    for sTBackgroundFunctionName in sTBackgroundFunctionNames:
        sTBackgroundFunction = sTBackgroundFunctions[sTBackgroundFunctionName]
        fitResult = sTHistogramNorm.Fit(sTBackgroundFunction, "M0NEI", "", inputArguments.sTMin, inputArguments.sTMax)
        fitStatus = int(fitResult)
        print("For background function \"{functionName}\", fit status = {fitStatus}".format(functionName=sTBackgroundFunctionName, fitStatus=fitStatus))
        if (fitStatus != 0 and fitStatus != 4000):
            print("WARNING: fit failed with status {fitStatus}".format(fitStatus=fitStatus))
        print("chiSq / ndf = {chiSq} / {ndf} = ".format(chiSq=sTBackgroundFunction.GetChisquare(), ndf=sTBackgroundFunction.GetNDF()))
        print("=====================")
        fitResultsOutputFile.write('{name}    {value0}    {error0}    {value1}    {error1}'.format(name=sTBackgroundFunctionName, value0=sTBackgroundFunction.GetParameter(0), error0=sTBackgroundFunction.GetParError(0), value1=sTBackgroundFunction.GetParameter(1), error1=sTBackgroundFunction.GetParError(1)))

    histogramScalesFile = open('analysis/normRatios_%s.txt'%(inputFilesSuffix), 'r')
    for histogramScalesLine in histogramScalesFile:
        nJets = int((histogramScalesLine.split())[0])
        normScale = float((histogramScalesLine.split())[1])
        print("nJets = {nJets}: normalization scale = {normScale}".format(nJets=nJets, normScale=normScale))
        hist_name = 'h_st_' + str(nJets) + 'Jets'
        plotHistogramWithFit(sTHistograms[hist_name], normScale, nJets)

    inputFile.Close()

#_____ Call main() ______#
if __name__ == '__main__':
    main()
