#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, ROOT, tmROOTUtils, array, pdb, math
from tmProgressBar import tmProgressBar

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Generate simple plot of ST distributions and scan over different norm regions to see in which one ST scaling is a good approximation.')
inputArgumentsParser.add_argument('--inputFilePath', required=True, help='Path to input file.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="STScaling", help='Output directory.',type=str)
inputArgumentsParser.add_argument('--outputFilePrefix', required=True, help='Name of output file.',type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--targetSTNorm', default=1150., help='Value of ST at which to normalize all histograms.',type=float)
inputArgumentsParser.add_argument('--nJetsMin', default=2, help='Min nJets bin.',type=int)
inputArgumentsParser.add_argument('--nJetsMax', default=6, help='Max nJets bin.',type=int)
inputArgumentsParser.add_argument('--nJetsNorm', default=2, help='Norm nJets bin.',type=int)
inputArgumentsParser.add_argument('--nScanPoints', default=51, help='Number of points for normalization region scan.',type=int)
inputArgumentsParser.add_argument('--normScanMin', default=900., help='Value of bin center at which to begin normalization region scan.',type=float)
inputArgumentsParser.add_argument('--normScanMax', default=1400., help='Value of bin center at which to end normalization region scan.',type=float)
inputArgumentsParser.add_argument('--target_STMax', default=2500., help='Max value of ST to target in normalization region scan.',type=float)
inputArgumentsParser.add_argument('--target_binWidth', default=100., help='Bin width to target in normalization region scan.',type=float)
inputArgumentsParser.add_argument('--maxSTForGlobalFit', default=1050., help='Max ST up to which to run global fit.',type=float)
inputArguments = inputArgumentsParser.parse_args()

histColors = {
    2: ROOT.kBlack,
    3: ROOT.kBlue,
    4: ROOT.kRed,
    5: ROOT.kGreen,
    6: ROOT.kViolet
}

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries)
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(3500.0) # Instead of infinity
n_STBins = len(STBoundaries) - 1
STRegionsAxis = ROOT.TAxis(n_STBins, array.array('d', STBoundaries))

# Load input TTrees into TChain
inputChain = ROOT.TChain("ggNtuplizer/EventTree")
inputChain.Add(inputArguments.inputFilePath)

nEvents = inputChain.GetEntries()
if (nEvents == 0): sys.exit("Number of available events is 0.")

STHistograms = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    STHistograms[nJetsBin] = ROOT.TH1F("h_STDistribution_{n}Jets".format(n=nJetsBin), "ST Distribution, nJets = {n};ST (GeV);A.U.".format(n = nJetsBin), n_STBins, array.array('d', STBoundaries))
    STHistograms[nJetsBin].Sumw2()

progressBar = tmProgressBar(nEvents)
progressBarUpdatePeriod = max(1, (nEvents//1000))
progressBar.initializeTimer()
STValues_nJetsMin = []
STValues_nJetsMax = []
for eventIndex in range(0,nEvents):
    treeStatus = inputChain.LoadTree(eventIndex)
    if treeStatus < 0:
        # sys.exit("Tree unreadable.")
        break
    evtStatus = inputChain.GetEntry(eventIndex)
    if evtStatus <= 0:
        # sys.exit("Event in tree unreadable.")
        continue
    if (eventIndex%progressBarUpdatePeriod == 0 or eventIndex == (nEvents - 1)): progressBar.updateBar(eventIndex/nEvents, eventIndex)

    nJetsBin = inputChain.b_nJets
    if (nJetsBin > inputArguments.nJetsMax): nJetsBin = inputArguments.nJetsMax
    if (nJetsBin < inputArguments.nJetsMin): sys.exit("Unexpected nJetsBin = {nJetsBin} at entry index = {entryIndex}".format(nJetsBin=nJetsBin, entryIndex=entryIndex))

    sT = inputChain.b_evtST
    STHistograms[nJetsBin].Fill(sT)
    if (nJetsBin == inputArguments.nJetsMin): STValues_nJetsMin.append(sT)
    if (nJetsBin == inputArguments.nJetsMax): STValues_nJetsMax.append(sT)
progressBar.terminate()

outputCanvas = ROOT.TCanvas("outputCanvas", "outputCanvas", 1024, 768)
upperPad = ROOT.TPad("upperPad", "upperPad", 0.0, 0.4, 1.0, 1.0)
lowerPad = ROOT.TPad("lowerPad", "lowerPad", 0.0, 0.0, 1.0, 0.4)
upperPad.SetBottomMargin(0.025)
lowerPad.SetTopMargin(0.025)
upperPad.Draw()
lowerPad.Draw()
upperPad.cd()
outputLegend = ROOT.TLegend(0.8, 0.65, 0.9, 0.9)
ROOT.gPad.SetLogy()
ROOT.gStyle.SetOptStat(0)
STHistograms[inputArguments.nJetsNorm].SetLineColor(histColors[inputArguments.nJetsNorm])
STHistograms[inputArguments.nJetsNorm].SetTitle("Comparison of ST Distributions")
STHistograms[inputArguments.nJetsNorm].Draw("HIST E1")
STHistograms[inputArguments.nJetsNorm].GetXaxis().SetLabelOffset(999);
STHistograms[inputArguments.nJetsNorm].GetXaxis().SetLabelSize(0);
legendEntry = outputLegend.AddEntry(STHistograms[inputArguments.nJetsNorm], "nJets = {n}".format(n = inputArguments.nJetsNorm))
legendEntry.SetTextColor(histColors[inputArguments.nJetsNorm])
legendEntry.SetLineColor(histColors[inputArguments.nJetsNorm])
nTargetEntries_normBin = STHistograms[inputArguments.nJetsNorm].GetBinContent(STHistograms[inputArguments.nJetsNorm].FindFixBin(inputArguments.targetSTNorm))
ratioHistograms = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if (nJetsBin == inputArguments.nJetsNorm): continue
    nEntries_normBin = STHistograms[nJetsBin].GetBinContent(STHistograms[nJetsBin].FindFixBin(inputArguments.targetSTNorm))
    STHistograms[nJetsBin].Scale(nTargetEntries_normBin/nEntries_normBin)
    ratioHistograms[nJetsBin] = ROOT.TH1F("h_STDistributionsRatio_{n}Jets".format(n=nJetsBin), "Ratio of ST Distributions, nJets = {n};ST (GeV);ratio".format(n = nJetsBin), n_STBins, array.array('d', STBoundaries))
    ratioHistograms[nJetsBin].Divide(STHistograms[nJetsBin], STHistograms[inputArguments.nJetsNorm])
    STHistograms[nJetsBin].SetLineColor(histColors[nJetsBin])
    STHistograms[nJetsBin].Draw("HIST E1 SAME")
    legendEntry = ROOT.TLegendEntry()
    if (nJetsBin == inputArguments.nJetsMax): legendEntry = outputLegend.AddEntry(STHistograms[nJetsBin], "nJets #geq {n}".format(n = inputArguments.nJetsMax))
    else: legendEntry = outputLegend.AddEntry(STHistograms[nJetsBin], "nJets = {n}".format(n = nJetsBin))
    legendEntry.SetTextColor(histColors[nJetsBin])
    legendEntry.SetLineColor(histColors[nJetsBin])
outputLegend.Draw()
lowerPad.cd()
isFirstToBeDrawn = True
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    if (nJetsBin == inputArguments.nJetsNorm): continue
    ratioHistograms[nJetsBin].SetLineColor(histColors[nJetsBin])
    if (isFirstToBeDrawn):
        ratioHistograms[nJetsBin].SetTitle("")
        ratioHistograms[nJetsBin].Draw("E1")
        ratioHistograms[nJetsBin].GetYaxis().SetRangeUser(0.0, 3.0)
        isFirstToBeDrawn = False
    else:
        ratioHistograms[nJetsBin].Draw("E1 SAME")
lineAt1 = ROOT.TLine(STBoundaries[0], 1.0, STBoundaries[-1], 1.0)
lineAt1.SetLineColor(histColors[inputArguments.nJetsNorm])
lineAt1.Draw()
outputCanvas.SaveAs("{oD}/{oFP}_STDistributions.png".format(oD=inputArguments.outputDirectory, oFP=inputArguments.outputFilePrefix))

def constantFunction(x, par):
    return (par[0])

def linearFunction(x, par):
    return (par[0] + x[0]*par[1])

def quadraticFunction(x, par):
    return (par[0] + x[0]*par[1] + x[0]*x[0]*par[2])

def globalFunction(x, par):
    if (x[0] >= par[4]): return (par[3])
    else: return (par[0] + x[0]*par[1] + x[0]*x[0]*par[2])

scalingCandidates = {}
maxJetDistributions = {}
minJetDistributions = {}
constantFitChi2PerNDF = ROOT.TGraph()
constantFitChi2PerNDF.SetName("constantFitChi2PerNDF")
constantFitChi2PerNDF.SetTitle(";lower edge of norm bin;#chi^{2}/n.d.f.")
linearFitChi2PerNDF = ROOT.TGraph()
linearFitChi2PerNDF.SetName("linearFitChi2PerNDF")
linearFitChi2PerNDF.SetTitle(";lower edge of norm bin;#chi^{2}/n.d.f.")
quadraticFitChi2PerNDF = ROOT.TGraph()
quadraticFitChi2PerNDF.SetName("quadraticFitChi2PerNDF")
quadraticFitChi2PerNDF.SetTitle(";lower edge of norm bin;#chi^{2}/n.d.f.")
ratioLinearToConstant = ROOT.TGraph()
ratioLinearToConstant.SetName("ratioLinearToConstant")
ratioLinearToConstant.SetTitle(";lower edge of norm bin;(#chi^{2}/n.d.f.)/(#chi^{2}/n.d.f.)_{constant fit}")
ratioQuadraticToConstant = ROOT.TGraph()
ratioQuadraticToConstant.SetName("ratioQuadraticToConstant")
ratioQuadraticToConstant.SetTitle(";lower edge of norm bin;(#chi^{2}/n.d.f.)/(#chi^{2}/n.d.f.)_{constant fit}")
ratioConstantToConstant = ROOT.TGraph() # lol
ratioConstantToConstant.SetName("ratioConstantToConstant")
ratioConstantToConstant.SetTitle(";lower edge of norm bin;(#chi^{2}/n.d.f.)/(#chi^{2}/n.d.f.)_{constant fit}")

globalFitChi2PerNDF = ROOT.TGraph()
globalFitChi2PerNDF.SetName("globalFitChi2PerNDF")
globalFitChi2PerNDF.SetTitle(";lower edge of norm bin;Best ST range")
for scalingHistogramIndex in range(0, inputArguments.nScanPoints):
    candidateNormMin = inputArguments.normScanMin + (scalingHistogramIndex)*(inputArguments.normScanMax - inputArguments.normScanMin)/(inputArguments.nScanPoints-1)
    candidateNormMax = candidateNormMin + inputArguments.target_binWidth
    nCandidateBins = int(0.5 + math.ceil((inputArguments.target_STMax - candidateNormMin)/inputArguments.target_binWidth))
    candidateSTMax = candidateNormMin + nCandidateBins*inputArguments.target_binWidth
    print("Generating fits for normMin: {n}, normMax: {x}, nBins: {nCB}, STMax: {cSTM}".format(n=candidateNormMin, x=candidateNormMax, nCB=nCandidateBins, cSTM=candidateSTMax))
    scalingCandidates[scalingHistogramIndex] = ROOT.TH1F("h_scalingCandidate_{i}".format(i=scalingHistogramIndex), "Distribution: max jets/min jets;ST (GeV);A.U.", nCandidateBins, candidateNormMin, candidateSTMax)
    scalingCandidates[scalingHistogramIndex].Sumw2()
    maxJetDistributions[scalingHistogramIndex] = ROOT.TH1F("h_maxJetDistribution_{i}".format(i=scalingHistogramIndex), "Distribution: max jets;ST (GeV);A.U.", nCandidateBins, candidateNormMin, candidateSTMax)
    maxJetDistributions[scalingHistogramIndex].Sumw2()
    minJetDistributions[scalingHistogramIndex] = ROOT.TH1F("h_minJetDistribution_{i}".format(i=scalingHistogramIndex), "Distribution: min jets;ST (GeV);A.U.", nCandidateBins, candidateNormMin, candidateSTMax)
    minJetDistributions[scalingHistogramIndex].Sumw2()

    # Fill min jet distributions
    for STValue in STValues_nJetsMin:
        if (STValue > candidateNormMin): minJetDistributions[scalingHistogramIndex].Fill(STValue)
    # Fill max jet distributions
    for STValue in STValues_nJetsMax:
        if (STValue > candidateNormMin): maxJetDistributions[scalingHistogramIndex].Fill(STValue)

    # Get ratio histogram
    nEffectiveBins = 0
    for ratioBinIndex in range(1, 1+nCandidateBins):
        numerator = maxJetDistributions[scalingHistogramIndex].GetBinContent(ratioBinIndex)
        numeratorError = maxJetDistributions[scalingHistogramIndex].GetBinError(ratioBinIndex)
        denominator = minJetDistributions[scalingHistogramIndex].GetBinContent(ratioBinIndex)
        denominatorError = minJetDistributions[scalingHistogramIndex].GetBinError(ratioBinIndex)
        if ((numerator > 0) and (denominator > 0)):
            scalingCandidates[scalingHistogramIndex].SetBinContent(ratioBinIndex, numerator/denominator)
            scalingCandidates[scalingHistogramIndex].SetBinError(ratioBinIndex, (numerator/denominator)*(numeratorError/numerator + denominatorError/denominator))
            nEffectiveBins += 1
        else:
            scalingCandidates[scalingHistogramIndex].SetBinContent(ratioBinIndex, 0.)
            scalingCandidates[scalingHistogramIndex].SetBinError(ratioBinIndex, 0.)
        # print("At bin index = {i}, bin content = {c}, bin error = {e}".format(i = ratioBinIndex, c = scalingCandidates[scalingHistogramIndex].GetBinContent(ratioBinIndex), e = scalingCandidates[scalingHistogramIndex].GetBinError(ratioBinIndex)))

    constantCandidate = ROOT.TF1("constantCandidateFit_{i}".format(i=scalingHistogramIndex), constantFunction, candidateNormMin, candidateSTMax, 1)
    constantCandidate.SetParameter(0, 1.)
    constantFitResult = scalingCandidates[scalingHistogramIndex].Fit(constantCandidate, "IEMN0S")
    constantFitChi2 = constantFitResult.Chi2()
    constantFitNDF = constantFitResult.Ndf()
    if not(constantFitNDF == (nEffectiveBins - 1)): sys.exit("Error in understanding: nEffectiveBins = {n}, NDF for constant fit = {ndf}".format(n=nEffectiveBins, ndf=constantFitNDF))
    print("Constant fit result: p0 = {c} +/- {e}; chi2 = {chi2}, ndf = {ndf}".format(c = constantFitResult.Parameter(0), e = constantFitResult.ParError(0), chi2=constantFitChi2, ndf=constantFitNDF))
    constantFitChi2PerNDF.SetPoint(constantFitChi2PerNDF.GetN(), candidateNormMin, constantFitChi2/constantFitNDF)
    ratioConstantToConstant.SetPoint(ratioConstantToConstant.GetN(), candidateNormMin, 1.)
    
    linearCandidate = ROOT.TF1("linearCandidateFit_{i}".format(i=scalingHistogramIndex), linearFunction, candidateNormMin, candidateSTMax, 2)
    linearCandidate.SetParameter(0, 1.)
    linearCandidate.SetParameter(1, 0.)
    linearFitResult = scalingCandidates[scalingHistogramIndex].Fit(linearCandidate, "IEMN0S")
    linearFitChi2 = linearFitResult.Chi2()
    linearFitNDF = linearFitResult.Ndf()
    if not(linearFitNDF == (nEffectiveBins - 2)): sys.exit("Error in understanding: nEffectiveBins = {n}, NDF for linear fit = {ndf}".format(n=nEffectiveBins, ndf=linearFitNDF))
    print("Linear fit result: p0 = {c} +/- {e}; p1 = {c1} +/- {e1}; chi2 = {chi2}, ndf = {ndf}".format(c = linearFitResult.Parameter(0), e = linearFitResult.ParError(0), c1 = linearFitResult.Parameter(1), e1 = linearFitResult.ParError(1), chi2 = linearFitChi2, ndf = linearFitNDF))
    linearFitChi2PerNDF.SetPoint(linearFitChi2PerNDF.GetN(), candidateNormMin, linearFitChi2/linearFitNDF)
    ratioLinearToConstant.SetPoint(ratioLinearToConstant.GetN(), candidateNormMin, (linearFitChi2*constantFitNDF)/(constantFitChi2*linearFitNDF))

    quadraticCandidate = ROOT.TF1("quadraticCandidateFit_{i}".format(i=scalingHistogramIndex), quadraticFunction, candidateNormMin, candidateSTMax, 3)
    quadraticCandidate.SetParameter(0, 1.)
    quadraticCandidate.SetParameter(1, 0.)
    quadraticCandidate.SetParameter(2, 0.)
    quadraticFitResult = scalingCandidates[scalingHistogramIndex].Fit(quadraticCandidate, "IEMN0S")
    quadraticFitChi2 = quadraticFitResult.Chi2()
    quadraticFitNDF = quadraticFitResult.Ndf()
    if not(quadraticFitNDF == (nEffectiveBins - 3)): sys.exit("Error in understanding: nEffectiveBins = {n}, NDF for quadratic fit = {ndf}".format(n=nEffectiveBins, ndf=quadraticFitNDF))
    print("Quadratic fit result: p0 = {c} +/- {e}; p1 = {c1} +/- {e1}; p2 = {c2} +/- {e2}; chi2 = {chi2}, ndf = {ndf}".format(c = quadraticFitResult.Parameter(0), e = quadraticFitResult.ParError(0), c1 = quadraticFitResult.Parameter(1), e1 = quadraticFitResult.ParError(1), c2 = quadraticFitResult.Parameter(2), e2 = quadraticFitResult.ParError(2), chi2 = quadraticFitChi2, ndf = quadraticFitNDF))
    quadraticFitChi2PerNDF.SetPoint(quadraticFitChi2PerNDF.GetN(), candidateNormMin, quadraticFitChi2/quadraticFitNDF)
    ratioQuadraticToConstant.SetPoint(ratioQuadraticToConstant.GetN(), candidateNormMin, (quadraticFitChi2*constantFitNDF)/(constantFitChi2*quadraticFitNDF))

    if (candidateNormMin <= inputArguments.maxSTForGlobalFit):
        globalCandidate = ROOT.TF1("globalCandidateFit_{i}".format(i=scalingHistogramIndex), globalFunction, candidateNormMin, candidateSTMax, 5)
        globalCandidate.SetParameter(0, 1.)
        globalCandidate.SetParameter(1, 0.)
        globalCandidate.SetParameter(2, 0.)
        globalCandidate.SetParameter(3, 1.)
        globalCandidate.SetParameter(4, 1000.)
        globalFitResult = scalingCandidates[scalingHistogramIndex].Fit(globalCandidate, "IEMN0S")
        globalFitChi2 = globalFitResult.Chi2()
        globalFitNDF = globalFitResult.Ndf()
        if not(globalFitNDF == (nEffectiveBins - 5)): sys.exit("Error in understanding: nEffectiveBins = {n}, NDF for global fit = {ndf}".format(n=nEffectiveBins, ndf=globalFitNDF))
        print("Global fit result: p0 = {c} +/- {e}; p1 = {c1} +/- {e1}; p2 = {c2} +/- {e2}; asymptotic constant = {c3} +/- {e3}; turn-on = {c4} +/- {e4}; chi2 = {chi2}, ndf = {ndf}".format(c = globalFitResult.Parameter(0), e = globalFitResult.ParError(0), c1 = globalFitResult.Parameter(1), e1 = globalFitResult.ParError(1), c2 = globalFitResult.Parameter(2), e2 = globalFitResult.ParError(2), c3 = globalFitResult.Parameter(3), e3 = globalFitResult.ParError(3), c4 = globalFitResult.Parameter(4), e4 = globalFitResult.ParError(4), chi2 = globalFitChi2, ndf = globalFitNDF))
        globalFitChi2PerNDF.SetPoint(globalFitChi2PerNDF.GetN(), candidateNormMin, globalFitResult.Parameter(4))

constantFitChi2PerNDF.SetLineColor(ROOT.kBlack)
ratioConstantToConstant.SetLineColor(ROOT.kBlack)
linearFitChi2PerNDF.SetLineColor(ROOT.kBlue)
ratioLinearToConstant.SetLineColor(ROOT.kBlue)
quadraticFitChi2PerNDF.SetLineColor(ROOT.kRed)
ratioQuadraticToConstant.SetLineColor(ROOT.kRed)
tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[constantFitChi2PerNDF, linearFitChi2PerNDF, quadraticFitChi2PerNDF], canvasName="c_chi2s", outputDocumentName="{oD}/{oFP}_chi2Values".format(oD=inputArguments.outputDirectory, oFP=inputArguments.outputFilePrefix), customOptStat=0)
tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[ratioConstantToConstant, ratioLinearToConstant, ratioQuadraticToConstant], canvasName="c_chi2Ratios", outputDocumentName="{oD}/{oFP}_chi2Ratios".format(oD=inputArguments.outputDirectory, oFP=inputArguments.outputFilePrefix), customOptStat=0, customYRange=[0., 2.])
tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[globalFitChi2PerNDF], canvasName="c_globalFits", outputDocumentName="{oD}/{oFP}_globalFits".format(oD=inputArguments.outputDirectory, oFP=inputArguments.outputFilePrefix), customOptStat=0)

print("All done!")
