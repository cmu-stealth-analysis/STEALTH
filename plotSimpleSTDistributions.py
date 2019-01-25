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
# inputArgumentsParser.add_argument('--maxSTForGlobalFit', default=1050., help='Max ST up to which to run global fit.',type=float)
inputArgumentsParser.add_argument('--nominalRho', default=1.25, help='Nominal value of rho to use to construct the pdf for the nJets = 2 and nJets = 6 data.',type=float)
inputArgumentsParser.add_argument('--usePDF', action='store_true', help="If this flag is set, then, for nJets >= 4, the number of events in a given bin is obtained from the PDF estimate instead of from the raw data.")
inputArgumentsParser.add_argument('--usePDFDenominator', action='store_true', help="If this flag is set, then, for nJets = 2, the number of events in a given bin is obtained from the PDF estimate instead of from the raw data.")
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

rooVar_ST = ROOT.RooRealVar("rooVar_ST", "rooVar_ST", inputArguments.normScanMin, inputArguments.target_STMax + inputArguments.target_binWidth, "GeV")
STValues = {}
STTrees = {}
STArrays = {}
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    STValues[nJetsBin] = []
    STTrees[nJetsBin] = ROOT.TTree("tree_{n}Jets".format(n=nJetsBin), "tree_{n}Jets".format(n=nJetsBin))
    STArrays[nJetsBin] = array.array('f', [0.])
    STTrees[nJetsBin].Branch('rooVar_ST', STArrays[nJetsBin], 'rooVar_ST/F')

progressBar = tmProgressBar(nEvents)
progressBarUpdatePeriod = max(1, (nEvents//1000))
progressBar.initializeTimer()
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

    ST = inputChain.b_evtST
    STHistograms[nJetsBin].Fill(ST)
    (STArrays[nJetsBin])[0] = ST
    STTrees[nJetsBin].Fill()
    STValues[nJetsBin].append(ST)
progressBar.terminate()

rooDataSets = {}
rooKernel_PDF_Estimators = {}
STFrame = rooVar_ST.frame(inputArguments.normScanMin, inputArguments.target_STMax + inputArguments.target_binWidth, int(0.5 + math.ceil((inputArguments.target_STMax - inputArguments.normScanMin)/inputArguments.target_binWidth)))
rooVar_ST.setRange("globalNormRange", inputArguments.normScanMax + inputArguments.target_binWidth, inputArguments.target_STMax)
for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    rooDataSets[nJetsBin] = ROOT.RooDataSet("rooDataSets_{n}Jets".format(n=nJetsBin), "rooDataSets_{n}Jets".format(n=nJetsBin), STTrees[nJetsBin], ROOT.RooArgSet(rooVar_ST))
    rooKernel_PDF_Estimators[nJetsBin] = ROOT.RooKeysPdf("kernelEstimate_{n}Jets".format(n=nJetsBin), "kernelEstimate_{n}Jets".format(n=nJetsBin), rooVar_ST, rooDataSets[nJetsBin], ROOT.RooKeysPdf.MirrorLeft, inputArguments.nominalRho)
    rooKernel_PDF_Estimators[nJetsBin].plotOn(STFrame, ROOT.RooFit.NormRange("globalNormRange"))
tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [STFrame], canvasName = "c_STKernels_linearScale", outputDocumentName="{oD}/{oFP}_STKernels_linearScale".format(oD=inputArguments.outputDirectory, oFP=inputArguments.outputFilePrefix))
tmROOTUtils.plotObjectsOnCanvas(listOfObjects = [STFrame], canvasName = "c_STKernels", outputDocumentName="{oD}/{oFP}_STKernels".format(oD=inputArguments.outputDirectory, oFP=inputArguments.outputFilePrefix), enableLogY = True)

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

# def globalFunction(x, par):
#     if (x[0] >= par[4]): return (par[3])
#     else: return (par[0] + x[0]*par[1] + x[0]*x[0]*par[2])

def resetSTRange():
    rooVar_ST.setRange(inputArguments.normScanMin, inputArguments.target_STMax + inputArguments.target_binWidth)

def getNormalizedIntegralOfPDFInNamedRange(inputRooPDF, normRangeName, targetRangeName):
    resetSTRange()
    integralObject_normRange = inputRooPDF.createIntegral(ROOT.RooArgSet(rooVar_ST), normRangeName)
    integralObject_targetRange = inputRooPDF.createIntegral(ROOT.RooArgSet(rooVar_ST), targetRangeName)
    normalizedIntegral = integralObject_targetRange.getVal() / integralObject_normRange.getVal()
    resetSTRange()
    return normalizedIntegral

scalingCandidates = {}
scalingCandidates_PDFBased = {}
jetDistributions = {}
constantFitsChi2PerNDF = {}
linearFitsChi2PerNDF = {}
quadraticFitsChi2PerNDF = {}

for nJetsBin in range(inputArguments.nJetsMin, 1 + inputArguments.nJetsMax):
    scalingCandidates[nJetsBin] = {}
    scalingCandidates_PDFBased[nJetsBin] = {}
    jetDistributions[nJetsBin] = {}
    constantFitsChi2PerNDF[nJetsBin] = ROOT.TGraph()
    constantFitsChi2PerNDF[nJetsBin].SetName("constantFitChi2PerNDF_{n}Jets".format(n=nJetsBin))
    constantFitsChi2PerNDF[nJetsBin].SetTitle(";lower edge of norm bin;#chi^{2}/n.d.f.")
    linearFitsChi2PerNDF[nJetsBin] = ROOT.TGraph()
    linearFitsChi2PerNDF[nJetsBin].SetName("linearFitChi2PerNDF_{n}Jets".format(n=nJetsBin))
    linearFitsChi2PerNDF[nJetsBin].SetTitle(";lower edge of norm bin;#chi^{2}/n.d.f.")
    quadraticFitsChi2PerNDF[nJetsBin] = ROOT.TGraph()
    quadraticFitsChi2PerNDF[nJetsBin].SetName("quadraticFitChi2PerNDF_{n}Jets".format(n=nJetsBin))
    quadraticFitsChi2PerNDF[nJetsBin].SetTitle(";lower edge of norm bin;#chi^{2}/n.d.f.")
    # ratioLinearToConstant = ROOT.TGraph()
    # ratioLinearToConstant.SetName("ratioLinearToConstant")
    # ratioLinearToConstant.SetTitle(";lower edge of norm bin;(#chi^{2}/n.d.f.)/(#chi^{2}/n.d.f.)_{constant fit}")
    # ratioQuadraticToConstant = ROOT.TGraph()
    # ratioQuadraticToConstant.SetName("ratioQuadraticToConstant")
    # ratioQuadraticToConstant.SetTitle(";lower edge of norm bin;(#chi^{2}/n.d.f.)/(#chi^{2}/n.d.f.)_{constant fit}")
    # ratioQuadraticToLinear = ROOT.TGraph()
    # ratioQuadraticToLinear.SetName("ratioQuadraticToLinear")
    # ratioQuadraticToLinear.SetTitle(";lower edge of norm bin;(#chi^{2}/n.d.f.)/(#chi^{2}/n.d.f.)_{constant fit}")
    # ratioConstantToConstant = ROOT.TGraph() # lol
    # ratioConstantToConstant.SetName("ratioConstantToConstant")
    # ratioConstantToConstant.SetTitle(";lower edge of norm bin;(#chi^{2}/n.d.f.)/(#chi^{2}/n.d.f.)_{constant fit}")
    # globalFitChi2PerNDF = ROOT.TGraph()
    # globalFitChi2PerNDF.SetName("globalFitChi2PerNDF")
    # globalFitChi2PerNDF.SetTitle(";lower edge of norm bin;Best ST range")
    for scalingHistogramIndex in range(0, inputArguments.nScanPoints):
        candidateNormMin = inputArguments.normScanMin + (scalingHistogramIndex)*(inputArguments.normScanMax - inputArguments.normScanMin)/(inputArguments.nScanPoints-1)
        candidateNormMax = candidateNormMin + inputArguments.target_binWidth
        nCandidateBins = int(0.5 + math.ceil((inputArguments.target_STMax - candidateNormMin)/inputArguments.target_binWidth))
        candidateSTMax = candidateNormMin + nCandidateBins*inputArguments.target_binWidth
        print("Generating fits for normMin: {n}, normMax: {x}, nBins: {nCB}, STMax: {cSTM}".format(n=candidateNormMin, x=candidateNormMax, nCB=nCandidateBins, cSTM=candidateSTMax))

        scalingCandidates[nJetsBin][scalingHistogramIndex] = ROOT.TH1F("h_scalingCandidate_{n}Jets_{i}".format(n=nJetsBin, i=scalingHistogramIndex), "Distribution: {n} jets/min jets;ST (GeV);A.U.".format(n=nJetsBin), nCandidateBins, candidateNormMin, candidateSTMax)
        scalingCandidates[nJetsBin][scalingHistogramIndex].Sumw2()
        scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex] = ROOT.TH1F("h_scalingCandidate_{n}Jets_PDFBased_{i}".format(n=nJetsBin, i=scalingHistogramIndex), "Distribution: {n} jets/min jets;ST (GeV);A.U.".format(n=nJetsBin), nCandidateBins, candidateNormMin, candidateSTMax)
        scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex].Sumw2()
        jetDistributions[nJetsBin][scalingHistogramIndex] = ROOT.TH1F("h_jetDistribution_{i}_{n}Jets".format(i=scalingHistogramIndex, n=nJetsBin), "Distribution: {n} jets;ST (GeV);A.U.".format(n=nJetsBin), nCandidateBins, candidateNormMin, candidateSTMax)
        jetDistributions[nJetsBin][scalingHistogramIndex].Sumw2()

        # Fill jet distributions
        for STValue in STValues[nJetsBin]:
            if (STValue > candidateNormMin): jetDistributions[nJetsBin][scalingHistogramIndex].Fill(STValue)

        # Normalization bin
        if (nJetsBin <= 3): continue
        normError = (jetDistributions[nJetsBin][scalingHistogramIndex].GetBinError(1)/jetDistributions[nJetsBin][scalingHistogramIndex].GetBinContent(1)) + (jetDistributions[inputArguments.nJetsNorm][scalingHistogramIndex].GetBinError(1)/jetDistributions[inputArguments.nJetsNorm][scalingHistogramIndex].GetBinContent(1))
        scalingCandidates[nJetsBin][scalingHistogramIndex].SetBinContent(1, 1.)
        scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex].SetBinContent(1, 1.)
        scalingCandidates[nJetsBin][scalingHistogramIndex].SetBinError(1, normError)
        scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex].SetBinError(1, normError)
        resetSTRange()
        rooVar_ST.setRange("scanPoint{i}_bin1".format(i=scalingHistogramIndex), (scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex]).GetXaxis().GetBinLowEdge(1), (scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex]).GetXaxis().GetBinUpEdge(1))
        resetSTRange()
        # Get ratio histogram
        nEffectiveBins = 1 # Normalization bin already set
        for ratioBinIndex in range(2, 1+nCandidateBins):
            resetSTRange()
            rooVar_ST.setRange("scanPoint{i}_bin{j}".format(i=scalingHistogramIndex, j=ratioBinIndex), (scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex]).GetXaxis().GetBinLowEdge(ratioBinIndex), (scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex]).GetXaxis().GetBinUpEdge(ratioBinIndex))
            resetSTRange()
            numerator = jetDistributions[nJetsBin][scalingHistogramIndex].GetBinContent(ratioBinIndex)/jetDistributions[nJetsBin][scalingHistogramIndex].GetBinContent(1)
            numeratorError = jetDistributions[nJetsBin][scalingHistogramIndex].GetBinError(ratioBinIndex)/jetDistributions[nJetsBin][scalingHistogramIndex].GetBinContent(1)
            denominator = jetDistributions[inputArguments.nJetsNorm][scalingHistogramIndex].GetBinContent(ratioBinIndex)/jetDistributions[inputArguments.nJetsNorm][scalingHistogramIndex].GetBinContent(1)
            denominatorError = jetDistributions[inputArguments.nJetsNorm][scalingHistogramIndex].GetBinError(ratioBinIndex)/jetDistributions[inputArguments.nJetsNorm][scalingHistogramIndex].GetBinContent(1)
            numerator_fromPDF = getNormalizedIntegralOfPDFInNamedRange(inputRooPDF=rooKernel_PDF_Estimators[nJetsBin], normRangeName="scanPoint{i}_bin1".format(i=scalingHistogramIndex), targetRangeName="scanPoint{i}_bin{j}".format(i=scalingHistogramIndex, j=ratioBinIndex))
            denominator_fromPDF = getNormalizedIntegralOfPDFInNamedRange(inputRooPDF=rooKernel_PDF_Estimators[inputArguments.nJetsNorm], normRangeName="scanPoint{i}_bin1".format(i=scalingHistogramIndex), targetRangeName="scanPoint{i}_bin{j}".format(i=scalingHistogramIndex, j=ratioBinIndex))
            denominator_toUse = denominator
            if (inputArguments.usePDFDenominator): denominator_toUse = denominator_fromPDF
            if ((numerator > 0) and (denominator > 0)):
                scalingCandidates[nJetsBin][scalingHistogramIndex].SetBinContent(ratioBinIndex, numerator/denominator_toUse)
                scalingCandidates[nJetsBin][scalingHistogramIndex].SetBinError(ratioBinIndex, (numerator/denominator_toUse)*(numeratorError/numerator + denominatorError/denominator_toUse))
                scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex].SetBinContent(ratioBinIndex, numerator_fromPDF/denominator_fromPDF)
                scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex].SetBinError(ratioBinIndex, (numerator_fromPDF/denominator_fromPDF)*(numeratorError/numerator_fromPDF + denominatorError/denominator_fromPDF))
                nEffectiveBins += 1
            else:
                scalingCandidates[nJetsBin][scalingHistogramIndex].SetBinContent(ratioBinIndex, 0.)
                scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex].SetBinContent(ratioBinIndex, 0.)
                scalingCandidates[nJetsBin][scalingHistogramIndex].SetBinError(ratioBinIndex, 0.)
                scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex].SetBinError(ratioBinIndex, 0.)
            # print("At bin index = {i}, bin content = {c}, bin error = {e}".format(i = ratioBinIndex, c = scalingCandidates[nJetsBin][scalingHistogramIndex].GetBinContent(ratioBinIndex), e = scalingCandidates[nJetsBin][scalingHistogramIndex].GetBinError(ratioBinIndex)))

        fitInput = scalingCandidates[nJetsBin][scalingHistogramIndex]
        if (inputArguments.usePDF): fitInput = scalingCandidates_PDFBased[nJetsBin][scalingHistogramIndex]

        # constantCandidate = ROOT.TF1("constantCandidateFit_{i}".format(i=scalingHistogramIndex), constantFunction, candidateNormMin, candidateSTMax, 1)
        # constantCandidate.SetParameter(0, 1.)
        # constantFitResult = fitInput.Fit(constantCandidate, "IEMS")
        # constantFitChi2 = constantFitResult.Chi2()
        # constantFitNDF = constantFitResult.Ndf()
        # if not(constantFitNDF == (nEffectiveBins - 1)): sys.exit("Error in understanding: nEffectiveBins = {n}, NDF for constant fit = {ndf}".format(n=nEffectiveBins, ndf=constantFitNDF))
        # print("Constant fit result: p0 = {c} +/- {e}; chi2 = {chi2}, ndf = {ndf}".format(c = constantFitResult.Parameter(0), e = constantFitResult.ParError(0), chi2=constantFitChi2, ndf=constantFitNDF))
        # constantFitChi2PerNDF.SetPoint(constantFitChi2PerNDF.GetN(), candidateNormMin, constantFitChi2/constantFitNDF)
        # ratioConstantToConstant.SetPoint(ratioConstantToConstant.GetN(), candidateNormMin, 1.)

        linearCandidate = ROOT.TF1("linearCandidateFit_{n}Jets_{i}".format(n=nJetsBin, i=scalingHistogramIndex), linearFunction, candidateNormMin, candidateSTMax, 2)
        linearCandidate.SetParameter(0, 1.)
        linearCandidate.SetParameter(1, 0.)
        linearFitResult = fitInput.Fit(linearCandidate, "IEMS")
        # linearFitChi2 = linearFitResult.Chi2()
        # linearFitNDF = linearFitResult.Ndf()
        # if not(linearFitNDF == (nEffectiveBins - 2)): sys.exit("Error in understanding: nEffectiveBins = {n}, NDF for linear fit = {ndf}".format(n=nEffectiveBins, ndf=linearFitNDF))
        # print("Linear fit result: p0 = {c} +/- {e}; p1 = {c1} +/- {e1}; chi2 = {chi2}, ndf = {ndf}".format(c = linearFitResult.Parameter(0), e = linearFitResult.ParError(0), c1 = linearFitResult.Parameter(1), e1 = linearFitResult.ParError(1), chi2 = linearFitChi2, ndf = linearFitNDF))
        # linearFitChi2PerNDF.SetPoint(linearFitChi2PerNDF.GetN(), candidateNormMin, linearFitChi2/linearFitNDF)
        # ratioLinearToConstant.SetPoint(ratioLinearToConstant.GetN(), candidateNormMin, (linearFitChi2*constantFitNDF)/(constantFitChi2*linearFitNDF))

        # quadraticCandidate = ROOT.TF1("quadraticCandidateFit_{i}".format(i=scalingHistogramIndex), quadraticFunction, candidateNormMin, candidateSTMax, 3)
        # quadraticCandidate.SetParameter(0, 1.)
        # quadraticCandidate.SetParameter(1, 0.)
        # quadraticCandidate.SetParameter(2, 0.)
        # quadraticFitResult = fitInput.Fit(quadraticCandidate, "IEMS")
        # quadraticFitChi2 = quadraticFitResult.Chi2()
        # quadraticFitNDF = quadraticFitResult.Ndf()
        # if not(quadraticFitNDF == (nEffectiveBins - 3)): sys.exit("Error in understanding: nEffectiveBins = {n}, NDF for quadratic fit = {ndf}".format(n=nEffectiveBins, ndf=quadraticFitNDF))
        # print("Quadratic fit result: p0 = {c} +/- {e}; p1 = {c1} +/- {e1}; p2 = {c2} +/- {e2}; chi2 = {chi2}, ndf = {ndf}".format(c = quadraticFitResult.Parameter(0), e = quadraticFitResult.ParError(0), c1 = quadraticFitResult.Parameter(1), e1 = quadraticFitResult.ParError(1), c2 = quadraticFitResult.Parameter(2), e2 = quadraticFitResult.ParError(2), chi2 = quadraticFitChi2, ndf = quadraticFitNDF))
        # quadraticFitChi2PerNDF.SetPoint(quadraticFitChi2PerNDF.GetN(), candidateNormMin, quadraticFitChi2/quadraticFitNDF)
        # ratioQuadraticToConstant.SetPoint(ratioQuadraticToConstant.GetN(), candidateNormMin, (quadraticFitChi2*constantFitNDF)/(constantFitChi2*quadraticFitNDF))
        # ratioQuadraticToLinear.SetPoint(ratioQuadraticToLinear.GetN(), candidateNormMin, (quadraticFitChi2*linearFitNDF)/(linearFitChi2*quadraticFitNDF))

        lineAt1 = ROOT.TLine(candidateNormMin, 1.0, candidateSTMax, 1.0)

        tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[fitInput, # constantCandidate,
                                                       linearCandidate, # quadraticCandidate,
                                                       lineAt1], canvasName="c_{n}Jets_step{i}_fits".format(n=nJetsBin, i=scalingHistogramIndex), outputDocumentName="{oD}/{oFP}_step{i}_{n}Jets_fits".format(oD=inputArguments.outputDirectory, oFP=inputArguments.outputFilePrefix, i=scalingHistogramIndex, n=nJetsBin), enableLogX = False, enableLogY = False, enableLogZ = False, customXRange=None, customYRange=None, customZRange=None)

        # if (candidateNormMin <= inputArguments.maxSTForGlobalFit):
        #     globalCandidate = ROOT.TF1("globalCandidateFit_{i}".format(i=scalingHistogramIndex), globalFunction, candidateNormMin, candidateSTMax, 5)
        #     globalCandidate.SetParameter(0, 1.)
        #     globalCandidate.SetParameter(1, 0.)
        #     globalCandidate.SetParameter(2, 0.)
        #     globalCandidate.SetParameter(3, 1.)
        #     globalCandidate.SetParameter(4, 1000.)
        #     globalFitResult = scalingCandidates[nJetsBin][scalingHistogramIndex].Fit(globalCandidate, "IEMN0S")
        #     globalFitChi2 = globalFitResult.Chi2()
        #     globalFitNDF = globalFitResult.Ndf()
        #     if not(globalFitNDF == (nEffectiveBins - 5)): sys.exit("Error in understanding: nEffectiveBins = {n}, NDF for global fit = {ndf}".format(n=nEffectiveBins, ndf=globalFitNDF))
        #     print("Global fit result: p0 = {c} +/- {e}; p1 = {c1} +/- {e1}; p2 = {c2} +/- {e2}; asymptotic constant = {c3} +/- {e3}; turn-on = {c4} +/- {e4}; chi2 = {chi2}, ndf = {ndf}".format(c = globalFitResult.Parameter(0), e = globalFitResult.ParError(0), c1 = globalFitResult.Parameter(1), e1 = globalFitResult.ParError(1), c2 = globalFitResult.Parameter(2), e2 = globalFitResult.ParError(2), c3 = globalFitResult.Parameter(3), e3 = globalFitResult.ParError(3), c4 = globalFitResult.Parameter(4), e4 = globalFitResult.ParError(4), chi2 = globalFitChi2, ndf = globalFitNDF))
        #     globalFitChi2PerNDF.SetPoint(globalFitChi2PerNDF.GetN(), candidateNormMin, globalFitResult.Parameter(4))

    # constantFitChi2PerNDF.SetLineColor(ROOT.kBlack)
    # ratioConstantToConstant.SetLineColor(ROOT.kBlack)
    # linearFitChi2PerNDF.SetLineColor(ROOT.kBlue)
    # ratioLinearToConstant.SetLineColor(ROOT.kBlue)
    # quadraticFitChi2PerNDF.SetLineColor(ROOT.kRed)
    # ratioQuadraticToConstant.SetLineColor(ROOT.kRed)
    # tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[constantFitChi2PerNDF, linearFitChi2PerNDF, quadraticFitChi2PerNDF], canvasName="c_chi2s", outputDocumentName="{oD}/{oFP}_chi2Values".format(oD=inputArguments.outputDirectory, oFP=inputArguments.outputFilePrefix), customOptStat=0)
    # tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[ratioConstantToConstant, ratioLinearToConstant, ratioQuadraticToConstant], canvasName="c_chi2Ratios", outputDocumentName="{oD}/{oFP}_chi2Ratios".format(oD=inputArguments.outputDirectory, oFP=inputArguments.outputFilePrefix), customOptStat=0, customYRange=[0., 2.])
    # tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[ratioConstantToConstant, ratioQuadraticToLinear], canvasName="c_chi2Ratios_quadraticToLinear", outputDocumentName="{oD}/{oFP}_chi2Ratios_quadraticToLinear".format(oD=inputArguments.outputDirectory, oFP=inputArguments.outputFilePrefix), customOptStat=0, customYRange=[0., 2.])
    # # tmROOTUtils.plotObjectsOnCanvas(listOfObjects=[globalFitChi2PerNDF], canvasName="c_globalFits", outputDocumentName="{oD}/{oFP}_globalFits".format(oD=inputArguments.outputDirectory, oFP=inputArguments.outputFilePrefix), customOptStat=0)

print("All done!")
