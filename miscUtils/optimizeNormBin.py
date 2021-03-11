#!/usr/bin/env python

from __future__ import print_function, division
import os, sys, argparse, pdb, math, json, subprocess, array, scipy.stats
import tmProgressBar, tmGeneralUtils
import stealthEnv, ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

STBoundariesSourceFile = "STRegionBoundaries_normOptimization.dat"

colors = {
    2: ROOT.kBlack,
    3: ROOT.kBlue+2,
    4: ROOT.kRed+1,
    5: ROOT.kGreen+3,
    6: ROOT.kViolet,
}

selection = "singlemedium"
identifier = "MC_GJet17"
year = "2017"
sourceFilePath  = stealthEnv.analysisRoot + "/STDistributions_singlephoton/distributions_{y}_{s}_{i}.root".format(y=year, i=identifier, s=selection)
outputDirectory = stealthEnv.analysisRoot + "/normBinOptimization_singlephoton"
STNormMax = 1450.0
STNormTarget = 1250.0

# selection = "signal_loose"
# identifier = "MC_GJet"
# year = "all"
# sourceFilePath  = stealthEnv.analysisRoot + "/STDistributions_doublephoton/distributions_{y}_{s}_{i}.root".format(y=year, i=identifier, s=selection)
# outputDirectory = stealthEnv.analysisRoot + "/normBinOptimization_doublephoton"
# STNormMax = 1450.0
# STNormTarget = 1250.0

if not(os.path.isdir(outputDirectory)): subprocess.check_call("mkdir -p {oD}".format(oD=outputDirectory), shell=True, executable="/bin/bash")

STRegionBoundariesFileObject = open(STBoundariesSourceFile)
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
STBoundaries.append(3500)
nSTSignalBins = len(STBoundaries) - 2 # First two lines are lower and upper boundaries for the normalization bin
n_STBins = len(STBoundaries) - 1

sourceFile = ROOT.TFile.Open(sourceFilePath, "READ")
if ((sourceFile.IsOpen() == ROOT.kFALSE) or (sourceFile.IsZombie())): sys.exit("ERROR: unable to open file with name {n}".format(n=sourceFilePath))

print("Getting ST datasets for source: {sF}".format(sF=sourceFilePath))
dataSets = {}
STDistributions = {}
for nJetsBin in range(2, 7):
    dataSets[nJetsBin] = ROOT.RooDataSet()
    sourceFile.GetObject("dataSet_{n}Jets".format(n=nJetsBin), dataSets[nJetsBin])
    dataSets[nJetsBin].SetName("dataSet_{n}Jets".format(n=nJetsBin))
    STDistributions[nJetsBin] = ROOT.TH1F()
    sourceFile.GetObject("h_ST_{n}JetsBin".format(n=nJetsBin), STDistributions[nJetsBin])
    STDistributions[nJetsBin].SetName("h_ST_{n}JetsBin".format(n=nJetsBin))

ratioHistograms = {}
isZeroBin = {}
for nJetsBin in range(3, 7):
    ratioHistograms[nJetsBin] = ROOT.TH1F("h_ST_ratio_{n}JetsTo2Jets".format(n=nJetsBin), "ST distribution: {n} Jets/2 Jets;ST".format(n=nJetsBin), n_STBins, array.array('d', STBoundaries))
    isZeroBin[nJetsBin] = {}
    for binCounter in range(1, 1 + ratioHistograms[nJetsBin].GetXaxis().GetNbins()):
        numerator = STDistributions[nJetsBin].GetBinContent(binCounter)
        denominator = STDistributions[2].GetBinContent(binCounter)
        ratio = 1.0
        ratioError = 1.0
        if ((denominator > 0) and (numerator > 0)):
            numeratorError = STDistributions[nJetsBin].GetBinError(binCounter)
            denominatorError = STDistributions[2].GetBinError(binCounter)
            ratio = numerator/denominator
            ratioError = ratio*math.sqrt(pow(numeratorError/numerator, 2) + pow(denominatorError/denominator, 2))
            isZeroBin[nJetsBin][binCounter] = False
        else:
            isZeroBin[nJetsBin][binCounter] = True
        ratioHistograms[nJetsBin].SetBinContent(binCounter, ratio)
        ratioHistograms[nJetsBin].SetBinError(binCounter, ratioError)

# kernels = {}
# shape_cumulatives = {}
# rooSTVar.setBins(nSTBins)
# STFrame = rooSTVar.frame(ROOT.RooFit.Bins(nSTBins), ROOT.RooFit.Title("kernels, cumulative".format(n=nJetsBin)))

# outputCanvas = ROOT.TCanvas("c_cumulativeShapes", "c_cumulativeShapes", 1024, 768)
# for nJetsBin in range(2, 7):
#     dataHist = ROOT.RooDataHist("dataHist_{n}Jets".format(n=nJetsBin), "data, {n} Jets".format(n=nJetsBin), ROOT.RooArgSet(rooSTVar), dataSets[nJetsBin])
#     kernels[nJetsBin] = ROOT.RooKeysPdf("kernelEstimate_{n}Jets".format(n=nJetsBin), "kernelEstimate_{n}Jets".format(n=nJetsBin), rooSTVar, dataSets[nJetsBin], ROOT.RooKeysPdf.MirrorLeft, 1.8)
#     shape_cumulatives[nJetsBin] = kernels[nJetsBin].createCdf(ROOT.RooArgSet(rooSTVar))
#     shape_cumulatives[nJetsBin].plotOn(STFrame, ROOT.RooFit.LineColor(colors[nJetsBin]))
# STFrame.Draw()
# outputCanvas.SaveAs("{oD}/STCumulativeShapes_{y}_{i}_{s}.pdf".format(oD=outputDirectory, y=year, i=identifier, s=selection))

def get_fTest_prob(chi2_1, chi2_2, ndf_1, ndf_2, printVerbose=False):
    d1 = (ndf_1-ndf_2)
    d2 = ndf_2
    N = (chi2_1-chi2_2)/d1
    D = chi2_2/d2
    fStat = N/D
    fProb = scipy.stats.f.cdf(fStat, d1, d2)
    expectedFStat = scipy.stats.distributions.f.isf(0.05, d1, d2)
    if printVerbose:
        print("chi2_1 = {chi2_1}, chi2_2 = {chi2_2}".format(chi2_1=chi2_1, chi2_2=chi2_2))
        print("ndf_1 = {ndf_1}, ndf_2 = {ndf_2}".format(ndf_1=ndf_1, ndf_2=ndf_2))
        print("d1, d2 = %d, %d"%(d1, d2))
        print("N, D = %f, %f"%(N, D))
        print("    f(%i,%i) = %f (expected at 95%%: %f)"%(d1,d2,fStat,expectedFStat))
        print("f.cdf(%i,%i) = %3.0f%%"%(d1,d2,100*fProb))
    return fProb

normBinIndex = 5
STNormTargetBin = STDistributions[2].GetXaxis().FindFixBin(STNormTarget)
fitTypes = ["const", "lin", "quad", "constrained_lin"]
fitFunctions = {
    "const": "pol0",
    "lin": "pol1",
    "quad": "pol2",
    "constrained_lin": "1.0+[0]*(x-{STNT})/1000.".format(STNT=STNormTarget)
}
chiSqPerNDFGraphTitles = {
    "const": "#chi^{2}/NDF, constant fit",
    "lin": "#chi^{2}/NDF, linear fit",
    "quad": "#chi^{2}/NDF, quadratic fit",
    "constrained_lin": "#chi^{2}/NDF, constrained linear fit"
}

chiSqPerNDFGraphs = {}
for fitType in fitTypes:
    chiSqPerNDFGraphs[fitType] = {}
    for nJetsBin in range(3, 7):
        chiSqPerNDFGraphs[fitType][nJetsBin] = ROOT.TGraph()
        chiSqPerNDFGraphs[fitType][nJetsBin].SetName("chiSqPerNDFs_fitType{t}_{n}JetsTo2Jets".format(t=fitType, n=nJetsBin))

fitParametersList = []
fTestProbValues = {"const_vs_lin": {}, "lin_vs_quad": {}, "constrained_lin_vs_lin": {}}
while True:
    STNorm = STDistributions[2].GetXaxis().GetBinCenter(normBinIndex)
    saveRatios = False
    if (STDistributions[2].GetXaxis().FindFixBin(STNorm) == STNormTargetBin): saveRatios = True
    if (STNorm > STNormMax): break
    print ("Trying ST norm: {n}".format(n=STNorm))

    histogramCopy_2Jets = ROOT.TH1F("h_ST_2JetsBin_normBin{i}".format(i=normBinIndex), "ST distribution: 2 Jets;ST", n_STBins-normBinIndex, array.array('d', STBoundaries[normBinIndex:]))
    histogramCopy_2Jets.Sumw2()
    for binCounter in range(1, 1 + histogramCopy_2Jets.GetXaxis().GetNbins()):
        binCenter = histogramCopy_2Jets.GetXaxis().GetBinCenter(binCounter)
        original_binIndex = STDistributions[2].FindFixBin(binCenter)
        histogramCopy_2Jets.SetBinContent(binCounter, STDistributions[2].GetBinContent(original_binIndex))
        histogramCopy_2Jets.SetBinError(binCounter, STDistributions[2].GetBinError(original_binIndex))
    # chisqPerNDF_avg = 0.
    for nJetsBin in range(3, 7):
        histogramCopy = ROOT.TH1F("h_ST_{n}JetsBin_normBin{i}".format(n=nJetsBin, i=normBinIndex), "ST distribution: {n} Jets;ST".format(n=nJetsBin), n_STBins-normBinIndex, array.array('d', STBoundaries[normBinIndex:]))
        histogramCopy.Sumw2()
        ratioGraph = ROOT.TGraphErrors()
        ratioGraph.SetName("ratioGraph_{n}JetsTo2Jets_normBin{i}".format(n=nJetsBin, i=normBinIndex))
        ratioGraphPoints = []
        for binCounter in range(1, 1 + histogramCopy.GetXaxis().GetNbins()):
            binCenter = histogramCopy.GetXaxis().GetBinCenter(binCounter)
            original_binIndex = STDistributions[nJetsBin].FindFixBin(binCenter)
            histogramCopy.SetBinContent(binCounter, STDistributions[nJetsBin].GetBinContent(original_binIndex))
            histogramCopy.SetBinError(binCounter, STDistributions[nJetsBin].GetBinError(original_binIndex))
            numerator = STDistributions[nJetsBin].GetBinContent(original_binIndex)/STDistributions[nJetsBin].GetBinContent(normBinIndex)
            denominator = STDistributions[2].GetBinContent(original_binIndex)/STDistributions[2].GetBinContent(normBinIndex)
            if ((numerator > 0) and (denominator > 0)):
                fracError_numerator_thisBin = STDistributions[nJetsBin].GetBinError(original_binIndex)/STDistributions[nJetsBin].GetBinContent(original_binIndex)
                fracError_numerator_normBin = STDistributions[nJetsBin].GetBinError(normBinIndex)/STDistributions[nJetsBin].GetBinContent(normBinIndex)

                fracError_denominator_thisBin = STDistributions[2].GetBinError(original_binIndex)/STDistributions[2].GetBinContent(original_binIndex)
                fracError_denominator_normBin = STDistributions[2].GetBinError(normBinIndex)/STDistributions[2].GetBinContent(normBinIndex)

                ratio = numerator/denominator
                ratioError = ratio*math.sqrt(pow(fracError_numerator_thisBin, 2) + pow(fracError_numerator_normBin, 2) + pow(fracError_denominator_thisBin, 2) + pow(fracError_denominator_normBin, 2))
                ratioGraphIndex = ratioGraph.GetN()
                ratioGraph.SetPoint(ratioGraphIndex, binCenter, ratio)
                ratioGraph.SetPointError(ratioGraphIndex, 0.5*(histogramCopy.GetXaxis().GetBinUpEdge(binCounter) - histogramCopy.GetXaxis().GetBinLowEdge(binCounter)), ratioError)
                ratioGraphPoints.append(tuple([binCenter, ratio, 0.5*(histogramCopy.GetXaxis().GetBinUpEdge(binCounter) - histogramCopy.GetXaxis().GetBinLowEdge(binCounter)), ratioError]))

        fits = {}
        fitResults = {}
        for fitType in fitTypes:
            fits[fitType] = ROOT.TF1("fit_type{t}_{n}JetsBin_normBin{i}".format(t=fitType, n=nJetsBin, i=normBinIndex), fitFunctions[fitType], histogramCopy.GetXaxis().GetBinLowEdge(2), histogramCopy.GetXaxis().GetBinUpEdge(histogramCopy.GetXaxis().GetNbins()))
            if (fitType == "constrained_lin"):
                fits[fitType].SetParName(0, "slopeDividedBy1000")
                fits[fitType].SetParameter(0, 0.)
                fits[fitType].SetParLimits(0, -0.5, 5.)
                fits[fitType].SetLineColor(colors[nJetsBin])
                fits[fitType].SetLineStyle(ROOT.kDashed)
            fitResults[fitType] = ratioGraph.Fit("fit_type{t}_{n}JetsBin_normBin{i}".format(t=fitType, n=nJetsBin, i=normBinIndex), "EX0QREMS+")
            try:
                chiSqPerNDF = fitResults[fitType].Chi2()/fitResults[fitType].Ndf()
                chiSqPerNDFGraphs[fitType][nJetsBin].SetPoint(chiSqPerNDFGraphs[fitType][nJetsBin].GetN(), STNorm, chiSqPerNDF)
            except ZeroDivisionError:
                print("WARNING: Zero division error for nJetsBin: {n}, fitType: {t}, STNorm: {STN}".format(n=nJetsBin, t=fitType, STN=STNorm))
                continue
            if (saveRatios and (fitType == "constrained_lin")):
                bestFitSlopeDividedBy1000 = fitResults[fitType].Parameter(0)
                fitParametersList.append(tuple(["float", "bestFitSlopeDividedBy1000_{n}Jets".format(n=nJetsBin), bestFitSlopeDividedBy1000]))
                bestFitSlopeErrorDividedBy1000 = fitResults[fitType].ParError(0)
                fitParametersList.append(tuple(["float", "bestFitSlopeErrorDividedBy1000_{n}Jets".format(n=nJetsBin), bestFitSlopeErrorDividedBy1000]))

        if saveRatios:
            for fitIndex in range(len(fitTypes)):
                fitType = fitTypes[fitIndex]
                outputCanvas = ROOT.TCanvas("c_residuals_{t}_{n}JetsBin".format(t=fitType, n=nJetsBin), "c_residuals_{t}_{n}JetsBin".format(t=fitType, n=nJetsBin), 1024, 256)
                residualsGraph = ROOT.TGraphErrors()
                residualsGraph.SetName("residuals_{t}_{n}JetsBin".format(t=fitType, n=nJetsBin))
                for index in range(0, ratioGraph.GetN()):
                    (binCenter, ratioValue, binWidth, ratioError) = ratioGraphPoints[index]
                    fitFunctionValue = 1.
                    if (index > 0): fitFunctionValue = fits[fitType].Eval(binCenter)
                    residual = ratioValue - fitFunctionValue
                    residualsGraph.SetPoint(index, binCenter, residual)
                    residualsGraph.SetPointError(index, binWidth, ratioError)
                residualsGraph.Draw("AP")
                lineAt0 = ROOT.TLine(histogramCopy.GetXaxis().GetBinLowEdge(2), 0., histogramCopy.GetXaxis().GetBinUpEdge(histogramCopy.GetXaxis().GetNbins()), 0.)
                lineAt0.SetLineStyle(ROOT.kSolid)
                lineAt0.Draw()    
                residualsGraph.SetTitle("Residuals for fitType: {t}".format(t=fitType))
                ROOT.gPad.Update()
                outputCanvas.SaveAs("{oD}/residuals_targetNorm_fitType_{t}_{n}JetsBin_{y}_{i}_{s}.pdf".format(oD=outputDirectory, t=fitType, n=nJetsBin, y=year, i=identifier, s=selection))

            outputCanvas = ROOT.TCanvas("c_ratios_{n}JetsBin".format(n=nJetsBin), "c_ratios_{n}JetsBin".format(n=nJetsBin), 1024, 768)
            ratioGraph.SetLineColor(colors[nJetsBin])
            ratioGraph.Draw("AP")
            outputCanvas.SaveAs("{oD}/ratios_targetNorm_{n}JetsBin_{y}_{i}_{s}.pdf".format(oD=outputDirectory, n=nJetsBin, y=year, i=identifier, s=selection))
            print("Getting f-test values for nJetsBin = {n}".format(n=nJetsBin))
            fTestProbValues["const_vs_lin"][nJetsBin] = get_fTest_prob(chi2_1=fitResults["const"].Chi2(), chi2_2=fitResults["lin"].Chi2(), ndf_1=fitResults["const"].Ndf(), ndf_2=fitResults["lin"].Ndf())
            fTestProbValues["lin_vs_quad"][nJetsBin] = get_fTest_prob(chi2_1=fitResults["lin"].Chi2(), chi2_2=fitResults["quad"].Chi2(), ndf_1=fitResults["lin"].Ndf(), ndf_2=fitResults["quad"].Ndf(), printVerbose=True)
            fTestProbValues["constrained_lin_vs_lin"][nJetsBin] = get_fTest_prob(chi2_1=fitResults["constrained_lin"].Chi2(), chi2_2=fitResults["lin"].Chi2(), ndf_1=fitResults["constrained_lin"].Ndf(), ndf_2=fitResults["lin"].Ndf())

        # comparisonType = "UU"
        # if getMCWeights: comparisonType = "WW"
        # comparisonType = "WW"
        # chi2 = histogramCopy.Chi2Test(histogramCopy_2Jets, "{cT} P CHI2/NDF".format(cT=comparisonType))
        # chisqPerNDF_avg += chi2/4 # four nJets bins
    # chisqPerNDFGraph.SetPoint(chisqPerNDFGraph.GetN(), STNorm, chisqPerNDF_avg)
    normBinIndex += 1
# chisqPerNDFGraph.SetMarkerStyle(ROOT.kCircle)
# chisqPerNDFGraph.SetMarkerSize(0.5)
# chisqPerNDFGraph.Draw("AP")
# chisqPerNDFGraph.GetXaxis().SetTitle("ST norm")
# chisqPerNDFGraph.GetYaxis().SetTitle("#chi^{2}/ndf")
# outputCanvas.Update()
# outputCanvas.SaveAs("{oD}/chiSqPerNDF_{y}_{i}_{s}.pdf".format(oD=outputDirectory, y=year, i=identifier, s=selection))
configurationParametersOutputFilePath = "{oD}/fitParameters_{y}_{i}_{s}.dat".format(oD=outputDirectory, y=year, i=identifier, s=selection)
tmGeneralUtils.writeConfigurationParametersToFile(configurationParametersList=fitParametersList, outputFilePath=configurationParametersOutputFilePath)
print("Configuration parameters written to file: {cPOFP}".format(cPOFP=configurationParametersOutputFilePath))

# Print f-test prob values in a LaTeX-formatted table
print("f-test prob values")
tmGeneralUtils.prettyPrintDictionary(inputDict=fTestProbValues)
print("f-test prob values (formatted):")
print("\\begin{tabular}{|l|c|c|c|}")
print("  \\hline")
print("  f-prob & const\\_vs\\_lin & lin\\_vs\\_quad & fixed\\_lin\\_vs\\_lin \\\\ \\hline")
print("  nJets = 3 & {fcl:.3f} & {flq:.3f} & {fcll:.3f} \\\\ \\hline".format(fcl=fTestProbValues["const_vs_lin"][3], flq=fTestProbValues["lin_vs_quad"][3], fcll=fTestProbValues["constrained_lin_vs_lin"][3]))
print("  nJets = 4 & {fcl:.3f} & {flq:.3f} & {fcll:.3f} \\\\ \\hline".format(fcl=fTestProbValues["const_vs_lin"][4], flq=fTestProbValues["lin_vs_quad"][4], fcll=fTestProbValues["constrained_lin_vs_lin"][4]))
print("  nJets = 5 & {fcl:.3f} & {flq:.3f} & {fcll:.3f} \\\\ \\hline".format(fcl=fTestProbValues["const_vs_lin"][5], flq=fTestProbValues["lin_vs_quad"][5], fcll=fTestProbValues["constrained_lin_vs_lin"][5]))
print("  nJets $\\geq$ 6 & {fcl:.3f} & {flq:.3f} & {fcll:.3f} \\\\ \\hline".format(fcl=fTestProbValues["const_vs_lin"][6], flq=fTestProbValues["lin_vs_quad"][6], fcll=fTestProbValues["constrained_lin_vs_lin"][6]))
print("\\end{tabular}")

# Print best fit slopes in LaTeX-formatted table
print("slope values (formatted):")
print("\\begin{tabular}{|l|c|}")
print("  \\hline")
print("  nJets Bin & Best-fit slope \\\\ \\hline")
print("  nJets = 3 & {s:.3f} $\\pm$ {se:.3f} \\\\ \\hline".format(s=fitParametersList[0][2], se=fitParametersList[1][2]))
print("  nJets = 4 & {s:.3f} $\\pm$ {se:.3f} \\\\ \\hline".format(s=fitParametersList[2][2], se=fitParametersList[3][2]))
print("  nJets = 5 & {s:.3f} $\\pm$ {se:.3f} \\\\ \\hline".format(s=fitParametersList[4][2], se=fitParametersList[5][2]))
print("  nJets $\\geq$ 6 & {s:.3f} $\\pm$ {se:.3f} \\\\ \\hline".format(s=fitParametersList[6][2], se=fitParametersList[7][2]))
print("\\end{tabular}")

for fitType in fitTypes:
    outputCanvas = ROOT.TCanvas("c_chiSqPerNDFs_fitType{t}".format(t=fitType), "c_chiSqPerNDFs_fitType{t}".format(t=fitType), 1024, 768)
    chiSqPerNDFsMultigraph = ROOT.TMultiGraph()
    for nJetsBin in range(3, 7):
        chiSqPerNDFGraphs[fitType][nJetsBin].SetLineColor(colors[nJetsBin])
        chiSqPerNDFGraphs[fitType][nJetsBin].SetMarkerStyle(ROOT.kCircle)
        chiSqPerNDFGraphs[fitType][nJetsBin].SetMarkerSize(0.75)
        chiSqPerNDFsMultigraph.Add(chiSqPerNDFGraphs[fitType][nJetsBin])
    chiSqPerNDFsMultigraph.Draw("APL")
    chiSqPerNDFsMultigraph.GetXaxis().SetTitle("ST norm")
    chiSqPerNDFsMultigraph.GetYaxis().SetTitle("#chi^{2}/NDF")
    chiSqPerNDFsMultigraph.SetTitle(chiSqPerNDFGraphTitles[fitType])
    outputCanvas.Update()
    outputCanvas.SaveAs("{oD}/chiSqPerNDFs_fitType_{t}_{y}_{i}_{s}.pdf".format(oD=outputDirectory, t=fitType, y=year, i=identifier, s=selection))

outputCanvas = ROOT.TCanvas("c_STDistributions", "c_STDistributions", 1024, 768)
legend = ROOT.TLegend(0.7, 0.6, 0.9, 0.9)
STDistributions[2].SetLineColor(colors[2])
STDistributions[2].Draw()
legendEntry = legend.AddEntry(STDistributions[2], "nJets = {nJetsBin}".format(nJetsBin=2))
legendEntry.SetMarkerColor(colors[2])
legendEntry.SetLineColor(colors[2])
legendEntry.SetTextColor(colors[2])
for nJetsBin in range(3, 7):
    STDistributions[nJetsBin].Scale(STDistributions[2].GetBinContent(STDistributions[2].FindFixBin(STNormTarget))/STDistributions[nJetsBin].GetBinContent(STDistributions[nJetsBin].FindFixBin(STNormTarget)))
    STDistributions[nJetsBin].SetLineColor(colors[nJetsBin])
    STDistributions[nJetsBin].Draw("SAME")
    legendText = "nJets = {nJetsBin}".format(nJetsBin=nJetsBin)
    if (nJetsBin == 6): legendText = "nJets #geq 6"
    legendEntry = legend.AddEntry(STDistributions[nJetsBin], legendText)
    legendEntry.SetMarkerColor(colors[nJetsBin])
    legendEntry.SetLineColor(colors[nJetsBin])
    legendEntry.SetTextColor(colors[nJetsBin])
legend.Draw()
ROOT.gPad.SetLogy()
ROOT.gStyle.SetOptStat(0)
outputCanvas.Update()
outputCanvas.SaveAs("{oD}/STDistributions_{y}_{i}_{s}.pdf".format(oD=outputDirectory, y=year, i=identifier, s=selection))

print("All done!")
