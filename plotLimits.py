#!/usr/bin/env python

from __future__ import print_function, division

import argparse, pdb, sys, math, array, os, subprocess
import ROOT, tmROOTUtils, tmGeneralUtils, tdrstyle, CMS_lumi, MCTemplateReader, stealthEnv, commonFunctions

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)
SIGNAL_CONTAMINATION_THRESHOLD = 0.1

inputArgumentsParser = argparse.ArgumentParser(description='Store expected and observed limits on signal strength and cross-section.')
inputArgumentsParser.add_argument('--crossSectionsFile', required=True, help='Path to dat file that contains cross-sections as a function of eventProgenitor mass, to use while weighting events.',type=str)
inputArgumentsParser.add_argument('--MCTemplatePath', required=True, help='Path to MC template.', type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArgumentsParser.add_argument('--eventProgenitor', required=True, help="Type of stealth sample. Two possible values: \"squark\" or \"gluino\".", type=str)
inputArgumentsParser.add_argument('--combineResultsDirectory', required=True, help='EOS path at which combine tool results can be found.',type=str)
inputArgumentsParser.add_argument('--combineOutputPrefix', default="fullChain", help='Prefix of Higgs combine results.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_rawOutput', default="limits", help='Output directory in which to store raw outputs.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_plots', default="publicationPlots", help='Output directory in which to store plots.',type=str)
inputArgumentsParser.add_argument('--outputSuffix', default="fullChain", help='Suffix to append to all results.',type=str)
inputArgumentsParser.add_argument('--maxAllowedRatio', default=10., help='Max allowed ratio for deviation between expected and observed limits.',type=float)
inputArgumentsParser.add_argument('--minNeutralinoMass', default=-1., help='Min value of the neutralino mass to plot.',type=float)
inputArgumentsParser.add_argument('--eventProgenitorMassOffset', default=-1., help='Min value of the event progenitor mass to plot is obtained by adding this offset to the template.',type=float)
inputArgumentsParser.add_argument('--minMassDifference', default=-1., help='Min difference between the masses of the event progenitor and neutralino.',type=float)
inputArgumentsParser.add_argument('--contour_signalStrength', default=1., help='Signal strength at which to obtain the contours.',type=float)
inputArgumentsParser.add_argument('--selectionsList', default="signal,loose_signal", help="Comma-separated list of selections, used to extract names of rate parameters.", type=str)
inputArgumentsParser.add_argument('--signalContaminationSource_signal', required=True, help="Path to source file for signal contamination histograms, signal selection.", type=str)
inputArgumentsParser.add_argument('--signalContaminationSource_signal_loose', required=True, help="Path to source file for signal contamination histograms, loose signal selection.", type=str)
inputArgumentsParser.add_argument('--signalContaminationMonitor_source_folder_eos', required=True, help="Path to EOS source file for signal contamination monitoring data files.", type=str)
inputArgumentsParser.add_argument('--plotObserved', action='store_true', help="If this flag is set, then the observed limits are plotted in addition to the expected limits.")
inputArguments = inputArgumentsParser.parse_args()

string_gluino = "#tilde{g}"
string_mass_gluino = "m_{" + string_gluino + "}"
string_squark = "#tilde{q}"
string_mass_squark = "m_{" + string_squark + "}"
string_neutralino = "#tilde{#chi}_{1}^{0}"
string_mass_neutralino = "m_{" + string_neutralino + "}"
string_singlino = "#tilde{S}"
string_mass_singlino = "m_{" + string_singlino + "}"
string_singlet = "S"
string_mass_singlet = "m_{" + string_singlet + "}"
string_gravitino = "#tilde{G}"
string_mass_gravitino = "m_{" + string_gravitino + "}"
string_photon = "#gamma"
string_eventProgenitor = None
if (inputArguments.eventProgenitor == "gluino"):
    string_eventProgenitor = string_gluino
elif (inputArguments.eventProgenitor == "squark"):
    string_eventProgenitor = string_squark
string_mass_eventProgenitor = "m_{" + string_eventProgenitor + "}"

if (inputArguments.eventProgenitor == "gluino"):
    decayChain = "pp#rightarrow" + string_gluino + string_gluino + ", " + string_gluino + "#rightarrow" + string_squark + "q, " + string_squark + "#rightarrow" + string_neutralino + "q, " + string_neutralino + "#rightarrow" + string_photon + string_singlino + ", " + string_singlino + "#rightarrow" + string_singlet + string_gravitino + ", " + string_singlet + "#rightarrowgg"
else:
    decayChain = "pp#rightarrow" + string_squark + string_squark + ", " + string_squark + "#rightarrow" + string_neutralino + "q, " + string_neutralino + "#rightarrow" + string_photon + string_singlino + ", " + string_singlino + "#rightarrow" + string_singlet + string_gravitino + ", " + string_singlet + "#rightarrowgg"
decayChain_supplementaryInfo1 = "(" + string_mass_singlino + " = 100 GeV, " + string_mass_singlet + " = 90 GeV, " + string_mass_gravitino + " = 0)"
decayChain_supplementaryInfo2 = "NNLO + NNLL exclusion"

selectionsToUse = []
for selection in ((inputArguments.selectionsList).strip()).split(","):
    if not(selection in ["signal", "signal_loose", "control"]): sys.exit("ERROR: Unrecognized region to use: {r}".format(r=selection))
    selectionsToUse.append(selection)

signalContaminationSourceFilePaths = {
    "signal": inputArguments.signalContaminationSource_signal,
    "signal_loose": inputArguments.signalContaminationSource_signal_loose
}
signalContaminationMonitoredQuantityLabels = ["fractionalSignalCorrection", "signalCorrectionOverBackground", "signalCorrectionSignificance", "signalCorrectionNormTermsOverFull"]
diagonal_down_shift = int(0.5 + inputArguments.minMassDifference)

STRegionBoundariesFileObject = open(inputArguments.inputFile_STRegionBoundaries, 'r')
nSTBoundaries = 0
STBoundaries = []
for STBoundaryString in STRegionBoundariesFileObject:
    if (STBoundaryString.strip()):
        nSTBoundaries += 1
        STBoundary = float(STBoundaryString.strip())
        STBoundaries.append(STBoundary)
nSTSignalBins = nSTBoundaries - 2 + 1 # First two lines are for the normalization bin, last boundary is at 3500
print("Using {n} signal bins for ST.".format(n = nSTSignalBins))
STRegionBoundariesFileObject.close()
STRegionTitles = {}
for STRegionIndex in range(1, nSTBoundaries):
    STRegionTitles[STRegionIndex] = "{l:.1f} < ST < {h:.1f}".format(l=STBoundaries[STRegionIndex-1], h=STBoundaries[STRegionIndex])
STRegionTitles[nSTBoundaries] = "ST > {b:.1f}".format(b=STBoundaries[nSTBoundaries-1])

def formatContours(contoursList, lineStyle, lineWidth, lineColor):
    contoursListIteratorNext = ROOT.TIter(contoursList)
    counter = 0
    while True:
        counter += 1
        contour = contoursListIteratorNext()
        if not(contour): break
        contour.SetLineStyle(lineStyle)
        contour.SetLineWidth(lineWidth)
        contour.SetLineColor(lineColor)

def drawContoursForLegend(ndc_x, delta_ndcx, ndc_y, delta_ndcy, lineWidth_middle, lineStyle_middle, lineColor_middle, lineWidth_topBottom, lineStyle_topBottom, lineColor_topBottom, middleLine_fudgeFactor):
    line_middle = ROOT.TLine()
    line_middle.SetLineWidth(lineWidth_middle)
    line_middle.SetLineStyle(lineStyle_middle)
    line_middle.SetLineColor(lineColor_middle)
    line_topBottom = ROOT.TLine()
    line_topBottom.SetLineWidth(lineWidth_topBottom)
    line_topBottom.SetLineStyle(lineStyle_topBottom)
    line_topBottom.SetLineColor(lineColor_topBottom)
    
    line_middle.DrawLineNDC(ndc_x, ndc_y+middleLine_fudgeFactor*delta_ndcy, ndc_x+delta_ndcx, ndc_y+middleLine_fudgeFactor*delta_ndcy) # middle
    line_topBottom.DrawLineNDC(ndc_x, ndc_y-0.5*delta_ndcy, ndc_x+delta_ndcx, ndc_y-0.5*delta_ndcy) # bottom
    line_topBottom.DrawLineNDC(ndc_x, ndc_y+1.5*delta_ndcy, ndc_x+delta_ndcx, ndc_y+1.5*delta_ndcy) # top

def getRGB(color):
    colorObject = ROOT.gROOT.GetColor(color)
    outputDictionary = {
        "red": colorObject.GetRed(),
        "green": colorObject.GetGreen(),
        "blue": colorObject.GetBlue()
    }
    return outputDictionary

def get_max_signal_contamination(eventProgenitorMass, neutralinoMass, inputHistogramsList, printDebug=False):
    signal_contamination_values = []
    for inputHistogram in inputHistogramsList:
        signal_contamination = inputHistogram.GetBinContent(inputHistogram.FindFixBin(eventProgenitorMass, neutralinoMass))
        if printDebug: print("from " + inputHistogram.GetName() + ": {s:.5f}".format(s=signal_contamination))
        signal_contamination_values.append(signal_contamination)
    return max(signal_contamination_values)

def get_signal_contamination_monitored_quantities(eventProgenitorMassBin, neutralinoMassBin):
    # Step 1: Copy file to tmp area
    tmp_output_directory = stealthEnv.scratchArea + "/plotLimits"
    monitor_file_name = "{cop}_signal_contamination_monitor_eventProgenitorMassBin{gBI}_neutralinoMassBin{nBI}.txt".format(cop=inputArguments.combineOutputPrefix, gBI=eventProgenitorMassBin, nBI=neutralinoMassBin)
    if not(os.path.isdir(tmp_output_directory)): subprocess.check_call("mkdir -p {oD}".format(oD=tmp_output_directory), shell=True, executable="/bin/bash")
    subprocess.check_call("xrdcp --force --silent --nopbar --streams 15 {i}/{mfn} {tod}/{mfn}".format(i="{p}/{inputF}".format(p=stealthEnv.EOSPrefix, inputF=inputArguments.signalContaminationMonitor_source_folder_eos), mfn=monitor_file_name, tod=tmp_output_directory), shell=True, executable="/bin/bash")

    # Step 2: Fetch parameters
    monitored_quantities = {}
    monitor_file_contents = tmGeneralUtils.getConfigurationFromFile(inputFilePath="{tod}/{mfn}".format(tod=tmp_output_directory, mfn=monitor_file_name))
    for label in signalContaminationMonitoredQuantityLabels:
        monitored_quantities[label] = {}
        for selection in selectionsToUse:
            monitored_quantities[label][selection] = {}
            for nJetsBin in range(4, 7):
                monitored_quantities[label][selection][nJetsBin] = {}
                for STRegionIndex in range(2, 8):
                    monitored_quantities[label][selection][nJetsBin][STRegionIndex] = monitor_file_contents["{l}_{s}_STRegion{r}_{n}Jets".format(l=label, s=selection, r=STRegionIndex, n=nJetsBin)]

    # Step 3: Remove tmp file from scratch area
    subprocess.check_call("rm -f {tod}/{mfn}".format(tod=tmp_output_directory, mfn=monitor_file_name), shell=True, executable="/bin/bash")
    return monitored_quantities

crossSectionsInputFileObject = open(inputArguments.crossSectionsFile, 'r')
crossSectionsDictionary = {}
crossSectionsFractionalUncertaintyDictionary = {}
for line in crossSectionsInputFileObject:
    crossSectionsData = line.split()
    eventProgenitorMass = int(0.5 + float(crossSectionsData[0]))
    crossSection = float(crossSectionsData[1])
    crossSectionFractionalUncertainty = 0.01*float(crossSectionsData[2])
    crossSectionsDictionary[eventProgenitorMass] = crossSection
    crossSectionsFractionalUncertaintyDictionary[eventProgenitorMass] = crossSectionFractionalUncertainty
crossSectionsInputFileObject.close()

templateReader = MCTemplateReader.MCTemplateReader(inputArguments.MCTemplatePath)
minEventProgenitorMass = (templateReader.minEventProgenitorMass + inputArguments.eventProgenitorMassOffset)
maxEventProgenitorMass = templateReader.maxEventProgenitorMass
minNeutralinoMass = inputArguments.minNeutralinoMass
maxNeutralinoMass = templateReader.maxNeutralinoMass

signalContaminationSourceFileHandles = {}
signalContaminationHistograms_input = {}
signalContaminationHistograms_cleaned = {}
for selection in selectionsToUse:
    signalContaminationSourceFileHandles[selection] = ROOT.TFile.Open(signalContaminationSourceFilePaths[selection], "READ")
    if (((signalContaminationSourceFileHandles[selection]).IsOpen() == ROOT.kFALSE) or ((signalContaminationSourceFileHandles[selection]).IsZombie())): sys.exit("ERROR: unable to open file at path {p}".format(p=signalContaminationSourceFilePaths[selection]))
    signalContaminationHistograms_input[selection] = {}
    signalContaminationHistograms_cleaned[selection] = {}
    for nJetsBin in range(2, 7):
        signalContaminationHistograms_input[selection][nJetsBin] = {}
        signalContaminationHistograms_cleaned[selection][nJetsBin] = {}
        STRegionsToFetch = None
        if (nJetsBin == 2): STRegionsToFetch = range(1, 8)
        elif (nJetsBin == 3): STRegionsToFetch = []
        else: STRegionsToFetch = [1]
        for STRegionIndex in STRegionsToFetch:
            localSignalBinLabel = "STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin)
            signalContaminationHistograms_input[selection][nJetsBin][STRegionIndex] = ROOT.TH2F()
            (signalContaminationSourceFileHandles[selection]).GetObject("h_signalContamination_{sBL}".format(sBL=localSignalBinLabel), signalContaminationHistograms_input[selection][nJetsBin][STRegionIndex])
            signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex] = ROOT.TH2F("h_signalContamination_cleaned_{sBL}".format(sBL=localSignalBinLabel),
                                                                                                  "h_signalContamination_cleaned_{sBL}".format(sBL=localSignalBinLabel),
                                                                                                  (signalContaminationHistograms_input[selection][nJetsBin][STRegionIndex]).GetXaxis().GetNbins(),
                                                                                                  (signalContaminationHistograms_input[selection][nJetsBin][STRegionIndex]).GetXaxis().GetXmin(),
                                                                                                  (signalContaminationHistograms_input[selection][nJetsBin][STRegionIndex]).GetXaxis().GetXmax(),
                                                                                                  (signalContaminationHistograms_input[selection][nJetsBin][STRegionIndex]).GetYaxis().GetNbins(),
                                                                                                  (signalContaminationHistograms_input[selection][nJetsBin][STRegionIndex]).GetYaxis().GetXmin(),
                                                                                                  (signalContaminationHistograms_input[selection][nJetsBin][STRegionIndex]).GetYaxis().GetXmax())
            (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).SetTitle((signalContaminationHistograms_input[selection][nJetsBin][STRegionIndex]).GetTitle())
signal_contamination_monitor_histograms = {}
for label in signalContaminationMonitoredQuantityLabels:
    signal_contamination_monitor_histograms[label] = {}
    for selection in selectionsToUse:
        signal_contamination_monitor_histograms[label][selection] = {}
        for nJetsBin in range(4, 7):
            signal_contamination_monitor_histograms[label][selection][nJetsBin] = {}
            for STRegionIndex in range(2, 8):
                signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex] = ROOT.TH2F("h_{l}_{s}_STRegion{r}_{n}JetsBin".format(l=label, s=selection, r=STRegionIndex, n=nJetsBin),
                                                                                                               "{l}, selection {s}, ST region {r}, {n} jets bin".format(l=label, s=selection, r=STRegionIndex, n=nJetsBin),
                                                                                                               templateReader.nEventProgenitorMassBins,
                                                                                                               templateReader.minEventProgenitorMass,
                                                                                                               templateReader.maxEventProgenitorMass,
                                                                                                               templateReader.nNeutralinoMassBins,
                                                                                                               templateReader.minNeutralinoMass,
                                                                                                               templateReader.maxNeutralinoMass)

limitsScanExpected=ROOT.TGraph2D()
limitsScanExpectedOneSigmaDown=ROOT.TGraph2D()
limitsScanExpectedOneSigmaUp=ROOT.TGraph2D()
crossSectionScanExpected = ROOT.TGraph2D()
limitsScanExpected.SetName("limitsScanExpected")
limitsScanExpectedOneSigmaDown.SetName("limitsScanExpectedOneSigmaDown")
limitsScanExpectedOneSigmaUp.SetName("limitsScanExpectedOneSigmaUp")
crossSectionScanExpected.SetName("crossSectionScanExpected")

limitsScanObserved=ROOT.TGraph2D()
limitsScanObservedOneSigmaDown=ROOT.TGraph2D()
limitsScanObservedOneSigmaUp=ROOT.TGraph2D()
crossSectionScanObserved = ROOT.TGraph2D()
limitsScanObserved.SetName("limitsScanObserved")
limitsScanObservedOneSigmaDown.SetName("limitsScanObservedOneSigmaDown")
limitsScanObservedOneSigmaUp.SetName("limitsScanObservedOneSigmaUp")
crossSectionScanObserved.SetName("crossSectionScanObserved")

signalStrengthScan = ROOT.TGraph2D()
signalStrengthScan.SetName("signalStrengthScan")

signalInjection_bestFitSignalStrengthScan = ROOT.TGraph2D()
signalInjection_bestFitSignalStrengthScan.SetName("signalInjection_bestFitSignalStrengthScan")

METCorrelationStudy_limitsRatioScan = ROOT.TGraph2D()
METCorrelationStudy_limitsRatioScan.SetName("METCorrelationStudy_limitsRatioScan")

# rateParamNames = []
# rateParamBestFitScans = {}
# abbreviated_selectionNames = {
#     "signal": "s",
#     "signal_loose": "l",
#     "control": "c"
# }
# for selection in selectionsToUse:
#     rateParamBestFitScans[selection] = {}
#     for rateParamType in ["const", "slope"]:
#         rateParamBestFitScans[selection][rateParamType] = {}
#         for nJetsBin in range(4, 7):
#             rateParamNames.append("{t}_{s}_{n}Jets".format(t=rateParamType, s=abbreviated_selectionNames[selection], n=nJetsBin))
#             rateParamBestFitScans[selection][rateParamType][nJetsBin] = ROOT.TGraph2D()
#             rateParamBestFitScans[selection][rateParamType][nJetsBin].SetName("rateParamBestFitScan_{r}_{t}_{n}Jets".format(r=selection, t=rateParamType, n=nJetsBin))

expectedCrossSectionLimits = []
observedCrossSectionLimits = []

maxValue_crossSectionScanExpected = -1
minValue_crossSectionScanExpected = -1
maxValue_crossSectionScanObserved = -1
minValue_crossSectionScanObserved = -1

def passesSanityCheck(observedUpperLimits, expectedUpperLimit):
    for observedUpperLimit in observedUpperLimits:
        ratio = observedUpperLimit/expectedUpperLimit
        if ((ratio > inputArguments.maxAllowedRatio) or (ratio < (1.0/inputArguments.maxAllowedRatio))):
            return False
    return True

unavailableBins = []
anomalousBinWarnings = []
minEventProgenitorMassBin = -1
maxEventProgenitorMassBin = -1
minNeutralinoMassBin = -1
maxNeutralinoMassBin = -1
for indexPair in templateReader.nextValidBin():
    eventProgenitorMassBin = indexPair[0]
    eventProgenitorMass = (templateReader.eventProgenitorMasses)[eventProgenitorMassBin]
    neutralinoMassBin = indexPair[1]
    neutralinoMass = (templateReader.neutralinoMasses)[neutralinoMassBin]
    if (neutralinoMass < inputArguments.minNeutralinoMass): continue
    if (eventProgenitorMass < minEventProgenitorMass): continue
    if ((eventProgenitorMass - neutralinoMass) < inputArguments.minMassDifference): continue

    crossSection = crossSectionsDictionary[int(0.5+eventProgenitorMass)]
    print("Analyzing bin at (eventProgenitorMassBin, neutralinoMassBin) = ({gMB}, {nMB}) ==> (eventProgenitorMass, neutralinoMass) = ({gM}, {nM})".format(gMB=eventProgenitorMassBin, gM=eventProgenitorMass, nMB=neutralinoMassBin, nM=neutralinoMass))

    # Check if signal contamination is in control at this mass point, and if so, fill "cleaned" signal contamination histograms
    inputHistogramsList = []
    for selection in selectionsToUse:
        for nJetsBin in range(2, 7):
            STRegionsToFetch = None
            if (nJetsBin == 3):
                STRegionsToFetch = []
            else:
                STRegionsToFetch = [1]
            for STRegionIndex in STRegionsToFetch:
                inputHistogramsList.append(signalContaminationHistograms_input[selection][nJetsBin][STRegionIndex])
    max_signal_contamination = get_max_signal_contamination(eventProgenitorMass, neutralinoMass, inputHistogramsList)
    if (max_signal_contamination > SIGNAL_CONTAMINATION_THRESHOLD):
        print("In relevant bins, max potential signal contamination = {s} is above threshold. Not using this bin for inference.".format(s=max_signal_contamination))
        continue
    print("Max potential signal contamination = {s} is below threshold.".format(s=max_signal_contamination))
    for selection in selectionsToUse:
        for nJetsBin in range(2, 7):
            STRegionsToFetch = None
            if (nJetsBin == 2): STRegionsToFetch = range(1, 8)
            elif (nJetsBin == 3): STRegionsToFetch = []
            else: STRegionsToFetch = [1]
            for STRegionIndex in STRegionsToFetch:
                (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).SetBinContent((signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).FindFixBin(eventProgenitorMass, neutralinoMass),
                                                                                                          (signalContaminationHistograms_input[selection][nJetsBin][STRegionIndex]).GetBinContent((signalContaminationHistograms_input[selection][nJetsBin][STRegionIndex]).FindFixBin(eventProgenitorMass, neutralinoMass)))

    # Get signal contamination monitored quantities
    signal_contamination_monitored_quantities = get_signal_contamination_monitored_quantities(eventProgenitorMassBin, neutralinoMassBin)
    for label in signalContaminationMonitoredQuantityLabels:
        for selection in selectionsToUse:
            for nJetsBin in range(4, 7):
                for STRegionIndex in range(2, 8):
                    if ((signal_contamination_monitored_quantities[label][selection][nJetsBin][STRegionIndex]) > 0.0):
                        (signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).SetBinContent((signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).FindFixBin(eventProgenitorMass, neutralinoMass), signal_contamination_monitored_quantities[label][selection][nJetsBin][STRegionIndex])

    if ((minEventProgenitorMassBin == -1) or (eventProgenitorMassBin < minEventProgenitorMassBin)): minEventProgenitorMassBin = eventProgenitorMassBin
    if ((maxEventProgenitorMassBin == -1) or (eventProgenitorMassBin > maxEventProgenitorMassBin)): maxEventProgenitorMassBin = eventProgenitorMassBin
    if ((minNeutralinoMassBin == -1) or (neutralinoMassBin < minNeutralinoMassBin)): minNeutralinoMassBin = neutralinoMassBin
    if ((maxNeutralinoMassBin == -1) or (neutralinoMassBin > maxNeutralinoMassBin)): maxNeutralinoMassBin = neutralinoMassBin

    # nominal
    expectedUpperLimit, expectedUpperLimitOneSigmaDown, expectedUpperLimitOneSigmaUp, observedUpperLimit = (-1, -1, -1, -1)
    try:
        expectedUpperLimit, expectedUpperLimitOneSigmaDown, expectedUpperLimitOneSigmaUp, observedUpperLimit = commonFunctions.get_expected_and_observed_limits_from_combine_output(combineOutputFilePath="{cRD}/higgsCombine_{cOP}_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin))
    except ValueError:
        unavailableBins.append((eventProgenitorMassBin, neutralinoMassBin))
        continue

    # cross section down
    observedUpperLimitOneSigmaDown = -1
    try:
        observedUpperLimitOneSigmaDown = commonFunctions.get_observed_limit_from_combine_output(combineOutputFilePath="{cRD}/higgsCombine_{cOP}_crossSectionsDown_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin))
    except ValueError:
        unavailableBins.append((eventProgenitorMassBin, neutralinoMassBin))
        continue

    # cross section up
    observedUpperLimitOneSigmaUp = -1
    try:
        observedUpperLimitOneSigmaUp = commonFunctions.get_observed_limit_from_combine_output(combineOutputFilePath="{cRD}/higgsCombine_{cOP}_crossSectionsUp_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin))
    except ValueError:
        unavailableBins.append((eventProgenitorMassBin, neutralinoMassBin))
        continue

    print("Limits: Observed: ({lobsdown}, {lobs}, {lobsup}); Expected: ({lexpdown}, {lexp}, {lexpup})".format(lobsdown=observedUpperLimitOneSigmaDown, lobs=observedUpperLimit, lobsup=observedUpperLimitOneSigmaUp, lexpdown=expectedUpperLimitOneSigmaDown, lexp=expectedUpperLimit, lexpup=expectedUpperLimitOneSigmaUp))
    if (inputArguments.plotObserved and not(passesSanityCheck(observedUpperLimits=[observedUpperLimit, observedUpperLimitOneSigmaUp, observedUpperLimitOneSigmaDown], expectedUpperLimit=expectedUpperLimit))):
        anomalousBinWarnings.append("WARNING: observed limits deviate too much from expected limits at eventProgenitorMass = {gM}, neutralinoMass={nM}".format(gM=eventProgenitorMass, nM=neutralinoMass))
    limitsScanExpected.SetPoint(limitsScanExpected.GetN(), eventProgenitorMass, neutralinoMass, expectedUpperLimit)
    limitsScanExpectedOneSigmaDown.SetPoint(limitsScanExpectedOneSigmaDown.GetN(), eventProgenitorMass, neutralinoMass, expectedUpperLimitOneSigmaDown)
    limitsScanExpectedOneSigmaUp.SetPoint(limitsScanExpectedOneSigmaUp.GetN(), eventProgenitorMass, neutralinoMass, expectedUpperLimitOneSigmaUp)
    crossSectionScanExpected.SetPoint(crossSectionScanExpected.GetN(), eventProgenitorMass, neutralinoMass, expectedUpperLimit*crossSection)
    expectedCrossSectionLimits.append(((eventProgenitorMass, neutralinoMass), expectedUpperLimit*crossSection))
    if ((minValue_crossSectionScanExpected == -1) or (expectedUpperLimit*crossSection < minValue_crossSectionScanExpected)): minValue_crossSectionScanExpected = expectedUpperLimit*crossSection
    if ((maxValue_crossSectionScanExpected == -1) or (expectedUpperLimit*crossSection > maxValue_crossSectionScanExpected)): maxValue_crossSectionScanExpected = expectedUpperLimit*crossSection

    limitsScanObserved.SetPoint(limitsScanObserved.GetN(), eventProgenitorMass, neutralinoMass, observedUpperLimit)
    limitsScanObservedOneSigmaDown.SetPoint(limitsScanObservedOneSigmaDown.GetN(), eventProgenitorMass, neutralinoMass, observedUpperLimitOneSigmaDown)
    limitsScanObservedOneSigmaUp.SetPoint(limitsScanObservedOneSigmaUp.GetN(), eventProgenitorMass, neutralinoMass, observedUpperLimitOneSigmaUp)
    crossSectionScanObserved.SetPoint(crossSectionScanObserved.GetN(), eventProgenitorMass, neutralinoMass, observedUpperLimit*crossSection)
    observedCrossSectionLimits.append(((eventProgenitorMass, neutralinoMass), observedUpperLimit*crossSection))
    if ((minValue_crossSectionScanObserved == -1) or (observedUpperLimit*crossSection < minValue_crossSectionScanObserved)): minValue_crossSectionScanObserved = observedUpperLimit*crossSection
    if ((maxValue_crossSectionScanObserved == -1) or (observedUpperLimit*crossSection > maxValue_crossSectionScanObserved)): maxValue_crossSectionScanObserved = observedUpperLimit*crossSection

    if (inputArguments.plotObserved): signalStrengthScan.SetPoint(signalStrengthScan.GetN(), eventProgenitorMass, neutralinoMass, observedUpperLimit)
    else: signalStrengthScan.SetPoint(signalStrengthScan.GetN(), eventProgenitorMass, neutralinoMass, expectedUpperLimit)

    try:
        signal_strength_best_fit = commonFunctions.get_best_fits_from_MultiDim_fitResult(multiDimFitResultFilePath="{cRD}/multidimfit_WITH_ADDED_SIGNAL_{cOP}_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin), parameter_names=["r"])
        signalInjection_bestFitSignalStrengthScan.SetPoint(signalInjection_bestFitSignalStrengthScan.GetN(), eventProgenitorMass, neutralinoMass, signal_strength_best_fit["r"])
        print("Best-fit signal strength from the multidim output for the signal-injected model: {v:.3f}".format(v=signal_strength_best_fit["r"]))
    except ValueError:
        anomalousBinWarnings.append("WARNING: best fit signal strength not available at eventProgenitorMass = {gM}, neutralinoMass={nM}".format(gM=eventProgenitorMass, nM=neutralinoMass))
    try:
        expectedUpperLimit_with_MET_uncertainties_uncorrelated = (commonFunctions.get_expected_and_observed_limits_from_combine_output(combineOutputFilePath="{cRD}/higgsCombine_{cOP}_METUncUncorrelated_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}.AsymptoticLimits.mH120.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin)))[0]
        METCorrelationStudy_limitsRatioScan.SetPoint(METCorrelationStudy_limitsRatioScan.GetN(), eventProgenitorMass, neutralinoMass, expectedUpperLimit_with_MET_uncertainties_uncorrelated/expectedUpperLimit)
        print("Ratio of limits with uncorrelated vs correlated MET uncertainties: {v:.3f}".format(v=expectedUpperLimit_with_MET_uncertainties_uncorrelated/expectedUpperLimit))
    except:
        anomalousBinWarnings.append("WARNING: MET correlation study limits not available at eventProgenitorMass = {gM}, neutralinoMass={nM}".format(gM=eventProgenitorMass, nM=neutralinoMass))

    # print("Now fetching best fit for rate params from multidim output...")
    # try:
    #     rateParam_bestFits = commonFunctions.get_best_fits_from_MultiDim_fitResult(multiDimFitResultFilePath="{cRD}/multidimfit_{cOP}_eventProgenitorMassBin{gMB}_neutralinoMassBin{nMB}.root".format(cRD=inputArguments.combineResultsDirectory, cOP=inputArguments.combineOutputPrefix, gMB=eventProgenitorMassBin, nMB=neutralinoMassBin), parameter_names=rateParamNames)
    #     print("rateParam_bestFits: {rPbF}".format(rPbF=rateParam_bestFits))

    #     for selection in selectionsToUse:
    #         for rateParamType in ["const", "slope"]:
    #             for nJetsBin in range(4, 7):
    #                 rateParamBestFitScans[selection][rateParamType][nJetsBin].SetPoint(rateParamBestFitScans[selection][rateParamType][nJetsBin].GetN(), eventProgenitorMass, neutralinoMass, rateParam_bestFits["{t}_{s}_{n}Jets".format(t=rateParamType, s=abbreviated_selectionNames[selection], n=nJetsBin)])

    # except ValueError:
    #     print("Unable to fetch rate params for this bin.")

for unavailableBin in unavailableBins:
    eventProgenitorMassBin = unavailableBin[0]
    eventProgenitorMass = (templateReader.eventProgenitorMasses)[eventProgenitorMassBin]
    neutralinoMassBin = unavailableBin[1]
    neutralinoMass = (templateReader.neutralinoMasses)[neutralinoMassBin]
    print("WARNING: Limits not available for (eventProgenitorMassBin, neutralinoMassBin) = ({ePMB}, {nMB}) ==> (eventProgenitorMass, neutralinoMass) = ({ePM}, {nM})".format(ePMB=eventProgenitorMassBin, ePM=eventProgenitorMass, nMB=neutralinoMassBin, nM=neutralinoMass))

for anomalousBinWarning in anomalousBinWarnings:
    print(anomalousBinWarning)

del templateReader

outputExpectedCrossSectionsFile=open("{oD}/expectedCrossSections_{s}.txt".format(oD=inputArguments.outputDirectory_rawOutput, s=inputArguments.outputSuffix), 'w')
outputExpectedCrossSectionsFile.write("{gMTitle:<19}{nMTitle:<19}{eXSTitle}\n".format(gMTitle="eventProgenitor mass", nMTitle="neutralino mass", eXSTitle="Expected limits on cross section (pb)"))
for expectedCrossSectionLimit in expectedCrossSectionLimits:
    outputExpectedCrossSectionsFile.write("{gM:<19.1f}{nM:<19.1f}{eXS:.3e}\n".format(gM=expectedCrossSectionLimit[0][0], nM=expectedCrossSectionLimit[0][1], eXS=expectedCrossSectionLimit[1]))
outputExpectedCrossSectionsFile.close()

if (inputArguments.plotObserved):
    outputObservedCrossSectionsFile=open("{oD}/observedCrossSections_{s}.txt".format(oD=inputArguments.outputDirectory_rawOutput, s=inputArguments.outputSuffix), 'w')
    outputObservedCrossSectionsFile.write("{gMTitle:<19}{nMTitle:<19}{eXSTitle}\n".format(gMTitle="eventProgenitor mass", nMTitle="neutralino mass", eXSTitle="Observed limits on cross section (pb)"))
    for observedCrossSectionLimit in observedCrossSectionLimits:
        outputObservedCrossSectionsFile.write("{gM:<19.1f}{nM:<19.1f}{oXS:.3e}\n".format(gM=observedCrossSectionLimit[0][0], nM=observedCrossSectionLimit[0][1], oXS=observedCrossSectionLimit[1]))
    outputObservedCrossSectionsFile.close()

listOf2DScans = [limitsScanExpected, limitsScanExpectedOneSigmaDown, limitsScanExpectedOneSigmaUp, crossSectionScanExpected, limitsScanObserved, limitsScanObservedOneSigmaDown, limitsScanObservedOneSigmaUp, crossSectionScanObserved, signalStrengthScan, signalInjection_bestFitSignalStrengthScan, METCorrelationStudy_limitsRatioScan]
# for selection in selectionsToUse:
#     for rateParamType in ["const", "slope"]:
#         for nJetsBin in range(4, 7):
#             listOf2DScans.append(rateParamBestFitScans[selection][rateParamType][nJetsBin])
for scan2D in listOf2DScans:
    scan2D.SetNpx(8*(1 + maxEventProgenitorMassBin - minEventProgenitorMassBin))
    scan2D.SetNpy(2*(1 + maxNeutralinoMassBin - minNeutralinoMassBin))

histogramExpectedLimits = limitsScanExpected.GetHistogram()
histogramExpectedLimits.SetName("histogramExpectedLimits")
histogramExpectedLimitsOneSigmaDown = limitsScanExpectedOneSigmaDown.GetHistogram()
histogramExpectedLimitsOneSigmaDown.SetName("histogramExpectedLimitsOneSigmaDown")
histogramExpectedLimitsOneSigmaUp = limitsScanExpectedOneSigmaUp.GetHistogram()
histogramExpectedLimitsOneSigmaUp.SetName("histogramExpectedLimitsOneSigmaUp")
histogramCrossSectionScanExpected = crossSectionScanExpected.GetHistogram()
histogramCrossSectionScanExpected.SetName("histogramCrossSectionScanExpected")

histogramObservedLimits = limitsScanObserved.GetHistogram()
histogramObservedLimits.SetName("histogramObservedLimits")
histogramObservedLimitsOneSigmaDown = limitsScanObservedOneSigmaDown.GetHistogram()
histogramObservedLimitsOneSigmaDown.SetName("histogramObservedLimitsOneSigmaDown")
histogramObservedLimitsOneSigmaUp = limitsScanObservedOneSigmaUp.GetHistogram()
histogramObservedLimitsOneSigmaUp.SetName("histogramObservedLimitsOneSigmaUp")
histogramCrossSectionScanObserved = crossSectionScanObserved.GetHistogram()
histogramCrossSectionScanObserved.SetName("histogramCrossSectionScanObserved")

expectedLimitContours = limitsScanExpected.GetContourList(inputArguments.contour_signalStrength)
expectedLimitContours.SetName("expectedLimitContours")
expectedLimitContoursOneSigmaDown = limitsScanExpectedOneSigmaDown.GetContourList(inputArguments.contour_signalStrength)
expectedLimitContoursOneSigmaDown.SetName("expectedLimitContoursOneSigmaDown")
expectedLimitContoursOneSigmaUp = limitsScanExpectedOneSigmaUp.GetContourList(inputArguments.contour_signalStrength)
expectedLimitContoursOneSigmaUp.SetName("expectedLimitContoursOneSigmaUp")

histogramSignalStrengthScan = signalStrengthScan.GetHistogram()
histogramSignalStrengthScan.SetName("histogramSignalStrengthScan")

histogram_signalInjection_bestFitSignalStrengthScan = signalInjection_bestFitSignalStrengthScan.GetHistogram()
histogram_signalInjection_bestFitSignalStrengthScan.SetName("histogram_signalInjection_bestFitSignalStrengthScan")

histogram_METCorrelationStudy_limitsRatioScan = METCorrelationStudy_limitsRatioScan.GetHistogram()
histogram_METCorrelationStudy_limitsRatioScan.SetName("histogram_METCorrelationStudy_limitsRatioScan")

# histogram_rateParamBestFitScans = {}
# for selection in selectionsToUse:
#     histogram_rateParamBestFitScans[selection] = {}
#     for rateParamType in ["const", "slope"]:
#         histogram_rateParamBestFitScans[selection][rateParamType] = {}
#         for nJetsBin in range(4, 7):
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin] = rateParamBestFitScans[selection][rateParamType][nJetsBin].GetHistogram()
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].SetName("best-fit rate parameter: {s} selection, type {t}, {n} Jets".format(s=selection, t=rateParamType, n=nJetsBin))

if (inputArguments.plotObserved):
    observedLimitContours = limitsScanObserved.GetContourList(inputArguments.contour_signalStrength)
    observedLimitContours.SetName("observedLimitContours")
    observedLimitContoursOneSigmaDown = limitsScanObservedOneSigmaDown.GetContourList(inputArguments.contour_signalStrength)
    observedLimitContoursOneSigmaDown.SetName("observedLimitContoursOneSigmaDown")
    observedLimitContoursOneSigmaUp = limitsScanObservedOneSigmaUp.GetContourList(inputArguments.contour_signalStrength)
    observedLimitContoursOneSigmaUp.SetName("observedLimitContoursOneSigmaUp")

color_expectedContours = ROOT.kRed
width_expectedContours_middle = 5
width_expectedContours_topBottom = 2
style_expectedContours_middle = 2
style_expectedContours_topBottom = 2
color_observedContours = ROOT.kBlack
width_observedContours_middle = 5
width_observedContours_topBottom = 2
style_observedContours_middle = 1
style_observedContours_topBottom = 1

formatContours(expectedLimitContours, style_expectedContours_middle, width_expectedContours_middle, color_expectedContours)
formatContours(expectedLimitContoursOneSigmaDown, style_expectedContours_topBottom, width_expectedContours_topBottom, color_expectedContours)
formatContours(expectedLimitContoursOneSigmaUp, style_expectedContours_topBottom, width_expectedContours_topBottom, color_expectedContours)
if (inputArguments.plotObserved):
    formatContours(observedLimitContours, style_observedContours_middle, width_observedContours_middle, color_observedContours)
    formatContours(observedLimitContoursOneSigmaDown, style_observedContours_topBottom, width_observedContours_topBottom, color_observedContours)
    formatContours(observedLimitContoursOneSigmaUp, style_observedContours_topBottom, width_observedContours_topBottom, color_observedContours)

# CMS_lumi.writeExtraText = False
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.lumi_13TeV = "138 fb^{-1}"
CMS_lumi.relPosX    = 0.15

H_ref = 600
W_ref = 800
W = W_ref
H  = H_ref
T = 0.08*H_ref
B = 0.12*H_ref
L = 0.12*W_ref
R = 0.04*W_ref

tdrstyle.setTDRStyle()
canvas = ROOT.TCanvas("c_{s}_observedLimits".format(s=inputArguments.outputSuffix), "c_{s}_observedLimits".format(s=inputArguments.outputSuffix), 50, 50, W, H)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin( L/W )
canvas.SetRightMargin( R/W )
canvas.SetTopMargin( T/H )
canvas.SetBottomMargin( B/H )
canvas.SetTickx(0)
canvas.SetTicky(0)
canvas.Draw()

ROOT.gPad.SetLogz()

# ROOT.gStyle.SetPalette(ROOT.kBird)
colorMin_RGB = getRGB(ROOT.kBlue+1)
colorMid_RGB = getRGB(ROOT.kYellow)
colorMax_RGB = getRGB(ROOT.kRed)
paletteRed = array.array('d', [colorMin_RGB["red"], colorMid_RGB["red"], colorMax_RGB["red"]])
paletteGreen = array.array('d', [colorMin_RGB["green"], colorMid_RGB["green"], colorMax_RGB["green"]])
paletteBlue = array.array('d', [colorMin_RGB["blue"], colorMid_RGB["blue"], colorMax_RGB["blue"]])
paletteStops = None
if (inputArguments.plotObserved):
    paletteStops = array.array('d', [0., math.log10(3.)/math.log10(maxValue_crossSectionScanObserved/minValue_crossSectionScanObserved), 1.]) # "mid" color is set to 5 times the min cross-section -- note that the axis is log-scaled
else:
    paletteStops = array.array('d', [0., math.log10(3.)/math.log10(maxValue_crossSectionScanExpected/minValue_crossSectionScanExpected), 1.]) # "mid" color is set to 5 times the min cross-section -- note that the axis is log-scaled

ROOT.TColor.CreateGradientColorTable(len(paletteStops), paletteStops, paletteRed, paletteGreen, paletteBlue, 999)
ROOT.gStyle.SetNumberContours(999)

ROOT.gPad.SetRightMargin(0.2)
ROOT.gPad.SetLeftMargin(0.15)
commonTitleSize = 0.046
if (inputArguments.plotObserved):
    histogramCrossSectionScanObserved.GetXaxis().SetTitle(string_mass_eventProgenitor + "(GeV)")
    histogramCrossSectionScanObserved.GetXaxis().SetTitleSize(commonTitleSize)
    histogramCrossSectionScanObserved.GetYaxis().SetTitle(string_mass_neutralino + "(GeV)")
    histogramCrossSectionScanObserved.GetYaxis().SetTitleOffset(1.)
    histogramCrossSectionScanObserved.GetYaxis().SetTitleSize(commonTitleSize)
    histogramCrossSectionScanObserved.GetZaxis().SetTitle("95% CL upper limit on cross-section (pb)")
    histogramCrossSectionScanObserved.GetZaxis().SetTitleOffset(1.)
    histogramCrossSectionScanObserved.GetZaxis().SetTitleSize(0.046)
    histogramCrossSectionScanObserved.Draw("colz")
    ROOT.gPad.Update()
    histogramCrossSectionScanObserved.GetXaxis().SetRangeUser(minEventProgenitorMass, maxEventProgenitorMass)
    histogramCrossSectionScanObserved.GetZaxis().SetRangeUser(minValue_crossSectionScanObserved, maxValue_crossSectionScanObserved)
    ROOT.gPad.Update()
else:
    histogramCrossSectionScanExpected.GetXaxis().SetTitle(string_mass_eventProgenitor + "(GeV)")
    histogramCrossSectionScanExpected.GetXaxis().SetTitleSize(commonTitleSize)
    histogramCrossSectionScanExpected.GetYaxis().SetTitle(string_mass_neutralino + "(GeV)")
    histogramCrossSectionScanExpected.GetYaxis().SetTitleOffset(1.)
    histogramCrossSectionScanExpected.GetYaxis().SetTitleSize(commonTitleSize)
    histogramCrossSectionScanExpected.GetZaxis().SetTitle("95% CL upper limit on cross-section (pb)")
    histogramCrossSectionScanExpected.GetZaxis().SetTitleOffset(1.)
    histogramCrossSectionScanExpected.GetZaxis().SetTitleSize(0.046)
    histogramCrossSectionScanExpected.Draw("colz")
    ROOT.gPad.Update()
    histogramCrossSectionScanExpected.GetXaxis().SetRangeUser(minEventProgenitorMass, maxEventProgenitorMass)
    histogramCrossSectionScanExpected.GetZaxis().SetRangeUser(minValue_crossSectionScanExpected, maxValue_crossSectionScanExpected)
    ROOT.gPad.Update()

contoursToDraw = [expectedLimitContours, expectedLimitContoursOneSigmaDown, expectedLimitContoursOneSigmaUp]
if (inputArguments.plotObserved):
    contoursToDraw.extend([observedLimitContours, observedLimitContoursOneSigmaDown, observedLimitContoursOneSigmaUp])

for contoursList in contoursToDraw:
    contoursList.Draw("SAME")

# line_eventProgenitorEqualsNeutralinoMass = ROOT.TLine(minEventProgenitorMass, minEventProgenitorMass, maxEventProgenitorMass, maxEventProgenitorMass)
# line_eventProgenitorEqualsNeutralinoMass.SetLineStyle(7)
# line_eventProgenitorEqualsNeutralinoMass.SetLineColor(ROOT.kBlack)
# line_eventProgenitorEqualsNeutralinoMass.SetLineWidth(3)
# line_eventProgenitorEqualsNeutralinoMass.Draw()
line_eventProgenitorEqualsNeutralinoMassShiftedDown = ROOT.TLine(minEventProgenitorMass, minEventProgenitorMass-diagonal_down_shift, maxEventProgenitorMass, maxEventProgenitorMass-diagonal_down_shift)
line_eventProgenitorEqualsNeutralinoMassShiftedDown.SetLineStyle(7)
line_eventProgenitorEqualsNeutralinoMassShiftedDown.SetLineColor(ROOT.kBlack)
line_eventProgenitorEqualsNeutralinoMassShiftedDown.SetLineWidth(3)
line_eventProgenitorEqualsNeutralinoMassShiftedDown.Draw()
ROOT.gPad.Update()

commonOffset = 0.178

latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(12)
latex.SetTextColor(ROOT.kBlack)

latex.SetTextSize(0.033)
latex.DrawLatexNDC(commonOffset, 0.89, decayChain)

latex.SetTextSize(0.025)
latex.DrawLatexNDC(commonOffset, 0.845, decayChain_supplementaryInfo1)

latex.SetTextSize(0.033)
latex.DrawLatexNDC(commonOffset, 0.815, decayChain_supplementaryInfo2)

drawContoursForLegend(commonOffset, 0.03, 0.765, 0.01, width_expectedContours_middle, style_expectedContours_middle, color_expectedContours, width_expectedContours_topBottom, style_expectedContours_topBottom, color_expectedContours, 0.59)
latex.SetTextSize(0.03)
latex.SetTextColor(color_expectedContours)
latex.DrawLatexNDC(commonOffset+0.04, 0.765, "Expected #pm 1#sigma_{experiment}")

if inputArguments.plotObserved:
    drawContoursForLegend(commonOffset, 0.03, 0.722, 0.01, width_observedContours_middle, style_observedContours_middle, color_observedContours, width_observedContours_topBottom, style_observedContours_topBottom, color_observedContours, 0.5775)
    latex.SetTextSize(0.03)
    latex.SetTextColor(color_observedContours)
    latex.DrawLatexNDC(commonOffset+0.04, 0.722, "Observed #pm 1#sigma_{theory}")

# latex.SetTextAlign(22)
# latex.SetTextColor(ROOT.kBlack)
# latex.SetTextSize(0.04)
# latex.SetTextAngle(tmROOTUtils.getTLineAngleInDegrees(ROOT.gPad, line_eventProgenitorEqualsNeutralinoMass))
# latex.DrawLatex(minEventProgenitorMass + 185.0, minEventProgenitorMass + 265.0, string_mass_neutralino + " = " + string_mass_eventProgenitor)

latex.SetTextAlign(22)
latex.SetTextColor(ROOT.kBlack)
latex.SetTextSize(0.032)
latex.SetTextAngle(tmROOTUtils.getTLineAngleInDegrees(ROOT.gPad, line_eventProgenitorEqualsNeutralinoMassShiftedDown))
latex.DrawLatex(minEventProgenitorMass + 200.0, minEventProgenitorMass + 175.0, string_mass_neutralino + " = " + string_mass_eventProgenitor + " - {s} GeV".format(s=diagonal_down_shift))

CMS_lumi.CMS_lumi(canvas, 4, 0)

ROOT.gPad.Update()
ROOT.gPad.RedrawAxis()
frame = ROOT.gPad.GetFrame()
frame.Draw()
canvas.Update()
if (inputArguments.plotObserved):
    canvas.SaveAs("{oD}/{s}_observedLimits.pdf".format(oD=inputArguments.outputDirectory_plots, s=inputArguments.outputSuffix))
    sys.exit(0) # if "plotObserved, none of the remaining plots are needed."
else:
    canvas.SaveAs("{oD}/{s}_expectedLimits.pdf".format(oD=inputArguments.outputDirectory_plots, s=inputArguments.outputSuffix))

# Save cleaned signal contamination plots
for selection in selectionsToUse:
    for nJetsBin in range(2, 7):
        STRegionsToFetch = None
        if (nJetsBin == 2): STRegionsToFetch = range(1, 8)
        elif (nJetsBin == 3): STRegionsToFetch = []
        else: STRegionsToFetch = [1]
        for STRegionIndex in STRegionsToFetch:
            outputCanvas = ROOT.TCanvas("signalContamination_STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin), "signalContamination_STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin), 1024, 1024)
            outputCanvas.SetFillColor(0)
            outputCanvas.SetBorderMode(0)
            outputCanvas.SetFrameFillStyle(0)
            outputCanvas.SetFrameBorderMode(0)
            # outputCanvas.SetLeftMargin( L/W )
            # outputCanvas.SetRightMargin( R/W )
            # outputCanvas.SetTopMargin( T/H )
            # outputCanvas.SetBottomMargin( B/H )
            outputCanvas.SetLeftMargin(0.12)
            outputCanvas.SetRightMargin(0.15)
            outputCanvas.SetTopMargin(0.1)
            outputCanvas.SetBottomMargin(0.1)
            outputCanvas.SetTickx(0)
            outputCanvas.SetTicky(0)
            outputCanvas.Draw()
            ROOT.gPad.SetLogz()
            ROOT.gStyle.SetPaintTextFormat(".1e")
            ROOT.gStyle.SetPalette(ROOT.kBird)
            (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).Draw("COLZ TEXT25")
            ROOT.gPad.Update()
            (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).GetXaxis().SetTitle(string_mass_eventProgenitor + "(GeV)")
            (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).GetXaxis().SetTitleSize(commonTitleSize)
            (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).GetYaxis().SetTitle(string_mass_neutralino + "(GeV)")
            (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).GetYaxis().SetTitleOffset(1.)
            (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).GetYaxis().SetTitleSize(commonTitleSize)
            (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).GetZaxis().SetTitle("Potential signal contamination, S/B")
            # (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).GetZaxis().SetTitleOffset(0.3)
            (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).GetZaxis().SetTitleSize(commonTitleSize)
            # (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).GetZaxis().SetRangeUser(0.00005, 0.2)
            (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).GetXaxis().SetRangeUser(minEventProgenitorMass, maxEventProgenitorMass)
            ROOT.gPad.Update()
            # line_eventProgenitorEqualsNeutralinoMass_forSC = ROOT.TLine(minEventProgenitorMass, minEventProgenitorMass, maxEventProgenitorMass, maxEventProgenitorMass)
            # line_eventProgenitorEqualsNeutralinoMass_forSC.SetLineStyle(7)
            # line_eventProgenitorEqualsNeutralinoMass_forSC.SetLineColor(ROOT.kBlack)
            # line_eventProgenitorEqualsNeutralinoMass_forSC.SetLineWidth(3)
            # line_eventProgenitorEqualsNeutralinoMass_forSC.Draw()
            line_eventProgenitorEqualsNeutralinoMassShiftedDown_forSC = ROOT.TLine(minEventProgenitorMass, minEventProgenitorMass-diagonal_down_shift, maxEventProgenitorMass, maxEventProgenitorMass-diagonal_down_shift)
            line_eventProgenitorEqualsNeutralinoMassShiftedDown_forSC.SetLineStyle(7)
            line_eventProgenitorEqualsNeutralinoMassShiftedDown_forSC.SetLineColor(ROOT.kBlack)
            line_eventProgenitorEqualsNeutralinoMassShiftedDown_forSC.SetLineWidth(3)
            line_eventProgenitorEqualsNeutralinoMassShiftedDown_forSC.Draw()
            ROOT.gPad.Update()
            # latex.SetTextAlign(22)
            # latex.SetTextColor(ROOT.kBlack)
            # latex.SetTextSize(0.04)
            # latex.SetTextAngle(tmROOTUtils.getTLineAngleInDegrees(ROOT.gPad, line_eventProgenitorEqualsNeutralinoMass_forSC))
            # latex.DrawLatex(minEventProgenitorMass + 185.0, minEventProgenitorMass + 265.0, string_mass_neutralino + " = " + string_mass_eventProgenitor)
            latex.SetTextAlign(22)
            latex.SetTextColor(ROOT.kBlack)
            latex.SetTextSize(0.032)
            latex.SetTextAngle(tmROOTUtils.getTLineAngleInDegrees(ROOT.gPad, line_eventProgenitorEqualsNeutralinoMassShiftedDown_forSC))
            latex.DrawLatex(minEventProgenitorMass + 200.0, minEventProgenitorMass + 175.0, string_mass_neutralino + " = " + string_mass_eventProgenitor + " - {s} GeV".format(s=diagonal_down_shift))
            latex.SetTextSize(0.04)
            latex.SetTextAngle(0.)
            latex.SetTextAlign(21)
            latex.DrawLatexNDC(0.5, 0.92, (signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).GetTitle())
            ROOT.gPad.Update()
            # print("Histogram title: {t}".format(t=(signalContaminationHistograms_cleaned[selection][nJetsBin][STRegionIndex]).GetTitle()))
            outputCanvas.SaveAs("{o}/signalContamination_cleaned_{s1}_{s2}_STRegion{r}_{n}Jets.pdf".format(o=inputArguments.outputDirectory_rawOutput, s1=selection, s2=inputArguments.outputSuffix, r=STRegionIndex, n=nJetsBin))

for selection in selectionsToUse:
    signalContaminationSourceFileHandles[selection].Close()

for label in signalContaminationMonitoredQuantityLabels:
    for selection in selectionsToUse:
        for nJetsBin in range(4, 7):
            for STRegionIndex in range(2, 8):
                outputCanvas = ROOT.TCanvas("signalContamination_STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin), "signalContamination_STRegion{r}_{n}Jets".format(r=STRegionIndex, n=nJetsBin), 1024, 1024)
                outputCanvas.SetFillColor(0)
                outputCanvas.SetBorderMode(0)
                outputCanvas.SetFrameFillStyle(0)
                outputCanvas.SetFrameBorderMode(0)
                # outputCanvas.SetLeftMargin( L/W )
                # outputCanvas.SetRightMargin( R/W )
                # outputCanvas.SetTopMargin( T/H )
                # outputCanvas.SetBottomMargin( B/H )
                outputCanvas.SetLeftMargin(0.12)
                outputCanvas.SetRightMargin(0.15)
                outputCanvas.SetTopMargin(0.1)
                outputCanvas.SetBottomMargin(0.1)
                outputCanvas.SetTickx(0)
                outputCanvas.SetTicky(0)
                outputCanvas.Draw()
                ROOT.gPad.SetLogz()
                ROOT.gStyle.SetPaintTextFormat(".1e")
                ROOT.gStyle.SetPalette(ROOT.kBird)
                (signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).Draw("COLZ TEXT25")
                ROOT.gPad.Update()
                (signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).GetXaxis().SetTitle(string_mass_eventProgenitor + "(GeV)")
                (signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).GetXaxis().SetTitleSize(commonTitleSize)
                (signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).GetYaxis().SetTitle(string_mass_neutralino + "(GeV)")
                (signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).GetYaxis().SetTitleOffset(1.)
                (signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).GetYaxis().SetTitleSize(commonTitleSize)
                # (signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).GetZaxis().SetTitleOffset(0.3)
                (signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).GetZaxis().SetTitleSize(commonTitleSize)
                # (signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).GetZaxis().SetRangeUser(0.00005, 0.2)
                (signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).GetXaxis().SetRangeUser(minEventProgenitorMass, maxEventProgenitorMass)
                (signal_contamination_monitor_histograms[label][selection][nJetsBin][STRegionIndex]).GetYaxis().SetRangeUser(minNeutralinoMass, maxNeutralinoMass)
                latex.SetTextSize(0.04)
                latex.SetTextAngle(0.)
                latex.SetTextAlign(21)
                nJetsLabel = "{n} Jets".format(n=nJetsBin)
                if (nJetsBin == 6): nJetsLabel = "#geq 6 Jets"
                latex.DrawLatexNDC(0.5, 0.92, "{ST}, {j}".format(ST=STRegionTitles[STRegionIndex], j=nJetsLabel))
                ROOT.gPad.Update()
                outputCanvas.SaveAs("{o}/{l}_{s1}_{s2}_STRegion{r}_{n}Jets.pdf".format(o=inputArguments.outputDirectory_rawOutput, l=label, s1=selection, s2=inputArguments.outputSuffix, r=STRegionIndex, n=nJetsBin))

signalStrengthCanvas = ROOT.TCanvas("c_{s}_signalStrengthScan".format(s=inputArguments.outputSuffix), "c_{s}_signalStrengthScan".format(s=inputArguments.outputSuffix), 50, 50, W, H)
signalStrengthCanvas.SetFillColor(0)
signalStrengthCanvas.SetBorderMode(0)
signalStrengthCanvas.SetFrameFillStyle(0)
signalStrengthCanvas.SetFrameBorderMode(0)
signalStrengthCanvas.SetLeftMargin( L/W )
signalStrengthCanvas.SetRightMargin( R/W )
signalStrengthCanvas.SetTopMargin( T/H )
signalStrengthCanvas.SetBottomMargin( B/H )
signalStrengthCanvas.SetTickx(0)
signalStrengthCanvas.SetTicky(0)
ROOT.gPad.SetRightMargin(0.2)
ROOT.gPad.SetLeftMargin(0.15)
signalStrengthCanvas.Draw()
histogramSignalStrengthScan.GetXaxis().SetTitle(string_mass_eventProgenitor + "(GeV)")
histogramSignalStrengthScan.GetXaxis().SetTitleSize(commonTitleSize)
histogramSignalStrengthScan.GetXaxis().SetRangeUser(minEventProgenitorMass, maxEventProgenitorMass)
histogramSignalStrengthScan.GetYaxis().SetTitle(string_mass_neutralino + "(GeV)")
histogramSignalStrengthScan.GetYaxis().SetTitleOffset(1.)
histogramSignalStrengthScan.GetYaxis().SetTitleSize(commonTitleSize)
if (inputArguments.plotObserved): histogramSignalStrengthScan.GetZaxis().SetTitle("Upper limit on observed signal strength.")
else: histogramSignalStrengthScan.GetZaxis().SetTitle("Upper limit on expected signal strength.")
histogramSignalStrengthScan.GetZaxis().SetTitleOffset(1.)
histogramSignalStrengthScan.GetZaxis().SetTitleSize(0.046)
ROOT.gPad.SetLogz()
histogramSignalStrengthScan.Draw("colz")
signalStrengthCanvas.SaveAs("{oD}/{s}_signalStrength.pdf".format(oD=inputArguments.outputDirectory_plots, s=inputArguments.outputSuffix))

specialPlotNames = ["signalInjection", "METCorr"]
specialPlots_histograms = {
    "signalInjection": histogram_signalInjection_bestFitSignalStrengthScan,
    "METCorr": histogram_METCorrelationStudy_limitsRatioScan
}
specialPlots_histograms_titles = {
    "signalInjection": "Best-fit signal strength for (signal + background) model.",
    "METCorr": "expected limits (uncorrelated)/expected limits (correlated)."
}
specialPlots_histograms_zlimits = {
    "signalInjection": tuple([0.995, 1.005]),
    "METCorr": tuple([0.975, 1.025])
}
specialPlots_histograms_targetFiles = {
    "signalInjection": "{s}_injectedSignalModel_bestFitSignalStrength".format(s=inputArguments.outputSuffix),
    "METCorr": "{s}_METUncCorrelationStudy".format(s=inputArguments.outputSuffix)
}

for specialPlotName in specialPlotNames:
    histogram_signalInjection_bestFitSignalStrengthCanvas = ROOT.TCanvas("c_{tF}".format(tF=specialPlots_histograms_targetFiles[specialPlotName]), "c_{tF}".format(tF=specialPlots_histograms_targetFiles[specialPlotName]), 50, 50, W, H)
    histogram_signalInjection_bestFitSignalStrengthCanvas.SetFillColor(0)
    histogram_signalInjection_bestFitSignalStrengthCanvas.SetBorderMode(0)
    histogram_signalInjection_bestFitSignalStrengthCanvas.SetFrameFillStyle(0)
    histogram_signalInjection_bestFitSignalStrengthCanvas.SetFrameBorderMode(0)
    histogram_signalInjection_bestFitSignalStrengthCanvas.SetLeftMargin( L/W )
    histogram_signalInjection_bestFitSignalStrengthCanvas.SetRightMargin( R/W )
    histogram_signalInjection_bestFitSignalStrengthCanvas.SetTopMargin( T/H )
    histogram_signalInjection_bestFitSignalStrengthCanvas.SetBottomMargin( B/H )
    histogram_signalInjection_bestFitSignalStrengthCanvas.SetTickx(0)
    histogram_signalInjection_bestFitSignalStrengthCanvas.SetTicky(0)
    ROOT.gPad.SetRightMargin(0.2)
    ROOT.gPad.SetLeftMargin(0.15)
    histogram_signalInjection_bestFitSignalStrengthCanvas.Draw()
    histogram_signalInjection_bestFitSignalStrengthScan.GetXaxis().SetTitle(string_mass_eventProgenitor + "(GeV)")
    histogram_signalInjection_bestFitSignalStrengthScan.GetXaxis().SetTitleSize(commonTitleSize)
    histogram_signalInjection_bestFitSignalStrengthScan.GetXaxis().SetRangeUser(minEventProgenitorMass, maxEventProgenitorMass)
    histogram_signalInjection_bestFitSignalStrengthScan.GetYaxis().SetTitle(string_mass_neutralino + "(GeV)")
    histogram_signalInjection_bestFitSignalStrengthScan.GetYaxis().SetTitleOffset(1.)
    histogram_signalInjection_bestFitSignalStrengthScan.GetYaxis().SetTitleSize(commonTitleSize)
    histogram_signalInjection_bestFitSignalStrengthScan.GetZaxis().SetTitle(specialPlots_histograms_titles[specialPlotName])
    histogram_signalInjection_bestFitSignalStrengthScan.GetZaxis().SetTitleOffset(1.)
    histogram_signalInjection_bestFitSignalStrengthScan.GetZaxis().SetTitleSize(0.046)
    histogram_signalInjection_bestFitSignalStrengthScan.Draw("colz")
    histogram_signalInjection_bestFitSignalStrengthScan.GetZaxis().SetRangeUser((specialPlots_histograms_zlimits[specialPlotName])[0], (specialPlots_histograms_zlimits[specialPlotName])[1])
    for contoursList in contoursToDraw:
        contoursList.Draw("SAME")
    histogram_signalInjection_bestFitSignalStrengthCanvas.Update()
    histogram_signalInjection_bestFitSignalStrengthCanvas.SaveAs("{oD}/{tF}.pdf".format(oD=inputArguments.outputDirectory_plots, tF=specialPlots_histograms_targetFiles[specialPlotName]))

# paletteStops = array.array('d', [0., 1., 5.]) # New palette for rate params
# ROOT.TColor.CreateGradientColorTable(len(paletteStops), paletteStops, paletteRed, paletteGreen, paletteBlue, 999)
# for selection in selectionsToUse:
#     for rateParamType in ["const", "slope"]:
#         for nJetsBin in range(4, 7):
#             bestFitRateParamCanvas = ROOT.TCanvas("c_{s}_{sel}_bestFitRateParam_type_{t}_{n}Jets".format(s=inputArguments.outputSuffix, sel=selection, t=rateParamType, n=nJetsBin), "c_{s}_{sel}_bestFitRateParam_type_{t}_{n}Jets".format(s=inputArguments.outputSuffix, sel=selection, t=rateParamType, n=nJetsBin), 50, 50, W, H)
#             bestFitRateParamCanvas.SetFillColor(0)
#             bestFitRateParamCanvas.SetBorderMode(0)
#             bestFitRateParamCanvas.SetFrameFillStyle(0)
#             bestFitRateParamCanvas.SetFrameBorderMode(0)
#             bestFitRateParamCanvas.SetLeftMargin( L/W )
#             bestFitRateParamCanvas.SetRightMargin( R/W )
#             bestFitRateParamCanvas.SetTopMargin( T/H )
#             bestFitRateParamCanvas.SetBottomMargin( B/H )
#             bestFitRateParamCanvas.SetTickx(0)
#             bestFitRateParamCanvas.SetTicky(0)
#             ROOT.gPad.SetRightMargin(0.2)
#             ROOT.gPad.SetLeftMargin(0.15)
#             bestFitRateParamCanvas.Draw()
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].GetXaxis().SetTitle(string_mass_eventProgenitor + "(GeV)")
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].GetXaxis().SetTitleSize(commonTitleSize)
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].GetXaxis().SetRangeUser(minEventProgenitorMass, maxEventProgenitorMass)
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].GetYaxis().SetTitle(string_mass_neutralino + "(GeV)")
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].GetYaxis().SetTitleOffset(1.)
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].GetYaxis().SetTitleSize(commonTitleSize)
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].GetZaxis().SetTitle("Best-fit value.")
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].GetZaxis().SetTitleOffset(1.)
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].GetZaxis().SetTitleSize(0.046)
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].SetMaximum(5.0)
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].SetMinimum(0.0)
#             nJetsLabel = None
#             if (nJetsBin < 6):
#                 nJetsLabel = "{n} Jets".format(n=nJetsBin)
#             else:
#                 nJetsLabel = "#geq 6 Jets"
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].SetTitle("Best fit, term: {t}, {jL}, selection {s}".format(t=rateParamType, jL=nJetsLabel, s=selection))
#             # ROOT.gPad.SetLogz()
#             histogram_rateParamBestFitScans[selection][rateParamType][nJetsBin].Draw("colz0")
#             bestFitRateParamCanvas.SaveAs("{oD}/{s}_{sel}_bestFitRateParam_type_{t}_{n}Jets.pdf".format(oD=inputArguments.outputDirectory_plots, s=inputArguments.outputSuffix, sel=selection, t=rateParamType, n=nJetsBin))

outputFileName = "{oD}/limits_{suffix}.root".format(oD=inputArguments.outputDirectory_rawOutput, suffix=inputArguments.outputSuffix)
outputFile=ROOT.TFile.Open(outputFileName, "RECREATE")
tObjectsToSave = [limitsScanExpected, limitsScanExpectedOneSigmaUp, limitsScanExpectedOneSigmaDown, crossSectionScanExpected, histogramExpectedLimits, histogramExpectedLimitsOneSigmaDown, histogramExpectedLimitsOneSigmaUp, histogramCrossSectionScanExpected, expectedLimitContours, expectedLimitContoursOneSigmaDown, expectedLimitContoursOneSigmaUp, histogramSignalStrengthScan, histogram_signalInjection_bestFitSignalStrengthScan, canvas, signalStrengthCanvas, histogram_signalInjection_bestFitSignalStrengthCanvas]
if (inputArguments.plotObserved):
    tObjectsToSave.extend([limitsScanObserved, limitsScanObservedOneSigmaUp, limitsScanObservedOneSigmaDown, crossSectionScanObserved, histogramObservedLimits, histogramObservedLimitsOneSigmaDown, histogramObservedLimitsOneSigmaUp, histogramCrossSectionScanObserved, observedLimitContours, observedLimitContoursOneSigmaDown, observedLimitContoursOneSigmaUp])
for tObject in tObjectsToSave:
    outputFile.WriteTObject(tObject)

outputFile.Close()
print("All done!")
