from __future__ import print_function, division

import argparse, ROOT, sys

QUANTILE_TOLERANCE=0.001

def get_nEvts_from_file(fileName, printDebug=False):
    ggIn = ROOT.TChain("ggNtuplizer/EventTree")
    ggIn.Add(fileName)
    nEntries = ggIn.GetEntries()
    if (printDebug): print("Entries in file \"{f}\": {e}".format(f=fileName, e=nEntries))
    return (nEntries)

def get_nEvts_from_fileList(inputFilesList, printDebug=False):
    listOfInputFiles = []
    inputFileNamesFileObject = open(inputFilesList, 'r')
    totalNEvents = 0
    for inputFileName in inputFileNamesFileObject:
        nEntriesInFile = get_nEvts_from_file(fileName=inputFileName.strip(), printDebug=printDebug)
        if (printDebug): print("Found {n} entries in file {f}.".format(n=nEntriesInFile, f=inputFileName))
        totalNEvents += nEntriesInFile
    inputFileNamesFileObject.close()
    return totalNEvents

def get_nEvts_from_fileList_check(inputFilesList, printDebug=False):
    listOfInputFiles = []
    inputFileNamesFileObject = open(inputFilesList, 'r')
    for inputFileName in inputFileNamesFileObject:
        listOfInputFiles.append(inputFileName.strip())
    inputFileNamesFileObject.close()

    # Load input TTrees into TChain
    ggIn = ROOT.TChain("ggNtuplizer/EventTree")

    for inputFile in listOfInputFiles:
        if (printDebug): print("Adding file {f} to chain.".format(f=inputFileName))
        ggIn.Add(inputFile)

    totalNEntries = ggIn.GetEntries()
    return totalNEntries

def get_number_of_lines_in_file(inputFilePath):
    fileObject = open(inputFilePath, 'r')
    nLines = len(fileObject.readlines())
    fileObject.close()
    return nLines

def read_limits_from_combine_output(combineOutputFilePath):
    combineOutputFile = None
    try:
        combineOutputFile = ROOT.TFile.Open(combineOutputFilePath, "READ")
        if ((combineOutputFile.IsZombie() == ROOT.kTRUE) or not(combineOutputFile.IsOpen() == ROOT.kTRUE)):
            # sys.exit("Error in opening file: {cOFP}".format(cOFP=combinOutputFilePath))
            raise ValueError
    except:
        raise ValueError
    limitTree = ROOT.TTree()
    combineOutputFile.GetObject("limit", limitTree)
    nEntriesFound = limitTree.GetEntries()
    if not(nEntriesFound == 6):
        raise ValueError
    limitsList = []
    for entryIndex in range(0, nEntriesFound):
        nBytesRead = limitTree.GetEntry(entryIndex)
        if (nBytesRead <= 0): sys.exit("ERROR in reading limit tree at index = {i}".format(i=entryIndex))
        upperLimit = limitTree.limit
        quantile = limitTree.quantileExpected
        limitsList.append((quantile, upperLimit))
    combineOutputFile.Close()
    return limitsList

def get_limit_with_quantile(limitsList, targetQuantile):
    for limitQuantileLine in limitsList:
        quantile = limitQuantileLine[0]
        if (abs((quantile/targetQuantile) - 1.0) < 0.001):
            return limitQuantileLine[1]
    sys.exit("ERROR: Unable to find limit at target quantile = {q}. limitsList: {lL}".format(q=targetQuantile, lL=limitsList))

def get_expected_and_observed_limits_from_combine_output(combineOutputFilePath):
    limitsList = read_limits_from_combine_output(combineOutputFilePath)
    expectedUpperLimit = get_limit_with_quantile(limitsList, 0.5)
    expectedUpperLimitOneSigmaDown = get_limit_with_quantile(limitsList, 0.16)
    expectedUpperLimitOneSigmaUp = get_limit_with_quantile(limitsList, 0.84)
    observedUpperLimit = get_limit_with_quantile(limitsList, -1.0)
    return (expectedUpperLimit, expectedUpperLimitOneSigmaDown, expectedUpperLimitOneSigmaUp, observedUpperLimit)

def get_observed_limit_from_combine_output(combineOutputFilePath):
    return (get_expected_and_observed_limits_from_combine_output(combineOutputFilePath)[3])

def get_best_fit_from_MultiDim_output(multiDimOutputFilePath):
    inputFile=ROOT.TFile.Open("{mDOFP}".format(mDOFP=multiDimOutputFilePath), "READ")
    if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)):
        sys.exit("Error in opening file: {mDOFP}".format(mDOFP=multiDimOutputFilePath))
    limitTree = ROOT.TTree()
    inputFile.GetObject("limit", limitTree)
    nEntriesFound = limitTree.GetEntries()
    if not(nEntriesFound == 1):
        inputFile.Close()
        raise ValueError
        # sys.exit("Error: multidim fit output not in expected format.")
    limitTree.GetEntry(0)
    bestFitValue = limitTree.r
    inputFile.Close()
    return bestFitValue

def get_best_fit_rateParams_from_MultiDim_fitResult(multiDimFitResultFilePath, paramNames):
    inputFile=ROOT.TFile.Open("{mDFRFP}".format(mDFRFP=multiDimFitResultFilePath), "READ")
    if not(inputFile): raise ValueError
    if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)):
        sys.exit("Error in opening file: {mDFRFP}".format(mDFRFP=multiDimOutputFilePath))
    outputDict = {}
    fitResult = ROOT.RooFitResult()
    inputFile.GetObject("fit_mdf", fitResult)
    for paramName in paramNames:
        try:
            outputDict[paramName] = fitResult.floatParsFinal().find(paramName).getVal()
        except:
            raise ValueError
    inputFile.Close()
    return outputDict
