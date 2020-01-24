from __future__ import print_function, division

import argparse, ROOT, sys

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

def get_expected_and_observed_limits_from_combine_output(combineOutputFilePath):
    combineOutputFile=ROOT.TFile.Open(combineOutputFilePath, "READ")
    if ((combineOutputFile.IsZombie() == ROOT.kTRUE) or not(combineOutputFile.IsOpen() == ROOT.kTRUE)):
        sys.exit("Error in opening file: {cOFP}".format(cOFP=combinOutputFilePath))
    limitTree = ROOT.TTree()
    combineOutputFile.GetObject("limit", limitTree)
    nEntriesFound = limitTree.GetEntries()
    if not(nEntriesFound == 6):
        raise ValueError
    limitTree.GetEntry(2)
    expectedUpperLimit = limitTree.limit
    limitTree.GetEntry(1)
    expectedUpperLimitOneSigmaDown = limitTree.limit
    limitTree.GetEntry(3)
    expectedUpperLimitOneSigmaUp = limitTree.limit
    limitTree.GetEntry(5)
    observedUpperLimit = limitTree.limit
    combineOutputFile.Close()
    return (expectedUpperLimit, expectedUpperLimitOneSigmaDown, expectedUpperLimitOneSigmaUp, observedUpperLimit)

def get_observed_limit_from_combine_output(combineOutputFilePath):
    combineOutputFile=ROOT.TFile.Open(combineOutputFilePath, "READ")
    if ((combineOutputFile.IsZombie() == ROOT.kTRUE) or not(combineOutputFile.IsOpen() == ROOT.kTRUE)):
        sys.exit("Error in opening file: {cOFP}".format(cOFP=combinOutputFilePath))
    limitTree = ROOT.TTree()
    combineOutputFile.GetObject("limit", limitTree)
    nEntriesFound = limitTree.GetEntries()
    if not(nEntriesFound == 6):
        raise ValueError
    limitTree.GetEntry(5)
    observedUpperLimit = limitTree.limit
    combineOutputFile.Close()
    return observedUpperLimit

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
