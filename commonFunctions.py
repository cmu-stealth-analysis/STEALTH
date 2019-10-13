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
