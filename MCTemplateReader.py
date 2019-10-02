from __future__ import print_function, division

import sys

import ROOT

class MCTemplateReader:
    nEventsFractionThreshold = 0.5
    gluinoMasses = {} # map from gluino mass bin index to gluino mass
    neutralinoMasses = {} # map from neutralino mass bin index to neutralino mass

    def __init__(self, templateSourceFilePath, templateName):
        self.templateFile = ROOT.TFile.Open(templateSourceFilePath, "READ")
        if (self.templateFile is None):
            sys.exit("ERROR: templateFile is None.")
        self.h_template = ROOT.TH2F()
        (self.templateFile).GetObject(templateName, self.h_template)
        if (self.h_template is None):
            sys.exit("ERROR: h_template is None.")
        self.nGluinoMassBins = ((self.h_template).GetXaxis()).GetNbins()
        self.minGluinoMass = ((self.h_template).GetXaxis()).GetXmin()
        self.maxGluinoMass = ((self.h_template).GetXaxis()).GetXmax()
        self.nNeutralinoMassBins = ((self.h_template).GetYaxis()).GetNbins()
        self.minNeutralinoMass = ((self.h_template).GetYaxis()).GetXmin()
        self.maxNeutralinoMass = ((self.h_template).GetYaxis()).GetXmax()
        self.gluinoMasses = {}
        self.neutralinoMasses = {}

        self.maxNEvents = -1
        for gluinoBinIndex in range(1, 1+self.nGluinoMassBins):
            xBinCenter = ((self.h_template).GetXaxis()).GetBinCenter(gluinoBinIndex)
            gluinoMass = int(0.5+xBinCenter)
            (self.gluinoMasses)[gluinoBinIndex] = gluinoMass
            for neutralinoBinIndex in range(1, self.nNeutralinoMassBins):
                yBinCenter = ((self.h_template).GetYaxis()).GetBinCenter(neutralinoBinIndex)
                neutralinoMass = int(0.5+yBinCenter)
                (self.neutralinoMasses)[neutralinoBinIndex] = neutralinoMass
                binContent = (self.h_template).GetBinContent(gluinoBinIndex, neutralinoBinIndex)
                print("At (gluinoMass, neutralinoMass) = ({gM}, {nM}), templateContents: {c}".format(gM=gluinoMass, nM=neutralinoMass, c=binContent))
                if ((self.maxNEvents < 0) || (self.maxNEvents < binContent)) self.maxNEvents = binContent

    def getTotalNEvents(self, gluinoBinIndex, neutralinoBinIndex):
        return ((self.h_template).GetBinContent(gluinoBinIndex, neutralinoBinIndex))

    def isValidBin(self, gluinoBinIndex, neutralinoBinIndex):
        nEvents = self.getTotalNEvents(gluinoBinIndex, neutralinoBinIndex)
        return (nEvents > (self.nEventsFractionThreshold)*(self.maxNEvents))

    def nextValidBin(self):
        for gluinoBinIndex in range(1, 1+self.nGluinoMassBins):
            for neutralinoBinIndex in range(1, self.nNeutralinoMassBins):
                if (self.isValidBin(gluinoBinIndex, neutralinoBinIndex)): yield (gluinoBinIndex, neutralinoBinIndex)

    def test(self):
        for indexPair in self.nextValidBin():
            gluinoBinIndex = indexPair[0]
            gluinoMass = (self.gluinoMasses)[gluinoBinIndex]
            neutralinoBinIndex = indexPair[1]
            neutralinoMass = (self.neutralinoMasses)[neutralinoBinIndex]
            print("Found valid bin at (gluinoMass, neutralinoMass) = ({gM}, {nM}); templateContents: {c}".format(gM=gluinoMass, nM=neutralinoMass, c=self.getTotalNEvents(gluinoBinIndex, neutralinoBinIndex)))

    def __del__(self):
        (self.templateFile).Close()
