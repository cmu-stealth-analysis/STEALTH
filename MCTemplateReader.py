from __future__ import print_function, division

import sys

import ROOT

class MCTemplateReader:
    nEventsFractionThreshold = 0.1

    def __init__(self, templateSourceFilePath):
        templateFile = ROOT.TFile.Open(templateSourceFilePath, "READ")
        if (templateFile is None):
            sys.exit("ERROR: templateFile is None.")
        h_template = ROOT.TH2F()
        templateFile.GetObject("h_masses", h_template)
        if (h_template is None):
            sys.exit("ERROR: h_template is None.")
        self.nGluinoMassBins = (h_template.GetXaxis()).GetNbins()
        self.minGluinoMass = (h_template.GetXaxis()).GetXmin()
        self.maxGluinoMass = (h_template.GetXaxis()).GetXmax()
        self.nNeutralinoMassBins = (h_template.GetYaxis()).GetNbins()
        self.minNeutralinoMass = (h_template.GetYaxis()).GetXmin()
        self.maxNeutralinoMass = (h_template.GetYaxis()).GetXmax()
        self.gluinoMasses = {} # map from gluino mass bin index to gluino mass 
        self.neutralinoMasses = {} # map from neutralino mass bin index to neutralino mass
        self.generated_nEvents = {} # first index: gluino mass bin, second index: neutralino mass bin; stores the number of events in a particular bin, useful for events weights
        self.maxNEvents = -1

        for gluinoBinIndex in range(1, 1+self.nGluinoMassBins):
            xBinCenter = (h_template.GetXaxis()).GetBinCenter(gluinoBinIndex)
            gluinoMass = int(0.5+xBinCenter)
            (self.gluinoMasses)[gluinoBinIndex] = gluinoMass
            (self.generated_nEvents)[gluinoBinIndex] = {}
            for neutralinoBinIndex in range(1, 1+self.nNeutralinoMassBins):
                yBinCenter = (h_template.GetYaxis()).GetBinCenter(neutralinoBinIndex)
                neutralinoMass = int(0.5+yBinCenter)
                (self.neutralinoMasses)[neutralinoBinIndex] = neutralinoMass
                binContent = h_template.GetBinContent(gluinoBinIndex, neutralinoBinIndex)
                # print("At (gluinoMass, neutralinoMass) = ({gM}, {nM}), templateContents: {c}".format(gM=gluinoMass, nM=neutralinoMass, c=binContent))
                ((self.generated_nEvents)[gluinoBinIndex])[neutralinoBinIndex] = binContent
                if ((self.maxNEvents < 0) or (self.maxNEvents < binContent)): self.maxNEvents = binContent
        templateFile.Close()

    def getTotalNEvents(self, gluinoBinIndex, neutralinoBinIndex):
        return ((self.generated_nEvents)[gluinoBinIndex][neutralinoBinIndex])

    def isValidBin(self, gluinoBinIndex, neutralinoBinIndex):
        nEvents = self.getTotalNEvents(gluinoBinIndex, neutralinoBinIndex)
        return (nEvents > (self.nEventsFractionThreshold)*(self.maxNEvents))

    def nextValidBin(self):
        for gluinoBinIndex in range(1, 1+self.nGluinoMassBins):
            for neutralinoBinIndex in range(1, 1+self.nNeutralinoMassBins):
                if (self.isValidBin(gluinoBinIndex, neutralinoBinIndex)): yield (gluinoBinIndex, neutralinoBinIndex)

    def test(self):
        for indexPair in self.nextValidBin():
            gluinoBinIndex = indexPair[0]
            gluinoMass = (self.gluinoMasses)[gluinoBinIndex]
            neutralinoBinIndex = indexPair[1]
            neutralinoMass = (self.neutralinoMasses)[neutralinoBinIndex]
            print("Found valid bin at (gluinoMass, neutralinoMass) = ({gM}, {nM}); templateContents: {c}".format(gM=gluinoMass, nM=neutralinoMass, c=self.getTotalNEvents(gluinoBinIndex, neutralinoBinIndex)))
