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
        self.nEventProgenitorMassBins = (h_template.GetXaxis()).GetNbins()
        self.minEventProgenitorMass = (h_template.GetXaxis()).GetXmin()
        self.maxEventProgenitorMass = (h_template.GetXaxis()).GetXmax()
        self.nNeutralinoMassBins = (h_template.GetYaxis()).GetNbins()
        self.minNeutralinoMass = (h_template.GetYaxis()).GetXmin()
        self.maxNeutralinoMass = (h_template.GetYaxis()).GetXmax()
        self.eventProgenitorMasses = {} # map from eventProgenitor mass bin index to eventProgenitor mass
        self.eventProgenitorMassBins = {} # map from eventProgenitor mass bin index to tuple containing lower and upper edges of eventProgenitor mass bin
        self.neutralinoMasses = {} # map from neutralino mass bin index to neutralino mass
        self.neutralinoMassBins = {} # map from neutralino mass bin index to tuple containing lower and upper edges of neutralino mass bin
        self.generated_nEvents = {} # first index: eventProgenitor mass bin, second index: neutralino mass bin; stores the number of events in a particular bin, useful for events weights
        self.maxNEvents = -1

        for eventProgenitorBinIndex in range(1, 1+self.nEventProgenitorMassBins):
            xBinCenter = (h_template.GetXaxis()).GetBinCenter(eventProgenitorBinIndex)
            xBinHi     = (h_template.GetXaxis()).GetBinUpEdge(eventProgenitorBinIndex)
            xBinLo     = (h_template.GetXaxis()).GetBinLowEdge(eventProgenitorBinIndex)
            (self.eventProgenitorMasses)[eventProgenitorBinIndex] = xBinCenter
            (self.eventProgenitorMassBins)[eventProgenitorBinIndex] = (xBinLo, xBinHi)
            (self.generated_nEvents)[eventProgenitorBinIndex] = {}
            for neutralinoBinIndex in range(1, 1+self.nNeutralinoMassBins):
                yBinCenter = (h_template.GetYaxis()).GetBinCenter(neutralinoBinIndex)
                yBinHi     = (h_template.GetYaxis()).GetBinUpEdge(neutralinoBinIndex)
                yBinLo     = (h_template.GetYaxis()).GetBinLowEdge(neutralinoBinIndex)
                (self.neutralinoMasses)[neutralinoBinIndex] = yBinCenter
                (self.neutralinoMassBins)[neutralinoBinIndex] = (yBinLo, yBinHi)
                binContent = h_template.GetBinContent(eventProgenitorBinIndex, neutralinoBinIndex)
                # print("At (eventProgenitorMass, neutralinoMass) = ({gM}, {nM}), templateContents: {c}".format(gM=eventProgenitorMass, nM=neutralinoMass, c=binContent))
                ((self.generated_nEvents)[eventProgenitorBinIndex])[neutralinoBinIndex] = binContent
                if ((self.maxNEvents < 0) or (self.maxNEvents < binContent)): self.maxNEvents = binContent
        templateFile.Close()

    def getTotalNEvents(self, eventProgenitorBinIndex, neutralinoBinIndex):
        return ((self.generated_nEvents)[eventProgenitorBinIndex][neutralinoBinIndex])

    def isValidBin(self, eventProgenitorBinIndex, neutralinoBinIndex):
        nEvents = self.getTotalNEvents(eventProgenitorBinIndex, neutralinoBinIndex)
        return (nEvents > (self.nEventsFractionThreshold)*(self.maxNEvents))

    def nextValidBin(self):
        for eventProgenitorBinIndex in range(1, 1+self.nEventProgenitorMassBins):
            for neutralinoBinIndex in range(1, 1+self.nNeutralinoMassBins):
                if (self.isValidBin(eventProgenitorBinIndex, neutralinoBinIndex)): yield (eventProgenitorBinIndex, neutralinoBinIndex)

    def test(self):
        for indexPair in self.nextValidBin():
            eventProgenitorBinIndex = indexPair[0]
            eventProgenitorMass = (self.eventProgenitorMasses)[eventProgenitorBinIndex]
            neutralinoBinIndex = indexPair[1]
            neutralinoMass = (self.neutralinoMasses)[neutralinoBinIndex]
            print("Found valid bin at (eventProgenitorMass, neutralinoMass) = ({gM}, {nM}); templateContents: {c}".format(gM=eventProgenitorMass, nM=neutralinoMass, c=self.getTotalNEvents(eventProgenitorBinIndex, neutralinoBinIndex)))
