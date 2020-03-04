#!/usr/bin/env python

from __future__ import print_function, division

import sys
import ROOT
import stealthEnv, MCTemplateReader

tDesignationsForProgenitor = ["t5", "t6"]
templateSources = {
    "t5": "{eP}/store/group/lpcsusystealth/MCGeneratedMasses/MCGeneratedMasses/MCGeneratedMasses_stealth_t5Wg_savedObjects.root".format(eP=stealthEnv.EOSPrefix),
    "t6": "{eP}/store/group/lpcsusystealth/MCGeneratedMasses/MCGeneratedMasses/MCGeneratedMasses_stealth_t6Wg_savedObjects.root".format(eP=stealthEnv.EOSPrefix)
}

for tDesignation in tDesignationsForProgenitor:
    print("Dumping contents of histogram in: {s}".format(s=templateSources[tDesignation]))
    templateReader = MCTemplateReader.MCTemplateReader(templateSources[tDesignation])
    for indexPair in templateReader.nextValidBin():
        eventProgenitorBinIndex = indexPair[0]
        eventProgenitorMass = (templateReader.eventProgenitorMasses)[eventProgenitorBinIndex]
        neutralinoBinIndex = indexPair[1]
        neutralinoMass = (templateReader.neutralinoMasses)[neutralinoBinIndex]
        print("(eventProgenitorMassBin, neutralinoMassBin)=({ePMB},{nMB}) is valid, corresponds to (eventProgenitorMass, neutralinoMass) = ({ePM}, {nM})".format(ePMB=eventProgenitorBinIndex, nMB=neutralinoBinIndex, ePM=eventProgenitorMass, nM=neutralinoMass))
