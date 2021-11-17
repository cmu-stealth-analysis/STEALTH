#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, subprocess, argparse

inputArgumentsParser = argparse.ArgumentParser(description='Print detailed information about a selected event.')
inputArgumentsParser.add_argument('--inputPath_with_prefix', required=True, action='append', help='Path to input n-tuple.', type=str)
inputArgumentsParser.add_argument('--ID_runLumiEvent', required=True, action='append', help='Event details for event to print, in format \"r:l:e\", where r, l, and e are respectively the run number, lumi number, and event ID.', type=str)
inputArguments = inputArgumentsParser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

events_to_print = []
for input_ID_runLumiEvent_raw in inputArguments.ID_runLumiEvent:
    input_ID_runLumiEvent_raw_split = input_ID_runLumiEvent_raw.split(":")
    if not(len(input_ID_runLumiEvent_raw_split) == 3): sys.exit("ERROR: Unexpected formatting for this value of ID_runLumiEvent: {v}".format(v=input_ID_runLumiEvent_raw))
    run = int(input_ID_runLumiEvent_raw_split[0])
    lumi = int(input_ID_runLumiEvent_raw_split[1])
    event_id = int(input_ID_runLumiEvent_raw_split[2])
    events_to_print.append((run, lumi, event_id))

inputChain = ROOT.TChain("ggNtuplizer/EventTree")
inputChain.SetMaxTreeSize(100000000000) # 1 TB
for inputPath in inputArguments.inputPath_with_prefix:
    inputChain.Add(inputPath)
nEntries = inputChain.GetEntries()

for eventIndex in range(0, nEntries):
    treeStatus = inputChain.LoadTree(eventIndex)
    if (treeStatus < 0):
        break
    evtStatus = inputChain.GetEntry(eventIndex)
    if (evtStatus <= 0):
        continue
    print_detailed_event_info = False
    for run_to_print, lumi_to_print, event_id_to_print in events_to_print:
        if (inputChain.event == event_id_to_print):
            if (inputChain.lumis == lumi_to_print):
                if (inputChain.run == run_to_print):
                    print_detailed_event_info = True
                    break
    if print_detailed_event_info:
        print("Printing detailed info for run: {r}, lumi: {l}, event_id: {e}".format(r=inputChain.run, l=inputChain.lumis, e=inputChain.event))
        print("genHT: {g}".format(g=inputChain.genHT))
        print("gen_pho1_et: {g}".format(g=inputChain.genPho1))
        print("gen_pho2_et: {g}".format(g=inputChain.genPho2))
        print("Total number of LHE-level particles: {n}".format(n=inputChain.nLHE))
        for LHE_index in range(inputChain.nLHE):
            print("At LHE_index: {i}, PDGID: {pdgid}".format(i=LHE_index, pdgid=inputChain.lhePID[LHE_index]))
        print("Total number of MC particles: {n}".format(n=inputChain.nMC))
        for MC_index in range(inputChain.nMC):
            if (inputChain.mcEt[MC_index] < 25.): continue
            status_flag_raw = inputChain.mcStatusFlag[MC_index]
            status_flag_fromHardProcessFinalState = ((status_flag_raw & 1) == 1)
            status_flag_isPromptFinalState = (((status_flag_raw >> 1) & 1) == 1)
            status_flag_isHardProcess = (((status_flag_raw >> 2) & 1) == 1)
            print("At index: {i}, found PDGID: {pdgid}, eta: {eta}, phi: {phi}, et: {et}, isPromptFinalState: {ipfs}, mom PDGID: {mompdgid}".format(i=MC_index, pdgid=inputChain.mcPID[MC_index], eta=inputChain.mcEta[MC_index], phi=inputChain.mcPhi[MC_index], et=inputChain.mcEt[MC_index], ipfs=status_flag_isPromptFinalState, mompdgid=inputChain.mcMomPID[MC_index]))
        print("Total number of photons: {n}".format(n=inputChain.nPho))
        for photon_index in range(inputChain.nPho):
            if (inputChain.phoEt[photon_index] < 25.): continue
            print("At index: {i}, eta: {eta}, phi: {phi}, pt: {pt}".format(i=photon_index, eta=inputChain.phoEta[photon_index], phi=inputChain.phoPhi[photon_index], pt=inputChain.phoEt[photon_index]))
        print("Selected photons info:")
        index_leading = inputChain.b_photonIndex_leading
        index_subLeading = inputChain.b_photonIndex_subLeading
        print("Leading photon: (index, pt, eta, phi) = ({i}, {pt}, {eta}, {phi})".format(i=index_leading, pt=inputChain.phoEt[index_leading], eta=inputChain.phoEta[index_leading], phi=inputChain.phoPhi[index_leading]))
        print("Subleading photon: (index, pt, eta, phi) = ({i}, {pt}, {eta}, {phi})".format(i=index_subLeading, pt=inputChain.phoEt[index_subLeading], eta=inputChain.phoEta[index_subLeading], phi=inputChain.phoPhi[index_subLeading]))

print("All done!")
