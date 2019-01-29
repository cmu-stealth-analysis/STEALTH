#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, os, sys, argparse, tmJDLInterface
inputArgumentsParser = argparse.ArgumentParser(description='Run Higgs combine tool on already generated datacards.')
inputArgumentsParser.add_argument('--dataCardsDirectory', default="analysis/dataCards", help='Path to directory containing already generated datacards.',type=str)
inputArgumentsParser.add_argument('--dataCardsPrefix', default="", help='Data cards prefix.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="root://cmseos.fnal.gov//store/user/lpcsusystealth/combineToolOutputs", help='EOS path on which to store combine tool outputs.',type=str)
inputArgumentsParser.add_argument('--MCTemplate', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--minGluinoMass', default=-1., help='Minimum gluino mass on which to run.', type=float)
inputArgumentsParser.add_argument('--isDryRun', action='store_true', help="Do not submit the actual jobs: instead, only print the shell command that would have been called.")
inputArguments = inputArgumentsParser.parse_args()

currentWorkingDirectory = os.getcwd()
# Make sure the CMSSW source tarball is the latest version
print("Updating CMSSW source tarball...")
updateCommand = "cd ~/private/stealth/cmssw && ./uploadTarball.sh && cd {cWD}".format(cWD=currentWorkingDirectory)
os.system(updateCommand)
# Copy event selection helper script into the working directory
copyCommand = "cp -u combineToolHelper.sh condor_working_directory/."
os.system(copyCommand)

generatedMCTemplate = ROOT.TFile(inputArguments.MCTemplate)
h_MCTemplate = generatedMCTemplate.Get("h_susyMasses_template")
for gluinoMassBin in range(1, 1+h_MCTemplate.GetXaxis().GetNbins()):
    gluinoMass = h_MCTemplate.GetXaxis().GetBinCenter(gluinoMassBin)
    if (inputArguments.minGluinoMass > 0 and gluinoMass < inputArguments.minGluinoMass): continue
    for neutralinoMassBin in range(1, 1+h_MCTemplate.GetYaxis().GetNbins()):
        if not(int(0.5 + h_MCTemplate.GetBinContent(gluinoMassBin, neutralinoMassBin)) == 1): continue
        neutralinoMass = h_MCTemplate.GetYaxis().GetBinCenter(neutralinoMassBin)
        print("gluino mass: {gM}, neutralino mass: {nM}".format(gM=gluinoMass, nM=neutralinoMass))
        dataCardPathsPrefix = "/uscms/home/tmudholk/private/stealth/STEALTH/{dCD}/{dCP}_dataCard_gluinoMassBin{gMB}_neutralinoMassBin{nMB}".format(dCD=inputArguments.dataCardsDirectory, dCP=inputArguments.dataCardsPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin)
        limitsConvergenceCheckScriptPath = "/uscms/home/tmudholk/private/stealth/STEALTH/checkLimitsConvergence.py"
        filesToTransfer = ["{dCPP}.txt".format(dCPP=dataCardPathsPrefix), "{dCPP}_crossSectionsDown.txt".format(dCPP=dataCardPathsPrefix), "{dCPP}_crossSectionsUp.txt".format(dCPP=dataCardPathsPrefix), limitsConvergenceCheckScriptPath]
        processIdentifier = "combineJob_{prefix}_gluinoMassBin{gMB}_neutralinoMassBin{nMB}".format(prefix=inputArguments.dataCardsPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin)
        jdlInterface = tmJDLInterface.tmJDLInterface(processName=processIdentifier, scriptPath="combineToolHelper.sh", outputDirectoryRelativePath="condor_working_directory")
        jdlInterface.addFilesToTransferFromList(filesToTransfer)
        # Arguments for script:
        jdlInterface.addScriptArgument("{oD}".format(oD=inputArguments.outputDirectory)) # Argument 1: output directory
        jdlInterface.addScriptArgument("{dCP}".format(dCP=inputArguments.dataCardsPrefix)) # Argument 2: data cards prefix
        jdlInterface.addScriptArgument("{gMB}".format(gMB=gluinoMassBin)) # Argument 3: gluino mass bin index
        jdlInterface.addScriptArgument("{nMB}".format(nMB=neutralinoMassBin)) # Argument 4: neutralino mass bin index
        # Write JDL
        jdlInterface.writeToFile()
        submissionCommand = "cd condor_working_directory && condor_submit {pI}.jdl && cd ..".format(pI=processIdentifier)
        print ("Generated command: {sC}".format(sC=submissionCommand))
        if (inputArguments.isDryRun):
            print("Not submitting due to dryRun flag.")
        else:
            os.system(submissionCommand)
            print ("Submitted.")
generatedMCTemplate.Close()
