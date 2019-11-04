#!/usr/bin/env python

from __future__ import print_function, division

import ROOT, os, sys, argparse, tmJDLInterface, MCTemplateReader, stealthEnv

inputArgumentsParser = argparse.ArgumentParser(description='Run Higgs combine tool on already generated datacards.')
inputArgumentsParser.add_argument('--dataCardsDirectory', default="analysis/dataCards", help='Path to directory containing already generated datacards.',type=str)
inputArgumentsParser.add_argument('--dataCardsPrefix', default="", help='Data cards prefix.',type=str)
inputArgumentsParser.add_argument('--outputDirectory', default="root://cmseos.fnal.gov//store/user/lpcsusystealth/combineToolOutputs", help='EOS path on which to store combine tool outputs.',type=str)
inputArgumentsParser.add_argument('--MCTemplatePath', default="plot_susyMasses_template.root", help='Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.', type=str)
inputArgumentsParser.add_argument('--minGluinoMass', default=-1., help='Minimum gluino mass on which to run.', type=float)
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--isDryRun', action='store_true', help="Do not submit the actual jobs: instead, only print the shell command that would have been called.")
inputArguments = inputArgumentsParser.parse_args()

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)

# Make sure the CMSSW source tarball is the latest version
print("Updating CMSSW source tarball...")
updateCommand = "cd {sCB}/.. && ./uploadTarball.sh && cd {sR}".format(sCB=stealthEnv.stealthCMSSWBase, sR=stealthEnv.stealthRoot)
os.system(updateCommand)
# Copy event selection helper script into the working directory
copyCommand = "mkdir -p {cWAR}/combine{oI} && cd {sR} && cp -u combineToolHelper.sh {cWAR}/combine{oI}/.".format(sR=stealthEnv.stealthRoot, cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier)
os.system(copyCommand)

templateReader = MCTemplateReader.MCTemplateReader(inputArguments.MCTemplatePath)
for indexPair in templateReader.nextValidBin():
    gluinoMassBin = indexPair[0]
    gluinoMass = (templateReader.gluinoMasses)[gluinoMassBin]
    neutralinoMassBin = indexPair[1]
    neutralinoMass = (templateReader.neutralinoMasses)[neutralinoMassBin]
    print("gluino mass: {gM}, neutralino mass: {nM}".format(gM=gluinoMass, nM=neutralinoMass))
    if ((inputArguments.dataCardsDirectory)[0] == "/"): # inputArguments.dataCardsDirectory is likely an absolute path, not a relative path
        dataCardPathsPrefix = "{dCD}/{dCP}_dataCard_gluinoMassBin{gMB}_neutralinoMassBin{nMB}".format(dCD=inputArguments.dataCardsDirectory, dCP=inputArguments.dataCardsPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin)
    else:
        dataCardPathsPrefix = "{sR}/{dCD}/{dCP}_dataCard_gluinoMassBin{gMB}_neutralinoMassBin{nMB}".format(sR=stealthEnv.stealthRoot, dCD=inputArguments.dataCardsDirectory, dCP=inputArguments.dataCardsPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin)
    limitsConvergenceCheckScriptPath = "{sR}/checkLimitsConvergence.py".format(sR=stealthEnv.stealthRoot)
    filesToTransfer = ["{dCPP}.txt".format(dCPP=dataCardPathsPrefix), "{dCPP}_crossSectionsDown.txt".format(dCPP=dataCardPathsPrefix), "{dCPP}_crossSectionsUp.txt".format(dCPP=dataCardPathsPrefix), limitsConvergenceCheckScriptPath]
    processIdentifier = "combineJob_{prefix}_gluinoMassBin{gMB}_neutralinoMassBin{nMB}".format(prefix=inputArguments.dataCardsPrefix, gMB=gluinoMassBin, nMB=neutralinoMassBin)
    jdlInterface = tmJDLInterface.tmJDLInterface(processName=processIdentifier, scriptPath="combineToolHelper.sh", outputDirectoryRelativePath="{cWAR}/combine{oI}".format(cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier))  # works even if "outputDirectoryRelativePath" is an absolute path
    jdlInterface.addFilesToTransferFromList(filesToTransfer)
    # Arguments for script:
    jdlInterface.addScriptArgument("{oD}".format(oD=inputArguments.outputDirectory)) # Argument 1: output directory
    jdlInterface.addScriptArgument("{dCP}".format(dCP=inputArguments.dataCardsPrefix)) # Argument 2: data cards prefix
    jdlInterface.addScriptArgument("{gMB}".format(gMB=gluinoMassBin)) # Argument 3: gluino mass bin index
    jdlInterface.addScriptArgument("{nMB}".format(nMB=neutralinoMassBin)) # Argument 4: neutralino mass bin index
    # Write JDL
    jdlInterface.writeToFile()
    submissionCommand = "cd {cWAR}/combine{oI}/ && condor_submit {pI}.jdl && cd {sR}".format(pI=processIdentifier, cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier, sR=stealthEnv.stealthRoot)
    print ("Generated command: {sC}".format(sC=submissionCommand))
    if (inputArguments.isDryRun):
        print("Not submitting due to dryRun flag.")
    else:
        os.system(submissionCommand)
        print ("Submitted.")
