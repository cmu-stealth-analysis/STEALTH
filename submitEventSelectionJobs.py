#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, re, ROOT, tmJDLInterface, stealthEnv, commonFunctions

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Submit jobs for final event selection.')
inputArgumentsParser.add_argument('--selectionsToRun', default="data,MC,MC_GJet,MC_QCD", help="Comma-separated list of selections to run. Allowed: \"data\", \"data_singlephoton\", \"data_jetHT\" \"MC\", \"MC_EMEnrichedQCD\", \"MC_GJet\", \"MC_GJet_singlephoton\", \"MC_QCD\", \"MC_QCD_singlephoton\", or \"MC_hgg\". For MC selections, disable HLT photon trigger and enable additional MC selection. Default is \"data,MC,MC_GJet,MC_QCD\".", type=str)
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.", type=str)
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_selections', default="{sER}/selections/DoublePhoton".format(sER=stealthEnv.stealthEOSRoot), help='Output directory name in which to store event selections.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_statistics', default="{sER}/statistics/DoublePhoton".format(sER=stealthEnv.stealthEOSRoot), help='Output directory name in which to store statistics histograms.',type=str)
inputArgumentsParser.add_argument('--disableJetSelection', action='store_true', help="Disable jet selection.")
inputArgumentsParser.add_argument('--invertElectronVeto', action='store_true', help="Invert electron veto.")
inputArgumentsParser.add_argument('--isProductionRun', action='store_true', help="By default, this script does not submit the actual jobs and instead only prints the shell command that would have been called. Passing this switch will execute the commands.")
inputArgumentsParser.add_argument('--preserveLogs', action='store_true', help="By default, this script moves all event selection logs to the archives. This switch will keep the logs where they are (but they may be overwritten).")
inputArgumentsParser.add_argument('--preserveInputFileLists', action='store_true', help="By default, this script regenerates the input file lists to be fed to the statistics and selection merging scripts. This switch will preserve the input file lists.")
inputArguments = inputArgumentsParser.parse_args()

def execute_in_env(commandToRun, printDebug=False):
    env_setup_command = "bash -c \"cd {sR} && source setupEnv.sh".format(sR=stealthEnv.stealthRoot)
    runInEnv = "{e_s_c} && set -x && {c} && set +x\"".format(e_s_c=env_setup_command, c=commandToRun)
    if (printDebug):
        print("About to execute command:")
        print("{c}".format(c=runInEnv))
    os.system(runInEnv)

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)

if not(inputArguments.preserveLogs):
    os.system("set -x && mkdir -p {sA}/logs && rsync --quiet --progress -a {cWAR}/selection{oI}/ {sA}/logs/ && rm -rf {cWAR}/selection{oI}/* && set +x".format(cWAR=stealthEnv.condorWorkAreaRoot, sA=stealthEnv.stealthArchives, oI=optional_identifier))

os.system("mkdir -p {cWAR}/selection{oI}".format(cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier))

selectionTypesToRun = []
for inputSelectionToRun in (inputArguments.selectionsToRun.split(",")):
    if (inputSelectionToRun == "data"):
        selectionTypesToRun.append("data")
    elif (inputSelectionToRun == "data_singlephoton"):
        selectionTypesToRun.append("data_singlephoton")
    elif (inputSelectionToRun == "data_jetHT"):
        selectionTypesToRun.append("data_jetHT")
    elif (inputSelectionToRun == "MC"):
        selectionTypesToRun.append("MC_stealth_t5")
        selectionTypesToRun.append("MC_stealth_t6")
    elif (inputSelectionToRun == "MC_EMEnrichedQCD"):
        selectionTypesToRun.append("MC_EMEnrichedQCD")
    elif (inputSelectionToRun == "MC_GJet"):
        # selectionTypesToRun.append("MC_GJet")
        selectionTypesToRun.append("MC_GJet1")
        selectionTypesToRun.append("MC_GJet2")
        selectionTypesToRun.append("MC_GJet3")
        selectionTypesToRun.append("MC_GJet4")
        # selectionTypesToRun.append("MC_GJet5")
    elif (inputSelectionToRun == "MC_GJet_singlephoton"):
        # selectionTypesToRun.append("MC_GJet_singlephoton")
        selectionTypesToRun.append("MC_GJet_singlephoton1")
        selectionTypesToRun.append("MC_GJet_singlephoton2")
        selectionTypesToRun.append("MC_GJet_singlephoton3")
        selectionTypesToRun.append("MC_GJet_singlephoton4")
        # selectionTypesToRun.append("MC_GJet_singlephoton5")
    elif (inputSelectionToRun == "MC_QCD"):
        # selectionTypesToRun.append("MC_QCD")
        selectionTypesToRun.append("MC_QCD1")
        selectionTypesToRun.append("MC_QCD2")
        selectionTypesToRun.append("MC_QCD3")
        selectionTypesToRun.append("MC_QCD4")
        selectionTypesToRun.append("MC_QCD5")
        selectionTypesToRun.append("MC_QCD6")
    elif (inputSelectionToRun == "MC_QCD_singlephoton"):
        # selectionTypesToRun.append("MC_QCD_singlephoton")
        selectionTypesToRun.append("MC_QCD_singlephoton1")
        selectionTypesToRun.append("MC_QCD_singlephoton2")
        selectionTypesToRun.append("MC_QCD_singlephoton3")
        selectionTypesToRun.append("MC_QCD_singlephoton4")
        selectionTypesToRun.append("MC_QCD_singlephoton5")
        selectionTypesToRun.append("MC_QCD_singlephoton6")
    elif (inputSelectionToRun == "MC_hgg"):
        selectionTypesToRun.append("MC_hgg")
    else:
        sys.exit("ERROR: invalid value for argument \"selectionsToRun\": {v}".format(v=inputSelectionToRun))

yearsToRun = []
if (inputArguments.year == "2016"):
    yearsToRun.append(2016)
elif (inputArguments.year == "2017"):
    yearsToRun.append(2017)
elif (inputArguments.year == "2018"):
    yearsToRun.append(2018)
elif (inputArguments.year == "all"):
    yearsToRun.append(2016)
    yearsToRun.append(2017)
    yearsToRun.append(2018)
else:
    sys.exit("ERROR: invalid value for argument \"year\": {v}".format(v=inputArguments.year))

fileLists = {
    # "MC_hgg": {
    #     2017: "fileLists/inputFileList_MC_Fall17_hgg.txt"
    # },
    # "MC_stealth_t6": {
    #     2017: "fileLists/inputFileList_MC_Fall17_stealth_t6Wg.txt"
    # },
    "MC_stealth_t5": {
        2016: "fileLists/inputFileList_MC_Fall17_stealth_t5Wg.txt",
        2017: "fileLists/inputFileList_MC_Fall17_stealth_t5Wg.txt",
        2018: "fileLists/inputFileList_MC_Fall17_stealth_t5Wg.txt"
    },
    "MC_stealth_t6": {
        2016: "fileLists/inputFileList_MC_Fall17_stealth_t6Wg.txt",
        2017: "fileLists/inputFileList_MC_Fall17_stealth_t6Wg.txt",
        2018: "fileLists/inputFileList_MC_Fall17_stealth_t6Wg.txt"
    },
    "MC_EMEnrichedQCD": {
        2016: "fileLists/inputFileList_MC_Fall17_MC_DoubleEMEnrichedQCD.txt",
        2017: "fileLists/inputFileList_MC_Fall17_MC_DoubleEMEnrichedQCD.txt",
        2018: "fileLists/inputFileList_MC_Fall17_MC_DoubleEMEnrichedQCD.txt"
    },
    # "MC_GJet": {
    #     2016: "fileLists/inputFileList_MC_Summer16_GJet.txt"
    # },
    # "MC_GJet_singlephoton": {
    #     2016: "fileLists/inputFileList_MC_Summer16_GJet.txt"
    # },
    # "MC_GJet1": {
    #     2016: "fileLists/inputFileList_MC_Summer16_GJet1.txt"
    # },
    # "MC_GJet_singlephoton1": {
    #     2016: "fileLists/inputFileList_MC_Summer16_GJet1.txt"
    # },
    # "MC_GJet2": {
    #     2016: "fileLists/inputFileList_MC_Summer16_GJet2.txt"
    # },
    # "MC_GJet_singlephoton2": {
    #     2016: "fileLists/inputFileList_MC_Summer16_GJet2.txt"
    # },
    # "MC_GJet3": {
    #     2016: "fileLists/inputFileList_MC_Summer16_GJet3.txt"
    # },
    # "MC_GJet_singlephoton3": {
    #     2016: "fileLists/inputFileList_MC_Summer16_GJet3.txt"
    # },
    # "MC_GJet4": {
    #     2016: "fileLists/inputFileList_MC_Summer16_GJet4.txt"
    # },
    # "MC_GJet_singlephoton4": {
    #     2016: "fileLists/inputFileList_MC_Summer16_GJet4.txt"
    # },
    # "MC_GJet5": {
    #     2016: "fileLists/inputFileList_MC_Summer16_GJet5.txt"
    # },
    # "MC_GJet_singlephoton5": {
    #     2016: "fileLists/inputFileList_MC_Summer16_GJet5.txt"
    # },
    "MC_GJet1": {
        2017: "fileLists/inputFileList_MC_Fall17_GJet1.txt"
    },
    "MC_GJet_singlephoton1": {
        2017: "fileLists/inputFileList_MC_Fall17_GJet1.txt"
    },
    "MC_GJet2": {
        2017: "fileLists/inputFileList_MC_Fall17_GJet2.txt"
    },
    "MC_GJet_singlephoton2": {
        2017: "fileLists/inputFileList_MC_Fall17_GJet2.txt"
    },
    "MC_GJet3": {
        2017: "fileLists/inputFileList_MC_Fall17_GJet3.txt"
    },
    "MC_GJet_singlephoton3": {
        2017: "fileLists/inputFileList_MC_Fall17_GJet3.txt"
    },
    "MC_GJet4": {
        2017: "fileLists/inputFileList_MC_Fall17_GJet4.txt"
    },
    "MC_GJet_singlephoton4": {
        2017: "fileLists/inputFileList_MC_Fall17_GJet4.txt"
    },
    "MC_QCD": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD.txt"
    },
    "MC_QCD_singlephoton": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD.txt"
    },
    "MC_QCD1": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD1.txt"
    },
    "MC_QCD_singlephoton1": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD1.txt"
    },
    "MC_QCD2": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD2.txt"
    },
    "MC_QCD_singlephoton2": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD2.txt"
    },
    "MC_QCD3": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD3.txt"
    },
    "MC_QCD_singlephoton3": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD3.txt"
    },
    "MC_QCD4": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD4.txt"
    },
    "MC_QCD_singlephoton4": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD4.txt"
    },
    "MC_QCD5": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD5.txt"
    },
    "MC_QCD_singlephoton5": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD5.txt"
    },
    "MC_QCD6": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD6.txt"
    },
    "MC_QCD_singlephoton6": {
        2017: "fileLists/inputFileList_MC_Fall17_MC_QCD6.txt"
    },
    "MC_hgg": {
        2016: "fileLists/inputFileList_MC_Fall17_hgg.txt",
        2017: "fileLists/inputFileList_MC_Fall17_hgg.txt",
        2018: "fileLists/inputFileList_MC_Fall17_hgg.txt"
    },
    "data": {
        2016: "fileLists/inputFileList_data_DoubleEG_2016_ntuplizedOct2019.txt",
        2017: "fileLists/inputFileList_data_DoubleEG_2017_ntuplizedOct2019.txt",
        2018: "fileLists/inputFileList_data_EGamma_2018_ntuplizedOct2019.txt"
    },
    "data_singlephoton": {
        2016: "fileLists/inputFileList_data_JetHT_2016_ntuplizedDec2019.txt",
        2017: "fileLists/inputFileList_data_JetHT_2017_ntuplizedDec2019.txt",
        2018: "fileLists/inputFileList_data_JetHT_2018_ntuplizedDec2019.txt"
    },
    "data_jetHT": {
        2016: "fileLists/inputFileList_data_JetHT_2016_ntuplizedDec2019.txt",
        2017: "fileLists/inputFileList_data_JetHT_2017_ntuplizedDec2019.txt",
        2018: "fileLists/inputFileList_data_JetHT_2018_ntuplizedDec2019.txt"
    }
}

target_nFilesPerJob = {
    "MC_stealth_t5": {
        2016: 25,
        2017: 25,
        2018: 25
    },
    "MC_stealth_t6": {
        2016: 25,
        2017: 25,
        2018: 25
    },
    "MC_EMEnrichedQCD": {
        2016: 100,
        2017: 100,
        2018: 100
    },
    "MC_GJet": {
        2017: 100
    },
    "MC_GJet_singlephoton": {
        2017: 20
    },
    "MC_GJet1": {
        2017: 200
    },
    "MC_GJet_singlephoton1": {
        2017: 40
    },
    "MC_GJet2": {
        2017: 200
    },
    "MC_GJet_singlephoton2": {
        2017: 40
    },
    "MC_GJet3": {
        2017: 100
    },
    "MC_GJet_singlephoton3": {
        2017: 20
    },
    "MC_GJet4": {
        2017: 50
    },
    "MC_GJet_singlephoton4": {
        2017: 50
    },
    "MC_QCD": {
        2017: 100
    },
    "MC_QCD_singlephoton": {
        2017: 20
    },
    "MC_QCD1": {
        2017: 200
    },
    "MC_QCD_singlephoton1": {
        2017: 40
    },
    "MC_QCD2": {
        2017: 200
    },
    "MC_QCD_singlephoton2": {
        2017: 40
    },
    "MC_QCD3": {
        2017: 200
    },
    "MC_QCD_singlephoton3": {
        2017: 40
    },
    "MC_QCD4": {
        2017: 100
    },
    "MC_QCD_singlephoton4": {
        2017: 20
    },
    "MC_QCD5": {
        2017: 50
    },
    "MC_QCD_singlephoton5": {
        2017: 10
    },
    "MC_QCD6": {
        2017: 25
    },
    "MC_QCD_singlephoton6": {
        2017: 5
    },
    "MC_hgg": {
        2016: 1,
        2017: 1,
        2018: 1
    },
    "data": {
        2016: 150,
        2017: 150,
        2018: 200
    },
    "data_singlephoton": {
        2016: 150,
        2017: 150,
        2018: 200
    },
    "data_jetHT": {
        2016: 150,
        2017: 150,
        2018: 200
    }
}

execute_in_env("eos {eP} mkdir -p {oD}{oI}".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier), printDebug=True)
execute_in_env("eos {eP} mkdir -p {oD}{oI}".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_statistics, oI=optional_identifier), printDebug=True)

# Make sure the tarballs to transfer are up to date
updateCommand = "cd {tUP} && ./update_tmUtilsTarball.sh && cd {sR} && ./update_eventSelectionTarball.sh && cd {sR}".format(tUP=stealthEnv.tmUtilsParent, sR=stealthEnv.stealthRoot)
os.system(updateCommand)
# Copy event selection helper script into the working directory
copyCommand = "cd {sR} && cp -u eventSelectionHelper.sh {cWAR}/selection{oI}/.".format(sR=stealthEnv.stealthRoot, cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier)
os.system(copyCommand)

disableJetSelectionString = "false"
allJetsString = ""
if (inputArguments.disableJetSelection):
    disableJetSelectionString = "true"
    allJetsString = "_allJets"

invertElectronVetoString = "false"
electronVetoString = ""
if (inputArguments.invertElectronVeto):
    invertElectronVetoString = "true"
    electronVetoString = "_invertElectronVeto"

for selectionType in selectionTypesToRun:
    for year in yearsToRun:
        if ((bool(re.match(r"^MC_GJet[0-9]*$", selectionType))) or (bool(re.match(r"^MC_GJet_singlephoton[0-9]*$", selectionType)))):
            if (year != 2017): continue # The only reason we need these is to calculate scaling systematics
        if (((bool(re.match(r"^MC_QCD[0-9]*$", selectionType))) or (bool(re.match(r"^MC_QCD_singlephoton[0-9]*$", selectionType))))
            or (selectionType == "MC_EMEnrichedQCD")):
            if (year != 2017): continue # The only reason we need these is to calculate ID efficiencies
        if not(inputArguments.preserveInputFileLists):
            os.system("cd {sR} && rm fileLists/inputFileList_selections_{t}{aJS}{eVS}_{y}{oI}_*.txt && rm fileLists/inputFileList_statistics_{t}{aJS}{eVS}_{y}{oI}.txt".format(oI=optional_identifier, t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, sR=stealthEnv.stealthRoot))
        inputPathsFile = fileLists[selectionType][year]
        nFilesPerJob = target_nFilesPerJob[selectionType][year]
        print("Submitting jobs for year={y}, selection type={t}".format(y=year, t=selectionType))

        total_nLines = commonFunctions.get_number_of_lines_in_file(inputFilePath=inputPathsFile)
        print("Total available nLines: {n}".format(n=total_nLines))

        if not(total_nLines > 0):
            os.system("rm -f submitEventSelectionJobs.lock")
            sys.exit("ERROR: Found 0 lines in input path {p}.".format(p=inputPathsFile))

        filesToTransfer = ["{xP}".format(xP=stealthEnv.x509Proxy), "{tUP}/tmUtils.tar.gz".format(tUP=stealthEnv.tmUtilsParent), "{tUP}/extract_tmUtilsTarball.sh".format(tUP=stealthEnv.tmUtilsParent), "{sR}/eventSelection.tar.gz".format(sR=stealthEnv.stealthRoot), "{sR}/extract_eventSelectionTarball.sh".format(sR=stealthEnv.stealthRoot), "{sR}/{iPF}".format(sR=stealthEnv.stealthRoot, iPF=inputPathsFile), "{sR}/STRegionBoundaries.dat".format(sR=stealthEnv.stealthRoot)]
        formatted_iPF = (inputPathsFile.split("/"))[-1]

        startLine = 1
        endLine = 0
        while endLine < total_nLines:
            endLine = startLine + nFilesPerJob - 1
            isLastIteration = (endLine >= total_nLines)
            if isLastIteration:
                endLine = total_nLines
            processIdentifier = "selectionJob_{sT}_{y}_begin_{sL}_end_{eL}".format(sT=selectionType, y=year, sL=startLine, eL=endLine)
            jdlInterface = tmJDLInterface.tmJDLInterface(processName=processIdentifier, scriptPath="eventSelectionHelper.sh", outputDirectoryRelativePath="{cWAR}/selection{oI}".format(cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier)) # works even if "outputDirectoryRelativePath" is an absolute path
            jdlInterface.addFilesToTransferFromList(filesToTransfer)
            # Arguments for script:
            # Note: it seems simpler and certainly more readable to just include the "=" signs with the argument names, but I'm not sure whether that will be understood correctly by Condor
            jdlInterface.addScriptArgument("{iPF}".format(iPF=formatted_iPF)) # Argument 1: inputPathsFile
            jdlInterface.addScriptArgument("{sT}".format(sT=selectionType)) # Argument 2: selectionType
            jdlInterface.addScriptArgument("{dJS}".format(dJS=disableJetSelectionString)) # Argument 3: disableJetSelection
            jdlInterface.addScriptArgument("{sL}".format(sL=startLine)) # Argument 4: lineNumberStartInclusive
            jdlInterface.addScriptArgument("{eL}".format(eL=endLine)) # Argument 5: lineNumberEndInclusive
            jdlInterface.addScriptArgument("{y}".format(y=year)) # Argument 6: year
            jdlInterface.addScriptArgument("{iEVS}".format(iEVS=invertElectronVetoString)) # Argument 7: invertElectronVeto

            # Other arguments:
            jdlInterface.addScriptArgument("{eP}".format(eP=stealthEnv.EOSPrefix)) # Argument 8: EOS prefix
            jdlInterface.addScriptArgument("{oD}{oI}".format(oD=inputArguments.outputDirectory_selections, oI=optional_identifier)) # Argument 9: selections output folder path
            jdlInterface.addScriptArgument("{oD}{oI}".format(oD=inputArguments.outputDirectory_statistics, oI=optional_identifier)) # Argument 10: statistics output folder path

            if (stealthEnv.habitat == "lxplus"):
                jdlInterface.setFlavor("tomorrow")
            # Write JDL
            jdlInterface.writeToFile()

            submissionCommand = "cd {cWAR}/selection{oI}/ && condor_submit {pI}.jdl && cd {sR}".format(pI=processIdentifier, cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier, sR=stealthEnv.stealthRoot)
            print ("Generated command: {sC}".format(sC=submissionCommand))
            if (inputArguments.isProductionRun):
                os.system(submissionCommand)
                print ("Submitted.")
            else:
                print("Not submitting because isProductionRun flag was not set.")
            if not(inputArguments.preserveInputFileLists):
                if ((selectionType == "data_singlephoton") or
                    (bool(re.match(r"^MC_GJet_singlephoton[0-9]*$", selectionType))) or
                    (bool(re.match(r"^MC_QCD_singlephoton[0-9]*$", selectionType)))):
                    os.system("echo \"{eP}/{oD}{oI}/selection_{t}{aJS}{eVS}_{y}_control_singlemedium_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}{aJS}{eVS}_{y}{oI}_control_singlemedium.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, b=startLine, e=endLine))
                    os.system("echo \"{eP}/{oD}{oI}/selection_{t}{aJS}{eVS}_{y}_control_singleloose_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}{aJS}{eVS}_{y}{oI}_control_singleloose.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, b=startLine, e=endLine))
                    os.system("echo \"{eP}/{oD}{oI}/selection_{t}{aJS}{eVS}_{y}_control_singlefake_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}{aJS}{eVS}_{y}{oI}_control_singlefake.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, b=startLine, e=endLine))
                elif (not(selectionType == "data_jetHT")):
                    os.system("echo \"{eP}/{oD}{oI}/selection_{t}{aJS}{eVS}_{y}_signal_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}{aJS}{eVS}_{y}{oI}_signal.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, b=startLine, e=endLine))
                    os.system("echo \"{eP}/{oD}{oI}/selection_{t}{aJS}{eVS}_{y}_signal_loose_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}{aJS}{eVS}_{y}{oI}_signal_loose.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, b=startLine, e=endLine))
                    os.system("echo \"{eP}/{oD}{oI}/selection_{t}{aJS}{eVS}_{y}_control_fakefake_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}{aJS}{eVS}_{y}{oI}_control_fakefake.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, b=startLine, e=endLine))
                os.system("echo \"{eP}/{oD}{oI}/statistics_{t}{aJS}{eVS}_{y}_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_statistics_{t}{aJS}{eVS}_{y}{oI}.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_statistics, oI=optional_identifier, t=selectionType, aJS=allJetsString, eVS=electronVetoString, y=year, b=startLine, e=endLine))
            if isLastIteration: break
            startLine = 1+endLine
            if (startLine > total_nLines): break

os.system("rm -f submitEventSelectionJobs.lock")
