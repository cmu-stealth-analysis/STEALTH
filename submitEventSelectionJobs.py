#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, re, json, math, subprocess
import ROOT
import tmJDLInterface, tmEOSUtils
import stealthEnv, commonFunctions

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Submit jobs for final event selection.')
inputArgumentsParser.add_argument('--selectionsToRun', default="data,MC,MC_DiPhotonJets,MC_GJetHT16,MC_GJetHT17,MC_GJetHT18,MC_HighHTQCD16,MC_HighHTQCD17,MC_HighHTQCD18,data_singlephoton,MC_DiPhotonJets_singlephoton,MC_GJetHT16_singlephoton,MC_GJetHT17_singlephoton,MC_GJetHT18_singlephoton,MC_HighHTQCD16_singlephoton,MC_HighHTQCD17_singlephoton,MC_HighHTQCD18_singlephoton", help="Comma-separated list of selections to run. Allowed: \"data\", \"data_singlephoton\", \"data_jetHT\", \"MC\", \"MC_DiPhotonJets\", \"MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)(16|17|18)(|_singlephoton)\", or \"MC_hgg\".", type=str)
inputArgumentsParser.add_argument('--year', default="all", help="Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.", type=str)
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_selections', default="{sER}/selections/DoublePhoton".format(sER=stealthEnv.stealthEOSRoot), help='Output directory name in which to store event selections.',type=str)
inputArgumentsParser.add_argument('--outputDirectory_statistics', default="{sER}/statistics/DoublePhoton".format(sER=stealthEnv.stealthEOSRoot), help='Output directory name in which to store statistics histograms.',type=str)
inputArgumentsParser.add_argument('--disablePhotonSelection', action='store_true', help="Disable photon selection.")
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
    os.system("set -x && rm -rf {cWAR}/selection{oI} && set +x".format(cWAR=stealthEnv.condorWorkAreaRoot, sA=stealthEnv.stealthArchives, oI=optional_identifier))

os.system("mkdir -p {cWAR}/selection{oI}".format(cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier))

subprocess.check_call("for WGTS_FNAME in `eos {ep} ls {ser}/MCWeights/*.json`; do echo ${{WGTS_FNAME}} && xrdcp --silent --nopbar --force --path {ep}/{ser}/MCWeights/${{WGTS_FNAME}} xSecLumiInfo/${{WGTS_FNAME}}; done".format(ep=stealthEnv.EOSPrefix, ser=stealthEnv.stealthEOSRoot), executable="/bin/bash", shell=True)

n_subsamples = {
    "MC_EMEnrichedGJetPt16": 3,
    "MC_EMEnrichedGJetPt17": 3,
    "MC_EMEnrichedGJetPt18": 3,
    "MC_HighHTQCD16": 7,
    "MC_HighHTQCD17": 8,
    "MC_HighHTQCD18": 8,
    "MC_GJetHT16": 5,
    "MC_GJetHT17": 4,
    "MC_GJetHT18": 4
}

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
    elif (inputSelectionToRun == "MC_hgg"):
        selectionTypesToRun.append("MC_hgg")
    elif ((inputSelectionToRun == "MC_DiPhotonJets") or (inputSelectionToRun == "MC_DiPhotonJets_singlephoton")):
        selectionTypesToRun.append(inputSelectionToRun)
    else:
        MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)(|_singlephoton)$", inputSelectionToRun)
        if MCBKGMatch:
            full_match = MCBKGMatch.group(0)
            dataset_id = MCBKGMatch.group(1)
            year_last_two_digits_str = MCBKGMatch.group(2)
            singlephoton_match = MCBKGMatch.group(3)
            for index_subsample in range(1, 1+n_subsamples["MC_" + dataset_id + year_last_two_digits_str]):
                selectionTypesToRun.append("{m}_{i}".format(m=full_match, i=index_subsample))
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
    "MC_stealth_t5": {
        2016: ("fileLists/inputFileList_MC_Fall17_stealth_t5Wg.txt", "{ep}/{ser}/MCWeights/PUWeights_stealth_t5_2016.root".format(ep=stealthEnv.EOSPrefix, ser=stealthEnv.stealthEOSRoot)),
        2017: ("fileLists/inputFileList_MC_Fall17_stealth_t5Wg.txt", "{ep}/{ser}/MCWeights/PUWeights_stealth_t5_2017.root".format(ep=stealthEnv.EOSPrefix, ser=stealthEnv.stealthEOSRoot)),
        2018: ("fileLists/inputFileList_MC_Fall17_stealth_t5Wg.txt", "{ep}/{ser}/MCWeights/PUWeights_stealth_t5_2018.root".format(ep=stealthEnv.EOSPrefix, ser=stealthEnv.stealthEOSRoot))
    },
    "MC_stealth_t6": {
        2016: ("fileLists/inputFileList_MC_Fall17_stealth_t6Wg.txt", "{ep}/{ser}/MCWeights/PUWeights_stealth_t6_2016.root".format(ep=stealthEnv.EOSPrefix, ser=stealthEnv.stealthEOSRoot)),
        2017: ("fileLists/inputFileList_MC_Fall17_stealth_t6Wg.txt", "{ep}/{ser}/MCWeights/PUWeights_stealth_t6_2017.root".format(ep=stealthEnv.EOSPrefix, ser=stealthEnv.stealthEOSRoot)),
        2018: ("fileLists/inputFileList_MC_Fall17_stealth_t6Wg.txt", "{ep}/{ser}/MCWeights/PUWeights_stealth_t6_2018.root".format(ep=stealthEnv.EOSPrefix, ser=stealthEnv.stealthEOSRoot))
    },
    "MC_hgg": {
        2016: "fileLists/inputFileList_MC_Summer16_hgg.txt",
        2017: "fileLists/inputFileList_MC_Fall17_hgg.txt",
        2018: "fileLists/inputFileList_MC_Autumn18_hgg.txt"
    },
    "data": {
        2016: "fileLists/inputFileList_data_DoubleEG_2016_ntuplizedOct2019.txt",
        2017: "fileLists/inputFileList_data_DoubleEG_2017_ntuplizedOct2019.txt",
        2018: "fileLists/inputFileList_data_EGamma_2018_ntuplizedOct2019.txt"
    },
    "data_singlephoton": {
        2016: "fileLists/inputFileList_data_SinglePhoton_2016_ntuplizedFeb2021.txt",
        2017: "fileLists/inputFileList_data_SinglePhoton_2017_ntuplizedFeb2021.txt",
        2018: "fileLists/inputFileList_data_EGamma_2018_ntuplizedFeb2021.txt"
    },
    "data_jetHT": {
        2016: "fileLists/inputFileList_data_JetHT_2016_ntuplizedDec2019.txt",
        2017: "fileLists/inputFileList_data_JetHT_2017_ntuplizedDec2019.txt",
        2018: "fileLists/inputFileList_data_JetHT_2018_ntuplizedDec2019.txt"
    }
}
fileLists["MC_DiPhotonJets"] = {}
fileLists["MC_DiPhotonJets_singlephoton"] = {}
for year_last_two_digits in [16, 17, 18]:
    year = 2000 + year_last_two_digits
    fileLists["MC_DiPhotonJets"][year] = ("fileLists/inputFileList_MC_DiPhotonJets_{y}.txt".format(y=year), "{ep}/{ser}/MCWeights/PUWeights_DiPhotonJets_{y}.root".format(ep=stealthEnv.EOSPrefix, ser=stealthEnv.stealthEOSRoot, y=year), "xSecLumiInfo/xsec_DiPhotonJets_{y}.json".format(y=year), "xSecLumiInfo/sumMCWeights_DiPhotonJets_{y}.json".format(y=year))
    fileLists["MC_DiPhotonJets_singlephoton"][year] = ("fileLists/inputFileList_MC_DiPhotonJets_{y}.txt".format(y=year), "{ep}/{ser}/MCWeights/PUWeights_DiPhotonJets_{y}.root".format(ep=stealthEnv.EOSPrefix, ser=stealthEnv.stealthEOSRoot, y=year), "xSecLumiInfo/xsec_DiPhotonJets_{y}.json".format(y=year), "xSecLumiInfo/sumMCWeights_DiPhotonJets_{y}.json".format(y=year))

for year_last_two_digits in [16, 17, 18]:
    year = 2000 + year_last_two_digits
    for MCBKGDatasetID in ["EMEnrichedGJetPt", "HighHTQCD", "GJetHT"]:
        for index_subsample in range(1, 1+n_subsamples["MC_{did}{y2}".format(did=MCBKGDatasetID, y2=year_last_two_digits)]):
            fileLists["MC_{did}{y2}_{i}".format(did=MCBKGDatasetID, y2=year_last_two_digits, i=index_subsample)] = {
                year: ("fileLists/inputFileList_MC_{did}{i}_{y}.txt".format(did=MCBKGDatasetID, i=index_subsample, y=year), "{ep}/{ser}/MCWeights/PUWeights_{did}{i}_{y}.root".format(ep=stealthEnv.EOSPrefix, ser=stealthEnv.stealthEOSRoot, did=MCBKGDatasetID, i=index_subsample, y=year), "xSecLumiInfo/xsec_{did}_{y}_{i}.json".format(did=MCBKGDatasetID, y=year, i=index_subsample), "xSecLumiInfo/sumMCWeights_{did}{i}_{y}.json".format(did=MCBKGDatasetID, i=index_subsample, y=year))
            }
            fileLists["MC_{did}{y2}_singlephoton_{i}".format(did=MCBKGDatasetID, y2=year_last_two_digits, i=index_subsample)] = {
                year: ("fileLists/inputFileList_MC_{did}{i}_{y}.txt".format(did=MCBKGDatasetID, i=index_subsample, y=year), "{ep}/{ser}/MCWeights/PUWeights_{did}{i}_{y}.root".format(ep=stealthEnv.EOSPrefix, ser=stealthEnv.stealthEOSRoot, did=MCBKGDatasetID, i=index_subsample, y=year), "xSecLumiInfo/xsec_{did}_{y}_{i}.json".format(did=MCBKGDatasetID, y=year, i=index_subsample), "xSecLumiInfo/sumMCWeights_{did}{i}_{y}.json".format(did=MCBKGDatasetID, i=index_subsample, y=year))
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
        2016: 30,
        2017: 30,
        2018: 20
    },
    "data_jetHT": {
        2016: 150,
        2017: 150,
        2018: 200
    }
}

target_nFilesPerJob["MC_DiPhotonJets"] = {}
target_nFilesPerJob["MC_DiPhotonJets_singlephoton"] = {}
for year_last_two_digits in [16, 17, 18]:
    year = 2000 + year_last_two_digits
    target_nFilesPerJob["MC_DiPhotonJets"][year] = 50
    target_nFilesPerJob["MC_DiPhotonJets_singlephoton"][year] = 10

for year_last_two_digits in [16, 17, 18]:
    year = 2000 + year_last_two_digits
    for MCBKGDataset in ["MC_EMEnrichedGJetPt", "MC_HighHTQCD", "MC_GJetHT"]:
        for index_subsample in range(1, 1+n_subsamples["{d}{y2}".format(d=MCBKGDataset, y2=year_last_two_digits)]):
            target_nFilesPerJob["{d}{y2}_{i}".format(d=MCBKGDataset, y2=year_last_two_digits, i=index_subsample)] = {year: 50}
            target_nFilesPerJob["{d}{y2}_singlephoton_{i}".format(d=MCBKGDataset, y2=year_last_two_digits, i=index_subsample)] = {year: 10}

def get_nFiles_per_job(selectionType, year):
    nFilesPerJob = target_nFilesPerJob[selectionType][year]
    if (inputArguments.disablePhotonSelection or inputArguments.disableJetSelection): nFilesPerJob = int(0.5 + max(1.0, nFilesPerJob/5.0))
    if ((selectionType == "MC_DiPhotonJets") and (year == 2016)): nFilesPerJob = int(0.5 + max(1.0, nFilesPerJob/5.0))
    MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)(|_singlephoton)_([0-9]*)$", selectionType)
    if MCBKGMatch:
        full_match = MCBKGMatch.group(0)
        dataset_id = MCBKGMatch.group(1)
        year_last_two_digits_str = MCBKGMatch.group(2)
        singlephoton_match = MCBKGMatch.group(3)
        if (((dataset_id == "HighHTQCD") or (dataset_id == "GJetHT")) and (year_last_two_digits_str == "16")):
            nFilesPerJob = int(0.5 + max(1.0, nFilesPerJob/5.0))
    return nFilesPerJob

execute_in_env("eos {eP} mkdir -p {oD}{oI}".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier), printDebug=True)
execute_in_env("eos {eP} mkdir -p {oD}{oI}".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_statistics, oI=optional_identifier), printDebug=True)

# Make sure the tarballs to transfer are up to date
updateCommand = "cd {tUP} && ./update_tmUtilsTarball.sh && cd {sR} && ./update_eventSelectionTarball.sh && cd {sR}".format(tUP=stealthEnv.tmUtilsParent, sR=stealthEnv.stealthRoot)
os.system(updateCommand)
# Copy event selection helper script into the working directory
copyCommand = "cd {sR} && cp -u eventSelectionHelper.sh {cWAR}/selection{oI}/.".format(sR=stealthEnv.stealthRoot, cWAR=stealthEnv.condorWorkAreaRoot, oI=optional_identifier)
os.system(copyCommand)

disablePhotonSelectionString = "false"
photonSelectionString = ""
if (inputArguments.disablePhotonSelection):
    disablePhotonSelectionString = "true"
    photonSelectionString = "_noPhotonSelection"

disableJetSelectionString = "false"
jetSelectionString = ""
if (inputArguments.disableJetSelection):
    disableJetSelectionString = "true"
    jetSelectionString = "_noJetSelection"

invertElectronVetoString = "false"
electronVetoString = ""
if (inputArguments.invertElectronVeto):
    invertElectronVetoString = "true"
    electronVetoString = "_invertElectronVeto"

overallIdentificationString = "{pSS}{jSS}{eVS}".format(pSS=photonSelectionString, jSS=jetSelectionString, eVS=electronVetoString)

for selectionType in selectionTypesToRun:
    for year in yearsToRun:
        MCBKGMatch = re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)(|_singlephoton)_([0-9]*)$", selectionType)
        if MCBKGMatch:
            year_last_two_digits_str = MCBKGMatch.group(2)
            year_MCBKGSample = 2000+int(0.5 + float(year_last_two_digits_str))
            if (year != year_MCBKGSample): continue
        fileListsInputPathsSource = fileLists[selectionType][year]
        inputPathsFile = None
        PUWeightsPath = "/dev/null"
        MCWeight = -1.0
        MCWeightPrecision = 6
        if isinstance(fileListsInputPathsSource, tuple):
            if (len(fileListsInputPathsSource) == 2):
                inputPathsFile, PUWeightsPath = fileListsInputPathsSource
            elif (len(fileListsInputPathsSource) == 4):
                inputPathsFile, PUWeightsPath, MCXSecInfoFile, MCSumWeightsFile = fileListsInputPathsSource
                cms_year_lumi = None
                xsec = None
                n_gen_events_raw_from_xsec_json = None
                n_gen_events_raw_from_sum_mcweights_json = None
                n_gen_events_weighted = None
                with open("xSecLumiInfo/lumi_run2.json", 'r') as lumi_json_file_handle:
                    lumi_values_raw_json = json.load(lumi_json_file_handle)
                    cms_year_lumi = lumi_values_raw_json[str(year)] # in inv pb
                with open(MCXSecInfoFile, 'r') as xsec_json_file_handle:
                    xsec_values_raw_json = json.load(xsec_json_file_handle)
                    xsec = xsec_values_raw_json["xsec"] # in pb
                    n_gen_events_raw_from_xsec_json = xsec_values_raw_json["nevents"]
                with open(MCSumWeightsFile, 'r') as sum_weights_json_file_handle:
                    sum_weights_raw_json = json.load(sum_weights_json_file_handle)
                    n_gen_events_raw_from_sum_mcweights_json = sum_weights_raw_json["total_nevts_raw"]
                    n_gen_events_weighted = sum_weights_raw_json["total_nevts_mc_weighted"]
                if not(n_gen_events_raw_from_xsec_json == n_gen_events_raw_from_sum_mcweights_json): sys.exit("ERROR: inconsistent number of events between the files {f1} and {f2}".format(f1=MCXSecInfoFile, f2=MCSumWeightsFile))
                MCWeight = (xsec*cms_year_lumi)/(1.0*n_gen_events_weighted)
                MCWeightPrecision = 6 + int(0.5 + max(0., math.log10(1.0/MCWeight)))
            else:
                sys.exit("ERROR: fileListsInputPathsSource in unexpected format: {f}".format(f=fileListsInputPathsSource))
        elif isinstance(fileListsInputPathsSource, basestring):
            inputPathsFile = fileListsInputPathsSource
        else:
            sys.exit("ERROR: fileListsInputPathsSource is neither a tuple nor a string. Its str representation is: {s}".format(s=str(fileListsInputPathsSource)))
        # nFilesPerJob = target_nFilesPerJob[selectionType][year]
        # if (inputArguments.disablePhotonSelection): nFilesPerJob = int(0.5 + max(1.0, nFilesPerJob/5.0))
        nFilesPerJob = get_nFiles_per_job(selectionType, year)
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
        isFirstIteration = True
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
            jdlInterface.addScriptArgument("{dPS}".format(dPS=disablePhotonSelectionString)) # Argument 3: disablePhotonSelection
            jdlInterface.addScriptArgument("{dJS}".format(dJS=disableJetSelectionString)) # Argument 4: disableJetSelection
            jdlInterface.addScriptArgument("{sL}".format(sL=startLine)) # Argument 5: lineNumberStartInclusive
            jdlInterface.addScriptArgument("{eL}".format(eL=endLine)) # Argument 6: lineNumberEndInclusive
            jdlInterface.addScriptArgument("{y}".format(y=year)) # Argument 7: year
            jdlInterface.addScriptArgument("{iEVS}".format(iEVS=invertElectronVetoString)) # Argument 8: invertElectronVeto
            jdlInterface.addScriptArgument(("{w:." + str(MCWeightPrecision)+ "f}").format(w=MCWeight)) # Argument 9: MCBackgroundWeight
            jdlInterface.addScriptArgument("{pwp}".format(pwp=PUWeightsPath)) # Argument 10: PUWeightsPathWithXRDPrefix

            # Other arguments:
            jdlInterface.addScriptArgument("{eP}".format(eP=stealthEnv.EOSPrefix)) # Argument 11: EOS prefix
            jdlInterface.addScriptArgument("{oD}{oI}".format(oD=inputArguments.outputDirectory_selections, oI=optional_identifier)) # Argument 12: selections output folder path
            jdlInterface.addScriptArgument("{oD}{oI}".format(oD=inputArguments.outputDirectory_statistics, oI=optional_identifier)) # Argument 13: statistics output folder path

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
                    (bool(re.match(r"^MC_DiPhotonJets_singlephoton$", selectionType))) or
                    (bool(re.match(r"^MC_(EMEnrichedGJetPt|HighHTQCD|GJetHT)([0-9]*)_singlephoton_([0-9]*)$", selectionType)))):
                    if (inputArguments.disablePhotonSelection):
                        if isFirstIteration:
                            os.system("rm -f fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_unified.txt && touch fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_unified.txt".format(t=selectionType, oIS=overallIdentificationString, y=year, oI=optional_identifier))
                        os.system("echo \"{eP}/{oD}{oI}/selection_{t}{oIS}_{y}_unified_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_unified.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, oIS=overallIdentificationString, y=year, b=startLine, e=endLine))
                    else:
                        if isFirstIteration:
                            for tempstring in ["medium", "loose", "fake"]:
                                os.system("rm -f fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_control_single{ts}.txt && touch fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_control_single{ts}.txt".format(t=selectionType, oIS=overallIdentificationString, y=year, oI=optional_identifier, ts=tempstring))
                        for tempstring in ["medium", "loose", "fake"]:
                            os.system("echo \"{eP}/{oD}{oI}/selection_{t}{oIS}_{y}_control_single{ts}_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_control_single{ts}.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, oIS=overallIdentificationString, y=year, b=startLine, e=endLine, ts=tempstring))
                elif (not(selectionType == "data_jetHT")):
                    if (inputArguments.disablePhotonSelection):
                        if isFirstIteration:
                            os.system("rm -f fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_unified.txt && touch fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_unified.txt".format(t=selectionType, oIS=overallIdentificationString, y=year, oI=optional_identifier))
                        os.system("echo \"{eP}/{oD}{oI}/selection_{t}{oIS}_{y}_unified_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_unified.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, oIS=overallIdentificationString, y=year, b=startLine, e=endLine))
                    else:
                        if isFirstIteration:
                            for tempstring in ["signal", "signal_loose", "control_fakefake"]:
                                os.system("rm -f fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_{ts}.txt && touch fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_{ts}.txt".format(t=selectionType, oIS=overallIdentificationString, y=year, oI=optional_identifier, ts=tempstring))
                        for tempstring in ["signal", "signal_loose", "control_fakefake"]:
                            os.system("echo \"{eP}/{oD}{oI}/selection_{t}{oIS}_{y}_{ts}_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_selections_{t}{oIS}_{y}{oI}_{ts}.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_selections, oI=optional_identifier, t=selectionType, oIS=overallIdentificationString, y=year, b=startLine, e=endLine, ts=tempstring))
                if isFirstIteration:
                    os.system("rm -f fileLists/inputFileList_statistics_{t}{oIS}_{y}{oI}.txt && touch fileLists/inputFileList_statistics_{t}{oIS}_{y}{oI}.txt".format(t=selectionType, oIS=overallIdentificationString, y=year, oI=optional_identifier))
                os.system("echo \"{eP}/{oD}{oI}/statistics_{t}{oIS}_{y}_begin_{b}_end_{e}.root\" >> fileLists/inputFileList_statistics_{t}{oIS}_{y}{oI}.txt".format(eP=stealthEnv.EOSPrefix, oD=inputArguments.outputDirectory_statistics, oI=optional_identifier, t=selectionType, oIS=overallIdentificationString, y=year, b=startLine, e=endLine))
            isFirstIteration = False
            if isLastIteration: break
            startLine = 1+endLine
            if (startLine > total_nLines): break

os.system("rm -f submitEventSelectionJobs.lock")
