#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, re, subprocess, time

# Environment variables
stealthRoot = os.getenv("STEALTH_ROOT")
stealthEOSRoot = os.getenv("STEALTH_EOS_ROOT")
stealthArchives = os.getenv("STEALTH_ARCHIVES")
EOSPrefix = os.getenv("EOSPREFIX")
tmUtilsParent = os.getenv("TM_UTILS_PARENT")
hostname = os.getenv("HOSTNAME")
x509Proxy = os.getenv("X509_USER_PROXY")
habitat = ""
if ("lxplus" in hostname):
    habitat = "lxplus"
elif ("fnal" in hostname):
    habitat = "fnal"
else:
    sys.exit("ERROR: Unrecognized hostname: {h}, seems to be neither lxplus nor fnal.".format(h=hostname))

print("Environment variables:")
print("stealthRoot={sR}".format(sR=stealthRoot))
print("stealthEOSRoot={sER}".format(sER=stealthEOSRoot))
print("stealthArchives={sA}".format(sA=stealthArchives))
print("EOSPrefix={eP}".format(eP=EOSPrefix))
print("tmUtilsParent={tUP}".format(tUP=tmUtilsParent))
print("hostname={hN}".format(hN=hostname))
print("x509Proxy={xP}".format(xP=x509Proxy))

# Register command line options
inputArgumentsParser = argparse.ArgumentParser(description='Get signal contamination for the control region.')
inputArgumentsParser.add_argument('--controlSelection', default="combined", help="Control region to use. Can be \"fakefake\", \"mediumfake\", or \"combined\".", type=str)
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the output selection and statistics folders carry this suffix.',type=str)
inputArguments = inputArgumentsParser.parse_args()

def execute_in_env(commandToRun, printDebug=False):
    env_setup_command = "bash -c \"cd {sR} && source setupEnv.sh".format(sR=stealthRoot)
    runInEnv = "{e_s_c} && set -x && {c} && set +x\"".format(e_s_c=env_setup_command, c=commandToRun)
    if (printDebug):
        print("About to execute command:")
        print("{c}".format(c=runInEnv))
    os.system(runInEnv)

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)

control_pattern = "*"
if not(inputArguments.controlSelection == "combined"):
    control_pattern = inputArguments.controlSelection

patterns = {
    "data": "{eP}/{sER}/selections/combined_DoublePhoton{oI}/merged_selection_data_2017_control_{cP}.root".format(eP=EOSPrefix, sER=stealthEOSRoot, oI=optional_identifier, cP=control_pattern),
    "MC": "{eP}/{sER}/selections/combined_DoublePhoton{oI}/merged_selection_MC_stealth_t5_2017_control_{cP}.root".format(eP=EOSPrefix, sER=stealthEOSRoot, oI=optional_identifier, cP=control_pattern)
}

execute_in_env("mkdir -p signalContamination/{dataEventHistograms,dataSystematics,MCEventHistograms,MCSystematics,signalContamination}")

# Step 1: Build data event histograms
print("Analysing data...")
execute_in_env("./getDataEventHistogramsAndSystematics.py --inputFilesList \"{iFL}\" --outputDirectory_eventHistograms \"signalContamination/dataEventHistograms/\" --outputDirectory_dataSystematics \"signalContamination/dataSystematics/\" --outputPrefix control_{cS}{oI}".format(iFL=patterns["data"], cS=inputArguments.controlSelection, oI=optional_identifier))

# Step 2: Build MC event histograms
print("Analyzing MC...")
execute_in_env("./getMCSystematics/bin/getEventHistograms inputMCPathMain={iFL} integratedLuminosityMain=41900.0 outputDirectory=signalContamination/MCEventHistograms/ outputPrefix=control_{cS}{oI}".format(iFL=patterns["MC"], cS=inputArguments.controlSelection, oI=optional_identifier))

# Step 3: Combine outputs of steps 1 and 2 to get signal contamination
print("Getting signal contamination...")
execute_in_env("./getMCSystematics/bin/getMCUncertainties inputNEventsFile=signalContamination/dataSystematics/control_{cS}{oI}_observedEventCounters.dat inputPath=signalContamination/MCEventHistograms/control_{cS}{oI}_savedObjects.root outputPrefix=control_{cS}{oI} outputDirectory=signalContamination/MCSystematics outputDirectory_signalContamination=signalContamination/signalContamination unrestrictedSignalContamination=true minGluinoMass=975.0 nGluinoMassBins=16".format(cS=inputArguments.controlSelection, oI=optional_identifier))
