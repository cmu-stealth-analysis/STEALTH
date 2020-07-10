from __future__ import print_function, division

import os, sys

print("Importing environment variables...")
stealthRoot = os.getenv("STEALTH_ROOT")
stealthEOSRoot = os.getenv("STEALTH_EOS_ROOT")
stealthArchives = os.getenv("STEALTH_ARCHIVES")
stealthCMSSWBase = os.getenv("STEALTH_CMSSW_BASE")
EOSPrefix = os.getenv("EOSPREFIX")
tmUtilsParent = os.getenv("TM_UTILS_PARENT")
hostname = os.getenv("HOSTNAME")
x509Proxy = os.getenv("X509_USER_PROXY")
condorWorkAreaRoot = os.getenv("CONDORWORKAREAROOT")
analysisRoot = os.getenv("ANALYSISROOT")
habitat = ""
if ("lxplus" in hostname):
    habitat = "lxplus"
elif ("fnal" in hostname):
    habitat = "fnal"
else:
    sys.exit("ERROR: Unrecognized hostname: {h}, seems to be neither lxplus nor fnal.".format(h=hostname))

print("stealthRoot: {sR}".format(sR=stealthRoot))
print("stealthEOSRoot: {sER}".format(sER=stealthEOSRoot))
print("stealthArchives: {sA}".format(sA=stealthArchives))
print("stealthCMSSWBase: {sCB}".format(sCB=stealthCMSSWBase))
print("EOSPrefix: {eP}".format(eP=EOSPrefix))
print("tmUtilsParent: {tUP}".format(tUP=tmUtilsParent))
print("hostname: {hN}".format(hN=hostname))
print("x509Proxy: {xP}".format(xP=x509Proxy))
print("condorWorkAreaRoot: {cWAR}".format(cWAR=condorWorkAreaRoot))
print("analysisRoot: {aR}".format(aR=analysisRoot))
print("Setting habitat = {h}".format(h=habitat))

def get_execution_command(commandToRun):
    env_setup_command = "bash -c \"cd {sR} && source setupEnv.sh".format(sR=stealthRoot)
    return "{e_s_c} && set -x && {c} && set +x\"".format(e_s_c=env_setup_command, c=commandToRun)

def execute_in_env(commandToRun, isDryRun=False, functionToCallIfCommandExitsWithError=None):
    executionCommand = get_execution_command(commandToRun)
    if (isDryRun):
        print("Dry-run, not executing:")
        print("{c}".format(c=executionCommand))
    else:
        print("Executing:")
        print("{c}".format(c=executionCommand))
        returnCode = os.system(executionCommand)
        if (returnCode != 0):
            if not(functionToCallIfCommandExitsWithError is None):
                if not(callable(functionToCallIfCommandExitsWithError)): sys.exit("ERROR in execute_in_env: command exited with error and unable to call functionToCallIfCommandExitsWithError")
                else: functionToCallIfCommandExitsWithError()
            sys.exit("ERROR in execute_in_env: command \"{c}\" returned status {rC}".format(c=executionCommand, rC=returnCode))
