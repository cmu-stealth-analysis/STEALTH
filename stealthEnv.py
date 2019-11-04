from __future__ import print_function, division

import os

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
