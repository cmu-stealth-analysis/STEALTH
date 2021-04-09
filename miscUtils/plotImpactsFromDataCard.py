#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, subprocess
import stealthEnv

inputArgumentsParser = argparse.ArgumentParser(description='Plot impacts of each nuisance parameter given a datacard.')
inputArgumentsParser.add_argument('--eos_path_to_txt_datacard', required=True, help='Path to datacard, code uses xrdcp to copy the datacard to a local file. Probably starts with: {s}'.format(s=stealthEnv.EOSPrefix),type=str)
inputArgumentsParser.add_argument('--outputFolder', default="/uscms_data/d3/tmudholk/analysisAreas/impacts", help='Path to folder in which to store output files.',type=str)
inputArgumentsParser.add_argument('--identifier', default=None, help='Identifier for output files.',type=str)
inputArguments = inputArgumentsParser.parse_args()

id_string = ""
if (inputArguments.identifier is not None): id_string = "_{i}".format(i=inputArguments.identifier)
name_argument = ""
if (inputArguments.identifier is not None): name_argument = " --name {i}".format(i=inputArguments.identifier)

if not(os.path.isdir(inputArguments.outputFolder)): subprocess.check_call("mkdir -p {oF}".format(oF=inputArguments.outputFolder), shell=True, executable="/bin/bash")

# Example workflow from Michael:
# $ text2workspace.py bkgfit_limit.txt -m 125
# $ combineTool.py -M Impacts -m 125 -d bkgfit_limit.root --doInitialFit --robustFit 1 --expectSignal 1 -t -1
# (-t -1 means it will use the Asimov dataset for the expected)
# $ for file in higgsCombine*.root; do mv "$file" "${file/123456.root/root}"; done
# $ combineTool.py -M Impacts -m 125 -d bkgfit_limit.root --robustFit 1 --doFits --parallel 10 --expectSignal 1 -t -1
# $ for file in higgsCombine*123456.root; do mv "$file" "${file/123456.root/root}"; done
# $ combineTool.py -M Impacts -m 125 -d bkgfit_limit.root -o impacts.json
# $ plotImpacts.py -i impacts.json -o impacts

# Instructions in the combine documentation:
# https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#nuisance-parameter-impacts

# Step 1: Copy datacard locally
stealthEnv.execute_in_env(commandToRun="cd {oF} && xrdcp --verbose --force --path {p} dataCard{i}.txt".format(oF=inputArguments.outputFolder, p=inputArguments.eos_path_to_txt_datacard, i=id_string))

# Step 2: Create Combine workspace from datacard
stealthEnv.execute_in_env(commandToRun="cd {oF} && text2workspace.py dataCard{i}.txt -m 125".format(oF=inputArguments.outputFolder, p=inputArguments.eos_path_to_txt_datacard, i=id_string))
workspacePath = "{oF}/dataCard{i}.root".format(oF=inputArguments.outputFolder, i=id_string)
if not(os.path.exists(workspacePath)): sys.exit("ERROR: expected to find file at location {wP}.root, but found none.".format(wP=workspacePath))

# Step 3: do initial fit to parameters of interest (just signal strength in our case) with --doInitialFit
stealthEnv.execute_in_env(commandToRun="cd {oF} && combineTool.py{na} -M Impacts -d dataCard{i}.root -m 125 --doInitialFit --robustFit 1 --expectSignal 0".format(oF=inputArguments.outputFolder, na=name_argument, i=id_string))

# Step 4: fit for each nuisance parameter with --doFits
stealthEnv.execute_in_env(commandToRun="cd {oF} && combineTool.py{na} -M Impacts -d dataCard{i}.root -m 125 --robustFit 1 --doFits --expectSignal 0 --parallel 8".format(oF=inputArguments.outputFolder, na=name_argument, i=id_string))

# Step 5: Collate outputs and write impacts into a json file
stealthEnv.execute_in_env(commandToRun="cd {oF} && combineTool.py{na} -M Impacts -d dataCard{i}.root -m 125 --expectSignal 0 -o impacts{i}.json".format(oF=inputArguments.outputFolder, na=name_argument, i=id_string))

# Step 6: Make the impact plots
stealthEnv.execute_in_env(commandToRun="cd {oF} && plotImpacts.py --per-page 18 --label-size 0.04 --left-margin 0.55 --height 300 -i impacts{i}.json -o impacts{i}".format(oF=inputArguments.outputFolder, i=id_string))
