#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, subprocess, math
import ROOT
import stealthEnv, commonFunctions

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

inputArgumentsParser = argparse.ArgumentParser(description='Plot impacts of each nuisance parameter given a datacard.')
inputArgumentsParser.add_argument('--outputFolder', required=True, help='Path to folder in which to store output files.',type=str)
inputArgumentsParser.add_argument('--datacardTemplateParentFolderWithPrefix', required=True, help='Path to EOS folder (including xrootd prefix) from which to fetch datacard.',type=str)
inputArgumentsParser.add_argument('--datacardTemplateFileName', required=True, help='Name of datacard.',type=str)
inputArgumentsParser.add_argument('--identifier', required=True, help='Human-readable ID for the output.',type=str)
inputArguments = inputArgumentsParser.parse_args()

# If "identifier" contains anything other than a letter from the alphabet, it can't be used in a TeX macro
if not((inputArguments.identifier).isalpha()): sys.exit("ERROR: argument \"identifier\" can only contain letters.")

signal_strengths = [0, 1]
output_folders = {
    0: "{oF}/{i}/asimov_signal0".format(oF=inputArguments.outputFolder, i=inputArguments.identifier),
    1: "{oF}/{i}/asimov_signal1".format(oF=inputArguments.outputFolder, i=inputArguments.identifier)
}
output_signal_strength_identifiers = {
    0: "signal_zero",
    1: "signal_injected"
}
output_signal_strength_identifiers_tex = {
    0: "SignalZero",
    1: "SignalInjected"
}

for signal_strength in signal_strengths:
    # Make sure output folder exists
    if not(os.path.isdir(output_folders[signal_strength])): subprocess.check_call("mkdir -p {oF}".format(oF=output_folders[signal_strength]), shell=True, executable="/bin/bash")

    # Step 1: Copy datacard locally
    stealthEnv.execute_in_env(commandToRun="cd {oF} && xrdcp --nopbar --force --path {iF}/{f} {f}".format(oF=output_folders[signal_strength], iF=inputArguments.datacardTemplateParentFolderWithPrefix, f=inputArguments.datacardTemplateFileName))
    # Need to add the line "shapes * * FAKE" to "fake" a shape analysis, because that is required by the FitDiagnostics tool
    stealthEnv.execute_in_env(commandToRun="cd {oF} && sed -i '/number of nuisance parameters/a shapes * * FAKE' {f}".format(oF=output_folders[signal_strength], f=inputArguments.datacardTemplateFileName))

    # Step 2: Run FitDiagnostics
    stealthEnv.execute_in_env(commandToRun="cd {oF} && combine -M FitDiagnostics -t -1 --expectSignal {s} --robustFit 1 --rMin -10 -d {f}".format(oF=output_folders[signal_strength], s=signal_strength, f=inputArguments.datacardTemplateFileName))

    # Step 3: Get best-fit value of r and write it out to a TeX file
    r_bestfit, r_error_lo, r_error_hi = commonFunctions.get_r_from_fit_diagnostics_output(input_file_path="{oF}/fitDiagnostics.root".format(oF=output_folders[signal_strength]), printDebug=False)
    output_tex_file_r_bestfit = open("{oF}/r_best_fit_{i}_{si}.tex".format(oF=output_folders[signal_strength], i=inputArguments.identifier, si=output_signal_strength_identifiers[signal_strength]), 'w')
    output_tex_file_r_bestfit.write("\\providecommand{{\\BestFitR{i}{si}}}{{Best fit $r$: {r_bestfit:.4G}~~-{r_error_lo:.3f}/+{r_error_hi:.3f}~~(68\\% CL)}}\n".format(i=inputArguments.identifier, si=output_signal_strength_identifiers_tex[signal_strength], r_bestfit=r_bestfit, r_error_lo=r_error_lo, r_error_hi=r_error_hi))
    output_tex_file_r_bestfit.close()

    # Step 4: Get pulls of nuisances using the standard "diffNuisances.py" script
    stealthEnv.execute_in_env(commandToRun="cd {oF} && python ${{CMSSW_BASE}}/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a fitDiagnostics.root > {oF}/diffNuisances_{i}_{si}_raw.txt 2>&1".format(oF=output_folders[signal_strength], i=inputArguments.identifier, si=output_signal_strength_identifiers[signal_strength]))

    # Step 5: Write out output of diffNuisances verbatim into a TeX file
    commonFunctions.write_diffNuisances_output_into_tex_file(input_file_path="{oF}/diffNuisances_{i}_{si}_raw.txt".format(oF=output_folders[signal_strength], i=inputArguments.identifier, si=output_signal_strength_identifiers[signal_strength]), output_file_path="{oF}/diffNuisances_{i}_{si}.tex".format(oF=output_folders[signal_strength], i=inputArguments.identifier, si=output_signal_strength_identifiers[signal_strength]))

    # Step 6: Create Combine workspace from datacard
    stealthEnv.execute_in_env(commandToRun="cd {oF} && text2workspace.py {f} -m 125".format(oF=output_folders[signal_strength], f=inputArguments.datacardTemplateFileName))
    workspace_path = "{oF}/{f}".format(oF=output_folders[signal_strength], f=(inputArguments.datacardTemplateFileName).replace(".txt", ".root"))
    print("Produced workspace at: {w}".format(w=workspace_path))
    if not(os.path.exists(workspace_path)): sys.exit("ERROR: expected to find file at location {w}, but found none.".format(w=workspace_path))

    # Step 7: do initial fit to parameters of interest (just signal strength in our case) with --doInitialFit
    stealthEnv.execute_in_env(commandToRun="cd {oF} && combineTool.py -M Impacts -d {w} -m 125 --doInitialFit --robustFit 1 --expectSignal {s} -t -1 --rMin -10".format(oF=output_folders[signal_strength], w=workspace_path, s=signal_strength))

    # Step 8: fit for each nuisance parameter with --doFits
    # stealthEnv.execute_in_env(commandToRun="cd {oF} && combineTool.py -M Impacts -d {w} -m 125 --robustFit 1 --doFits --expectSignal {s} -t -1 --rMin -10 --parallel 12".format(oF=output_folders[signal_strength], w=workspace_path, s=signal_strength))
    stealthEnv.execute_in_env(commandToRun="cd {oF} && combineTool.py -M Impacts -d {w} -m 125 --robustFit 1 --doFits --expectSignal {s} -t -1 --rMin -10".format(oF=output_folders[signal_strength], w=workspace_path, s=signal_strength))

    # Step 9: Collate outputs and write impacts into a json file
    stealthEnv.execute_in_env(commandToRun="cd {oF} && combineTool.py -M Impacts -d {w} -m 125 --expectSignal {s} -t -1 --rMin -10 -o impacts.json".format(oF=output_folders[signal_strength], w=workspace_path, s=signal_strength))

    # Step 10: Make the impact plots
    stealthEnv.execute_in_env(commandToRun="cd {oF} && plotImpacts.py -i impacts.json -o impacts_{i}_{si} --label-size 0.04 --left-margin 0.55 --height 500 --per-page 21".format(oF=output_folders[signal_strength], i=inputArguments.identifier, si=output_signal_strength_identifiers[signal_strength]))

    # Step 11: Rerun FitDiagnostics with extra plots
    stealthEnv.execute_in_env(commandToRun="cd {oF} && combine -M FitDiagnostics -t -1 --expectSignal {s} --robustFit 1 --rMin -10 --saveWithUncertainties --saveOverallShapes --numToysForShapes 200 --plots -d {f}".format(oF=output_folders[signal_strength], s=signal_strength, f=inputArguments.datacardTemplateFileName))

    # Step 12: Save high-res versions of 2D correlation plots, and print important values
    commonFunctions.print_and_save_high_res_correlations(input_file_path="{oF}/fitDiagnostics.root".format(oF=output_folders[signal_strength]), output_folder=output_folders[signal_strength], suffix="{i}_{si}".format(i=inputArguments.identifier, si=output_signal_strength_identifiers[signal_strength]), list_correlations_to_save=["correlation_b", "correlation_s", "correlation_bins_b", "correlation_bins_s"])

print("All done!")
