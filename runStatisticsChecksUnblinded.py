#!/usr/bin/env python

from __future__ import print_function, division

import os, sys, argparse, subprocess, math, array
import ROOT
import tmHEPDataInterface
import stealthEnv, commonFunctions

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

inputArgumentsParser = argparse.ArgumentParser(description='Plot impacts of each nuisance parameter given a datacard.')
inputArgumentsParser.add_argument('--outputFolder', required=True, help='Path to folder in which to store output files.',type=str)
inputArgumentsParser.add_argument('--datacardTemplateParentFolderWithPrefix', required=True, help='Path to EOS folder (including xrootd prefix) from which to fetch datacard.',type=str)
inputArgumentsParser.add_argument('--datacardTemplateFileName', required=True, help='Name of datacard.',type=str)
inputArgumentsParser.add_argument('--identifier', required=True, help='Human-readable ID for the output.',type=str)
inputArgumentsParser.add_argument('--inputFile_STRegionBoundaries', default="STRegionBoundaries.dat", help='Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.', type=str)
inputArguments = inputArgumentsParser.parse_args()

# If "identifier" contains anything other than a letter from the alphabet, it can't be used in a TeX macro
if not((inputArguments.identifier).isalpha()): sys.exit("ERROR: argument \"identifier\" can only contain letters.")

output_folder = "{oF}/{i}".format(oF=inputArguments.outputFolder, i=inputArguments.identifier)
# Make sure output folder exists
if not(os.path.isdir(output_folder)): subprocess.check_call("mkdir -p {oF}".format(oF=output_folder), shell=True, executable="/bin/bash")

# Step 1: Copy datacard locally
stealthEnv.execute_in_env(commandToRun="cd {oF} && xrdcp --nopbar --force --path {iF}/{f} {f}".format(oF=output_folder, iF=inputArguments.datacardTemplateParentFolderWithPrefix, f=inputArguments.datacardTemplateFileName))
# Need to add the line "shapes * * FAKE" to "fake" a shape analysis, because that is required by the FitDiagnostics tool
stealthEnv.execute_in_env(commandToRun="cd {oF} && sed -i '/number of nuisance parameters/a shapes * * FAKE' {f}".format(oF=output_folder, f=inputArguments.datacardTemplateFileName))

# Step 2: Run FitDiagnostics
stealthEnv.execute_in_env(commandToRun="cd {oF} && combine -M FitDiagnostics --robustFit 1 --rMin -10 --expectSignal 0 --saveShapes --saveWithUncertainties --plots -d {f}".format(oF=output_folder, f=inputArguments.datacardTemplateFileName))

# Step 3: Get best-fit value of r and write it out to a TeX file
r_bestfit, r_error_lo, r_error_hi = commonFunctions.get_r_from_fit_diagnostics_output(input_file_path="{oF}/fitDiagnostics.root".format(oF=output_folder), printDebug=False)
output_tex_file_r_bestfit = open("{oF}/r_best_fit_{i}.tex".format(oF=output_folder, i=inputArguments.identifier), 'w')
output_tex_file_r_bestfit.write("\\providecommand{{\\BestFitR{i}}}{{Best fit $r$: {r_bestfit:.4G}~~-{r_error_lo:.3f}/+{r_error_hi:.3f}~~(68\\% CL)}}\n".format(i=inputArguments.identifier, r_bestfit=r_bestfit, r_error_lo=r_error_lo, r_error_hi=r_error_hi))
output_tex_file_r_bestfit.close()

# Step 4: Get pulls of nuisances using the standard "diffNuisances.py" script
stealthEnv.execute_in_env(commandToRun="cd {oF} && python ${{CMSSW_BASE}}/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a fitDiagnostics.root > {oF}/diffNuisances_{i}_raw.txt 2>&1".format(oF=output_folder, i=inputArguments.identifier))

# Step 5: Write out output of diffNuisances verbatim into a TeX file
commonFunctions.write_diffNuisances_output_into_tex_file(input_file_path="{oF}/diffNuisances_{i}_raw.txt".format(oF=output_folder, i=inputArguments.identifier), output_file_path="{oF}/diffNuisances_{i}.tex".format(oF=output_folder, i=inputArguments.identifier))

# Step 6: Create Combine workspace from datacard
stealthEnv.execute_in_env(commandToRun="cd {oF} && text2workspace.py {f} -m 125".format(oF=output_folder, f=inputArguments.datacardTemplateFileName))
workspace_path = "{oF}/{f}".format(oF=output_folder, f=(inputArguments.datacardTemplateFileName).replace(".txt", ".root"))
print("Produced workspace at: {w}".format(w=workspace_path))
if not(os.path.exists(workspace_path)): sys.exit("ERROR: expected to find file at location {w}, but found none.".format(w=workspace_path))

# Step 7: do initial fit to parameters of interest (just signal strength in our case) with --doInitialFit
stealthEnv.execute_in_env(commandToRun="cd {oF} && combineTool.py -M Impacts -d {w} -m 125 --doInitialFit --robustFit 1 --expectSignal 0 --rMin -10".format(oF=output_folder, w=workspace_path))

# Step 8: fit for each nuisance parameter with --doFits
stealthEnv.execute_in_env(commandToRun="cd {oF} && combineTool.py -M Impacts -d {w} -m 125 --robustFit 1 --doFits --expectSignal 0 --rMin -10 --parallel 12".format(oF=output_folder, w=workspace_path))

# Step 9: Collate outputs and write impacts into a json file
stealthEnv.execute_in_env(commandToRun="cd {oF} && combineTool.py -M Impacts -d {w} -m 125 --expectSignal 0 --rMin -10 -o impacts.json".format(oF=output_folder, w=workspace_path))

# Step 10: Make the impact plots
stealthEnv.execute_in_env(commandToRun="cd {oF} && plotImpacts.py -i impacts.json -o impacts_{i} --label-size 0.04 --left-margin 0.55 --height 500 --per-page 21".format(oF=output_folder, i=inputArguments.identifier))

# Step 11: Rerun FitDiagnostics with extra plots
stealthEnv.execute_in_env(commandToRun="cd {oF} && combine -M FitDiagnostics --robustFit 1 --rMin -10 --expectSignal 0 --saveWithUncertainties --saveOverallShapes --numToysForShapes 200 --plots -d {f}".format(oF=output_folder, f=inputArguments.datacardTemplateFileName))

# Step 12: Save high-res versions of 2D correlation plots, and print important values
commonFunctions.print_and_save_high_res_correlations(input_file_path="{oF}/fitDiagnostics.root".format(oF=output_folder), output_folder=output_folder, suffix=inputArguments.identifier, list_correlations_to_save=["correlation_b", "correlation_s", "correlation_bins_b", "correlation_bins_s"])

# Step 13: get covariance matrix of the b-only fit and save it to a HEPData-formatted yaml
STRegionsAxis = None
n_STBins = 0
with open(inputArguments.inputFile_STRegionBoundaries, 'r') as STRegionBoundariesFileObject:
    STBoundaries = []
    for STBoundaryString in STRegionBoundariesFileObject:
        if (STBoundaryString.strip()):
            STBoundary = float(STBoundaryString.strip())
            STBoundaries.append(STBoundary)
    STBoundaries.append(3500.0) # Instead of infinity
    n_STBins = len(STBoundaries) - 1
    STRegionsAxis = ROOT.TAxis(n_STBins, array.array('d', STBoundaries))
input_file_handle = ROOT.TFile.Open("{oF}/fitDiagnostics.root".format(oF=output_folder), "READ")
if ((input_file_handle.IsZombie() == ROOT.kTRUE) or not(input_file_handle.IsOpen() == ROOT.kTRUE)): sys.exit("ERROR: unable to open file at location {oF}/fitDiagnostics.root".format(oF=output_folder))
for hist_source, out_path_for_hepdata_yaml in [
        ("shapes_prefit/overall_total_covar", "{oF}/covariance_bins_prefit.yaml".format(oF=output_folder)),
        ("shapes_fit_b/overall_total_covar", "{oF}/covariance_bins_b_only_fit.yaml".format(oF=output_folder))
]:
    input_th2d_covariance = ROOT.TH2D()
    input_file_handle.GetObject(hist_source, input_th2d_covariance)
    data_for_hepdata_yaml = {
        'ST Bin (axis 1)': {
            'units': 'GeV',
            'data': []
        },
        'NJets Bin (axis 1)': {
            'units': None,
            'data': []
        },
        'ST Bin (axis 2)': {
            'units': 'GeV',
            'data': []
        },
        'NJets Bin (axis 2)': {
            'units': None,
            'data': []
        },
        'Covariance': {
            'units': None,
            'data': []
        }
    }
    indep_vars_for_hepdata_yaml = ['ST Bin (axis 1)', 'NJets Bin (axis 1)', 'ST Bin (axis 2)', 'NJets Bin (axis 2)']
    dep_vars_for_hepdata_yaml = ['Covariance']
    # validate binning of the covariance histogram
    if not input_th2d_covariance.GetXaxis().GetNbins() == 18:
        raise ValueError
    if not input_th2d_covariance.GetYaxis().GetNbins() == 18:
        raise ValueError
    i_bin = 1
    st_reg_n_jets_to_i_bin = {}
    for STRegionIndex in range(2, 1+STRegionsAxis.GetNbins()):
        for nJetsBin in range(4, 7):
            expected_name = "sST{s}J{n}_0".format(s=STRegionIndex, n=nJetsBin)
            if not (input_th2d_covariance.GetXaxis().GetBinLabel(i_bin) == expected_name):
                print('Expected bin label: {e}, got: {g}'.format(e=expected_name, g=input_th2d_covariance.GetXaxis().GetBinLabel(i_bin)))
                raise ValueError
            if not (input_th2d_covariance.GetYaxis().GetBinLabel(i_bin) == expected_name):
                print('Expected bin label: {e}, got: {g}'.format(e=expected_name, g=input_th2d_covariance.GetYaxis().GetBinLabel(i_bin)))
                raise ValueError
            st_reg_n_jets_to_i_bin[(STRegionIndex, nJetsBin)] = i_bin
            i_bin += 1
    # fill in the covariance data
    for STRegionIndex1 in range(2, 1+STRegionsAxis.GetNbins()):
        STMin1 = STRegionsAxis.GetBinLowEdge(STRegionIndex1)
        STLabel1 = None
        if STRegionIndex1 == STRegionsAxis.GetNbins():
            STLabel1 = ('> {v:.1f}'.format(v=STMin1), [])
        else:
            STMax1 = STRegionsAxis.GetBinUpEdge(STRegionIndex1)
            STLabel1 = (STMin1, STMax1, [])
        for nJetsBin1 in range(4, 7):
            ix = st_reg_n_jets_to_i_bin[(STRegionIndex1, nJetsBin1)]
            for STRegionIndex2 in range(2, 1+STRegionsAxis.GetNbins()):
                STMin2 = STRegionsAxis.GetBinLowEdge(STRegionIndex2)
                STLabel2 = None
                if STRegionIndex2 == STRegionsAxis.GetNbins():
                    STLabel2 = ('> {v:.1f}'.format(v=STMin2), [])
                else:
                    STMax2 = STRegionsAxis.GetBinUpEdge(STRegionIndex2)
                    STLabel2 = (STMin2, STMax2, [])
                for nJetsBin2 in range(4, 7):
                    iy = st_reg_n_jets_to_i_bin[(STRegionIndex2, nJetsBin2)]
                    cov = input_th2d_covariance.GetBinContent(ix, iy)
                    data_for_hepdata_yaml['ST Bin (axis 1)']['data'].append(STLabel1)
                    data_for_hepdata_yaml['NJets Bin (axis 1)']['data'].append((nJetsBin1, []))
                    data_for_hepdata_yaml['ST Bin (axis 2)']['data'].append(STLabel2)
                    data_for_hepdata_yaml['NJets Bin (axis 2)']['data'].append((nJetsBin2, []))
                    data_for_hepdata_yaml['Covariance']['data'].append((cov, []))
    tmHEPDataInterface.save_to_yaml(data_for_hepdata_yaml,
                                    indep_vars_for_hepdata_yaml,
                                    dep_vars_for_hepdata_yaml,
                                    out_path_for_hepdata_yaml)
input_file_handle.Close()
print("All done!")
