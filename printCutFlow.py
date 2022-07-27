#!/usr/bin/env python

from __future__ import print_function, division

import subprocess, os, sys, argparse
import stealthEnv

inputArgumentsParser = argparse.ArgumentParser(description='Print cut-flow to output tex file.')
inputArgumentsParser.add_argument('--selection', choices=['data', 'MC_stealth_t5', 'MC_stealth_t6'], required=True, help='selection type: data or MC.',type=str)
inputArgumentsParser.add_argument('--optionalIdentifier', default="", help='If set, the input statistics files are read from a folder with this suffix.',type=str)
inputArguments = inputArgumentsParser.parse_args()

event_cuts = ["HLTSelection", "MCGenInformation", "overlap", "doublePhoton", "invariantMass", "NJets", "ST"]
cuts_masked = None
if (inputArguments.selection == 'data'):
    cuts_masked = set(["MCGenInformation", "overlap"])
elif ((inputArguments.selection == 'MC_stealth_t5') or (inputArguments.selection == 'MC_stealth_t6')):
    cuts_masked = set(["HLTSelection", "overlap"])
else:
    sys.exit('Unexpected selection type: {s}'.format(s=inputArguments.selection))

cut_descriptions = {
    'HLTSelection': r'diphoton HLT',
    'MCGenInformation': r'MC gen cuts',
    'overlap': r'overlap removal',
    'doublePhoton': r'$N_{\mathrm{medium}~\gamma} = 2$',
    'invariantMass': r'$m_{\gamma \gamma} \geq 90~\GeV$',
    'NJets': r'$N_{\mathrm{jets}} \geq 2$',
    'ST': r'$\st \geq 1000~\GeV$'
}
if ((inputArguments.selection == 'MC_stealth_t5') or (inputArguments.selection == 'MC_stealth_t6')):
    cut_descriptions['HLTSelection'] = r'HLT (emulated)'

def parseCutFlow(cut_flow_lines):
    cut_flow_raw = dict()
    if not(len(cut_flow_lines) == 2+len(event_cuts)):
        sys.exit("Unexpected number of lines produced by printCutFlow: {n}".format(n=len(cut_flow_lines)))
    cut_flow_raw['N_analyzed'] = int((cut_flow_lines[0]).strip())
    cut_flow_raw['N_selected'] = int((cut_flow_lines[1]).strip())
    cut_flow_lines = cut_flow_lines[2:]
    for cut_index, cut in enumerate(event_cuts):
        cuts_split = ((cut_flow_lines[cut_index]).strip()).split()
        if (not(len(cuts_split) == 4)):
            sys.exit('cutflow line {l} in unexpected format.'.format(l=cut_flow_lines[cut_index]))
        if (not(cuts_split[0] == cut)):
            sys.exit("Expected {c1}, got {c2}.".format(c1=cut, c2=cuts_split[0]))
        cut_flow_raw[cut] = dict()
        cut_flow_raw[cut]['N_passing_cut'] = int(cuts_split[1])
        cut_flow_raw[cut]['N_passing_all_cuts_upto'] = int(cuts_split[2])
        cut_flow_raw[cut]['N_passing_all_cuts_besides'] = int(cuts_split[3])
    return cut_flow_raw

def getFloatLatex(r):
    r_str = '{r:.5e}'.format(r=r)
    r_str_split = r_str.split('e')
    base = float(r_str_split[0])
    exponent = int(r_str_split[1])
    return r'{b:.3f} \times 10^{{{e}}}$'.format(b=base, e=exponent)

def getStringRepFraction(num, denom, forceFloatLatex):
    fraction_str = r'$\frac{' + str(num) + r'}{' + str(denom) + r'} = '
    ratio = num/denom
    value_str = r'{r:.2f}\%$'.format(r=100*ratio)
    if (ratio < 0.001):
        value_str = getFloatLatex(ratio)
    return fraction_str + value_str

def getStringRepsForCut(cut_flow_raw, cut):
    global_str = getStringRepFraction(cut_flow_raw[cut]['N_passing_cut'], cut_flow_raw['N_analyzed'], True)
    n_minus_1_str = getStringRepFraction(cut_flow_raw['N_selected'], cut_flow_raw[cut]['N_passing_all_cuts_besides'], False)
    cumul_denom = cut_flow_raw['N_analyzed']
    if not(cut == event_cuts[0]):
        cumul_denom = cut_flow_raw[event_cuts[event_cuts.index(cut)-1]]['N_passing_all_cuts_upto']
    cumul_str = getStringRepFraction(cut_flow_raw[cut]['N_passing_all_cuts_upto'], cumul_denom, False)
    return (global_str, n_minus_1_str, cumul_str)

def writeCutFlowToTex(output_handle, cut_flow_raw):
    output_handle.write(r'\begin{tabular}{|m{0.18\textwidth}|m{0.25\textwidth}|m{0.25\textwidth}|m{0.25\textwidth}|@{}m{0pt}@{}}' + '\n')
    output_handle.write(r'  \hline' + '\n')
    output_handle.write(r'  Event Cut & $E_{\mathrm{global}}$ & $E_{\mathrm{N}-1}$ & $E_{\mathrm{cumulative}}$' + r'\\ \hline' + '\n')
    for cut in event_cuts:
        if cut in cuts_masked:
            continue
        global_str, n_minus_1_str, cumul_str = getStringRepsForCut(cut_flow_raw, cut)
        output_handle.write('  {cd} & {g} & {n} & {c}'.format(cd=cut_descriptions[cut], g=global_str, n=n_minus_1_str, c=cumul_str) + r'&\\[10pt] \hline' + '\n')
    output_handle.write(r'\end{tabular}' + '\n')

optional_identifier = ""
if (inputArguments.optionalIdentifier != ""): optional_identifier = "_{oI}".format(oI=inputArguments.optionalIdentifier)

analysisOutputDirectory = "{aR}/analysis{oI}".format(aR=stealthEnv.analysisRoot, oI=optional_identifier)
output_folder = "{aOD}/cutFlow".format(aOD=analysisOutputDirectory)

if not(os.path.isdir("{oF}".format(oF=output_folder))): subprocess.check_call("mkdir -p {oF}".format(oF=output_folder), shell=True, executable="/bin/bash")

year_prefix = "{eP}/{sER}/selections/combined_DoublePhoton{oI}/merged_selection_{s}_".format(eP=stealthEnv.EOSPrefix, oI=optional_identifier, sER=stealthEnv.stealthEOSRoot, s=inputArguments.selection)
year_suffix = "_signal.root"
target_raw = "{of}/cut_flow_raw_{s}.txt".format(of=output_folder, s=inputArguments.selection)
subprocess.check_call("./eventSelection/bin/printCutFlow {yp}2016{ys} {yp}2017{ys} {yp}2018{ys} {t}".format(yp=year_prefix, ys=year_suffix, t=target_raw), shell=True, executable="/bin/bash")
cut_flow_raw = None
with open(target_raw) as cut_flow_raw_file_handle:
    cut_flow_raw = parseCutFlow(cut_flow_raw_file_handle.readlines())

target = "{of}/cut_flow_{s}.tex".format(of=output_folder, s=inputArguments.selection)
with open(target, 'w') as cut_flow_output_file_handle:
    writeCutFlowToTex(cut_flow_output_file_handle, cut_flow_raw)

print("All done!")
