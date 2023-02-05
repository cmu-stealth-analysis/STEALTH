from __future__ import print_function, division

import argparse, sys, math
import ROOT

QUANTILE_TOLERANCE=0.001
CORRELATION_INTERESTING_THRESHOLD = 0.1

def get_nEvts_from_file(fileName, printDebug=False):
    ggIn = ROOT.TChain("ggNtuplizer/EventTree")
    ggIn.Add(fileName)
    nEntries = ggIn.GetEntries()
    if (printDebug): print("Entries in file \"{f}\": {e}".format(f=fileName, e=nEntries))
    return (nEntries)

def get_nEvts_from_fileList(inputFilesList, printDebug=False):
    listOfInputFiles = []
    inputFileNamesFileObject = open(inputFilesList, 'r')
    totalNEvents = 0
    for inputFileName in inputFileNamesFileObject:
        nEntriesInFile = get_nEvts_from_file(fileName=inputFileName.strip(), printDebug=printDebug)
        if (printDebug): print("Found {n} entries in file {f}.".format(n=nEntriesInFile, f=inputFileName))
        totalNEvents += nEntriesInFile
    inputFileNamesFileObject.close()
    return totalNEvents

def get_nEvts_from_fileList_check(inputFilesList, printDebug=False):
    listOfInputFiles = []
    inputFileNamesFileObject = open(inputFilesList, 'r')
    for inputFileName in inputFileNamesFileObject:
        listOfInputFiles.append(inputFileName.strip())
    inputFileNamesFileObject.close()

    # Load input TTrees into TChain
    ggIn = ROOT.TChain("ggNtuplizer/EventTree")

    for inputFile in listOfInputFiles:
        if (printDebug): print("Adding file {f} to chain.".format(f=inputFileName))
        ggIn.Add(inputFile)

    totalNEntries = ggIn.GetEntries()
    return totalNEntries

def get_number_of_lines_in_file(inputFilePath):
    fileObject = open(inputFilePath, 'r')
    nLines = len(fileObject.readlines())
    fileObject.close()
    return nLines

def read_limits_from_combine_output(combineOutputFilePath):
    combineOutputFile = None
    try:
        combineOutputFile = ROOT.TFile.Open(combineOutputFilePath, "READ")
        if ((combineOutputFile.IsZombie() == ROOT.kTRUE) or not(combineOutputFile.IsOpen() == ROOT.kTRUE)):
            # sys.exit("Error in opening file: {cOFP}".format(cOFP=combinOutputFilePath))
            raise ValueError
    except:
        raise ValueError
    limitTree = ROOT.TTree()
    combineOutputFile.GetObject("limit", limitTree)
    nEntriesFound = limitTree.GetEntries()
    if not(nEntriesFound == 6):
        raise ValueError
    limitsList = []
    for entryIndex in range(0, nEntriesFound):
        nBytesRead = limitTree.GetEntry(entryIndex)
        if (nBytesRead <= 0): sys.exit("ERROR in reading limit tree at index = {i}".format(i=entryIndex))
        upperLimit = limitTree.limit
        quantile = limitTree.quantileExpected
        limitsList.append((quantile, upperLimit))
    combineOutputFile.Close()
    return limitsList

def get_limit_with_quantile(limitsList, targetQuantile):
    for limitQuantileLine in limitsList:
        quantile = limitQuantileLine[0]
        if (abs((quantile/targetQuantile) - 1.0) < 0.001):
            return limitQuantileLine[1]
    sys.exit("ERROR: Unable to find limit at target quantile = {q}. limitsList: {lL}".format(q=targetQuantile, lL=limitsList))

def get_expected_and_observed_limits_from_combine_output(combineOutputFilePath):
    limitsList = read_limits_from_combine_output(combineOutputFilePath)
    expectedUpperLimit = get_limit_with_quantile(limitsList, 0.5)
    expectedUpperLimitTwoSigmaDown = get_limit_with_quantile(limitsList, 0.025)
    expectedUpperLimitOneSigmaDown = get_limit_with_quantile(limitsList,  0.16)
    expectedUpperLimitOneSigmaUp   = get_limit_with_quantile(limitsList,  0.84)
    expectedUpperLimitTwoSigmaUp   = get_limit_with_quantile(limitsList, 0.975)
    observedUpperLimit = get_limit_with_quantile(limitsList, -1.0)
    return (expectedUpperLimit, expectedUpperLimitTwoSigmaDown, expectedUpperLimitOneSigmaDown, expectedUpperLimitOneSigmaUp, expectedUpperLimitTwoSigmaUp, observedUpperLimit)

def get_observed_limit_from_combine_output(combineOutputFilePath):
    return (get_expected_and_observed_limits_from_combine_output(combineOutputFilePath)[5])

def get_best_fit_from_MultiDim_output(multiDimOutputFilePath):
    inputFile=ROOT.TFile.Open("{mDOFP}".format(mDOFP=multiDimOutputFilePath), "READ")
    if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)):
        sys.exit("Error in opening file: {mDOFP}".format(mDOFP=multiDimOutputFilePath))
    limitTree = ROOT.TTree()
    inputFile.GetObject("limit", limitTree)
    nEntriesFound = limitTree.GetEntries()
    if not(nEntriesFound == 1):
        inputFile.Close()
        raise ValueError
        # sys.exit("Error: multidim fit output not in expected format.")
    limitTree.GetEntry(0)
    bestFitValue = limitTree.r
    inputFile.Close()
    return bestFitValue

def get_best_fits_from_MultiDim_fitResult(multiDimFitResultFilePath, parameter_names):
    inputFile=ROOT.TFile.Open("{mDFRFP}".format(mDFRFP=multiDimFitResultFilePath), "READ")
    if not(inputFile): raise ValueError
    if ((inputFile.IsZombie() == ROOT.kTRUE) or not(inputFile.IsOpen() == ROOT.kTRUE)):
        sys.exit("Error in opening file: {mDFRFP}".format(mDFRFP=multiDimOutputFilePath))
    outputDict = {}
    fitResult = ROOT.RooFitResult()
    inputFile.GetObject("fit_mdf", fitResult)
    for parameter_name in parameter_names:
        try:
            outputDict[parameter_name] = fitResult.floatParsFinal().find(parameter_name).getVal()
        except:
            raise ValueError
    inputFile.Close()
    return outputDict

def write_ten_times_expected_upper_limit_from_combine_output_to_file(combineOutputFilePath, outputFilePath):
    """
    This is used in the combine helper script in order to set rmax while calling MultiDimFit.
    """
    outputFileHandle = open(outputFilePath, 'w')
    try:
        expectedUpperLimit, unused1_, expectedUpperLimitOneSigmaDown, expectedUpperLimitOneSigmaUp, unused2_, observedUpperLimit = get_expected_and_observed_limits_from_combine_output(combineOutputFilePath)
        outputFileHandle.write("{l:.4f}\n".format(l=10.0*expectedUpperLimit))
    except:
        outputFileHandle.write("unavailable\n")
    outputFileHandle.close()

def get_r_from_fit_diagnostics_output(input_file_path, printDebug):
    input_file_handle = ROOT.TFile.Open(input_file_path, "READ")
    if ((input_file_handle.IsZombie() == ROOT.kTRUE) or not(input_file_handle.IsOpen() == ROOT.kTRUE)): sys.exit("ERROR: unable to open file at location {l}".format(l=input_file_path))
    input_folder_name = None
    input_tree_handle = ROOT.TTree()
    input_file_handle.GetObject("tree_fit_sb", input_tree_handle)
    n_tree_entries = input_tree_handle.GetEntries()
    if printDebug:
        print("Found tree \"tree_fit_sb\" with {n} entries.".format(n=n_tree_entries))
    if not(n_tree_entries == 1): sys.exit("Unable to get best-fit r from diagnostics output because it is in an unexpected format.")
    n_bytes_read = input_tree_handle.GetEntry(0)
    if (n_bytes_read <= 0): sys.exit("ERROR in reading tree \"tree_fit_sb\" at index = {i}".format(i=entryIndex))
    r_bestfit = float(input_tree_handle.r)
    r_error_lo = float(input_tree_handle.rLoErr)
    r_error_hi = float(input_tree_handle.rHiErr)
    if printDebug:
        print("r_bestfit: {r} + {rerror_up} - {rerror_down}".format(r=r_bestfit, rerror_up=r_error_hi, rerror_down = r_error_lo))
    input_file_handle.Close()
    return tuple(list([r_bestfit, r_error_lo, r_error_hi]))

def write_diffNuisances_output_into_tex_file(input_file_path, output_file_path, printDebug=False):
    input_file_handle = open(input_file_path, 'r')
    output_file_handle = open(output_file_path, 'w')
    output_file_handle.write("\n")
    output_file_handle.write("\\begin{scriptsize}\n")
    output_file_handle.write("\\begin{verbatim}\n")
    begin_reading = False
    for line in input_file_handle:
        line_contents = (line.strip()).split()
        if printDebug: print("Found line: {l}".format(l=line_contents))
        if ((len(line_contents) > 0) and (line_contents[0] == "name")): begin_reading = True
        if not(begin_reading): continue
        output_file_handle.write(line)
    output_file_handle.write("\\end{verbatim}\n")
    output_file_handle.write("\\end{scriptsize}\n")
    output_file_handle.close()
    input_file_handle.close()

def convert_covariance_to_correlation(input_th2d_covariance, target_name):
    output_th2d = input_th2d_covariance.Clone(target_name)
    output_th2d.SetTitle("Bin-by-bin correlations")
    output_th2d.GetXaxis().SetTitle("")
    output_th2d.GetYaxis().SetTitle("")
    for bin_index_x in range(1, 1+output_th2d.GetXaxis().GetNbins()):
        for bin_index_y in range(1, 1+output_th2d.GetYaxis().GetNbins()):
            cov_ij = input_th2d_covariance.GetBinContent(bin_index_x, bin_index_y)
            cov_ii = input_th2d_covariance.GetBinContent(bin_index_x, bin_index_x)
            cov_jj = input_th2d_covariance.GetBinContent(bin_index_y, bin_index_y)
            corr_ij = cov_ij/math.sqrt(cov_ii*cov_jj)
            output_th2d.SetBinContent(bin_index_x, bin_index_y, corr_ij)
    return output_th2d

def print_and_save_high_res_correlations(input_file_path, output_folder, suffix, list_correlations_to_save):
    input_file_handle = ROOT.TFile.Open(input_file_path, "READ")
    if ((input_file_handle.IsZombie() == ROOT.kTRUE) or not(input_file_handle.IsOpen() == ROOT.kTRUE)): sys.exit("ERROR: unable to open file at location {l}".format(l=input_file_path))
    th2s_to_fetch_info = {
        "correlation_b": ("correlation", "covariance_fit_b"),
        "correlation_s": ("correlation", "covariance_fit_s"),
        "correlation_bins_b": ("covariance", "shapes_fit_b/overall_total_covar"),
        "correlation_bins_s": ("covariance", "shapes_fit_s/overall_total_covar"),
    }
    for th2_to_fetch_name in list_correlations_to_save:
        correlations1D = ROOT.TH1D("hist1D_" + th2_to_fetch_name, "Correlation values", 2000, -1.0, 1.0)
        th2_to_fetch_info = th2s_to_fetch_info[th2_to_fetch_name]
        th2_to_fetch_type = th2_to_fetch_info[0]
        th2_to_fetch_path = th2_to_fetch_info[1]
        input_th2 = None
        if (th2_to_fetch_type == "correlation"):
            input_th2 = ROOT.TH2D()
            input_file_handle.GetObject(th2_to_fetch_path, input_th2)
        elif (th2_to_fetch_type == "covariance"):
            input_th2d_covariance = ROOT.TH2D()
            input_file_handle.GetObject(th2_to_fetch_path, input_th2d_covariance)
            input_th2 = convert_covariance_to_correlation(input_th2d_covariance, th2_to_fetch_name)
        interesting_correlations_output_file_handle = open("{o}/interesting_correlations_{n}_{s}.tex".format(o=output_folder, n=th2_to_fetch_name, s=suffix), 'w')
        interesting_correlations_output_file_handle.write("\n")
        interesting_correlations_output_file_handle.write("\\begin{scriptsize}\n")
        interesting_correlations_output_file_handle.write("\\begin{verbatim}\n")
        for bin_index_x in range(1, 1+input_th2.GetXaxis().GetNbins()):
            label_x = input_th2.GetXaxis().GetBinLabel(bin_index_x)
            for bin_index_y in range(1, 1+bin_index_x):
                label_y = input_th2.GetYaxis().GetBinLabel(bin_index_y)
                correlation = input_th2.GetBinContent(bin_index_x, bin_index_y)
                correlations1D.Fill(correlation)
                is_mundane = ((math.fabs(correlation) < CORRELATION_INTERESTING_THRESHOLD) or (math.fabs(correlation - 1.0) < CORRELATION_INTERESTING_THRESHOLD))
                if not(is_mundane): interesting_correlations_output_file_handle.write("({nx}, {ny}) --> {c:.3f}\n".format(nx=label_x, ny=label_y, c=correlation))
        interesting_correlations_output_file_handle.write("\\end{verbatim}\n")
        interesting_correlations_output_file_handle.write("\\end{scriptsize}\n")
        interesting_correlations_output_file_handle.close()
        output_canvas_correlations1D = ROOT.TCanvas("c_corr1D_" + th2_to_fetch_name, "c_corr1D_" + th2_to_fetch_name)
        ROOT.gStyle.SetOptStat(110010)
        correlations1D.Draw()
        output_canvas_correlations1D.SaveAs("{o}/correlations1D_{n}_{s}.pdf".format(o=output_folder, n=th2_to_fetch_name, s=suffix))
        output_canvas = ROOT.TCanvas("c_" + th2_to_fetch_name, "c_" + th2_to_fetch_name, 1920, 1920)
        input_th2.GetXaxis().SetLabelSize(0.01)
        input_th2.GetYaxis().SetLabelSize(0.01)
        ROOT.gPad.Update()
        ROOT.gStyle.SetOptStat(0)
        ROOT.gPad.Update()
        ROOT.gPad.SetBottomMargin(0.18)
        ROOT.gPad.SetLeftMargin(0.18)
        ROOT.gPad.SetTopMargin(0.1)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.Update()
        input_th2.Draw("colz")
        ROOT.gPad.Update()
        right_edge = ROOT.gPad.GetX2()
        bottom_edge = ROOT.gPad.GetY1()
        top_edge = ROOT.gPad.GetY2()
        palette_axis = input_th2.GetListOfFunctions().FindObject("palette")
        palette_axis.SetY1NDC(0.18)
        palette_axis.SetY2NDC(0.9)
        palette_axis.SetX1NDC(0.85)
        palette_axis.SetX2NDC(0.9)
        palette_axis.SetLabelSize(0.025)
        ROOT.gPad.Update()
        output_canvas.SaveAs("{o}/{n}_{s}_high_res.pdf".format(o=output_folder, n=th2_to_fetch_name, s=suffix))
    input_file_handle.Close()
